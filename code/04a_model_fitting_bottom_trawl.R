# remove.packages("sdmTMB")
# Needs sdmTMB main with profile_breakpt():
# remotes::install_github("pbs-assess/sdmTMB", dependencies = TRUE)
library(sdmTMB)
library(sp)
library(dplyr)
library(tidyr)
library(pals)
library(purrr)
library(gt)
library(openxlsx2)
library(here)

here()

#Load functions
source("code/helper_funs.R")

#Load data
files <- list.files(path = "data/processed_data/fish_filtered", pattern = ".rds", full.names=T)
dat <- map(files,readRDS)
dat <- bind_rows(dat)

#remove iphc
remove_iphc <- T
if(remove_iphc){
  dat <- filter(dat, survey!="iphc")
}

#Remove weird depths
dat <- filter(dat, depth>0)

#Remove oxygen outliers
dat <- filter(dat, O2_umolkg<1500)

dat <- dat  %>%
  drop_na(depth,year, mi1,mi2,mi3, X, Y, catch_weight)


#Combine NBS & EBS into a BS region?
combine_bs <- F
if(combine_bs){
  dat$region <- ifelse((dat$region=="ebs"|dat$region=="nbs"), "bs", dat$region)
}

#Remove NBS 
remove_nbs <- T
if(remove_nbs){
  dat <- filter(dat, region!="nbs")
}

#Remove Aleutian Islands
remove_ai <- T
if(remove_ai){
  dat <- filter(dat, region!="ai")
}

##Species to do--full list
species <- c(unique(dat$common_name))

#Remove catch outliers?
remove_outlier <- T

#Years
dat <- filter(dat, year>=2003)

##Clean up survey names
dat <- dat %>%
  mutate(across('survey', str_replace, 'EBS', 'afsc_bsai'),across('survey', str_replace, 'GOA', 'afsc_goa') )

#Spatio-temporal variation
spatio_temp <- F

#Filter depths?
filter_depth <- T

#Name of output folder
output_folder <- "output"

#Fit models 1-6? (Quadratic depth, not cubic depth)
quad_depth_m1_6 <- F # not used downstream (05/06 use models 7, 8, 13-15)

#Fit models 
cub_depth_m7_16 <- T

#Fit models 9-12 and 16? Not used downstream
fit_unused <- F

#depth & latitude info from species table
species_table <- read_excel("data/species_table.xlsx")
species_table$common_name <- tolower(species_table$common_name)

#Just in case, remove duplicates
dat <- unique(dat)

#fit models
for (h in 1:length(species)){
  this_species = species[h]
  print(this_species)
  dat.2.use <- dplyr::filter(dat, common_name == this_species)
  
  #Years for adding extra_time later
  min_year <- min(dat.2.use$year)
  max_year <- max(dat.2.use$year)
  
  #Rename column
  dat.2.use$catch <- dat.2.use$catch_weight
  
  #Get depth habitat
  if(filter_depth){
    closest_value <- species_table$depth[species_table$common_name==this_species]
    dat.2.use <- filter(dat.2.use, depth<closest_value)
  }
  ##Separate out regions for 
  #Make list to store filtered data
  subs <- list()
  regions <- unique(dat.2.use$region)
  for(i in 1:length(regions)){
    region.2.use <- regions[i]
    sub <- filter(dat.2.use, region==region.2.use)
    num <- filter(sub, catch>0)
    #Only include if there are more than 50 observations in the region
    if(nrow(num)>50){
      subs[[i]] <- sub
      names(subs)[[i]] <- region.2.use
    }
  }
  #Remove elements with no obs
  subs <- subs[(names(subs)=="cc"|names(subs)=="bc"|names(subs)=="ebs"|names(subs)=="goa")]
  
  #Add global dataset
  dat.global <- bind_rows(subs)
  if(length(unique(dat.global$region))>1){
  subs[[length(subs)+1]] <- dat.global
  }
  
  ##Add variables and clean up datasets
  dats <- list()
  for(j in 1:length(subs)){
    sub <- subs[[j]]
    #Quadratic and scaled temp
    sub$temp_scaled <- scale(sub$temperature_C)
    sub$temp_scaled2 <- sub$temp_scaled^2
    
    #Quadratic and scaled depth
    sub$log_depth_scaled <- scale(log(sub$depth))
    sub$log_depth_scaled2 <- sub$log_depth_scaled^2
    sub$log_depth_scaled3 <- sub$log_depth_scaled^3
    
    #Re-scale MIs
    sub$mi1_s <- scale(sub$mi1)
    sub$mi2_s <- scale(sub$mi2)
    sub$mi3_s <- scale(sub$mi3)
    sub$po2_s <- scale(sub$po2)
    
    #Remove catch outliers
    if(remove_outlier==T){
      #Remove outliers = catch > 10 sd above the mean
      sub$catch_s <- scale(sub$catch)
      sub <- dplyr::filter(sub, catch_s <=20)
    }
    dats[[j]] <- sub
  }
  
  #Create objects to test
  #List data names
  dat_names <- unique(dat.global$region)
  if(length(unique(dat.global$region))>1){
  dat_names[length(dat_names)+1] <- "coastwide"
  }
  
  for(i in 1:length(dats)){
    sub <- dats[[i]]
    sub <- unique(sub)
    print(dat_names[i])
    #Make mesh
    bnd <- INLA::inla.nonconvex.hull(cbind(sub$X, sub$Y), 
                                     convex = -0.05)
    inla_mesh <- INLA::inla.mesh.2d(
      boundary = bnd,
      max.edge = c(150, 1000),
      offset = -0.1, # default -0.1
      cutoff = 50,
      min.angle = 5 # default 21
    )
    spde <- make_mesh(sub, c("X", "Y"), mesh = inla_mesh)
    
    priors <- sdmTMBpriors(
      matern_s = pc_matern(
        range_gt = 50, range_prob = 0.05, #A value one expects the range is greater than with 1 - range_prob probability.
        sigma_lt = 25, sigma_prob = 0.05 #A value one expects the marginal SD (sigma_O or sigma_E internally) is less than with 1 - sigma_prob probability.
      ),
      matern_st = pc_matern(
        range_gt = 50, range_prob = 0.05,
        sigma_lt = 25, sigma_prob = 0.05
      ),
      #  ar1_rho = normal(0.7,0.1),
      #tweedie_p = normal(1.5,0.2)
    )
    # refactor to avoid identifiability errors
    sub$region <- as.factor(as.character(sub$region))
    sub$year <- as.factor(as.character(sub$year))
    
    saveRDS(sub, file = paste0("output/",output_folder, "/", this_species, "_", dat_names[i], "_dat.rds"))
    
    print(paste(this_species))
    if(quad_depth_m1_6){
    # Model 1: null
    print(paste(this_species))
    print("fitting m1")
    if(dat_names[i]!="coastwide"){
      formula =   "catch ~ -1 + year+log_depth_scaled+ log_depth_scaled2"
      st_type="iid"
    } else {
      formula = "catch ~ -1 + year+region+log_depth_scaled+ log_depth_scaled2"
      if(!spatio_temp){
        st_type <- "off"
        formula = "catch ~ -1 + year+ region+log_depth_scaled+ log_depth_scaled2"
      } else{
        st_type="iid"
        formula = "catch ~ -1 + year+region +log_depth_scaled+ log_depth_scaled2"
      }
    }
    start = Sys.time()
    m1 <- try(sdmTMB(
      formula = as.formula(formula),
      mesh = spde,
      time = "year",
      family = tweedie(link = "log"),
      data = sub,
      priors = priors,
      share_range = TRUE,
      spatial = "on",
      spatiotemporal = paste(st_type),
      control = sdmTMBcontrol(multiphase = FALSE, 
                              newton_loops = 3,
                              nlminb_loops = 2)
      # extra_time = (min_year:max_year)[which(min_year:max_year %in% unique(sub$year) == FALSE)]
    ))
    print( Sys.time() - start )
    saveRDS(m1, file = paste0("output/",output_folder, "/", this_species, "_", dat_names[i], "_model1.rds"))
    
    # Model 2: quadratic temp (uniform across regions)
    print("fitting m2")
    if(dat_names[i]!="coastwide"){
      formula =  "catch ~ -1+year+temp_scaled + temp_scaled2 + log_depth_scaled+ log_depth_scaled2"
      st_type="iid"
    } else {
      if(!spatio_temp){
        st_type <- "off"
        formula ="catch ~ -1 + region+year+temp_scaled + temp_scaled2 + log_depth_scaled+ log_depth_scaled2"
      } else{
        st_type="iid"
        formula ="catch ~ -1 + region+year+temp_scaled + temp_scaled2 + log_depth_scaled+ log_depth_scaled2"
      }
    }
    start = Sys.time()
    m2 <- try(sdmTMB(
      formula = as.formula(formula),
      mesh = spde,
      time = "year",
      family = tweedie(link = "log"),
      data = sub,
      priors = priors,
      share_range = TRUE,
      spatial = "on",
      spatiotemporal = paste(st_type),
      control = sdmTMBcontrol(multiphase = FALSE,
                              newton_loops = 3,
                              nlminb_loops = 2)
      # extra_time = (min_year:max_year)[which(min_year:max_year %in% unique(sub$year) == FALSE)]
    ))
    print( Sys.time() - start )
    saveRDS(m2, file = paste0("output/",output_folder, "/", this_species,  "_", dat_names[i],"_model2.rds"))
    
    # Model 3: breakpoint MI low
    print("fitting m3")
    if(dat_names[i]!="coastwide"){
      formula =   "catch ~ -1 +year+breakpt(mi1_s)+ log_depth_scaled+ log_depth_scaled2"
      st_type="iid"
    } else {
      if(!spatio_temp){
        st_type <- "off"
        formula = "catch ~ -1 + region+year +breakpt(mi1_s)+ log_depth_scaled+ log_depth_scaled2"
      } else{
        st_type="iid"
        formula = "catch ~ -1 + region +year+breakpt(mi1_s)+ log_depth_scaled+ log_depth_scaled2"
      }
    }
    start = Sys.time()
    m3 <-try(profile_breakpt(
      formula = as.formula(formula),
      mesh = spde,
      time = "year",
      family = tweedie(link = "log"),
      data = sub,
      priors = priors,
      share_range = TRUE,
      spatial = "on",
      spatiotemporal = paste(st_type),
      control = sdmTMBcontrol(multiphase = FALSE,
                              newton_loops = 3,
                              nlminb_loops = 2)
      #extra_time = (min_year:max_year)[which(min_year:max_year %in% unique(sub$year) == FALSE)]
    ))
    print( Sys.time() - start )
    saveRDS(m3, file = paste0("output/",output_folder, "/", this_species, "_", dat_names[i], "_model3.rds"))
    
    # Model 4: Breakpoint(mi median)
    print("fitting m4")
    start = Sys.time()
    if(dat_names[i]!="coastwide"){
      formula =   "catch ~ -1+year +breakpt(mi2_s)+ log_depth_scaled+ log_depth_scaled2"
      st_type="iid"
    } else {
      if(!spatio_temp){
        st_type <- "off"
        formula =  "catch ~ -1 + region+year +breakpt(mi2_s)+ log_depth_scaled+ log_depth_scaled2"
      } else{
        st_type="iid"
        formula =  "catch ~ -1 + region +year+breakpt(mi2_s)+ log_depth_scaled+ log_depth_scaled2"
      }
    }
    m4 <-try(profile_breakpt(
      formula = as.formula(formula),
      mesh = spde,
      time = "year",
      family = tweedie(link = "log"),
      data = sub,
      priors = priors,
      share_range = TRUE,
      spatial = "on",
      spatiotemporal = paste(st_type),
      control = sdmTMBcontrol(multiphase = FALSE,
                              newton_loops = 3,
                              nlminb_loops = 2)
      # extra_time = (min_year:max_year)[which(min_year:max_year %in% unique(sub$year) == FALSE)]
    ))
    print( Sys.time() - start )
    saveRDS(m4, file = paste0("output/",output_folder, "/",this_species, "_", dat_names[i], "_model4.rds"))
    
    # Model 5: breakpoint(mi high)
    print("fitting m5")
    start = Sys.time()
    if(dat_names[i]!="coastwide"){
      formula =   "catch ~ -1+year +breakpt(mi3_s)+ log_depth_scaled+ log_depth_scaled2"
      st_type="iid"
    } else {
      if(!spatio_temp){
        st_type <- "off"
        formula =  "catch ~ -1 + region+year +breakpt(mi3_s)+ log_depth_scaled+ log_depth_scaled2"
      } else{
        st_type="iid"
        formula =  "catch ~ -1 + region +year+breakpt(mi3_s)+ log_depth_scaled+ log_depth_scaled2"
      }
    }
    m5 <- try(profile_breakpt(
      formula = as.formula(formula),
      mesh = spde,
      time = "year",
      family = tweedie(link = "log"),
      data = sub,
      priors = priors,
      share_range = TRUE,
      spatial = "on",
      spatiotemporal = paste(st_type),
      control = sdmTMBcontrol(multiphase = FALSE,
                              newton_loops = 3,
                              nlminb_loops = 2)
      #extra_time = (min_year:max_year)[which(min_year:max_year %in% unique(sub$year) == FALSE)]
    ))
    print( Sys.time() - start )
    saveRDS(m5, file = paste0("output/",output_folder, "/", this_species, "_", dat_names[i], "_model5.rds"))
   
     # Model 6: breakpoint(mi high)
    print("fitting m6")
    start = Sys.time()
    if(dat_names[i]!="coastwide"){
      formula =   "catch ~ -1 +year+breakpt(po2_s)+ log_depth_scaled+ log_depth_scaled2"
      st_type="iid"
    } else {
      if(!spatio_temp){
        st_type <- "off"
        formula =  "catch ~ -1 + region+year +breakpt(po2_s)+ log_depth_scaled+ log_depth_scaled2"
      } else{
        st_type="iid"
        formula =  "catch ~ -1 + region +year+breakpt(po2_s)+ log_depth_scaled+ log_depth_scaled2"
      }
    }
    m6 <- try(profile_breakpt(
      formula = as.formula(formula),
      mesh = spde,
      time = "year",
      family = tweedie(link = "log"),
      data = sub,
      priors = priors,
      share_range = TRUE,
      spatial = "on",
      spatiotemporal = paste(st_type),
      control = sdmTMBcontrol(multiphase = FALSE,
                              newton_loops = 3,
                              nlminb_loops = 2)
      #extra_time = (min_year:max_year)[which(min_year:max_year %in% unique(sub$year) == FALSE)]
    ))
    print( Sys.time() - start )
    saveRDS(m6, file = paste0("output/",output_folder, "/", this_species, "_", dat_names[i], "_model6.rds"))
    }
    if(cub_depth_m7_16){
    print("fitting m7")
    if(dat_names[i]!="coastwide"){
      formula =   "catch ~ -1 + year+ log_depth_scaled+ log_depth_scaled2+log_depth_scaled3"
        st_type="iid"
    } else {
      formula = "catch ~ -1 + year+region+log_depth_scaled+ log_depth_scaled2+log_depth_scaled3"
      if(!spatio_temp){
        st_type <- "off"
        formula = "catch ~ -1 + year+region+log_depth_scaled+ log_depth_scaled2+log_depth_scaled3"
      } else{
        st_type="iid"
        formula = "catch ~ -1 + year+region+log_depth_scaled+ log_depth_scaled2+log_depth_scaled3"
      }
    }
    start = Sys.time()
    m7 <- try(sdmTMB(
      formula = as.formula(formula),
      mesh = spde,
      time = "year",
      family = tweedie(link = "log"),
      data = sub,
      priors = priors,
      share_range = TRUE,
      spatial = "on",
      spatiotemporal = paste(st_type),
      control = sdmTMBcontrol(multiphase = FALSE, 
                              newton_loops = 3,
                              nlminb_loops = 2)
     # extra_time = (min_year:max_year)[which(min_year:max_year %in% unique(sub$year) == FALSE)]
    ))
    print( Sys.time() - start )
    saveRDS(m7, file = paste0("output/",output_folder, "/", this_species, "_", dat_names[i], "_model7.rds"))
    
    # Model 8: quadratic temp (uniform across regions)
    print("fitting m8")
    if(dat_names[i]!="coastwide"){
      formula =  "catch ~ -1+year+temp_scaled + temp_scaled2 + log_depth_scaled+ log_depth_scaled2+log_depth_scaled3"
      st_type="iid"
    } else {
      if(!spatio_temp){
        st_type <- "off"
        formula ="catch ~ -1 + year+region+temp_scaled + temp_scaled2 + log_depth_scaled+ log_depth_scaled2+log_depth_scaled3"
      } else{
        st_type="iid"
        formula ="catch ~ -1 + year+region+temp_scaled + temp_scaled2 + log_depth_scaled+ log_depth_scaled2+log_depth_scaled3"
      }
    }
    start = Sys.time()
    m8 <- try(sdmTMB(
      formula = as.formula(formula),
      mesh = spde,
      time = "year",
      family = tweedie(link = "log"),
      data = sub,
      priors = priors,
      share_range = TRUE,
      spatial = "on",
      spatiotemporal = paste(st_type),
      control = sdmTMBcontrol(multiphase = FALSE,
                              newton_loops = 3,
                              nlminb_loops = 2)
     # extra_time = (min_year:max_year)[which(min_year:max_year %in% unique(sub$year) == FALSE)]
    ))
    print( Sys.time() - start )
    saveRDS(m8, file = paste0("output/",output_folder, "/", this_species,  "_", dat_names[i],"_model8.rds"))
    
    if(fit_unused){ # models 9-12
    # Model 9: breakpoint MI low
    print("fitting m9")
    if(dat_names[i]!="coastwide"){
      formula =   "catch ~ -1 +year+breakpt(mi1_s)+ log_depth_scaled+ log_depth_scaled2+log_depth_scaled3"
      st_type="iid"
    } else {
      if(!spatio_temp){
        st_type <- "off"
        formula = "catch ~ -1 + region+year +breakpt(mi1_s)+ log_depth_scaled+ log_depth_scaled2+log_depth_scaled3"
      } else{
        st_type="iid"
        formula = "catch ~ -1 + year+ region +breakpt(mi1_s)+ log_depth_scaled+ log_depth_scaled2+log_depth_scaled3"
      }
    }
    start = Sys.time()
    m9 <-try(profile_breakpt(
      formula = as.formula(formula),
      mesh = spde,
      time = "year",
      family = tweedie(link = "log"),
      data = sub,
      priors = priors,
      share_range = TRUE,
      spatial = "on",
      spatiotemporal = paste(st_type),
      control = sdmTMBcontrol(multiphase = FALSE,
                              newton_loops = 2)
      #extra_time = (min_year:max_year)[which(min_year:max_year %in% unique(sub$year) == FALSE)]
    ))
    print( Sys.time() - start )
    saveRDS(m9, file = paste0("output/",output_folder, "/", this_species, "_", dat_names[i], "_model9.rds"))
    
    # Model 10: Breakpoint(mi median)
    print("fitting m10")
    start = Sys.time()
    if(dat_names[i]!="coastwide"){
      formula =   "catch ~ -1 +year+breakpt(mi2_s)+ log_depth_scaled+ log_depth_scaled2+log_depth_scaled3"
      st_type="iid"
    } else {
      if(!spatio_temp){
        st_type <- "off"
        formula =  "catch ~ -1 + year+region +breakpt(mi2_s)+ log_depth_scaled+ log_depth_scaled2+log_depth_scaled3"
      } else{
        st_type="iid"
        formula =  "catch ~ -1 + year+ region +breakpt(mi2_s)+ log_depth_scaled+ log_depth_scaled2+log_depth_scaled3"
      }
    }
    m10 <-try(profile_breakpt(
      formula = as.formula(formula),
      mesh = spde,
      time = "year",
      family = tweedie(link = "log"),
      data = sub,
      priors = priors,
      share_range = TRUE,
      spatial = "on",
      spatiotemporal = paste(st_type),
      control = sdmTMBcontrol(multiphase = FALSE,
                              newton_loops = 2)
     # extra_time = (min_year:max_year)[which(min_year:max_year %in% unique(sub$year) == FALSE)]
    ))
    print( Sys.time() - start )
    saveRDS(m10, file = paste0("output/",output_folder, "/",this_species, "_", dat_names[i], "_model10.rds"))
    
    # Model 11: breakpoint(mi high)
    print("fitting m11")
    start = Sys.time()
    if(dat_names[i]!="coastwide"){
      formula =   "catch ~ -1+year +breakpt(mi3_s)+ log_depth_scaled+ log_depth_scaled2+log_depth_scaled3"
      st_type="iid"
    } else {
      if(!spatio_temp){
        st_type <- "off"
        formula =  "catch ~ -1 + year+region +breakpt(mi3_s)+ log_depth_scaled+ log_depth_scaled2+log_depth_scaled3"
      } else{
        st_type="iid"
        formula =  "catch ~ -1 + year+region+breakpt(mi3_s)+ log_depth_scaled+ log_depth_scaled2+log_depth_scaled3"
      }
    }
    m11 <- try(profile_breakpt(
      formula = as.formula(formula),
      mesh = spde,
      time = "year",
      family = tweedie(link = "log"),
      data = sub,
      priors = priors,
      share_range = TRUE,
      spatial = "on",
      spatiotemporal = paste(st_type),
      control = sdmTMBcontrol(multiphase = FALSE,
                              newton_loops = 3,
                              nlminb_loops = 2)
      #extra_time = (min_year:max_year)[which(min_year:max_year %in% unique(sub$year) == FALSE)]
    ))
    print( Sys.time() - start )
    saveRDS(m11, file = paste0("output/",output_folder, "/", this_species, "_", dat_names[i], "_model11.rds"))
    # Model 12: breakpoint(o2)
    print("fitting m12")
    start = Sys.time()
    if(dat_names[i]!="coastwide"){
      formula =   "catch ~ -1+year +breakpt(po2_s)+ log_depth_scaled+ log_depth_scaled2+log_depth_scaled3"
      st_type="iid"
    } else {
      if(!spatio_temp){
        st_type <- "off"
        formula =  "catch ~ -1 + year+region+breakpt(po2_s)+ log_depth_scaled+ log_depth_scaled2+log_depth_scaled3"
      } else{
        st_type="iid"
        formula =  "catch ~ -1 + year+ region +breakpt(po2_s)+ log_depth_scaled+ log_depth_scaled2+log_depth_scaled3"
      }
    }
    m12 <- try(profile_breakpt(
      formula = as.formula(formula),
      mesh = spde,
      time = "year",
      family = tweedie(link = "log"),
      data = sub,
      priors = priors,
      share_range = TRUE,
      spatial = "on",
      spatiotemporal = paste(st_type),
      control = sdmTMBcontrol(multiphase = FALSE,
                              newton_loops = 3,
                              nlminb_loops = 2)
     # extra_time = (min_year:max_year)[which(min_year:max_year %in% unique(sub$year) == FALSE)]
    ))
    print( Sys.time() - start )
    saveRDS(m12, file = paste0("output/",output_folder, "/", this_species, "_", dat_names[i], "_model12.rds"))
    }
  # Model 13: breakpoint MI low+quad temp
  print("fitting m13")
  if(dat_names[i]!="coastwide"){
    formula =   "catch ~ -1+ year +breakpt(mi1_s)+ +temp_scaled + temp_scaled2+log_depth_scaled+ log_depth_scaled2+log_depth_scaled3"
    st_type="iid"
  } else {
    if(!spatio_temp){
      st_type <- "off"
      formula = "catch ~ -1 + year+region+breakpt(mi1_s)+ +temp_scaled + temp_scaled2+log_depth_scaled+ log_depth_scaled2+log_depth_scaled3"
    } else{
      st_type="iid"
      formula = "catch ~ -1 + year+ region +breakpt(mi1_s)+ +temp_scaled + temp_scaled2+log_depth_scaled+ log_depth_scaled2+log_depth_scaled3"
    }
  }
  start = Sys.time()
  m13 <-try(profile_breakpt(
    formula = as.formula(formula),
    mesh = spde,
    time = "year",
    family = tweedie(link = "log"),
    data = sub,
    priors = priors,
    share_range = TRUE,
    spatial = "on",
    spatiotemporal = paste(st_type),
    control = sdmTMBcontrol(multiphase = FALSE,
                            newton_loops = 3,
                            nlminb_loops = 2)
    #extra_time = (min_year:max_year)[which(min_year:max_year %in% unique(sub$year) == FALSE)]
  ))
  print( Sys.time() - start )
  #Refits at fixed cutpoints for breakpoint uncertainty (used in 06)
  if(!inherits(m13, "try-error")) m13$mix <- try(breakpt_mix(m13))
  saveRDS(m13, file = paste0("output/",output_folder, "/", this_species, "_", dat_names[i], "_model13.rds"))
  
  # Model 14: Breakpoint(mi median)
  print("fitting m14")
  start = Sys.time()
  if(dat_names[i]!="coastwide"){
    formula =   "catch ~ -1 + year +breakpt(mi2_s)+temp_scaled + temp_scaled2+ log_depth_scaled+ log_depth_scaled2+log_depth_scaled3"
    st_type="iid"
  } else {
    if(!spatio_temp){
      st_type <- "off"
      formula =  "catch ~ -1 + year+region +breakpt(mi2_s)+temp_scaled + temp_scaled2+ log_depth_scaled+ log_depth_scaled2+log_depth_scaled3"
    } else{
      st_type="iid"
      formula =  "catch ~ -1 + year + region +breakpt(mi2_s)+temp_scaled + temp_scaled2+ log_depth_scaled+ log_depth_scaled2+log_depth_scaled3"
    }
  }
  m14 <-try(profile_breakpt(
    formula = as.formula(formula),
    mesh = spde,
    time = "year",
    family = tweedie(link = "log"),
    data = sub,
    priors = priors,
    share_range = TRUE,
    spatial = "on",
    spatiotemporal = paste(st_type),
    control = sdmTMBcontrol(multiphase = FALSE,
                            newton_loops = 3,
                            nlminb_loops = 2)
    # extra_time = (min_year:max_year)[which(min_year:max_year %in% unique(sub$year) == FALSE)]
  ))
  print( Sys.time() - start )
  #Refits at fixed cutpoints for breakpoint uncertainty (used in 06)
  if(!inherits(m14, "try-error")) m14$mix <- try(breakpt_mix(m14))
  saveRDS(m14, file = paste0("output/",output_folder, "/", this_species, "_", dat_names[i], "_model14.rds"))
  
  # Model 15: breakpoint(mi high)
  print("fitting m15")
  start = Sys.time()
  if(dat_names[i]!="coastwide"){
    formula =   "catch ~ -1 + year +breakpt(mi3_s)+ temp_scaled + temp_scaled2+log_depth_scaled+ log_depth_scaled2+log_depth_scaled3"
    st_type="iid"
  } else {
    if(!spatio_temp){
      st_type <- "off"
      formula =  "catch ~ -1 + year+region +breakpt(mi3_s)+temp_scaled + temp_scaled2+ log_depth_scaled+ log_depth_scaled2+log_depth_scaled3"
    } else{
      st_type="iid"
      formula =  "catch ~ -1 + year+ region +breakpt(mi3_s)+temp_scaled + temp_scaled2+ log_depth_scaled+ log_depth_scaled2+log_depth_scaled3"
    }
  }
  m15 <- try(profile_breakpt(
    formula = as.formula(formula),
    mesh = spde,
    time = "year",
    family = tweedie(link = "log"),
    data = sub,
    priors = priors,
    share_range = TRUE,
    spatial = "on",
    spatiotemporal = paste(st_type),
    control = sdmTMBcontrol(multiphase = FALSE,
                            newton_loops = 3,
                            nlminb_loops = 2)
    #extra_time = (min_year:max_year)[which(min_year:max_year %in% unique(sub$year) == FALSE)]
  ))
  print( Sys.time() - start )
  #Refits at fixed cutpoints for breakpoint uncertainty (used in 06)
  if(!inherits(m15, "try-error")) m15$mix <- try(breakpt_mix(m15))
  saveRDS(m15, file = paste0("output/",output_folder, "/", this_species, "_", dat_names[i], "_model15.rds"))
  if(fit_unused){ # model 16
  # Model 16: breakpoint(o2)
  print("fitting m16")
  start = Sys.time()
  if(dat_names[i]!="coastwide"){
    formula =   "catch ~ -1+ year +breakpt(po2_s)+ temp_scaled + temp_scaled2+log_depth_scaled+ log_depth_scaled2+log_depth_scaled3"
    st_type="iid"
  } else {
    if(!spatio_temp){
      st_type <- "off"
      formula =  "catch ~ -1 + year+region +breakpt(po2_s)+temp_scaled + temp_scaled2+ log_depth_scaled+ log_depth_scaled2+log_depth_scaled3"
    } else{
      st_type="iid"
      formula =  "catch ~ -1 + year+ region +breakpt(po2_s)+temp_scaled + temp_scaled2+ log_depth_scaled+ log_depth_scaled2+log_depth_scaled3"
    }
  }
  m16 <- try(profile_breakpt(
    formula = as.formula(formula),
    mesh = spde,
    time = "year",
    family = tweedie(link = "log"),
    data = sub,
    priors = priors,
    share_range = TRUE,
    spatial = "on",
    spatiotemporal = paste(st_type),
    control = sdmTMBcontrol(multiphase = FALSE,
                            newton_loops = 3,
                            nlminb_loops = 2)
    # extra_time = (min_year:max_year)[which(min_year:max_year %in% unique(sub$year) == FALSE)]
  ))
  print( Sys.time() - start )
  saveRDS(m16, file = paste0("output/",output_folder, "/", this_species, "_", dat_names[i], "_model16.rds"))
  }
  gc()
    }
}
}

