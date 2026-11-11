##Unit conversions and equations
calc_po2_sat <- function(salinity, temp, depth, oxygen, lat, long, umol_m3, ml_L) {
  # Input:       S = Salinity (pss-78)
  #              T = Temp (deg C) ! use potential temp
  #depth is in meters
  
  #Pena et al. ROMS and GLORYS are in mmol per m^3 (needs to be converted to umol per kg) (o2 from trawl data was in mL/L, so had to do extra conversions)
  gas_const = 8.31
  partial_molar_vol = 0.000032
  kelvin = 273.15
  boltz = 0.000086173324
  
  #this was dumb because it is actually equivalent to umol/kg
  #convert mmol to umol
  umol_m3 <- oxygen*1000
  #convert m3 to l
  umol_l <- umol_m3/1000
  #convert from molality (moles per volume) to molarity (moles per mass)
  #1 L of water = 1 kg of water, so no equation needed?? Right??
  o2_umolkg <- umol_l * 1/1
  
  SA = gsw_SA_from_SP(salinity,depth,long,lat) #absolute salinity for pot T calc
  pt = gsw_pt_from_t(SA,temp,depth) #potential temp at a particular depth
  #this is for if using data that has oxygen in ml/L
  #CT = gsw_CT_from_t(SA,temp,depth) #conservative temp
  #sigma0 = gsw_sigma0(SA,CT)
  # o2_umolkg = oxygen*44660/(sigma0+1000) 
  
  O2_Sat0 = gsw_O2sol_SP_pt(salinity,pt)
  
  #= o2satv2a(sal,pt) #uses practical salinity and potential temp - solubity at p =1 atm
  press = exp(depth*10000*partial_molar_vol/gas_const/(temp+kelvin))
  O2_satdepth = O2_Sat0*press
  
  #solubility at p=0
  sol0 = O2_Sat0/0.209
  sol_Dep = sol0*press
  po2 = o2_umolkg/sol_Dep
  po2 <- po2 * 101.325 # convert to kPa
  return(po2)
  
}

calc_mi <- function(Eo, Ao, W, n,po2, inv.temp, fancy) {
  if(fancy==T){
  mi = W^n*Ao*po2 *exp(Eo * inv.temp)
  } else {
    mi = po2 *exp(Eo * inv.temp)
  }
  return(mi)
}

calc_po2_crit <- function(inv.temp, taxa, mi, body_size, model, fancy) {
  W <- body_size
  inv.temp <- inv.temp
  ###Calculate Metabolic index 
  ##Species parameters from Tim's paper
  #Read table in
  mi_pars <- read.csv("data/taxa_table.csv")
  mi_pars$Group <- tolower(mi_pars$Group)
  mi_pars <- filter(mi_pars,Group==taxa)
  
  V <- mi_pars$logV
  n <- mi_pars$n
  Eo <- c(mi_pars$Eolow, mi_pars$Eo, mi_pars$Eohigh)
  Ao <- 1/exp(V)

  Ao <- 1/exp(V)
  n <- n
  Eo <-  if(model=="model3"|model=="model9") Eo[1] else if(model=="model4"|model=="model10") Eo[2] else Eo[3]
  if(fancy==T){
  po2 = mi/(W^n*Ao*exp(Eo * inv.temp))
  } else{
    po2 = mi/(exp(Eo * inv.temp)) 
  }
  return(po2)
}

# calc o2 solubility, relies on o2 in umol/kg
gsw_O2sol_SP_pt <- function(sal,pt) {
  x = sal
  pt68 = pt*1.00024
  y = log((298.15 - pt68)/(273.15 + pt68))
  
  a0 =  5.80871
  a1 =  3.20291
  a2 =  4.17887
  a3 =  5.10006
  a4 = -9.86643e-2
  a5 =  3.80369
  b0 = -7.01577e-3
  b1 = -7.70028e-3
  b2 = -1.13864e-2
  b3 = -9.51519e-3
  c0 = -2.75915e-7
  
  O2sol = exp(a0 + y*(a1 + y*(a2 + y*(a3 + y*(a4 + a5*y)))) + x*(b0 + y*(b1 + y*(b2 + b3*y)) + c0*x))
  return(O2sol)
}

back.convert <- function(x, mean_orig, sd_orig) {
  x* sd_orig+mean_orig
}

convert_class <- function(x) {
  for (i in 1:ncol(x)) x[,i] <- as(x[,i], Class = "matrix")
  return(x)
}

###DATA PROCESSING ###

print_species <- function(type){
  biomass <- combine_all(type)
  species <- biomass[1:2]
  return(species)
}

length_expand_nwfsc2 <- function(spc, sci_name){
  bio <- readRDS("data/fish_raw/NOAA/nwfsc_bio.rds")
  names(bio) = tolower(names(bio))
  bio$scientific_name <- tolower(bio$scientific_name)
  bio$common_name <- tolower(bio$common_name)
  bio$event_id = as.numeric(bio$trawl_id)
  bio$common_name <- ifelse(bio$common_name=="pacific spiny dogfish", "spiny dogfish", bio$common_name)

 bio_sub = dplyr::filter(bio, common_name == spc)
 #If there is weight or length data present in data, get haul mean of weight
 if(nrow(bio_sub)>0) {
   ##Add weights if only length data:
   # find length-weight relationship parameters for one species
   pars <- rfishbase::length_weight(
     #convert scientific name to first letter capitalized for fishbase
     species_list = stringr::str_to_sentence(sci_name, locale = "en"))
   #Get mean
   a <- mean(pars$a)
   b <- mean(pars$b)
   #For longspine thornyhead, which is not in database for some reason, got these values from fishbase direct page
   if(sci_name=="sebastolobus altivelis"){
     a <- 0.00912
     b=3.09
   }
   ##Assign a and b values if there are none available for the species in fishbase
   if(a=="NaN"){
     a <- 0.00675	
   }
   if(b=="NaN"){
     b <- 3
   }
   #Calculate weight (and convert from g to kg)
   bio_sub <- dplyr::mutate(bio_sub, weight = ifelse(is.na(weight), ((a*length_cm^b)*0.001), weight))
   bio_summ <- group_by(bio_sub, event_id) %>% summarize(haul_weight=median(weight))
   return(bio_summ)
 }
}

length_expand_nwfsc <- function(spc, sci_name) {
  # load, clean, and join data
  bio <- readRDS("data/fish_raw/NOAA/nwfsc_bio.rds")
  catch <- readRDS("data/fish_raw/NOAA/nwfsc_catch.rds")
  names(catch) = tolower(names(catch))
  names(bio) = tolower(names(bio))
  bio$scientific_name <- tolower(bio$scientific_name)
  bio$common_name <- tolower(bio$common_name)
  catch$common_name <- tolower(catch$common_name)
  
  bio$trawl_id = as.character(bio$trawl_id)
  catch$trawl_id=as.character(catch$trawl_id)
  
  #fix dogfish
  catch$common_name <- ifelse(catch$common_name=="pacific spiny dogfish", "spiny dogfish", catch$common_name)
  bio$common_name <- ifelse(bio$common_name=="pacific spiny dogfish", "spiny dogfish", bio$common_name)
  
  #haul$sampling_end_hhmmss = as.numeric(haul$sampling_end_hhmmss)
  #haul$sampling_start_hhmmss = as.numeric(haul$sampling_start_hhmmss)
  
  #Combine data
  dat = dplyr::left_join(filter(bio[,c("trawl_id", "year", "scientific_name", "common_name", "weight", "ageing_lab", "length_cm", "width_cm", "sex", "age")], !is.na(length_cm)), 
                         catch[,c("trawl_id","common_name", "subsample_count","area_swept_ha","longitude_dd", "latitude_dd","subsample_wt_kg","total_catch_numbers","total_catch_wt_kg","cpue_kg_km2")],
                         relationship = "many-to-many")
  
  
  # filter out species of interest from joined (catch/haul/bio) dataset
  dat_sub = dplyr::filter(dat, common_name == spc)
  
  # fit length-weight regression by year to predict fish weights that have lengths only.
  # note a rank-deficiency warning may indicate there is insufficient data for some year/sex combinations (likely for unsexed group)
  #Create one set of data with only year/survey combinations with at least some weight data:
  
  if(nrow(dat_sub)>0) {
    #Add column counting number of length observations per catch
    dat_sub$is_length <- ifelse(!is.na(dat_sub$length_cm), 1,0) 
    dat_sub <- group_by(dat_sub, trawl_id) %>% mutate(nlength=sum(is_length)) %>% ungroup()
    dat_sub$trawl_id <- as.numeric(dat_sub$trawl_id)
    
    fitted <-  filter(dat_sub, !is.na(length_cm)) %>%
      group_by(year) %>%
      mutate(sum = sum(weight, na.rm=T)) %>%
      filter(sum>0) %>%
      ungroup()
    fitted$weight <- ifelse(fitted$weight==0, NA, fitted$weight)
    #And one for those without any weight data:
    not_fitted <-  filter(dat_sub, !is.na(length_cm)) %>%
      group_by(year) %>%
      mutate(sum = sum(weight, na.rm=T)) %>%
      filter(sum==0) %>%
      ungroup()
    
    if(nrow(fitted)>0){
      fitted <-  
        group_nest(fitted, year) %>%
        mutate(
          model = purrr::map(data, ~ lm(log(weight) ~ log(length_cm), data = .x)),
          tidied = purrr::map(model, broom::tidy),
          augmented = purrr::map(model, broom::augment),
          predictions = purrr::map2(data, model, modelr::add_predictions)
        )
      
      # replace missing weights with predicted weights
      dat_pos2 = fitted %>%
        tidyr::unnest(predictions) %>%
        dplyr::select(-data, -model, -tidied, -augmented) %>%
        dplyr::mutate(weight = ifelse(is.na(weight), exp(pred), weight))
      
      dat_pos <- bind_rows(dat_pos2, not_fitted)
    }
    
    if(nrow(fitted)==0 & nrow(not_fitted)>0){
      dat_pos <- not_fitted
    }
    
    if(nrow(fitted)>0 | nrow(not_fitted)>0){
      #If there is no weight data available to get weight-length empirical interpolation, use fishbase to fill in
      # find length-weight relationship parameters for one species
      pars <- rfishbase::length_weight(
        #convert scientific name to first letter capitalized for fishbase
        species_list = stringr::str_to_sentence(sci_name, locale = "en"))
      #Get mean
      a <- mean(pars$a)
      b <- mean(pars$b)
      ##For this species, which is for some reason not in database
      if(sci_name=="sebastolobus altivelis"){
        a <- 0.00912
        b=3.09
      }
      ##Assign a and b values if there are none available for the species in fishbase
      if(a=="NaN"){
        a <- 0.00675	
      }
      if(b=="NaN"){
        b <- 3
      }
      #Make NAs where text
      dat_pos$weight <- ifelse(dat_pos$weight=="NaN", NA, dat_pos$weight)
      #Remove any zero weights
      dat_pos$weight <- ifelse(dat_pos$weight==0, NA, dat_pos$weight)
      #Calculate weight (and convert from g to kg)
      dat_pos <- dplyr::mutate(dat_pos, weight = ifelse(is.na(weight), ((a*length_cm^b)*0.001), weight))
      #Add column getting the mean individual weight 
      dat_test <- group_by(dat_pos, trawl_id) %>% mutate(haul_weight=mean(weight)) %>% ungroup()
      
      trawlids <- unique(dat_pos$trawl_id)
      if(length(trawlids!=0)){
        p <- data.frame(trawl_id = trawlids,
                        p1 = 0,
                        p2 = 0,
                        p3 = 0,
                        p4 = 0)
        
        sizethresholds <- quantile(dat_pos$weight, c(0.15, 0.5, 0.85, 1), na.rm = T)
        for (i in 1:length(trawlids)) {
          haul_sample<- dplyr::filter(dat_pos, trawl_id == trawlids[i])
          if(nrow(haul_sample) > 0 | var(haul_sample$weight >0)) {
            # fit kernel density to weight frequency
            smoothed_w <- KernSmooth::bkde(haul_sample$weight, range.x = c(min(dat_pos$weight), max(dat_pos$weight)), bandwidth = 2)
            # make sure smoother predicts positive or zero density
            smoothed_w$y[smoothed_w$y<0] <- 0
            # calculate proportion by biomass and by number
            p_w_byweight <- smoothed_w$y * smoothed_w$x / sum(smoothed_w$x*smoothed_w$y)
            
            p_w_byweight[p_w_byweight<0] <- 0
            #p_w_bynum[p_w_bynum<0] <- 0
            
            p1 <- sum(p_w_byweight[smoothed_w$x<=sizethresholds[1]])
            p2 <- sum(p_w_byweight[smoothed_w$x>sizethresholds[1] & smoothed_w$x <=sizethresholds[2]])
            p3 <- sum(p_w_byweight[smoothed_w$x>sizethresholds[2] & smoothed_w$x <=sizethresholds[3]])
            p4 <- sum(p_w_byweight[smoothed_w$x>sizethresholds[3]])
            
            p[i,2:5] <- c(p1, p2, p3, p4)
            
          }
          else {
            indx <- which(sizethresholds>haul_sample$weight)
            p[i, min(indx)+1] <- 1
          }
        }
        # add hauls with zero catch back in
        absent = filter(dat_sub, cpue_kg_km2 == 0)
        trawlids <- unique(absent$trawl_id)
        #  absent.df <- data.frame(trawl_id = trawlids,
        #      p1 = 0,
        #     p2 = 0,
        #     p3 = 0,
        #     p4 = 0)
        all_hauls <- p
        # all_hauls <- rbind(p, absent.df)
        all_hauls$trawl_id <- as.numeric(all_hauls$trawl_id)
        dat_sub$median_weight <- median(dat_sub$weight, na.rm=T)
        nlengths <- unique(dat_sub[,c("trawl_id","nlength", "median_weight")])
        meanweight <- unique(dat_test[,c("trawl_id","haul_weight")])
        all_hauls2 <- left_join(all_hauls, nlengths)
        all_hauls2 <- left_join(all_hauls2, meanweight)
        return(all_hauls2)
      }
    }
  }
  if(nrow(dat_sub)>0){
    if(nrow(fitted)==0 & nrow(not_fitted)==0){
      trawlids <- unique(dat_sub$trawl_id)
      absent.df <- data.frame(trawl_id = trawlids,
                              p1 = NA,
                              p2 = NA,
                              p3 = NA,
                              p4 = NA,
                              nlength=0,
                              haul_weight=NA, 
                              median_weight=NA)
      return(absent.df)
    }
  }
  if(nrow(dat_sub)>0){
    if(length(trawlids)==0){
      trawlids <- unique(dat_sub$trawl_id)
      absent.df <- data.frame(trawl_id = trawlids,
                              p1 = NA,
                              p2 = NA,
                              p3 = NA,
                              p4 = NA,
                              nlength=0,
                              haul_weight=NA)
      return(absent.df)
    }
  }
if(nrow(dat_sub)==0){
  return(warning("species not present in data"))
}
}


load_data_nwfsc <- function(spc,sci_name, dat.by.size, length) {
  dat <- readRDS("data/fish_raw/NOAA/nwfsc_catch.rds")
  names(dat) = tolower(names(dat))
  dat$common_name <- tolower(dat$common_name)
  dat$scientific_name <- tolower(dat$scientific_name)
  #Fix dogfish
  dat$common_name <- ifelse(dat$common_name=="pacific spiny dogfish", "spiny dogfish", dat$common_name)
  dat = dplyr::filter(dat, scientific_name == sci_name)
  if(is.data.frame(dat.by.size)){
    dat.by.size$trawl_id <- as.character(dat.by.size$trawl_id)
    dat <- left_join(dat, dat.by.size, by = "trawl_id")
    # remove tows where there was positive catch but no length measurements
  }
  if(!is.data.frame(dat.by.size)){
    dat$nlength <- NA
    dat$median_weight <- NA
    dat$haul_weight <- NA
    dat$p1 <- NA
    dat$p2 <- NA
    dat$p3 <- NA
    dat$p4 <- NA
  }
  # remove tows where there was positive catch but no length measurements
  if(length){
    dat <- dplyr::filter(dat, !is.na(p1))
  }
  # analyze or years and hauls with adequate oxygen and temperature data, within range of occurrence
  
  #O2 from trawl data is in ml/l 
  # just in case, remove any missing or nonsense values from sensors
  # dat <- dplyr::filter(dat, !is.na(o2), !is.na(sal), !is.na(temp), is.finite(sal))
  # dat <- calc_po2_mi(dat)
  # dat <- dplyr::filter(dat, !is.na(temp), !is.na(mi))
  
  # prepare data and models -------------------------------------------------
  dat$longitude <- dat$longitude_dd
  dat$latitude <- dat$latitude_dd
  dat$event_id <- dat$trawl_id
  dat$date <- as.POSIXct(as.Date(as.POSIXct("1970-01-01")+as.difftime(dat$date,units="days")))
  # get julian day
  dat$julian_day <- rep(NA, nrow(dat))
  for (i in 1:nrow(dat)){ 
    dat$julian_day[i] <- as.POSIXlt(dat$date[i], format = "%Y-%b-%d")$yday
  }
  dat <- dplyr::select(dat, event_id, common_name, project, vessel, tow, year, date, longitude_dd, latitude_dd, longitude, latitude, cpue_kg_km2,
                       depth_m, julian_day, nlength, median_weight, haul_weight, pass, p1, p2, p3, p4)
  
  # UTM transformation
  dat_ll = dat
  sp::coordinates(dat_ll) <- c("longitude_dd", "latitude_dd")
  sp::proj4string(dat_ll) <- sp::CRS("+proj=longlat +datum=WGS84")
  # convert to utm with spTransform
  dat_utm = sp::spTransform(dat_ll, 
                            sp::CRS("+proj=utm +zone=10 +datum=WGS84 +units=km"))
  # convert back from sp object to data frame
  dat = as.data.frame(dat_utm)
  dat = dplyr::rename(dat, X = coords.x1,
                      Y = coords.x2)
  dat$scientific_name <- sci_name
  dat$depth <- dat$depth_m
  dat$depth_m <- NULL
  dat$survey <- "nwfsc"
  return(dat)
}

length_expand_afsc2 <- function(spc, sci_name) {
  bio2 <-readRDS("data/fish_raw/NOAA/ak_bts_goa_ebs_nbs_indivero_all_levels.RDS")
  haul <- bio2$haul
  specimen <- bio2$specimen
  species <- bio2$species
  size <- bio2$size
  
  #make lowercase
  names(haul) <- tolower(names(haul))
  names(specimen) <- tolower(names(specimen))
  names(species) <- tolower(names(species))
  names(size) <- tolower(names(bio2$size))
  
  species$species_name <- tolower(species$species_name)
  
  #Combine size and species
  lengths <- dplyr::left_join(size, species)
  lengths <- dplyr::left_join(lengths, haul)
  
  #Expand to make separate row for each measurement
  lengths <- dplyr::filter(lengths, !is.na(frequency))
  lengths2 <- tidyr::uncount(lengths, weights=frequency)
  
  #Combine
  bio <- dplyr::left_join(specimen, species)
  bio <- dplyr::left_join(bio, haul)
  
  #Combine specimen and length data
  bio3 <- dplyr::bind_rows(bio, lengths)
  
  #Convert length from mm to cm
  bio3$length_cm <- bio3$length*0.1
  
  #Convert weight from g to kg
  bio3$weight <- bio3$weight*0.001
  
  #Event_id
  bio3$event_id <-bio3$hauljoin
  
  bio_sub = dplyr::filter(bio3, species_name == sci_name)
  
  #If there is weight or length data present in data, get haul mean of weight
  if(nrow(bio_sub)>0) {
    ##Add weights if only length data:
    # find length-weight relationship parameters for one species
    pars <- rfishbase::length_weight(
      #convert scientific name to first letter capitalized for fishbase
      species_list = stringr::str_to_sentence(sci_name, locale = "en"))
    #Get mean
    a <- mean(pars$a)
    b <- mean(pars$b)
    #For longspine thornyhead, which is not in database for some reason, got these values from fishbase direct page
    if(sci_name=="sebastolobus altivelis"){
      a <- 0.00912
      b=3.09
    }
    ##Assign a and b values if there are none available for the species in fishbase
    if(a=="NaN"){
      a <- 0.00675	
    }
    if(b=="NaN"){
      b <- 3
    }
    #Calculate weight (and convert from g to kg)
    bio_sub <- dplyr::mutate(bio_sub, weight = ifelse(is.na(weight), ((a*length_cm^b)*0.001), weight))
    bio_summ <- group_by(bio_sub, event_id) %>% summarize(haul_weight=median(weight))
    return(bio_summ)
  }
}
length_expand_afsc <- function(spc, sci_name) {
  # load, clean, and join data
  bio2 <-readRDS("data/fish_raw/NOAA/ak_bts_goa_ebs_nbs_indivero_all_levels.RDS")
  catch2 <- readRDS("data/fish_raw/NOAA/ak_bts_goa_ebs_nbs_indivero_cpue_zerofilled.RDS")
  
  #Isolate necessary parts of full data to get specimen weights/lengths per haul
  haul <- bio2$haul
  specimen <- bio2$specimen
  species <- bio2$species
  size <- bio2$size
  
  #make lowercase
  names(haul) <- tolower(names(haul))
  names(specimen) <- tolower(names(specimen))
  names(species) <- tolower(names(species))
  names(size) <- tolower(names(bio2$size))
  
  species$species_name <- tolower(species$species_name)
  
  #Combine size and species
  lengths <- dplyr::left_join(size, species)
  lengths <- dplyr::left_join(lengths, haul)
  
  #Expand to make separate row for each measurement
  lengths <- dplyr::filter(lengths, !is.na(frequency))
  lengths2 <- tidyr::uncount(lengths, weights=frequency)
  
  #Combine
  bio <- dplyr::left_join(specimen, species)
  bio <- dplyr::left_join(bio, haul)
  
  #Combine specimen and length data
  bio3 <- dplyr::bind_rows(bio, lengths)
  
  #Convert length from mm to cm
  bio3$length_cm <- bio3$length*0.1
  
  #Convert weight from g to kg
  bio3$weight <- bio3$weight*0.001
  
  #Combine catch data with species data
  names(catch2) <- tolower(names(catch2))
  catch <- dplyr::left_join(catch2, species)
  
  # filter out species of interest from joined (catch/haul/bio) dataset
  catch_sub = dplyr::filter(catch, species_name == sci_name)
  bio_sub = dplyr::filter(bio3, species_name == sci_name)
  
  #Select only necessary columns for joining
  catch4 <- catch_sub[,c("hauljoin", "survey", "year", "depth_m", "latitude_dd_start", "longitude_dd_start", "cpue_kgkm2", "species_name", "common_name")]
  bio4 <- dplyr::filter(bio_sub[,c("hauljoin", "performance", "species_name", "common_name", "length_cm", "weight", "sex", "age")], !is.na(length_cm))
  #Combine data
  dat <-dplyr::left_join(bio4, catch4, relationship = "many-to-many")
  dat <- dplyr::mutate(dat, trawl_id=hauljoin)
  #According to the codebook https://repository.library.noaa.gov/view/noaa/50147, 0 means Good performance, and the other numbers are for "Satisfactory, and then a "but"..."; negative numbers are Unsatisfactory
  #Dataset already includes only Good and Satisfactory hauls, Unsatisfactory are removed
  
  #If years=T in the function, this will subset data to 1999 onward to remove possibly funky data prior to 1999
  years <- F
  if(years){
    dat_sub = dplyr::filter(dat, year>1999)
  }
  
  if(!years){
    dat_sub = dat
  }
  
  #Add column counting number of length observations per catch
  dat_sub$is_length <- ifelse(!is.na(dat_sub$length_cm), 1,0) 
  dat_sub <- group_by(dat_sub, trawl_id) %>% mutate(nlength=sum(is_length)) %>% ungroup()
  dat_sub$trawl_id <- as.numeric(dat_sub$trawl_id)
  
  # fit length-weight regression by year to predict fish weights that have lengths only.
  # note a rank-deficiency warning may indicate there is insufficient data for some year/sex combinations (likely for unsexed group)
  
  if(nrow(dat_sub)>0){
    fitted = dat_sub
    # #If region=T in function, this will do the length-weight regression for each broad geographic region (EBS, NBS, and GOA)
    # if(region){
    #   # dplyr::select(trawl_id,year,
    #   #               subsample_wt_kg, total_catch_wt_kg, area_swept_ha_der, cpue_kg_km2,
    #   #               individual_tracking_id, sex, length_cm, weight) %>%
    #   
    # #Create one set of data with only year/survey combinations with at least some weight data:
    #   fitted <-  filter(fitted, !is.na(length_cm)) %>%
    #     group_by(year,survey) %>%
    #     mutate(sum = sum(weight, na.rm=T)) %>%
    #     filter(sum>0) %>%
    #     ungroup()
    #   
    #   #Remove erroneous zero weight observations because cause error in lm()
    #   fitted$weight <- ifelse(fitted$weight==0, NA, fitted$weight)
    #   
    #  #And one for those without any weight data:
    #   not_fitted <-  filter(dat_sub, !is.na(length_cm)) %>%
    #     group_by(year,survey) %>%
    #     mutate(sum = sum(weight, na.rm=T)) %>%
    #     filter(sum==0) %>%
    #     ungroup()
    #   
    # #Fit regression model for years with data
    #   
    #   fitted <-
    #   group_nest(fitted, year,survey) %>%
    #   mutate(
    #     model = purrr::map(data, ~ lm(log(weight) ~ log(length_cm), data = .x)),
    #     tidied = purrr::map(model, broom::tidy),
    #     augmented = purrr::map(model, broom::augment),
    #     predictions = purrr::map2(data, model, modelr::add_predictions)
    #   )
    # }
    
    #If region=F in function, this will do length-weight regression for all of Alaska combined
    # dplyr::select(trawl_id,year,
    #               subsample_wt_kg, total_catch_wt_kg, area_swept_ha_der, cpue_kg_km2,
    #               individual_tracking_id, sex, length_cm, weight) %>%
    
    #Create one set of data with only year/survey combinations with at least some weight data:
    fitted <-  filter(fitted, !is.na(length_cm)) %>%
      group_by(year) %>%
      mutate(sum = sum(weight, na.rm=T)) %>%
      filter(sum>0) %>%
      ungroup()
    #And one for those without any weight data:
    not_fitted <-  filter(dat_sub, !is.na(length_cm)) %>%
      group_by(year,survey) %>%
      mutate(sum = sum(weight, na.rm=T)) %>%
      filter(sum==0) %>%
      ungroup()
    
    if(nrow(fitted)>0){
      #Fit regression model for years with data
      #Remove weird zero weights
      fitted$weight <- ifelse(fitted$weight==0, NA, fitted$weight)
      fitted <-
        group_nest(fitted, year) %>%
        mutate(
          model = purrr::map(data, ~ lm(log(weight) ~ log(length_cm), data = .x)),
          tidied = purrr::map(model, broom::tidy),
          augmented = purrr::map(model, broom::augment),
          predictions = purrr::map2(data, model, modelr::add_predictions)
        )
      # replace missing weights with predicted weights for years with data
      dat_pos2 = fitted %>%
        tidyr::unnest(predictions) %>%
        dplyr::select(-data, -model, -tidied, -augmented) %>%
        try(dplyr::mutate(weight = ifelse(is.na(weight), exp(pred), weight)))
      
      #combine back with data for years without data
      dat_pos <- bind_rows(dat_pos2, not_fitted)
    }
    if(nrow(fitted)==0){
      dat_pos <- not_fitted
    }
    ##Checked that this works by checking number of rows--looks like it passes!
    
    #If there is no weight data available to get weight-length empirical interpolation, use fishbase to fill in
    # find length-weight relationship parameters for one species
    pars <- rfishbase::length_weight(
      #convert scientific name to first letter capitalized for fishbase
      species_list = stringr::str_to_sentence(sci_name, locale = "en"))
    #Get mean
    a <- mean(pars$a)
    b <- mean(pars$b)
    #For longspine thornyhead, which is not in database for some reason, got these values from fishbase direct page
    if(sci_name=="sebastolobus altivelis"){
      a <- 0.00912
      b=3.09
    }
    ##Assign a and b values if there are none available for the species in fishbase
    if(a=="NaN"){
      a <- 0.00675	
    }
    if(b=="NaN"){
      b <- 3
    }
    #Calculate weight (and convert from g to kg)
    dat_pos <- dplyr::mutate(dat_pos, weight = ifelse(is.na(weight), ((a*length_cm^b)*0.001), weight))
    #convert cm-g units
    #haul level weight
    dat_test <- group_by(dat_pos, trawl_id) %>% mutate(haul_weight=mean(weight)) %>% ungroup()
    
    #make column of trawl_id, which is called hauljoin originally in the AFSC data
    trawlids <- unique(dat_pos$trawl_id)
    if(length(trawlids!=0)){
      p <- data.frame(trawl_id = trawlids,
                      p1 = 0,
                      p2 = 0,
                      p3 = 0,
                      p4 = 0)
      
      sizethresholds <- quantile(dat_pos$weight, c(0.15, 0.5, 0.85, 1), na.rm = T)
      for (i in 1:length(trawlids)) {
        haul_sample<- dplyr::filter(dat_pos, trawl_id == trawlids[i])
        if(nrow(haul_sample) > 0 | var(haul_sample$weight >0)) {
          # fit kernel density to weight frequency
          smoothed_w <- KernSmooth::bkde(haul_sample$weight, range.x = c(min(dat_pos$weight), max(dat_pos$weight)), bandwidth = 2)
          # make sure smoother predicts positive or zero density
          smoothed_w$y[smoothed_w$y<0] <- 0
          # calculate proportion by biomass and by number
          p_w_byweight <- smoothed_w$y * smoothed_w$x / sum(smoothed_w$x*smoothed_w$y)
          
          
          p_w_byweight[p_w_byweight<0] <- 0
          #p_w_bynum[p_w_bynum<0] <- 0
          
          p1 <- sum(p_w_byweight[smoothed_w$x<=sizethresholds[1]])
          p2 <- sum(p_w_byweight[smoothed_w$x>sizethresholds[1] & smoothed_w$x <=sizethresholds[2]])
          p3 <- sum(p_w_byweight[smoothed_w$x>sizethresholds[2] & smoothed_w$x <=sizethresholds[3]])
          p4 <- sum(p_w_byweight[smoothed_w$x>sizethresholds[3]])
          
          p[i,2:5] <- c(p1, p2, p3, p4)
          
        }
        else {
          indx <- which(sizethresholds>haul_sample$weight)
          p[i, min(indx)+1] <- 1
        }
      }
      
      # add hauls with zero catch back in
      absent = filter(dat_sub, cpue_kgkm2 == 0)
      if(nrow(absent)>0){
        trawlids <- unique(absent$trawl_id)
        absent.df <- data.frame(trawl_id = trawlids,
                                p1 = 0,
                                p2 = 0,
                                p3 = 0,
                                p4 = 0)
        
        all_hauls <- rbind(p, absent.df)
      }
      if(nrow(absent)==0){
        all_hauls <- p
      }
      all_hauls$trawl_id <- as.numeric(all_hauls$trawl_id)
      dat_sub$median_weight <- median(dat_sub$weight, na.rm=T)
      nlengths <- unique(dat_sub[,c("trawl_id","nlength", "median_weight")])
      meanweight <- unique(dat_test[,c("trawl_id","haul_weight")])
      all_hauls2 <- left_join(all_hauls, nlengths)
      all_hauls2 <- left_join(all_hauls2, meanweight)
      return(all_hauls2)
    }
  }
  if(nrow(dat_sub)>0){
    if(nrow(fitted)==0 & nrow(not_fitted)==0){
      trawlids <- unique(dat_sub$trawl_id)
      absent.df <- data.frame(trawl_id = trawlids,
                              p1 = NA,
                              p2 = NA,
                              p3 = NA,
                              p4 = NA,
                              nlength=0,
                              haul_weight=NA, 
                              median_weight=NA)
      return(absent.df)
    }
  }
  #If there are no hauls at all with any length measurements, do this instead
  if(nrow(dat_sub)>0){
    trawlids <- unique(dat_sub$trawl_id)
    if(length(trawlids)==0){
      absent.df <- data.frame(trawl_id = trawlids,
                              p1 = NA,
                              p2 = NA,
                              p3 = NA,
                              p4 = NA,
                              nlength=0,
                              haul_weight=NA)
      return(absent.df)
    }
  }
if(nrow(dat_sub)==0){
  
  return(warning("species not present in data"))
}
  }


load_data_afsc <- function(sci_name, spc, dat.by.size, length=T) {
  bio2 <-readRDS("data/fish_raw/NOAA/ak_bts_goa_ebs_nbs_indivero_all_levels.RDS")
  dat <-readRDS("data/fish_raw/NOAA/ak_bts_goa_ebs_nbs_indivero_cpue_zerofilled.RDS")
  species <- bio2$species
  names(dat) = tolower(names(dat))
  names(species) =tolower(names(species))
  species$species_name <- tolower(species$species_name)
  species$common_name <- tolower(species$common_name)
  dat <- dplyr::left_join(dat, species)
  dat.by.size$trawl_id <- as.character(dat.by.size$trawl_id)
  dat$trawl_id <-as.character(dat$hauljoin)
  dat = dplyr::filter(dat, species_name ==sci_name)
  dat <- dplyr::left_join(dat, dat.by.size, by = "trawl_id")
  # remove tows where there was positive catch but no length measurements
  if(length){
    dat <- dplyr::filter(dat, !is.na(p1))
  }
  # analyze or years and hauls with adequate oxygen and temperature data, within range of occurrence
  
  # get julian day
  dat$julian_day <- rep(NA, nrow(dat))
  haul <- bio2$haul
  names(haul) <- tolower(names(haul))
  haul <- haul[,c("hauljoin", "start_time")]
  dat <- left_join(dat, haul)
  for (i in 1:nrow(dat)) dat$julian_day[i] <- as.POSIXlt(dat$start_time[i], format = "%Y-%b-%d")$yday
  
  #O2 from trawl data is in ml/l 
  # just in case, remove any missing or nonsense values from sensors
  # dat <- dplyr::filter(dat, !is.na(o2), !is.na(sal), !is.na(temp), is.finite(sal))
  # dat <- calc_po2_mi(dat)
  # dat <- dplyr::filter(dat, !is.na(temp), !is.na(mi))
  
  # prepare data and models -------------------------------------------------
  dat$longitude_dd <- dat$longitude_dd_start
  dat$latitude_dd <- dat$latitude_dd_start
  dat$longitude <- dat$longitude_dd
  dat$latitude <- dat$latitude_dd
  dat$scientific_name <- dat$species_name
  dat$cpue_kg_km2 <- dat$cpue_kgkm2
  dat$project <- dat$survey
  dat$event_id <- dat$trawl_id
  dat$date <- as.POSIXct(as.Date(dat$start_time, format = "%Y-%b-%d"))
  dat <- dplyr::select(dat, event_id, common_name, scientific_name, project, survey, year, date, longitude_dd, latitude_dd, longitude, latitude, cpue_kg_km2,
                       depth_m, julian_day, nlength, median_weight, haul_weight, p1, p2, p3, p4)
  
  
  # UTM transformation
  dat_ll = dat
  sp::coordinates(dat_ll) <- c("longitude_dd", "latitude_dd")
  sp::proj4string(dat_ll) <- sp::CRS("+proj=longlat +datum=WGS84")
  # convert to utm with spTransform
  dat_utm = sp::spTransform(dat_ll, 
                            sp::CRS("+proj=utm +zone=10 +datum=WGS84 +units=km"))
  # convert back from sp object to data frame
  dat = as.data.frame(dat_utm)
  dat = dplyr::rename(dat, X = coords.x1,
                      Y = coords.x2)
  dat$depth <- dat$depth_m
  dat$depth_m <- NULL
  return(dat)
}

length_expand_bc2 <- function(spc, sci_name) {
  # load, clean, and join data
  itis <- readRDS("data/fish_raw/BC/species-table.rds")
  bio2 <- readRDS("data/fish_raw/BC/pbs-bio-samples.rds")
  
  #Combine with species data
  bio <- dplyr::full_join(bio2, itis, by="itis", relationship="many-to-many")
  
  #Rename dogfish
  bio$common_name <- ifelse(bio$common_name=="north pacific spiny dogfish", "spiny dogfish", bio$common_name)
  
  #Put bio data in the same format as the NOAA bio data
  #Convert g to kg
  bio$weight <- bio$weight*0.001
  
  #rename columnns
  bio$length_cm <- bio$length
  
  #Rename missing species
  if(sci_name=="sebastes aleutianus"){
    bio$scientific_name <- ifelse(str_detect(bio$scientific_name, "sebastes aleutianus"), "sebastes aleutianus", bio$scientific_name)
    bio$common_name <- ifelse(str_detect(bio$common_name, "rougheye"), "rougheye rockfish", bio$common_name)
  }

  # filter out species of interest from joined (catch/haul/bio) dataset
  bio_sub = dplyr::filter(bio[,c("event_id", "length_cm", "weight", "scientific_name")], scientific_name==sci_name)
  
  #If there is weight or length data present in data, get haul mean of weight
  if(nrow(bio_sub)>0) {
   ##Add weights if only length data:
    # find length-weight relationship parameters for one species
    pars <- rfishbase::length_weight(
      #convert scientific name to first letter capitalized for fishbase
      species_list = stringr::str_to_sentence(sci_name, locale = "en"))
    #Get mean
    a <- mean(pars$a)
    b <- mean(pars$b)
    #For longspine thornyhead, which is not in database for some reason, got these values from fishbase direct page
    if(sci_name=="sebastolobus altivelis"){
      a <- 0.00912
      b=3.09
    }
    ##Assign a and b values if there are none available for the species in fishbase
    if(a=="NaN"){
      a <- 0.00675	
    }
    if(b=="NaN"){
      b <- 3
    }
    #Calculate weight (and convert from g to kg)
    bio_sub <- dplyr::mutate(bio_sub, weight = ifelse(is.na(weight), ((a*length_cm^b)*0.001), weight))
    bio_summ <- group_by(bio_sub, event_id) %>% summarize(haul_weight=median(weight))
    return(bio_summ)
  }
}
  

length_expand_bc <- function(spc, sci_name) {
  # load, clean, and join data
  itis <- readRDS("data/fish_raw/BC/species-table.rds")
  haul <- readRDS("data/fish_raw/BC/pbs-haul.rds")
  catch <- readRDS("data/fish_raw/BC/pbs-catch.rds")
  bio2 <- readRDS("data/fish_raw/BC/pbs-bio-samples.rds")
  
  catch$common_name <- catch$species_common_name
  catch$species_science_name <- NULL
  catch$species_common_name <- NULL
  
  #Rename dogfish
  catch$common_name <- ifelse(catch$common_name=="north pacific spiny dogfish", "spiny dogfish", catch$common_name)
  
  #Merge the official BC bio data and the official BC haul data to get metadata (from Sean, or here: https://open.canada.ca/data/en/dataset/86af7918-c2ab-4f1a-ba83-94c9cebb0e6c)
  bio <- dplyr::full_join(bio2, haul, by="event_id", relationship="many-to-many")
  
  #Combine with species data
  bio <- dplyr::full_join(bio, itis, by="itis", relationship="many-to-many")
  catch <- dplyr::full_join(catch, itis, relationship="many-to-many")

  #Put bio data in the same format as the NOAA bio data
  #Convert g to kg
  bio$weight <- bio$weight*0.001
  
  #rename columnns
  bio$length_cm <- bio$length
  
  #Clean catch data
  names(catch) = tolower(names(catch))
  
  #Merge with ITIS information to get scientific name
  catch <- left_join(catch,itis)
  
  #haul$sampling_end_hhmmss = as.numeric(haul$sampling_end_hhmmss)
  #haul$sampling_start_hhmmss = as.numeric(haul$sampling_start_hhmmss)
  
  #Combine catch data with haul data
  dat <- dplyr::left_join(catch, haul, relationship = "many-to-many")
  
  #Rename missing species
  if(sci_name=="sebastes aleutianus"){
    bio$scientific_name <- ifelse(str_detect(bio$scientific_name, "sebastes aleutianus"), "sebastes aleutianus", bio$scientific_name)
    bio$common_name <- ifelse(str_detect(bio$common_name, "rougheye"), "rougheye rockfish", bio$common_name)
  }
  #Combine bio/haul data with catch data
  dat <- dplyr::left_join(dat, filter(bio[,c("event_id", "age","length_cm", "weight", "scientific_name", "common_name")], !is.na(length_cm)), relationship = "many-to-many")
  
  # filter out species of interest from joined (catch/haul/bio) dataset
  dat_sub = dplyr::filter(dat, scientific_name==sci_name)

  dat_sub$event_id <- as.numeric(dat_sub$event_id)
  trawlids <- unique(dat_sub$event_id)
  
  if(nrow(dat_sub)>0) {
    #Add column counting number of length observations per catch
    dat_sub$is_length <- ifelse(!is.na(dat_sub$length_cm), 1,0) 
    dat_sub <- group_by(dat_sub, event_id) %>% mutate(nlength=sum(is_length)) %>% ungroup()
    
    # fit length-weight regression by year to predict fish weights that have lengths only.
    # note a rank-deficiency warning may indicate there is insufficient data for some year/sex combinations (likely for unsexed group)
    
    #Create one set of data with only year/survey combinations with at least some weight data:
    fitted <-  filter(dat_sub, !is.na(length_cm)) %>%
      group_by(year) %>%
      mutate(sum = sum(weight, na.rm=T)) %>%
      filter(sum>0) %>%
      ungroup()
    fitted$weight <- ifelse(fitted$weight==0, NA, fitted$weight)
    #And one for those without any weight data:
    not_fitted <-  filter(dat_sub, !is.na(length_cm)) %>%
      group_by(year) %>%
      mutate(sum = sum(weight, na.rm=T)) %>%
      filter(sum==0) %>%
      ungroup()
    
    #Fit regression model for years with data
    fitted <-
      group_nest(fitted, year) %>%
      mutate(
        model = purrr::map(data, ~ lm(log(weight) ~ log(length_cm), data = .x)),
        tidied = purrr::map(model, broom::tidy),
        augmented = purrr::map(model, broom::augment),
        predictions = purrr::map2(data, model, modelr::add_predictions)
      )
    
    # replace missing weights with predicted weights
    dat_pos2 = fitted %>%
      tidyr::unnest(predictions) %>%
      dplyr::select(-data, -model, -tidied, -augmented) %>%
      try(dplyr::mutate(weight = ifelse(is.na(weight), exp(pred), weight)))
    
    #combine back with data for years without data
    dat_pos <- bind_rows(dat_pos2, not_fitted)
    
    ##Checked that this works by checking number of rows--looks like it passes!
    
    #If there is no weight data available to get weight-length empirical interpolation, use fishbase to fill in
    # find length-weight relationship parameters for one species
    pars <- rfishbase::length_weight(
      #convert scientific name to first letter capitalized for fishbase
      species_list = stringr::str_to_sentence(sci_name, locale = "en"))
    #Get mean
    a <- mean(pars$a)
    b <- mean(pars$b)
    #For longspine thornyhead, which is not in database for some reason, got these values from fishbase direct page
    if(sci_name=="sebastolobus altivelis"){
      a <- 0.00912
      b=3.09
    }
    ##Assign a and b values if there are none available for the species in fishbase
    if(a=="NaN"){
      a <- 0.00675	
    }
    if(b=="NaN"){
      b <- 3
    }
    #Calculate weight (and convert from g to kg)
    dat_pos <- dplyr::mutate(dat_pos, weight = ifelse(is.na(weight), ((a*length_cm^b)*0.001), weight))
    dat_test <- group_by(dat_pos, event_id) %>% mutate(haul_weight=mean(weight)) %>% ungroup()
    
    trawlids <- unique(dat_pos$event_id)
    if(length(trawlids!=0)){
      p <- data.frame(event_id = trawlids,
                      p1 = 0,
                      p2 = 0,
                      p3 = 0,
                      p4 = 0)
      
      sizethresholds <- quantile(dat_pos$weight, c(0.15, 0.5, 0.85, 1), na.rm = T)
      for (i in 1:length(trawlids)) {
        haul_sample<- dplyr::filter(dat_pos, event_id == trawlids[i])
        if(nrow(haul_sample) > 0 | var(haul_sample$weight >0)) {
          # fit kernel density to weight frequency
          smoothed_w <- KernSmooth::bkde(haul_sample$weight, range.x = c(min(dat_pos$weight), max(dat_pos$weight)), bandwidth = 2)
          # make sure smoother predicts positive or zero density
          smoothed_w$y[smoothed_w$y<0] <- 0
          # calculate proportion by biomass and by number
          p_w_byweight <- smoothed_w$y * smoothed_w$x / sum(smoothed_w$x*smoothed_w$y)
          
          p_w_byweight[p_w_byweight<0] <- 0
          #p_w_bynum[p_w_bynum<0] <- 0
          
          p1 <- sum(p_w_byweight[smoothed_w$x<=sizethresholds[1]])
          p2 <- sum(p_w_byweight[smoothed_w$x>sizethresholds[1] & smoothed_w$x <=sizethresholds[2]])
          p3 <- sum(p_w_byweight[smoothed_w$x>sizethresholds[2] & smoothed_w$x <=sizethresholds[3]])
          p4 <- sum(p_w_byweight[smoothed_w$x>sizethresholds[3]])
          
          p[i,2:5] <- c(p1, p2, p3, p4)
          
        }
        else {
          indx <- which(sizethresholds>haul_sample$weight)
          p[i, min(indx)+1] <- 1
        }
      }
      
      # add hauls with zero catch back in
      absent = filter(dat_sub, catch_weight == 0)
      trawlids <- unique(absent$event_id)
      absent.df <- data.frame(event_id = trawlids,
                              p1 = 0,
                              p2 = 0,
                              p3 = 0,
                              p4 = 0)
      
      all_hauls <- rbind(p, absent.df)
      all_hauls$event_id <- as.numeric(all_hauls$event_id)
      dat_sub$median_weight <- median(dat_sub$weight, na.rm=T)
      nlengths <- unique(dat_sub[,c("event_id","nlength", "median_weight")])
      meanweight <- unique(dat_test[,c("event_id","haul_weight")])
      all_hauls2 <- left_join(all_hauls, nlengths)
      all_hauls2 <- left_join(all_hauls2, meanweight)
      return(all_hauls2)
    }
  }
  if(nrow(dat_sub)>0){
    if(nrow(fitted)==0 & nrow(not_fitted)==0){
      trawlids <- unique(dat_sub$event_id)
      absent.df <- data.frame(event_id = trawlids,
                              p1 = NA,
                              p2 = NA,
                              p3 = NA,
                              p4 = NA,
                              nlength=0,
                              haul_weight=NA, 
                              median_weight=NA)
      return(absent.df)
    }
  }
  if(nrow(dat_sub)>0){
    if(length(trawlids)==0){
      absent.df <- data.frame(event_id = trawlids,
                              p1 = NA,
                              p2 =NA,
                              p3 = NA,
                              p4 = NA,
                              nlength=0,
                              haul_weight=NA)
      return(absent.df)
    }
  }
  if(nrow(dat_sub)==0){
    
    return(warning("species not present in data"))
  }
  }

load_data_bc <- function(sci_name,dat.by.size, length=T, spc) {
  catch <- readRDS("data/fish_raw/BC/pbs-catch.rds")
  haul <- readRDS("data/fish_raw/BC/pbs-haul.rds")
  itis <- readRDS("data/fish_raw/BC/species-table.rds")
  #Rename dogfish
  catch$common_name <- ifelse(catch$species_common_name=="north pacific spiny dogfish", "spiny dogfish", catch$species_common_name)
  catch$species_common_name <- NULL
  catch <- left_join(catch,itis)
  
  #haul$sampling_end_hhmmss = as.numeric(haul$sampling_end_hhmmss)
  #haul$sampling_start_hhmmss = as.numeric(haul$sampling_start_hhmmss)
  
  #Combine catch data with haul data
  dat <- dplyr::left_join(catch, haul, relationship = "many-to-many")
  
  # dat.by.size$event_id <- as.character(dat.by.size$event_id)
  dat = dplyr::filter(dat, scientific_name == sci_name)
  dat <- left_join(dat, dat.by.size, by = "event_id")
  # remove tows where there was positive catch but no length measurements
  if(length){
    dat <- dplyr::filter(dat, !is.na(p1))
  }
  # analyze or years and hauls with adequate oxygen and temperature data, within range of occurrence
  
  # get julian day
  dat$julian_day <- rep(NA, nrow(dat))
  for (i in 1:nrow(dat)) dat$julian_day[i] <- as.POSIXlt(dat$date[i], format = "%Y-%b-%d")$yday
  
  
  #O2 from trawl data is in ml/l 
  # just in case, remove any missing or nonsense values from sensors
  # dat <- dplyr::filter(dat, !is.na(o2), !is.na(sal), !is.na(temp), is.finite(sal))
  # dat <- calc_po2_mi(dat)
  # dat <- dplyr::filter(dat, !is.na(temp), !is.na(mi))
  
  # prepare data and models -------------------------------------------------
  dat$cpue_kg_km2 <- dat$catch_weight
  dat$longitude_dd <- dat$lon_start
  dat$latitude_dd <- dat$lat_start
  dat$longitude <- dat$lon_start
  dat$latitude <- dat$lat_start
  dat$event_id <- as.character(dat$event_id)
  dat$year <- substr(dat$date, start=1, stop=4)
  dat$year <- as.integer(dat$year)
  dat$date <- as.POSIXct(as.Date(dat$date, format = "%Y-%b-%d"))
  dat$project <- dat$survey_name
  dat$salinity_psu <- dat$salinity_PSU
  dat$salinity_PSU <- NULL
  dat <- dplyr::select(dat, event_id, scientific_name, project, year, date, longitude_dd, latitude_dd, longitude, latitude, cpue_kg_km2,
                       depth_m, julian_day, nlength,median_weight, haul_weight, pass, p1, p2, p3, p4, temperature_C, do_mlpL, salinity_psu)
  dat <- filter(dat, !is.na(latitude_dd))
  
  
  # UTM transformation
  dat_ll = dat
  sp::coordinates(dat_ll) <- c("longitude_dd", "latitude_dd")
  sp::proj4string(dat_ll) <- sp::CRS("+proj=longlat +datum=WGS84")
  # convert to utm with spTransform
  dat_utm = sp::spTransform(dat_ll, 
                            sp::CRS("+proj=utm +zone=10 +datum=WGS84 +units=km"))
  # convert back from sp object to data frame
  dat = as.data.frame(dat_utm)
  dat = dplyr::rename(dat, X = coords.x1,
                      Y = coords.x2)
  dat$common_name <- spc
  dat$depth <- dat$depth_m
  dat$depth_m <- NULL
  dat$survey <- "dfo"
  
  #convert oxygen mg/L to umol_kg
  SA = gsw_SA_from_SP(dat$salinity_psu,dat$depth,dat$longitude,dat$latitude) #absolute salinity for pot T calc
  pt = gsw_pt_from_t(SA,dat$temperature_C,dat$depth) #potential temp at a particular depth
  CT = gsw_CT_from_t(SA,dat$temperature_C,dat$depth) #conservative temp
  dat$sigma0_kgm3 = gsw_sigma0(SA,CT)
  dat$O2_umolkg = dat$do_mlpL*44660/(dat$sigma0_kgm3+1000) 
  
  return(dat)
}

IPHC <- function (catch, adjustment, spc, sci_name) {
  #make lowercase
  colnames(catch) <- tolower(colnames(catch))
  colnames(adjustment) <- tolower(colnames(adjustment))
  
  if(spc!="pacific halibut"){
    catch$scientific_name <- tolower(catch$'scientific name')
    catch$common_name<- paste(spc)
    catch <- subset(catch, catch$scientific_name==sci_name)
    catch$cpue_count <- (catch$'number observed'/catch$hooksobserved)*catch$hooksretrieved
    catch$cpue_weight <- NA
    #Get CTD data and sample metadata to join
    dat <-  read_excel("data/fish_raw/IPHC/Set and Pacific halibut data.xlsx")
    colnames(dat) <- tolower(colnames(dat))
    dat <- dat[,c("year","stlkey","setno","date", 'midlat fished','midlon fished', "avgdepth (fm)", "temp c", "salinity psu", "oxygen_ml")]
    #Join
    catch <- left_join(catch, dat)
    ##Add zeros
    #All sets
    sets <- unique(dat)
    #Which are not in the catch data (positive catches)
    missing <- anti_join(sets, catch)
    #Add species columns back
    missing$common_name <- paste(spc)
    missing$scientific_name <- paste(sci_name)
    #Add weight and count
    missing$cpue_weight <- 0
    missing$cpue_count <- 0
    #Add to dataset
    catch <- bind_rows(catch, missing)
     }
  if(spc=="pacific halibut"){
  #Calculate CPUE as U32 and O32 count and weight divided by effective skates hauled
  catch$cpue_O32_count <- catch$`o32 pacific halibut count`/catch$`effective skates hauled`
  catch$cpue_U32_count <- catch$`u32 pacific halibut count`/catch$`effective skates hauled`
  catch$cpue_O32_weight <- catch$`o32 pacific halibut weight`/catch$`effective skates hauled`
  catch$cpue_U32_weight <- catch$`u32 pacific halibut weight`/catch$`effective skates hauled`
  #Sum weight classes together
  catch$cpue_count <- catch$cpue_O32_count+catch$cpue_U32_count
  catch$cpue_weight <- catch$cpue_O32_weight+catch$cpue_U32_weight
  catch$common_name <- "pacific halibut"
  catch$scientific_name <-"hippoglossus stenolepis"
  }
  #join
  adjustment$stlkey <- as.character(adjustment$stlkey)
  data <- left_join(catch, adjustment, by="stlkey")
  data$h.adj <- as.numeric(data$h.adj)
  
  #Calculate ajustment factor
  data$cpue_count <- data$cpue_count * data$h.adj
  data$cpue_weight <- data$cpue_weight * data$h.adj
  #Convert biomass to kg
  data$cpue_weight <- data$cpue_weight*0.453592
  #Extract columns of interest
  dat <- data[,c("year.x", "date.x", "common_name", "scientific_name", "midlat fished", "midlon fished", "avgdepth (fm)","temp c", "salinity psu", "oxygen_ml", "cpue_weight", "cpue_count")]
  #Re-name columns to match the NOAA and BC data
  colnames(dat) <- c("year", "date", "common_name", "scientific_name", "latitude", "longitude", "depth", "temperature_C", "salinity_psu", "do_mlpL", "cpue_weight", "cpue_count")

  #Add columns to match the NOAA and BC data
  dat$survey <- "iphc"
  
  ##Add region
  # load regional polygons
  regions.hull <- readRDS("data/processed_data/regions_hull.rds")
  #make dataframe an sf object
  dat <- drop_na(dat, latitude, longitude)
  dat_df <-  st_as_sf(dat, coords = c("longitude", "latitude"), crs = st_crs(4326))
  dat_df$latitude <- dat$latitude
  dat_df$longitude <- dat$longitude
  # cycle through all regions
  region_list <- c("ai", "bc", "cc", "ebs", "goa")
  dats <- list()
  for (i in 1:length(region_list)) {
    region <- region_list[i]
    poly <- regions.hull[i,2]
    # pull out observations within each region
    region_dat  <- st_filter(dat_df, poly)
    region_dat$region <- paste(region_list[i])
    dats[[i]] <- as.data.frame(region_dat)
  }
  #Bind back together
  dat <- bind_rows(dats)
  #Date and month in right format
  dat$month <- case_when(grepl("May",dat$date) ~5,
                         grepl("Jun",dat$date)  ~6,
                         grepl("Jul",dat$date)  ~7,
                         grepl("Aug",dat$date)  ~8,
                         grepl("Sep",dat$date)  ~9,
                         grepl("Oct",dat$date)  ~10)
  dat$day <- as.numeric(substr(dat$date, 1,2))
  dat$date <-  as.POSIXct(as.Date(with(dat,paste(year,month,day,sep="-")),"%Y-%m-%d"))
  dat$doy <- as.POSIXlt(dat$date, format = "%Y-%b-%d")$yday
  dat$year <- as.numeric(dat$year)
  dat$day <- NULL
  #convert fathoms to m
  dat$depth <- dat$depth*1.8288
  
  #convert oxygen mg/L to umol_kg
  SA = gsw_SA_from_SP(dat$salinity_psu,dat$depth,dat$longitude,dat$latitude) #absolute salinity for pot T calc
  pt = gsw_pt_from_t(SA,dat$temperature_C,dat$depth) #potential temp at a particular depth
  CT = gsw_CT_from_t(SA,dat$temperature_C,dat$depth) #conservative temp
  dat$sigma0_kgm3 = gsw_sigma0(SA,CT)
  dat$O2_umolkg = dat$do_mlpL*44660/(dat$sigma0_kgm3+1000) 
  
  #Convert coordinates
  dat <- subset(dat, !is.na(latitude))
  dat <- dat %>%
    st_as_sf(coords=c('longitude','latitude'),crs=4326,remove = F) %>%  
    st_transform(crs = "+proj=utm +zone=10 +datum=WGS84 +units=km") %>% 
    mutate(X=st_coordinates(.)[,1],Y=st_coordinates(.)[,2]) 
 
  return(dat)
}

combine_all <- function(type){
  #BC
  #BC raw data
  bio_qc <- read.csv("data/fish_raw/BC/QCS_biology.csv")
  bio_vi <- read.csv("data/fish_raw/BC/WCVI_biology.csv")
  bio_hs <- read.csv("data/fish_raw/BC/HS_biology.csv")
  bio_hg <- read.csv("data/fish_raw/BC/WCHG_biology.csv")
  itis  <- readRDS("data/fish_raw/BC/itis_bc.rds")
  
  haul <- readRDS("data/fish_raw/BC/pbs-haul.rds")
  catch <- readRDS("data/fish_raw/BC/pbs-catch.rds")
  
  haul_qc <- read.csv("data/fish_raw/BC/QCS_effort.csv")
  haul_vi <- read.csv("data/fish_raw/BC/WCVI_effort.csv")
  haul_hs <- read.csv("data/fish_raw/BC/HS_effort.csv")
  haul_hg <- read.csv("data/fish_raw/BC/WCHG_effort.csv")
  
  #Combine BC bio and haul data to combine together
  bio2 <- rbind(bio_hg, bio_hs, bio_qc, bio_vi)
  haul2 <- rbind(haul_hg, haul_hs, haul_qc, haul_vi)
  
  #Merge the official BC bio data and the official BC haul data to get metadata (from here: https://open.canada.ca/data/en/dataset/86af7918-c2ab-4f1a-ba83-94c9cebb0e6c)
  bio2$Set.number <- bio2$Tow.number
  bio <- dplyr::left_join(bio2, haul2, relationship="many-to-many")
  
  names(bio) = tolower(names(bio))
  names(haul) = tolower(names(haul))
  
  #Put bio data in the same format as the NOAA bio data
  bio$scientific_name <- tolower(bio$scientific.name)
  bio$common_name <- tolower(bio$english.common.name)
  bio$weight <- bio$weight..g.*0.001
  bio$length_cm <- bio$fork.length..mm.*0.1
  bio$length_cm <- ifelse(is.na(bio$length_cm), (bio$total.length..mm.*0.1), bio$length_cm)
  
  #Make column names consistent so can join to surveyjoin haul data to get event_id
  haul$year <- as.character(substr(haul$date, start=1, stop=4))
  bio$year <- as.character(bio$survey.year)
  haul$set.date <- substr(haul$date, start=1, stop=10)
  bio$lat_start <- bio$start.latitude
  bio$lat_end <- bio$end.latitude
  bio$lon_start <- bio$start.longitude
  bio$lon_end <- bio$end.longitude
  
  #Combine BC bio data and haul metadata with surveyjoin metadata
  bio3 <- dplyr::left_join(bio[,c("age", "sex", "length_cm", "weight", "scientific_name", "common_name", "lat_start", "lon_start", "lat_end", "lon_end", "set.date")], haul, relationship="many-to-many")
  bio <- bio3
  
  #Clean catch data
  names(catch) = tolower(names(catch))
  
  catch <- dplyr::left_join(catch,itis)
  
  #haul$sampling_end_hhmmss = as.numeric(haul$sampling_end_hhmmss)
  #haul$sampling_start_hhmmss = as.numeric(haul$sampling_start_hhmmss)
  
  #Combine catch data with haul data
  dat_bc <- dplyr::left_join(catch, haul, relationship = "many-to-many")
  #Combine bio/haul data with catch data
  #Only positive catches
  dat_bc <- subset(dat_bc, catch_weight>0)
  dat_bc$cpue_kg_km2 <- dat_bc$catch_weight
  #Biomass
  if(type=="hauls"){
    counts_bc <- aggregate(cpue_kg_km2~scientific_name, dat_bc, FUN=length)
  }
  if(type=="biomass"){
    counts_bc <- aggregate(cpue_kg_km2~scientific_name, dat_bc, FUN=sum)
  }
  
  #Add common name
  counts_bc <- unique(dplyr::left_join(counts_bc, bio[,c("scientific_name", "common_name")]))
  counts_bc <- dplyr::mutate(counts_bc, bc=cpue_kg_km2)
  counts_bc$cpue_kg_km2 <- NULL
  rm(list=setdiff(ls(), "counts_bc"))
  
  #NWFSC
  # load, clean, and join data
  bio <- readRDS("data/fish_raw/NOAA/nwfsc_bio.rds")
  load("data/fish_raw/NOAA/nwfsc_haul.rda")
  haul <- nwfsc_haul
  catch <- readRDS("data/fish_raw/NOAA/nwfsc_catch.rds")
  names(catch) = tolower(names(catch))
  names(bio) = tolower(names(bio))
  names(haul) = tolower(names(haul))
  bio$scientific_name <- tolower(bio$scientific_name)
  bio$common_name <- tolower(bio$common_name)
  catch$common_name <- tolower(catch$common_name)
  
  bio$trawl_id = as.character(bio$trawl_id)
  haul$trawl_id = as.character(haul$event_id)
  catch$trawl_id=as.character(catch$trawl_id)
  haul$year <- as.character(substr(haul$date, start=1, stop=4))
  bio$year <- as.character(bio$year)
  catch$date <- NULL
  bio$date <- NULL
  bio$year <- NULL
  
  #haul$sampling_end_hhmmss = as.numeric(haul$sampling_end_hhmmss)
  #haul$sampling_start_hhmmss = as.numeric(haul$sampling_start_hhmmss)
  
  #Combine data
  dat = dplyr::left_join(catch[,c("trawl_id","common_name", "subsample_count","area_swept_ha","longitude_dd", "latitude_dd",
                                  "subsample_wt_kg","total_catch_numbers","total_catch_wt_kg","cpue_kg_km2")], haul, relationship = "many-to-many") %>%
    dplyr::left_join(filter(bio[,c("trawl_id", "scientific_name", "common_name", "weight", "ageing_lab", "oto_id", "length_cm", "width_cm", "sex", "age")], !is.na(length_cm)), relationship = "many-to-many") %>%
    filter(performance == "Satisfactory")
  
  #Only positive catches
  dat_nw <- subset(dat, cpue_kg_km2>0)
  #Biomass
  if(type=="hauls"){
    counts_nw <- aggregate(cpue_kg_km2~scientific_name, dat_nw, FUN=length)
  }
  if(type=="biomass"){
    counts_nw <- aggregate(cpue_kg_km2~scientific_name, dat_nw, FUN=sum)
  }
  counts_nw <- unique(left_join(counts_nw, bio[,c("scientific_name", "common_name")]))
  counts_nw <- dplyr::mutate(counts_nw, nw=cpue_kg_km2)
  counts_nw$cpue_kg_km2 <- NULL
  
  rm(list=setdiff(ls(), c("counts_bc", "counts_nw")))
  
  ##Alaska
  bio2 <-readRDS("data/fish_raw/NOAA/ak_bts_goa_ebs_nbs_all_levels.RDS")
  catch2 <- readRDS("data/fish_raw/NOAA/ak_bts_goa_ebs_nbs_cpue_zerofilled.RDS")
  
  #Isolate necessary parts of full data to get specimen weights/lengths per haul
  haul <- bio2$haul
  specimen <- bio2$specimen
  species <- bio2$species
  size <- bio2$size
  
  #make lowercase
  names(haul) <- tolower(names(haul))
  names(specimen) <- tolower(names(specimen))
  names(species) <- tolower(names(species))
  names(size) <- tolower(names(bio2$size))
  
  species$species_name <- tolower(species$species_name)
  
  #Combine catch data with species data
  names(catch2) <- tolower(names(catch2))
  catch <- dplyr::left_join(catch2, species)
  
  #Select only necessary columns for joining
  catch4 <- catch[,c("hauljoin", "survey", "year", "depth_m", "latitude_dd_start", "longitude_dd_start", "cpue_kgkm2", "species_name", "common_name")]
  dat <- dplyr::filter(catch4, year>1998)
  #bio4 <- dplyr::filter(bio4, year>1998)
  
  #Combine data
  dat <- dplyr::mutate(dat, trawl_id=hauljoin)
  #According to the codebook https://repository.library.noaa.gov/view/noaa/50147, 0 means Good performance, and the other numbers are for "Satisfactory, and then a "but"..."; negative numbers are Unsatisfactory
  #Dataset already includes only Good and Satisfactory hauls, Unsatisfactory are removed
  
  #Only positive catches
  dat$cpue_kg_km2 <- dat$cpue_kgkm2
  dat_ak <- subset(dat, cpue_kg_km2>0)
  dat_ak$scientific_name <- dat_ak$species_name
  #Biomass
  if(type=="hauls"){
    counts_ak <- aggregate(cpue_kg_km2~scientific_name, dat_ak, FUN=length)
  }
  if(type=="biomass"){
    counts_ak <- aggregate(cpue_kg_km2~scientific_name, dat_ak, FUN=sum)
  }
  counts_ak <- dplyr::mutate(counts_ak, ak=cpue_kg_km2)
  counts_ak$cpue_kg_km2 <- NULL
  counts_ak <- dplyr::left_join(counts_ak, species[,c("common_name", "species_name")], by=c("scientific_name"="species_name"))
  counts_ak$common_name <- tolower(counts_ak$common_name)
  
  #counts_ak <- unique(left_join(counts_nw, bio4[,c("scientific_name", "common_name")]))
  
  #Combine
  counts <- full_join(counts_nw, counts_bc)
  counts <- full_join(counts, counts_ak)
  return(counts)
}

###LOAD DATA ###
prepare_data <- function(spc,sci_name,mi,iphc, file_name, taxa){
  dat.by.size <- try(length_expand_bc(spc, sci_name))
  gc()
  if(is.data.frame(dat.by.size)){
    dat3 <- try(load_data_bc(sci_name = sci_name, spc=spc, dat.by.size = dat.by.size, length=F))
    dat3 <- try(unique(dat3))
  }
  gc()
  rm(dat.by.size)
  
  dat.by.size <- try(length_expand_nwfsc(spc=spc, sci_name=sci_name))
  gc()
  if(is.data.frame(dat.by.size)){
    dat2 <- try(load_data_nwfsc(spc= spc, sci_name=sci_name, dat.by.size = dat.by.size, length=F))
    dat2 <- try(unique(dat2))
  }
  
  gc()
  rm(dat.by.size)
  
  dat.by.size <- try(length_expand_afsc(spc, sci_name))
  gc()
  if(is.data.frame(dat.by.size)){
    dat5 <- try(load_data_afsc(sci_name = sci_name, spc=spc, dat.by.size = dat.by.size, length=F))
    dat5 <- try(unique(dat5))
  }

  gc()
  
  #All regions present
  if(exists("dat3") & exists("dat2") & exists("dat5")){
    dat4 <- bind_rows(dat3, dat2, dat5)
  }
  #Only BC
  if(exists("dat3") & !exists("dat2") & !exists("dat5")){
    dat4 <- dat3
  }
  #Only AK
  if(!exists("dat3") & !exists("dat2") & exists("dat5")){
    dat4 <- dat5
  }
  
  #Only NWFSC
  if(!exists("dat3") & exists("dat2") & !exists("dat5")){
    dat4 <- dat2
    
  }
  #BC & NW
  if(exists("dat3") & exists("dat2") & !exists("dat5")){
    dat4 <- bind_rows(dat3, dat2)
  }
  #BC & AK
  if(exists("dat3") & !exists("dat2")& exists("dat5")){
    dat4 <- bind_rows(dat3, dat5)
  }
  #AK & NW
  if(!exists("dat3") & exists("dat2")& exists("dat5")){
    dat4 <- bind_rows(dat2, dat5)
  }
  gc()

  #  adjustment <- read_excel("~/Dropbox/choke species/code/choke-species-data/data/fish_raw/IPHC/iphc-2023-fiss-hadj-20231031.xlsx")
   # dat_IPHC <- IPHC(catch, adjustment)
  #  dat4 <- bind_rows(dat4, dat_IPHC)
#  }
  #Make a duplicate copy object
  dat <- dat4
  
  ###Add in situ data####
    #Isolate just NWFSC and EBS/NBS/GOA data for joining, because  dfo are already in the dataset
    insitu <- readRDS("data/processed_data/o2/insitu_combined.rds")
    insitu <- filter(insitu, survey!="iphc")
    insitu <- filter(insitu, survey!="dfo")
    insitu$event_id <- as.character(insitu$event_id)
    #Northwest
    dat6 <- filter(dat4, survey=="nwfsc")
    #Just columns of interest
    dat6 <- left_join(dat6[,c("event_id", "date", "year", "survey", "latitude", "longitude", "depth", "X", "Y", "cpue_kg_km2", "julian_day", "nlength", "median_weight", "haul_weight", "pass", "p1", "p2", "p3", "p4", "common_name", "scientific_name", "depth", "vessel", "tow")], insitu, by=c("date", "latitude", "longitude"))
    #Edit columns
    dat6$event_id <- dat6$event_id.x
    dat6$event_id.x <- NULL
    dat6$event_id.y <- NULL
    dat6$year <- dat6$year.x
    dat6$year.y <- NULL
    dat6$year.x <- NULL
    dat6$depth <- dat6$depth.x
    dat6$depth.x <- NULL
    dat6$depth.y <- NULL
    dat6$X <- dat6$X.x
    dat6$Y <- dat6$Y.x
    dat6$X.x <- NULL
    dat6$X.y <- NULL
    dat6$Y.x <- NULL
    dat6$Y.y <- NULL
    dat6$depth.1 <- NULL
    dat6$month <- NULL
    dat6$doy <- NULL
    dat6$survey <- dat6$survey.x
    dat6$survey.x <- NULL
    dat6$survey.y <- NULL
    dat6$time <- NULL
    
    dat7 <- dat6
    
    #GOA, EBS & NBS 
    dat6 <- filter(dat4, survey=="EBS"|survey=="NBS"|survey=="GOA")
    dat6 <- left_join(dat6[,c("event_id", "date", "year", "survey", "latitude", "longitude", "depth", "X", "Y", "cpue_kg_km2", "julian_day", "nlength", "median_weight", "haul_weight", "pass", "p1", "p2", "p3", "p4", "common_name", "scientific_name", "depth", "vessel", "tow")], insitu[,c("event_id", "temperature_C", "do_mlpL", "salinity_psu", "sigma0_kgm3", "O2_umolkg")], by=c("event_id"))
    
    #Edit columns
    dat6$depth.1 <- NULL

    #Recombine back with DFO data
    if(exists("dat3")){
      dat <- bind_rows(dat6, dat7,dat3)
    }
    
    if(!exists("dat3")){
      dat <- bind_rows(dat6, dat7)
    }
  
    #Some other columns
    dat$doy <- dat$julian_day
    dat$julian_day <- NULL
    dat$cpue_kg_km2_sub <- dat$cpue_kg_km2 * (dat$p2+dat$p3)
    dat$cpue_kg_km2_sub <- ifelse(dat$cpue_kg_km2==0, 0, dat$cpue_kg_km2_sub)
    dat$p1 <- ifelse(dat$cpue_kg_km2==0, 0, dat$p1)
    dat$p2 <- ifelse(dat$cpue_kg_km2==0, 0, dat$p2)
    dat$p3 <- ifelse(dat$cpue_kg_km2==0, 0, dat$p3)
    dat$p4 <- ifelse(dat$cpue_kg_km2==0, 0, dat$p4)
    dat$survey <- tolower(dat$survey)
    dat$survey <- ifelse(dat$survey=="nwfsc", "wcbts", dat$survey)
    dat$region<- case_when(dat$survey=="nbs"~"nbs",
                           dat$survey=="wcbts"~"cc",
                           dat$survey=="ebs"~"ebs",
                           dat$survey=="goa"~"goa",
                           dat$survey=="dfo"~"bc")
     
    #Combine with IPHC if needed
    if(iphc){
      if(spc=="pacific halibut"){
        catch <-  read_excel("data/fish_raw/IPHC/Set and Pacific halibut data.xlsx")
      }
      if(spc!="pacific halibut"){
        catch <-  read_excel("data/fish_raw/IPHC/Non_Pacific_halibut_data.xlsx")
      }
      adjustment <- read_excel("data/fish_raw/IPHC/iphc-2024-fiss-hadj.xlsx")
      dat_IPHC <- IPHC(catch, adjustment, spc, sci_name)
      dat <- bind_rows(dat, dat_IPHC)  
    }
    
    #Some other columns
    dat$log_depth_scaled <- scale(log(dat$depth))
    dat$log_depth_scaled2 <- with(dat, log_depth_scaled ^ 2)
    dat$month <- month(dat$date)    #Some other columns
    
    if(mi==T){
    ##Calculate metabolic index
    #Calculate pO2 from umol kg
    dat$po2 <- calc_po2_sat(salinity=dat$salinity_psu, temp=dat$temperature_C, depth=dat$depth, oxygen=dat$O2_umolkg, lat=dat$latitude, long=dat$longitude, umol_m3=T, ml_L=F)
    
    ##Calculate inverse temp
    kelvin = 273.15
    boltz = 0.000086173324
    tref <- 12
    dat$invtemp <- (1 / boltz)  * ( 1 / (dat$temperature_C + 273.15) - 1 / (tref + 273.15))
    
    ###Calculate Metabolic index 
    ##Species parameters from Tim's paper
    #Read table in
    mi_pars <- read.csv("data/taxa_table.csv")
    mi_pars$Group <- tolower(mi_pars$Group)
    mi_pars <- filter(mi_pars,Group==taxa)
    
    V <- mi_pars$logV
    n <- mi_pars$n
    Eo <- c(mi_pars$Eolow, mi_pars$Eo, mi_pars$Eohigh)
    Ao <- 1/exp(V)
    
    #Abbreviated equation
    #dat$mi1 <- dat$po2*exp(Eo1* dat$invtemp)
    #dat$mi2 <- dat$po2*exp(Eo2* dat$invtemp)
    #dat$mi2 <- dat$po2*exp(Eo2* dat$invtemp)
    
    ##If using full equation and weight
    #Could link to Tim's estimates
    #MI_pars <- readRDS("metabolic_index-main/MI_pars.rds")
    #Use average weight across hauls of the species if haul is missing the weight
    dat$haul_weight <- ifelse(is.na(dat$haul_weight), mean(dat$median_weight, na.rm=T), dat$haul_weight)
    #If no weights at all for the species, use an average mid-size weight from the sablefish paper (# this works for the adult class (0.5 - 6 kg).  for the large adult, adjust)
    dat$haul_weight <- ifelse(is.na(dat$haul_weight), 3, dat$haul_weight)
    
    #Calculate metabolic index with each Eo
    #Calculate Ao from V above (using median)
    
    dat$mi1 <- calc_mi(Eo[1], Ao, dat$haul_weight,  n, dat$po2, dat$invtemp)
    dat$mi2 <- calc_mi(Eo[2], Ao, dat$haul_weight, n, dat$po2, dat$invtemp)
    dat$mi3 <- calc_mi(Eo[3], Ao, dat$haul_weight, n, dat$po2, dat$invtemp)
    
    #For average body size from each haul
    #dat$mi5 <- calc_mi(Eo1, Ao, dat$haul_weight,  n, dat$po2, dat$invtemp)
    
    #dat$po2_s <- (scale(dat$po2))
    #dat$mi1_s <-(scale(dat$mi1))
    #dat$mi2_s <-(scale(dat$mi2))
   # dat$mi3_s <-(scale(dat$mi3))
    #dat$mi5_s <- scale(dat$mi5)
    }
    
    #Reorder columns
    if(iphc==F){
      dat <- select(dat, scientific_name, common_name, survey_name, survey, region, year, date, depth, longitude, latitude, catch_weight, catch_count,  salinity_psu, temperature_C, sigma0_kgm3,do_mlpL, O2_umolkg, invtemp, po2, mi1,mi2,mi3, X, Y, haul_weight, pass, vessel, tow)
    }
    if(iphc==T){
      dat <- select(dat, scientific_name, common_name, survey_name, survey,region, year, date,  depth, longitude, latitude, catch_weight, catch_count,cpue_weight, cpue_count, salinity_psu, temperature_C, sigma0_kgm3, do_mlpL, O2_umolkg, invtemp, po2,mi1,mi2,mi3,X, Y, haul_weight, pass, vessel, tow)
      
    }
#Remove duplicates
    dat <- unique(dat)
    #Remove geometry
    dat$geometry <- NULL
    saveRDS(dat, file = paste("data/processed_data/fish/dat_", file_name, ".rds", sep=""))
  try(return(dat))
}

##Function to get separate species
prepare_data2 <- function(full_data, spc, sci_name, taxa, iphc, file_name, fancy){
  #Filter species
  dat_sub<- filter(dat, common_name==spc)
  dat_sub$catch_weight <- dat_sub$catch_weight/dat_sub$effort
  dat_sub$catch_numbers <- dat_sub$catch_numbers/dat_sub$effort
  
  #If species is not in the surveyjoin data, use the complete separate catch data
  #NWFSC
  if(nrow(dat_sub)==0){
    catch <- readRDS("data/fish_raw/NOAA/nwfsc_catch.rds")
    names(catch) = tolower(names(catch))
    catch$common_name <- tolower(catch$common_name)
    catch$scientific_name <- tolower(catch$scientific_name)
    #Fix dogfish
    catch$common_name <- ifelse(catch$common_name=="pacific spiny dogfish", "spiny dogfish", catch$common_name)
    catch$common_name <- ifelse(grepl("rougheye", catch$common_name), "rougheye rockfish", catch$common_name)
    #Filter to species
    catch = dplyr::filter(catch, scientific_name == sci_name)
    if(nrow(catch)>0){
    catch$event_id <- as.numeric(catch$trawl_id)
    catch$latitude <- catch$latitude_dd
    catch$longitude <- catch$longitude_dd
    # UTM transformation
    dat_ll = catch
    sp::coordinates(dat_ll) <- c("longitude_dd", "latitude_dd")
    sp::proj4string(dat_ll) <- sp::CRS("+proj=longlat +datum=WGS84")
    # convert to utm with spTransform
    dat_utm = sp::spTransform(dat_ll, 
                              sp::CRS("+proj=utm +zone=10 +datum=WGS84 +units=km"))
    # convert back from sp object to data frame
    catch = as.data.frame(dat_utm)
    catch = dplyr::rename(catch, X = coords.x1,
                        Y = coords.x2)
    catch$scientific_name <- sci_name
    catch$depth <- catch$depth_m
    catch$survey <- "nwfsc"
    catch$region <- "cc"
    catch$survey_name <- catch$project
    #Convert to per hectare 
    catch$catch_weight <- catch$cpue_kg_km2/100
    catch$catch_numbers <- catch$total_catch_numbers/catch$area_swept_ha
    dat_sub1 <-dplyr::select(catch, event_id, common_name, scientific_name, survey_name, survey, region, year, catch_weight, catch_numbers, longitude, latitude, X, Y,
                             depth)
    rm(catch)
    gc()
    }  else {
      dat_sub1 <- NULL
      rm(catch)
    }

    ##BC
    catch <- readRDS("data/fish_raw/BC/pbs-catch.rds")
    haul <- readRDS("data/fish_raw/BC/pbs-haul.rds")
    itis <- readRDS("data/fish_raw/BC/species-table.rds")
    catch <- left_join(catch,itis)
    #Rename dogfish
    catch$common_name <- ifelse(catch$common_name=="north pacific spiny dogfish", "spiny dogfish", catch$common_name)
    catch$common_name <- ifelse(grepl("rougheye", catch$common_name), "rougheye rockfish", catch$common_name)
    
    #Combine catch data with haul data
    catch <- dplyr::left_join(catch, haul, relationship = "many-to-many")
    
    catch = dplyr::filter(catch, scientific_name == sci_name)
  
    if(nrow(catch)>0){
    catch$longitude <- catch$lon_start
    catch$latitude <- catch$lat_start

    # UTM transformation
   dat_ll = catch
    sp::coordinates(dat_ll) <- c("longitude", "latitude")
    sp::proj4string(dat_ll) <- sp::CRS("+proj=longlat +datum=WGS84")
    # convert to utm with spTransform
    dat_utm = sp::spTransform(dat_ll, 
                              sp::CRS("+proj=utm +zone=10 +datum=WGS84 +units=km"))
    # convert back from sp object to  dataframe
    catch = as.data.frame(dat_utm)
    catch= dplyr::rename(catch, X = coords.x1,
                        Y = coords.x2)
    catch$depth <- dat$depth_m
    catch$survey <- "dfo"
    catch$region <- "bc"
    catch$longitude <- catch$lon_start
    catch$latitude <- catch$lat_start
    catch$depth <- catch$depth_m
    dat_sub2 <- dplyr::select(catch, event_id, common_name, scientific_name, survey_name, survey, region, year, catch_weight, catch_numbers, longitude, latitude,
                           depth, X, Y)
    catch$catch_weight <- catch$catch_weight/catch$effort
    catch$catch_numbers <- catch$catch_numbers/catch$effort
    rm(catch)
    gc()
    
    }else {
      dat_sub2 <- NULL
      rm(catch)
    }

  #AFSC
  bio2 <-readRDS("data/fish_raw/NOAA/ak_bts_goa_ebs_nbs_indivero_all_levels.RDS")
  catch <-readRDS("data/fish_raw/NOAA/ak_bts_goa_ebs_nbs_indivero_cpue_zerofilled.RDS")
  species <- bio2$species
  names(catch) = tolower(names(catch))
  names(species) =tolower(names(species))
  species$species_name <- tolower(species$species_name)
  species$common_name <- tolower(species$common_name)
  catch <- dplyr::left_join(catch, species)
  catch$common_name <- ifelse(catch$common_name=="north pacific spiny dogfish", "spiny dogfish", catch$common_name)
  catch$common_name <- ifelse(grepl("rougheye", catch$common_name), "rougheye rockfish", catch$common_name)
  
  catch <- filter(catch, common_name==spc)
  
  if(nrow(catch)>0){
  catch$longitude <- catch$longitude_dd_start
  catch$latitude <- catch$latitude_dd_start
  catch$scientific_name <- catch$species_name
  catch$event_id <- as.numeric(catch$hauljoin)
  catch$depth <- catch$depth_m
  catch$catch_weight <- catch$cpue_kgkm2/100
  catch$catch_numbers <- catch$count/catch$area_swept_km2
  catch$catch_numbers <- catch$catch_numbers/100
  catch$region <- tolower(catch$survey)
  # UTM transformation
  dat_ll = catch
  sp::coordinates(dat_ll) <- c("longitude", "latitude")
  sp::proj4string(dat_ll) <- sp::CRS("+proj=longlat +datum=WGS84")
  # convert to utm with spTransform
  dat_utm = sp::spTransform(dat_ll, 
                            sp::CRS("+proj=utm +zone=10 +datum=WGS84 +units=km"))
  # convert back from sp object to data frame
  catch = as.data.frame(dat_utm)
  catch = dplyr::rename(catch, X = coords.x1,
                      Y = coords.x2)
  catch$latitude <- catch$latitude_dd_start
  catch$longitude <- catch$longitude_dd_start
  dat_sub3 <- dplyr::select(catch, event_id, common_name, scientific_name,survey, region, year, catch_weight, catch_numbers, longitude, latitude,
                            depth, X, Y)
  rm(catch)
  gc()
  } else {
    dat_sub3 <- NULL
  }
  #Combine which exist
  # List the data frames
  dfs <- list(df1 = dat_sub1, df2 = dat_sub2, df3 = dat_sub3)
  
  # Filter out the NULL values (non-existing data frames)
  dfs <- dfs[sapply(dfs, is.data.frame)]
  
  # Combine the remaining data frames
  if (length(dfs) > 0) {
    dat_sub <- bind_rows(dfs)
    #Combine with in situ data
    insitu <- readRDS("data/processed_data/o2/insitu_combined.rds")
    #Remove IPHC, because this is already in the IPHC fish catch data
    insitu <- filter(insitu, survey!="iphc")
    insitu <- unique(insitu)
    #Just columns of interest
    dat_sub <- left_join(dat_sub, insitu[,c("temperature_C", "do_mlpL", "salinity_psu", "sigma0_kgm3", "O2_umolkg", "event_id")], by="event_id")
    
  } else {
    dat_sub <- NULL
  }
  }
  
  ##Add IPHC if needed
  if(iphc){
    if(spc=="pacific halibut"){
      catch <-  read_excel("data/fish_raw/IPHC/Set and Pacific halibut data.xlsx")
    }
    if(spc!="pacific halibut"){
      catch <-  read_excel("data/fish_raw/IPHC/Non_Pacific_halibut_data.xlsx")
    }
    adjustment <- read_excel("data/fish_raw/IPHC/iphc-2024-fiss-hadj.xlsx")
    dat_IPHC <- IPHC(catch, adjustment, spc, sci_name)
    if(nrow(dat_sub)>0){
    dat_sub <- bind_rows(dat_sub, dat_IPHC)
    } else {
      dat_sub <- dat_IPHC
      dat_sub$event_id <- NA
    }
    
  }
  
  #Add in mean weight per haul & proportional catch by size
  df1 <- try(length_expand_bc2(spc, sci_name))
  df2 <- try(length_expand_nwfsc2(spc, sci_name))
  df3 <- try(length_expand_afsc2(spc, sci_name))
  
  #Combine
  # List the data frames
  dfs <- list(df1 = df1, df2 = df2, df3 = df3)
  
  # Filter out the NULL values (non-existing data frames)
  dfs <- dfs[sapply(dfs, is.data.frame)]
  
  # Combine the remaining data frames
  if (length(dfs) > 0) {
    combined_df <- bind_rows(dfs)
    dat_sub <- left_join(dat_sub, combined_df, by="event_id")
  } else {
    combined_df <- NULL
    dat_sub$haul_weight <- NA
  }
  
  #Fill in missing haul weights
  #Use median weight across hauls of the species if haul is missing the weight
  dat_sub$haul_weight <- ifelse(is.na(dat_sub$haul_weight), median(dat_sub$haul_weight, na.rm=T), dat_sub$haul_weight)
  #If no weights at all for the species, use an average mid-size weight from the sablefish paper (# this works for the adult class (0.5 - 6 kg).  for the large adult, adjust)
  dat_sub$haul_weight <- ifelse(is.na(dat_sub$haul_weight), 2, dat_sub$haul_weight)
  
  ###Calculate MI
  ##Calculate metabolic index
  #Calculate pO2 from umol kg
  dat_sub$po2 <- calc_po2_sat(salinity=dat_sub$salinity_psu, temp=dat_sub$temperature_C, depth=dat_sub$depth, oxygen=dat_sub$O2_umolkg, lat=dat_sub$latitude, long=dat_sub$longitude, umol_m3=T, ml_L=F)
  
  ##Calculate inverse temp
  kelvin = 273.15
  boltz = 0.000086173324
  tref <- 12
  dat_sub$invtemp <- (1 / boltz)  * ( 1 / (dat_sub$temperature_C + 273.15) - 1 / (tref + 273.15))
  
  ###Calculate Metabolic index 
  ##Species parameters from Tim's paper
  #Read table in
  mi_pars <- read.csv("data/taxa_table.csv")
  mi_pars$Group <- tolower(mi_pars$Group)
  mi_pars <- filter(mi_pars,Group==taxa)
  
  V <- mi_pars$logV
  n <- mi_pars$n
  Eo <- c(mi_pars$Eolow, mi_pars$Eo, mi_pars$Eohigh)
  Ao <- 1/exp(V)
  
  #Abbreviated equation
  #dat$mi1 <- dat$po2*exp(Eo1* dat$invtemp)
  #dat$mi2 <- dat$po2*exp(Eo2* dat$invtemp)
  #dat$mi2 <- dat$po2*exp(Eo2* dat$invtemp)
  
  #Calculate each Ao
  dat_sub$mi1 <- calc_mi(Eo[1], Ao, dat_sub$haul_weight,  n, dat_sub$po2, dat_sub$invtemp, fancy)
  dat_sub$mi2 <- calc_mi(Eo[2], Ao, dat_sub$haul_weight, n, dat_sub$po2, dat_sub$invtemp, fancy)
  dat_sub$mi3 <- calc_mi(Eo[3], Ao, dat_sub$haul_weight, n, dat_sub$po2, dat_sub$invtemp, fancy)
  dat_sub$geometry <- NULL
  dat_sub <- as.data.frame(dat_sub)
  saveRDS(dat_sub, file = paste("data/processed_data/fish3/dat_", file_name, ".rds", sep=""))
}

lookup_taxa_t <- function(taxa.name, t.range, w.2.use) {
  all.dat <- load_data()
  #do the usual filtering of data
  all.dat <- filter_dat(all.dat)
  wref <- 5
  tref <- 15
  kb <-  8.617333262145E-5
  sim_beta <- readRDS("analysis/taxa_sims.RDS")
  ltaxa.name <- tolower(taxa.name)
  taxa.name<- tolower(taxa.name)
  lookup.taxa <- ltaxa.name %in% tolower(sim_beta$Group)
  sim_beta$Group <- tolower(sim_beta$Group)
  if(!lookup.taxa){
    ltaxa.name <- "teleostei"
    lookup.taxa <- ltaxa.name %in% tolower(sim_beta$Group)
    taxa.name<- ltaxa.name
  }
  
  if(lookup.taxa) {
    options(warn = -1)
    sims <- dplyr::filter(sim_beta, Group == taxa.name)
    if (sims$level[1] == 1) taxa_dat <- dplyr::filter(all.dat, tolower(Class) == ltaxa.name)
    if (sims$level[1] == 2) taxa_dat <- dplyr::filter(all.dat, tolower(Order) == ltaxa.name)
    if (sims$level[1] == 3) taxa_dat <- dplyr::filter(all.dat, tolower(Family) == ltaxa.name)
  }
  
  Wmed <- median(w.2.use/ wref)
  
  pcrit_df <- tibble(Temp = NA, 
                     Pcritlower90 = NA,
                     Pcritupper90= NA,
                     Pcritlower50 = NA,
                     Pcritupper50= NA)
  
  for (i in seq_along(t.range)) {
    # calculate Pcrit for given temperature
    inv.temp <- (1 / kb) * (1 / kelvin(t.range[i]) - 1 / kelvin(tref) )
    pcrit <- sims$logV - sims$n * log(Wmed) - sims$Eo * inv.temp
    pcrit_90 <- get_x_range(x = exp(pcrit), p = 0.9)
    pcrit_50 = get_x_range(x = exp(pcrit), p = 0.5)
    pcrit_df <- pcrit_df %>% add_row(Temp = t.range[i],
                                     Pcritlower90 = pcrit_90[1], 
                                     Pcritupper90 = pcrit_90[2],
                                     Pcritlower50 = pcrit_50[1],
                                     Pcritupper50 = pcrit_50[2]
    )
    
  }
  # delete placeholder NA row
  pcrit_df <- na.omit(pcrit_df)
  # smooth lower and upper bounnds
  lower90mod = gam(Pcritlower90~ s(Temp), data = pcrit_df)
  pcrit_df$lower90s = predict.gam(lower90mod, data = pcrit_df)
  upper90mod = gam(Pcritupper90~ s(Temp), data = pcrit_df)
  pcrit_df$upper90s = predict.gam(upper90mod, data = pcrit_df)
  
  lower50mod = gam(Pcritlower50~ s(Temp), data = pcrit_df)
  pcrit_df$lower50s = predict.gam(lower50mod, data = pcrit_df)
  upper50mod = gam(Pcritupper50~ s(Temp), data = pcrit_df)
  pcrit_df$upper50s = predict.gam(upper50mod, data = pcrit_df)
  
  
  return(pcrit_df)
} 

# master function to solve for the density that produces a tail probability equal
# to target, and to give the limits of that range (the inner p percentile interval)
get_x_range <- function(x, p) {
  dens <- uniroot(f = fit_dens, interval = c(0.01, 1), x = x, p = p)
  smooth <- bkde (x)
  index <- which(smooth$y >=dens$root)
  xbounds = c(min(smooth$x[index]), max(smooth$x[index]))
  return(xbounds)
}
# Functions for calculating credibility interval
# function to return the probability in the tails below a theshold density
getprob <- function(dens, x) {
  smooth <- bkde (x)
  deltax <- diff(smooth$x)[1]
  # turn into probability
  probs <- smooth$y * deltax
  index <- which(smooth$y <dens)
  return(sum(probs[index]))
}
# function to get the difference between tail probability and target probability
fit_dens <- function(dens, x, p) {
  pstar <- getprob(dens, x)
  target_p <- 1 - p
  return(target_p - pstar)
}

# packages that might not be installed - install any not already installed
packages <- c("devtools", "rfishbase", "gsw")
install.packages(setdiff(packages, rownames(installed.packages())))
library(dplyr)
library(devtools)
library(gsw)
library(readxl)
library(stringr)
library(lubridate)
library(sf)

calc_po2_sat <- function(salinity, temp, depth, oxygen, lat, long, umol_m3, ml_L) {
  # Input:       S = Salinity (pss-78)
  #              T = Temp (deg C) ! use potential temp
  #depth is in meters
  
  #Pena et al. ROMS and GLORYS are in mmol per m^3 (needs to be converted to umol per kg) (o2 from trawl data was in mL/L, so had to do extra conversions)
  gas_const = 8.31
  partial_molar_vol = 0.000032
  kelvin = 273.15
  boltz = 0.000086173324
  
  #convert mmol to umol
  umol_m3 <- oxygen*1000
  #convert m3 to l
  umol_l <- umol_m3/1000
  #convert from molality (moles per volume) to molarity (moles per mass)
  #1 L of water = 1 kg of water, so no equation needed?? Right??
  o2_umolkg <- umol_l * 1/1
  
  SA = gsw_SA_from_SP(salinity,depth,long,lat) #absolute salinity for pot T calc
  pt = gsw_pt_from_t(SA,temp,depth) #potential temp at a particular depth
  #this is for if using data that has oxygen in ml/L
  #CT = gsw_CT_from_t(SA,temp,depth) #conservative temp
  #sigma0 = gsw_sigma0(SA,CT)
  # o2_umolkg = oxygen*44660/(sigma0+1000)
  
  O2_Sat0 = gsw_O2sol_SP_pt(salinity,pt)
  
  #= o2satv2a(sal,pt) #uses practical salinity and potential temp - solubity at p =1 atm
  press = exp(depth*10000*partial_molar_vol/gas_const/(temp+kelvin))
  O2_satdepth = O2_Sat0*press
  
  #solubility at p=0
  sol0 = O2_Sat0/0.209
  sol_Dep = sol0*press
  po2 = o2_umolkg/sol_Dep
  po2 <- po2 * 101.325 # convert to kPa
  return(po2)
  
}

# calc o2 solubility, relies on o2 in umol/kg
gsw_O2sol_SP_pt <- function(sal,pt) {
  x = sal
  pt68 = pt*1.00024
  y = log((298.15 - pt68)/(273.15 + pt68))
  
  a0 =  5.80871
  a1 =  3.20291
  a2 =  4.17887
  a3 =  5.10006
  a4 = -9.86643e-2
  a5 =  3.80369
  b0 = -7.01577e-3
  b1 = -7.70028e-3
  b2 = -1.13864e-2
  b3 = -9.51519e-3
  c0 = -2.75915e-7
  
  O2sol = exp(a0 + y*(a1 + y*(a2 + y*(a3 + y*(a4 + a5*y)))) + x*(b0 + y*(b1 + y*(b2 + b3*y)) + c0*x))
  return(O2sol)
}

# calc o2 solubility, relies on o2 in umol/kg
calc_o2_sol <- function(sal, temp, depth, long, lat) {
  SA = gsw_SA_from_SP(sal,depth,long,lat) #absolute salinity for pot T calc
  pt = gsw_pt_from_t(SA,temp,depth) #potential temp at a particular depth
  pt68 = pt*1.00024
  y = log((298.15 - pt68)/(273.15 + pt68))
  
  a0 =  5.80871
  a1 =  3.20291
  a2 =  4.17887
  a3 =  5.10006
  a4 = -9.86643e-2
  a5 =  3.80369
  b0 = -7.01577e-3
  b1 = -7.70028e-3
  b2 = -1.13864e-2
  b3 = -9.51519e-3
  c0 = -2.75915e-7
  
  O2sol = exp(a0 + y*(a1 + y*(a2 + y*(a3 + y*(a4 + a5*y)))) + sal*(b0 + y*(b1 + y*(b2 + b3*y)) + c0*sal))
  return(O2sol)
}

calc_sigma <- function(s, t, p) {
  constants <- list(B0 = 8.24493e-1, B1 = -4.0899e-3, B2 = 7.6438e-5, B3 = -8.2467e-7, B4 = 5.3875e-9, C0 = -5.72466e-3, C1 = 1.0227e-4, C2 = -1.6546e-6, D0 = 4.8314e-4, A0 = 999.842594, A1 = 6.793952e-2, A2 = -9.095290e-3, A3 = 1.001685e-4, A4 = -1.120083e-6, A5 = 6.536332e-9, FQ0 = 54.6746, FQ1 = -0.603459, FQ2 = 1.09987e-2, FQ3 = -6.1670e-5, G0 = 7.944e-2, G1 = 1.6483e-2, G2 = -5.3009e-4, i0 = 2.2838e-3, i1 = -1.0981e-5, i2 = -1.6078e-6, J0 =1.91075e-4, M0 = -9.9348e-7, M1 = 2.0816e-8, M2 = 9.1697e-10, E0 = 19652.21, E1 = 148.4206, E2 = -2.327105, E3 = 1.360477e-2, E4 = -5.155288e-5, H0 = 3.239908, H1 = 1.43713e-3, H2 = 1.16092e-4, H3 = -5.77905e-7, K0 = 8.50935e-5, K1 =-6.12293e-6, K2 = 5.2787e-8)
  
  list2env(constants, environment())
  constant_names <- names(constants)
  t2 = t*t
  t3 = t*t2
  t4 = t*t3
  t5 = t*t4
  #  if (s <= 0.0) s = 0.000001
  s32 = s^ 1.5
  p = p / 10.0 # convert decibars to bars */
  sigma = A0 + A1*t + A2*t2 + A3*t3 + A4*t4 + A5*t5 + (B0 + B1*t + B2*t2 + B3*t3 + B4*t4)*s + (C0 + C1*t + C2*t2)*s32 + D0*s*s
  kw = E0 + E1*t + E2*t2 + E3*t3 + E4*t4
  aw = H0 + H1*t + H2*t2 + H3*t3
  bw = K0 + K1*t + K2*t2
  k = kw + (FQ0 + FQ1*t + FQ2*t2 + FQ3*t3)*s + (G0 + G1*t + G2*t2)*s32 + (aw + (i0 + i1*t + i2*t2)*s + (J0*s32))*p + (bw + (M0 + M1*t + M2*t2)*s)*p*p
  val = 1 - p / k
  sigma = sigma / val - 1000.0
  rm(constant_names)
  return(sigma)
}

# convert ml /l to umol / kg
convert_o2 <- function(o2ml_l, sigma){
  (o2ml_l * 44.660)/((1000 + sigma)/1000)
}

rsq <- function(x, y) cor(x,y) * cor(x,y)

#Function to calculate breakpoint effect
breakpoint_calc <- function(x, b_slope, b_thresh){
  if (x < b_thresh) {
    pred = x  * b_slope
    #pred = (x-b_thresh) * b_slope
  } else {
    pred=b_thresh * b_slope
    #pred = 0
  }
  return(pred)
}

update_df <- function(species_df, taxa.name, w, pcrit_type = "smr", rep, ParentChild_gz) {
  # get parameter estimates
  taxa_estimates <- estimate_taxa(taxa.name, 
                                  w = w ,
                                  temperature = 10,
                                  method = pcrit_type,
                                  rep  = model.fit$rep, 
                                  ParentChild_gz = model.fit$ParentChild_gz)
  # use pmap to run calc_p_po2
  betas <- taxa_estimates$parameters[,1]
  var_covar <- taxa_estimates$var_covar
  species_df$p_pcrit <- pmap_dbl(list(po2 = species_df$po2, temp = species_df$temperature_C),
                                 ~calc_p_po2(po2 = ..1, temp = ..2, w = w,  pcrit_type = pcrit_type, betas = betas, var_covar = var_covar)
  )
  
  return(species_df)
}


## Profile-likelihood breakpoint helpers (sdmTMB::profile_breakpt()) ----

#sdmTMB fit from either sdmTMB() or profile_breakpt() output
get_fit <- function(x) {
  if (inherits(x, "breakpt_profile")) x$fit else x
}

#Refit at evenly spaced fixed cutpoints where the profile is within `within`
#NLL units of its minimum, weighting each fit by its likelihood. See the
#threshold-models article in sdmTMB ("To include uncertainty in the
#breakpoint..."). Each end is extended to where a straight line to the next
#profile fit out crosses `within`. Returns one row per refit with the cut, NLL,
#weight, and the estimated slope parameter and its SE (cutpoint fixed). With
#breakpt(positive = TRUE), that parameter is log(slope) (`log_scale` is TRUE).
#`envir` is where the original call's variables (formula, priors, etc.) live.
breakpt_refits <- function(pr, within = 4, n_cuts = 10, envir = parent.frame()) {
  prof <- pr$profile[!is.na(pr$profile$nll), ]
  prof <- prof[order(prof$cut), ]
  inside <- which(prof$dnll < within)
  lo <- min(inside)
  hi <- max(inside)
  cross <- function(i, j) {
    prof$cut[i] + (within - prof$dnll[i]) *
      (prof$cut[j] - prof$cut[i]) / (prof$dnll[j] - prof$dnll[i])
  }
  lwr <- if (lo > 1) cross(lo, lo - 1) else prof$cut[lo]
  upr <- if (hi < nrow(prof)) cross(hi, hi + 1) else prof$cut[hi]
  cuts <- seq(lwr, upr, length.out = n_cuts)

  fit <- pr$fit
  ctrl <- fit$control
  ctrl$start <- fit$tmb_obj$env$parList() # start from the best fit
  ctrl$map <- list(b_threshold = factor(c(1, NA))) # estimate slope, fix cut
  ctrl$multiphase <- FALSE
  call <- fit$call
  call$mesh <- fit$spde
  call$data <- fit$data

  out <- lapply(cuts, function(cut) {
    ctrl$start$b_threshold[2] <- cut
    call$control <- ctrl
    f <- try(suppressMessages(eval(call, envir)), silent = TRUE)
    if (inherits(f, "try-error") || is.null(f$sd_report)) {
      return(data.frame(cut = cut, nll = NA, par = NA, par_se = NA))
    }
    sdr <- f$sd_report
    i <- names(sdr$par.fixed) == "b_threshold"
    data.frame(
      cut = cut, nll = f$model$objective,
      par = unname(sdr$par.fixed[i]),
      par_se = unname(sqrt(diag(sdr$cov.fixed))[i])
    )
  })
  out <- do.call(rbind, out)
  out$log_scale <- identical(as.integer(fit$tmb_data$threshold_func), 3L)
  out$slope <- if (out$log_scale[1]) exp(out$par) else out$par
  ok <- is.finite(out$nll) & is.finite(out$par_se)
  #e.g., a slope near 0 with breakpt(positive = TRUE), where log(slope) runs away
  if (!any(ok)) stop("No fixed-cutpoint refit has a finite NLL and slope SE.")
  out$weight <- 0
  out$weight[ok] <- exp(min(out$nll[ok]) - out$nll[ok])
  out$weight <- out$weight / sum(out$weight)
  out
}

#n draws of the slope and cutpoint from profile_breakpt() output with $refits
#from breakpt_refits(). Draws are split among the fixed-cut refits by weight;
#within a refit, the slope is drawn from its estimate and SE (on the log scale
#for breakpt(positive = TRUE), so slopes stay positive).
#Returns a matrix with columns slope and cut (scaled covariate units).
#If breakpt_refits() failed, repeats the point estimate (no uncertainty).
breakpt_draws <- function(pr, n = 100) {
  refits <- pr$refits
  if (!is.data.frame(refits)) {
    warning("No breakpt_refits() results in `pr$refits`; using the point estimate.")
    return(cbind(slope = rep(pr$estimate$slope, n), cut = rep(pr$estimate$cut, n)))
  }
  k <- sample(nrow(refits), n, replace = TRUE, prob = refits$weight)
  slope <- rnorm(n, refits$par[k], refits$par_se[k])
  if (refits$log_scale[1]) slope <- exp(slope)
  cbind(slope = slope, cut = refits$cut[k])
}


#Save a fitted model with strip_sdmTMB() applied first (smaller .rds files).
#The joint precision matrix is dropped too, so predict(nsim = ) won't work on
#saved fits; 06 draws from the fixed-cut refits instead.
#Handles sdmTMB fits and profile_breakpt() results (a list holding the fit in
#`$fit`); try-errors and anything else are saved as is.
save_model <- function(x, file) {
  if (inherits(x, "sdmTMB")) {
    x <- sdmTMB::strip_sdmTMB(x, strip_joint_precision = TRUE)
  } else if (is.list(x) && inherits(x$fit, "sdmTMB")) {
    x$fit <- sdmTMB::strip_sdmTMB(x$fit, strip_joint_precision = TRUE)
  }
  saveRDS(x, file = file)
}
