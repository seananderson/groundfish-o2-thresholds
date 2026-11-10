library(readxl)
library(sdmTMB)
library(dplyr)
library(tidyr)
library(gt)
library(openxlsx2)
library(here)

here()

#Load functions
source("code/helper_funs.R")

#Output folder
output_folder <- "output"

#Depth structure
cubic <- T
quadratic <- F

#Models
models <- c("model1", "model2", "model3", "model4", "model5")
mi_models <- c("model3", "model4", "model5")
cubic <- T
if(cubic){
models <- c("model7", "model8", "model13", "model14","model15")
mi_models <- c("model13", "model14", "model15")
}
aic_table = as.data.frame(matrix(NA, 1, length(models)+5))
error_table <- as.data.frame(matrix(NA, 1, length(models)+5))
dat_names <- c("cc", "bc", "goa", "ebs", "coastwide", "cc _iphc", "bc _iphc", "goa _iphc", "ebs _iphc", "coastwide _iphc")
species <- read_excel("data/species_table.xlsx")
species$common_name <- tolower(species$common_name)
species$scientific_name <- tolower(species$scientific_name)
species <- unique(species$common_name)

for(i in 1:length(species)) {
  this_species = species[i]
  print(this_species)
  for(h in 1:length(dat_names)){
    temp <- matrix(NA, 1, length(models)+5)
    errors <- matrix(NA, 1, length(models)+5)
    this_dat <- dat_names[h]
    print(this_dat)
    dat <- try(readRDS(file = paste0("output/",output_folder, "/", this_species, "_", this_dat, "_dat.rds")))
    if(class(dat)=="try-error"){
      errors[1,1:length(models)] <- paste("no data")
    }
    if(class(dat)!="try-error"){
      for(j in 1:length(models)) {
        this_model <- models[j]
        print(this_model)
        fit <- try(readRDS(file = paste0("output/", output_folder, "/", this_species,"_",this_dat, "_", this_model,".rds")))
        if(class(fit)=="try-error"){
          errors[1,j] <- paste0(fit)
        }
        if(class(fit)!="try-error"){
          s <- try(sanity(get_fit(fit), silent=TRUE))
          if(class(s)!="try-error"){
            if((s$hessian_ok + s$eigen_values_ok +s$nlminb_ok)==3){
              temp[1,j] <- AIC(fit)
              if(this_model %in% mi_models){
                pars <- as.data.frame(tidy(get_fit(fit), effects="fixed", conf.int=T))
                #Cutpoint is fixed at its profile estimate, so it has no SE; check the slope
                thresh <- filter(pars,grepl("slope", term))
                sd_thresh <- thresh$std.error
                if(!is.finite(sd_thresh)){
                  temp[1,j] <- NA
                }
              }
            }
            if((s$hessian_ok + s$eigen_values_ok +s$nlminb_ok)< 3) {
            temp[1,j] <- NA
            errors[1,j] <- paste0(s$hessian_ok, ", ", s$eigen_values_ok, ", ", s$nlminb_ok)
          } 
          }
        }
      }
      mins <- apply(temp, 1, min, na.rm=T)
      for(c in 1:length(models)){
        temp[1,c] <- abs(temp[1,c]-mins[1])
      }
      
      temp[1,length(models)+1] <- this_species
      temp[1,length(models)+2] <- this_dat
      temp[1,length(models)+3]  <- nrow(filter(dat, catch>0))
      temp[1,length(models)+4]  <- length(unique(dat$year))
      temp[1,length(models)+5] <- length(unique(dat$region))
    
      #Add to table
      aic_table <- bind_rows(aic_table, as.data.frame(temp))
    }
    errors[1,length(models)+1] <- this_species
    errors[1,length(models)+2] <- this_dat
    error_table <- bind_rows(error_table, as.data.frame(errors))
  }
}

aic_table <- aic_table[2:nrow(aic_table),]
colnames(aic_table) <- c("model7", "model8",  "model13", "model14","model15", "species", "data type", "N obs", "N years", "N regions")
#How many in each model type?
counts <- aic_table %>%
  summarise(across(everything(), ~ sum(. == 0, na.rm=T)))


#Round
aic_table <- aic_table %>%
  mutate(across(1:length(models), as.numeric))
aic_table <- aic_table %>%
  mutate(across(1:length(models), round, 8))


#Replace NA with --
aic_table2 <- aic_table
aic_table2[is.na(aic_table2)] <- "--"
aic <- aic_table2 %>%
  gt() %>%
  # tab_style(
  #   style = cell_fill(color = "gray"),
  #   locations = cells_body(columns = model1, rows = model1==0)
  # ) %>%
  # tab_style(
  #   style = cell_fill(color = "lightcoral"),
  #   locations = cells_body(columns = model2, rows = model2==0)
  # ) %>%
  # tab_style(
  #   style = cell_fill(color = "thistle2"),
  #   locations = cells_body(columns = model3, rows = model3==0)
  # ) %>%
  # tab_style(
  #   style = cell_fill(color = "thistle"),
  #   locations = cells_body(columns = model4, rows = model4==0)
  # ) %>%
  # tab_style(
  #   style = cell_fill(color = "thistle3"),
  #   locations = cells_body(columns = model5, rows = model5==0)
  # ) %>%
  # tab_style(
  #   style = cell_fill(color = "lightblue"),
  #   locations = cells_body(columns = model6, rows = model6==0)
  # ) %>%
  tab_style(
    style = cell_fill(color = "gray"),
    locations = cells_body(columns = model7, rows = model7==0)
  ) %>%
  tab_style(
    style = cell_fill(color = "indianred"),
    locations = cells_body(columns = model8, rows = model8==0)
  ) %>%
  # tab_style(
  #   style = cell_fill(color = "thistle2"),
  #   locations = cells_body(columns = model9, rows = model9==0)
  # ) %>%
  # tab_style(
  #   style = cell_fill(color = "thistle"),
  #   locations = cells_body(columns = model10, rows = model10==0)
  # ) %>%
  # tab_style(
  #   style = cell_fill(color = "thistle3"),
  #   locations = cells_body(columns = model11, rows = model11==0)
  # ) %>%
  # tab_style(
  #   style = cell_fill(color = "lightblue"),
  #   locations = cells_body(columns = model12, rows = model12==0)
  # ) %>%
      tab_style(
        style = cell_fill(color = "thistle2"),
        locations = cells_body(columns = model13, rows = model13==0)
      ) %>%
      tab_style(
        style = cell_fill(color = "thistle"),
        locations = cells_body(columns = model14, rows = model14==0)
      ) %>%
      tab_style(
        style = cell_fill(color = "thistle3"),
        locations = cells_body(columns = model15, rows = model15==0)
      )

#save
write_xlsx(aic_table, file=paste0("output/", output_folder, "/aic_table.xlsx"))
saveRDS(aic_table, file=paste0("output/", output_folder, "/aic_table.rds"))
gtsave(aic, file=paste0("output/", output_folder, "/aic_table.html"))

#Evaluate errors
#Remove species with no data
error_table <- error_table[2:nrow(error_table),1:7]
colnames(error_table) <- c("model7", "model8",  "model13", "model14","model15", "species", "data type")


