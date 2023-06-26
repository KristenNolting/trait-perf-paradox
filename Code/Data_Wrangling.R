# Loading Protea trait data from Nolting et a. (2021). Cleaning code for present analyses.

# load libraries
library(dplyr)
library(tidyverse)

# clean working environment
rm(list = ls())

# Read in full dataframe (from Nolting et al. 2021)
full_data <- read.csv("Data/Individual_Data_for_Analysis.csv", na.strings = ".")


# Drop populations not used (these were not presented in Nolting et al. 2021 either)
# Create variable for Instantaneous WUE
# Select relevant columns
subset_data <- full_data %>%
  filter(Species != "LDSA") %>%
  filter(Site != "McGregor") %>% 
  filter(Site != "Jonaskop_2") %>%
  mutate(WUE_Instantaneous = Photo/Trmmol) %>%
  select(Species, Site, Ind_No, Age, Height_in, Canopy_area_in2, Total_Fruit_No, HV_m2_cm2, WD, 
         Ks, KL_Units, TotalAssim, BW_dry_avg, Photo, Trmmol, Cond, WUE_Intrinsic,
         Thick_cm, Area_cm2, LMA_g.cm2, LWR, LD_g.cm3, Length_mm_Top, Density_mm2_Top, 
         SPI_Top, WUE_Instantaneous) %>%
  drop_na()

# Dropping the Cederberg, Protea repens population (presented in Nolting et al. 2021, dropped here for simplicity)
subset_data$keep <- "yes"
subset_data$keep[subset_data$Site == "Cederberg" & subset_data$Species == "PRRE"] <- "no"

Protea_data <- subset_data %>%
  filter(keep == "yes") %>%
  droplevels()
# check
Protea_data[Protea_data$Site == "Cederberg" & Protea_data$Species == "PRRE",]

# This centers the trait data at zero and scales it by 1 standard deviation
# We use centered and scaled data in all analyses
Protea_data <- Protea_data %>%
  mutate(
    BW_scaled = scale(BW_dry_avg),
    WD_scaled = scale(WD),
    Leaf_Area_scaled = scale(Area_cm2),
    LMA_scaled = scale(LMA_g.cm2),
    LD_scaled = scale(LD_g.cm3),
    LWR_scaled = scale(LWR),
    Stom_L_scaled = scale(Length_mm_Top),
    Stom_D_scaled = scale(Density_mm2_Top),
    SPI_scaled = scale(SPI_Top),
    Ks_scaled = scale(Ks),
    LSC_scaled = scale(KL_Units),
    Photo_scaled = scale(Photo),
    Total_Assim_scaled = scale(TotalAssim),
    Cond_scaled = scale(Cond),
    WUE_Instan_scaled = scale(WUE_Instantaneous),
    Height_scaled = scale(Height_in),
    Canopy_scaled = scale(Canopy_area_in2)
  )

# select variables used in subsequent analyses
Protea_data <- Protea_data %>%
  select(Species, Site, BW_scaled, WD_scaled, Leaf_Area_scaled, LMA_scaled, LD_scaled,
         LWR_scaled, Stom_L_scaled, Stom_D_scaled, Ks_scaled, LSC_scaled, Photo_scaled,
         Total_Assim_scaled, Cond_scaled, WUE_Instan_scaled)

#check
dim(Protea_data)

# save as R object to load in associated analysis files
save(Protea_data, file = "Data/Protea_data.RData")

# clean up, as we'll work with Protea_data
rm(full_data)
rm(subset_data)
