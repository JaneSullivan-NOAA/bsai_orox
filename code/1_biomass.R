# Survey and random effects model biomass estimates
# Contact: jane.sullivan@noaa.gov
# Last updated: Sep 2020

# Set up ----

# Assessment year (most recent year with complete data set)
YEAR <- 2019
dat_path <- paste0("data/", YEAR) # directory where source data is contained
out_path <- paste0("results/", YEAR) # directory for results/output
dir.create(out_path)

libs <- c("tidyverse", "R2admb")
if(length(libs[which(libs %in% rownames(installed.packages()) == FALSE )]) > 0) {
  install.packages(libs[which(libs %in% rownames(installed.packages()) == FALSE)])}
lapply(libs, library, character.only = TRUE)

# Data ----

biom <- read_csv(paste0(dat_path, "/bsai_orox_biomass_", YEAR, ".csv"))

# Summary ----

# Stock structure
biom <- biom %>% mutate(full_complex = "All-OR",
                        complex_bySST = ifelse(species_code == 30020, "SST", "non-SST"),
                        complex_bySSTandDusky = ifelse(species_code == 30020, "SST",
                                                       ifelse(species_code == 30152, "Dusky", "non-SST-or-Dusky"))) 

biom %>%
  pivot_longer(cols = c("full_complex", "complex_bySST", "complex_bySSTandDusky"), 
               names_to = "complex", values_to = "species") %>% 
  group_by(complex, species, survey, area, year) %>% 
  summarize(biomass = sum(biomass),
            var = sum(var)) %>% 
  ungroup() %>% 
  mutate(cv = sqrt(var) / biomass) %>% 
  arrange(complex, year) %>% 
  View() #print(n = Inf)

ebs_slope %>% 
  left_join(spp %>% select(species_code, common_name)) %>% # species name
  mutate(complex_bySST = ifelse(species_code == 30020, "SST", "non-SST")) %>% 
  group_by(year, complex_bySST) %>% 
  summarize(biomass = sum(biomass),
            var = sum(var)) %>% 
  ungroup() %>% 
  mutate(cv = sqrt(var) / biomass) %>% 
  arrange(complex_bySST, year) %>% 
  print(n = Inf)
# AI non-SST and AI SST
# BS Slope non-SST and SST BS Slope
# BS Shelf all OR