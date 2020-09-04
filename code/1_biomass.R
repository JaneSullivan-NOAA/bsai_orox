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

full_biom <- read_csv(paste0(dat_path, "/bsai_orox_biomass_", YEAR, ".csv"))

# Summary ----

# Stock structure for assessment: splits by SST and non-SST, also AI survey
# split into Southern Bering Sea and AI
biom <- full_biom %>% 
  mutate(species = ifelse(species_code == 30020, "SST", "non-SST"),
         region = ifelse(survey == "AI" & area == "SBS", "SBS",
                         ifelse(survey == "AI" & grepl("AI", area), "AI",
                                survey))) %>% 
  group_by(species, region, year) %>% 
  summarize(biomass = sum(biomass),
            var = sum(var)) %>% 
  ungroup() %>% 
  # Currently years with 0 biomass are removed from RE model
  mutate(cv = ifelse(biomass == 0, 0, sqrt(var) / biomass)) %>% 
  arrange(species, region, year) 

full_estimates <- biom %>% 
  pivot_wider(id_cols = c("species", "year"), names_from = "region", 
              values_from = "biomass", values_fill = 0) 

full_cv <- biom %>% 
  pivot_wider(id_cols = c("species", "year"), names_from = "region", 
              values_from = "cv", values_fill = 0) 

unname(srv_est) # gets rid of column names
srv_est[srv_est == 0] <- "-9" # ADMB flag 
unname(srv_SE)
srv_SE[srv_SE == 0] <- "-9"

cat("# Model start and end years","\n",yrs,"\n",
    "# Number of survey indices fit (i.e., regions/depth strata)","\n",nregs,"\n",
    "# Number or process error parameters","\n",1,"\n",
    "# Process error index","\n",PEI,"\n",
    "# Number of surveys","\n",nobs,"\n",
    "# Survey years","\n",yrs_srv,"\n",
    "# Survey biomass","\n", sep=" ",file="re.dat")

write.table(srv_est, file = "re.dat", sep = " ", append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(paste0("# Survey biomass SE"), file = "re.dat", sep = " ", append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(srv_SE, file = "re.dat", sep = " ", append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE)


# other ways to split:
full_biom %>% 
  mutate(full_complex = "All-OR",
       complex_bySST = ifelse(species_code == 30020, "SST", "non-SST"),
       complex_bySSTandDusky = ifelse(species_code == 30020, "SST",
                                      ifelse(species_code == 30152, "Dusky", "non-SST-or-Dusky"))) %>% 
  pivot_longer(cols = c("full_complex", "complex_bySST", "complex_bySSTandDusky"), 
               names_to = "complex", values_to = "species") 
  
# AI non-SST and AI SST
# BS Slope non-SST and SST BS Slope
# BS Shelf all OR