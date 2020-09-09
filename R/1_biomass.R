# Survey and random effects model biomass estimates
# Contact: jane.sullivan@noaa.gov
# Last updated: Sep 2020

# Set up ----

# Assessment year (most recent year with complete data set)
YEAR <- 2019

libs <- c("tidyverse", "R2admb")
if(length(libs[which(libs %in% rownames(installed.packages()) == FALSE )]) > 0) {
  install.packages(libs[which(libs %in% rownames(installed.packages()) == FALSE)])}
lapply(libs, library, character.only = TRUE)

# Data ----

full_biom <- read_csv(paste0(dat_path, "/bsai_orox_biomass_", YEAR, ".csv"))

# Directory setup ----

root <- getwd()
dat_path <- paste0("data/", YEAR) # directory where source data is contained
out_path <- paste0("results/", YEAR) # directory for results/output
dir.create(out_path)
tpl_dir <- file.path("admb")
files <- list.files(file.path(tpl_dir)) # searchable file list
tpl <- file.path(tpl_dir, files[which(grepl(".tpl", files))]) # find .tpl
exe <- file.path(tpl_dir, files[which(grepl(".exe", files))]) # find .exe
std <- paste0("RE.std")

setwd(tpl_dir)

# If there's not an exe for the RE model, create one.
if(!file.exists("RE.exe")) {
  setup_admb()
  compile_admb("RE", verbose = TRUE)
}

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

input_biom <- biom %>% 
  pivot_wider(id_cols = c("species", "year"), names_from = "region", 
              values_from = "biomass", values_fill = 0) %>% 
  mutate(value = "mu") %>% 
  bind_rows(biom %>% 
              pivot_wider(id_cols = c("species", "year"), names_from = "region", 
                          values_from = "cv", values_fill = 0) %>% 
              mutate(value = "cv")) %>% 
  arrange(species, value, year)

# Run model ----

spp <- unique(full_estimates$species)
yrs <- unique(input_biom$year)
areas <- input_biom %>% select(-c(species, year, value)) %>% colnames
process_error <- rep(1, length(areas))

for(i in 1:length(spp)) {
  
  # Create new directory and copy over tpl and exe files
  setwd(root)
  spp_out <- file.path(paste0(out_path, "/biomass_", spp[i]))
  dir.create(spp_out)
  file.copy(from = tpl, to = spp_out, overwrite = TRUE)
  file.copy(from = exe, to = spp_out, overwrite = TRUE)
  setwd(spp_out)

  # Biomass estimates
  srv_est <- input_biom %>% 
    filter(species == spp[i] & value == "mu") %>% 
    select(-c(species, year, value)) %>% 
    as.matrix()
  srv_est <- unname(srv_est) # gets rid of column names
  srv_est[srv_est == 0] <- "-9" # ADMB flag
  
  # Biomass variability
  srv_cv <- input_biom %>% 
    filter(species == spp[i] & value == "cv") %>% 
    select(-c(species, year, value)) %>% 
    as.matrix()
  srv_cv <- unname(srv_cv) # gets rid of column names
  srv_cv[srv_cv == 0] <- "-9" # ADMB flag
  
  # Write dat file
  cat("# Model start and end years","\n", min(yrs), max(yrs), "\n",
      "# Number of survey indices fit (i.e., regions/depth strata)","\n", length(areas), "\n",
      "# Number or process error parameters","\n",1,"\n",
      "# Process error index","\n", process_error, "\n",
      "# Number of surveys","\n", length(yrs), "\n",
      "# Survey years","\n", yrs, "\n",
      "#", areas, "\n", 
      "# Survey biomass","\n", sep = " ", file = "RE.dat")
  
  write.table(srv_est, file = "RE.dat", sep = " ", append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE)
  write.table(paste0("# Survey biomass SE"), file = "RE.dat", sep = " ", append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE)
  write.table(srv_cv, file = "RE.dat", sep = " ", append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE)

  # remove pre-existing std file if it exists and run admb
  if (file.exists(std)) {file.remove(std, showWarnings = FALSE)}
  run_admb("RE", verbose = TRUE)
}

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