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

options(scipen = 999)

# ADMB helper fxns ----
source("R/admb_helper.R")

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

# Data ----

full_biom <- read_csv(paste0(dat_path, "/bsai_orox_biomass_", YEAR, ".csv"))

# Biomass summary ----

# Stock structure for assessment: splits by SST and non-SST, also AI survey
# split into Southern Bering Sea and AI
biom <- full_biom %>% 
  mutate(species = ifelse(species_code == 30020, "SST", "nonSST"),
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

# expand to include all years, species, regions for input to RE model
input_biom <- biom %>% 
  expand(year = min(biom$year):YEAR, species, region) %>% 
  left_join(biom)

input_biom <- input_biom %>% 
  pivot_wider(id_cols = c("species", "year"), names_from = "region", 
              values_from = "biomass", values_fill = 0) %>% 
  mutate(value = "mu") %>% 
  bind_rows(input_biom %>% 
              pivot_wider(id_cols = c("species", "year"), names_from = "region", 
                          values_from = "cv", values_fill = 0) %>% 
              mutate(value = "cv")) %>% 
  arrange(species, value, year) %>% 
  replace(is.na(.), 0)

# Create a look up table of species x area combinations to include in the RE
# model runs (remove areas where species are not sampled, e.g. SST in EBS_SHELF)
model_lkup <- input_biom %>% 
  filter(value == "mu") %>% 
  pivot_longer(cols = unique(biom$region), names_to = "area", values_to = "biomass") %>% 
  group_by(species, area) %>% 
  summarize(tot_biom = sum(biomass)) %>% 
  ungroup() %>% 
  filter(tot_biom > 0) %>% 
  select(species, area)

# Executable ----

# If there's not an exe for the RE model, create one. If you update the tpl,
# you'll need to recompile it manually.

setwd(tpl_dir)

if(!file.exists("RE.exe")) {
  setup_admb()
  compile_admb("RE", verbose = TRUE)
}

# Run model ----

spp <- unique(model_lkup$species)
yrs <- unique(input_biom$year)
process_error <- 1 # TODO: talk with PH about this assumption and alternative options

diagnostics <- list() # for model diagnostics (convergence, mgc)
re_output <- list() # random effect biomass estimate output
area_sd <- list() # sd of the estimated biomass sd generated, property of sdreport obj

for(i in 1:length(spp)) {
  
  name <- paste0("RE_", spp[i]) # will be tpl name
  areas <- model_lkup %>% filter(species == spp[i]) %>% pull(area) # areas in model
  error <- rep(process_error, length(areas))
  
  # Create new directory and copy over tpl and exe files. Rename them to
  # specific species/complex
  setwd(root)
  spp_out <- file.path(paste0(out_path, "/RE_", spp[i]))
  dir.create(spp_out)
  file.copy(from = tpl, to = spp_out, overwrite = TRUE)
  file.copy(from = exe, to = spp_out, overwrite = TRUE)
  setwd(spp_out)
  file.rename("RE.tpl", paste0(name, ".tpl"))
  file.rename("RE.exe", paste0(name, ".exe"))
  
  # Biomass estimates
  srv_est <- input_biom %>% 
    filter(species == spp[i] & value == "mu") %>% 
    select(-c(species, year, value)) %>% 
    select(areas) %>% 
    as.matrix()
  srv_est <- unname(srv_est) # gets rid of column names
  srv_est[srv_est == 0] <- "-9" # ADMB flag for NA or zero data
  
  # Biomass variability
  srv_cv <- input_biom %>% 
    filter(species == spp[i] & value == "cv") %>% 
    select(-c(species, year, value)) %>% 
    select(areas) %>% 
    as.matrix()
  srv_cv <- unname(srv_cv) # gets rid of column names
  srv_cv[srv_cv == 0] <- "-9" # ADMB flag for NA or zero data
  
  # Write dat file
  cat("# Model start and end years","\n", min(yrs), max(yrs), "\n",
      "# Number of survey indices fit (i.e., regions/depth strata)","\n", length(areas), "\n",
      "# Number or process error parameters","\n",1,"\n",
      "# Process error index","\n", error, "\n",
      "# Number of surveys","\n", length(yrs), "\n",
      "# Survey years","\n", yrs, "\n",
      "#", areas, "\n", 
      "# Survey biomass","\n", sep = " ", file = paste0(name, ".dat"))
  
  write.table(srv_est, file = paste0(name, ".dat"), sep = " ", append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE)
  write.table(paste0("# Survey biomass SE"), file = paste0(name, ".dat"), sep = " ", append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE)
  write.table(srv_cv, file = paste0(name, ".dat"), sep = " ", append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE)

  # remove pre-existing std file if it exists and run admb
  if (file.exists(paste0(name, ".std"))) {file.remove(paste0(name, ".std"), showWarnings = FALSE)}
  run_admb(name, verbose = TRUE)
  file.rename("rwout.rep", paste0(name, "_rwout.rep"))
  
  # full admb output
  df <- read_admb(tplfile = name, repfile = paste0(name, "_rwout"))
  
  # dataframe with species complex and year
  id <- tibble(species = spp[i], 
               year = df$yrs)
  
  # total biomass + 95% ci
  re_output[[i]] <- id %>% 
    mutate(region = "TOTAL",
           re_est = df$biom_TOT,
           re_lci = df$biom_TOT_LCI,
           re_uci = df$biom_TOT_UCI) %>% 
    # area-specific biomass estimates + 95% ci
    bind_rows(id %>% 
                bind_cols(setNames(as_tibble(df$biomA), areas)) %>% # colnames = areas where survey occurred
                pivot_longer(cols = areas, names_to = "region", values_to = "re_est") %>% 
                left_join(id %>% 
                            bind_cols(setNames(as_tibble(df$LCI), areas)) %>%  
                            pivot_longer(cols = areas, names_to = "region", values_to = "re_lci")) %>% 
                left_join(id %>% 
                            bind_cols(setNames(as_tibble(df$UCI), areas)) %>%  
                            pivot_longer(cols = areas, names_to = "region", values_to = "re_uci")))
  
  # Area-specific sd of biomass sd - used to calculate BS or AI wide RE estimates
  area_sd[[i]] <- id %>% 
    bind_cols(setNames(as_tibble(df$biomsd), areas)) %>% 
    pivot_longer(cols = areas, names_to = "region", values_to = "biomsd") %>% 
    left_join(id %>% 
                bind_cols(setNames(as_tibble(df$biomsd.sd), areas)) %>% 
                pivot_longer(cols = areas, names_to = "region", values_to = "biomsd_sd"))
    
  # diagnostics - just mgc and test for positive definite hessian
  diagnostics[[i]] <- tibble(species = spp[i],
                             logDetHess = df$fit$logDetHess,
                             mgc = df$fit$maxgrad) %>% 
    mutate(converged = mgc < 0.001 & logDetHess > 0)
}

re_output <- do.call(rbind, re_output)
area_sd <- do.call(rbind, area_sd)
(diagnostics <- do.call(rbind, diagnostics))

# Clean RE output ----

setwd(root)

# Survey biomass 95% confidence intervals (note that this is just replicating
# how RACE/GAP calculates them)
biom <- biom %>% 
  mutate(lci = biomass - 1.96 * sqrt(var),
         uci = biomass + 1.96 * sqrt(var),
         lci = ifelse(lci < 0, 0, lci),
         region = factor(region, levels = c("AI", "SBS", "EBS_SHELF", "EBS_SLOPE")))

re_total <- re_output %>% filter(region == "TOTAL")
re_area <- re_output %>% filter(region != "TOTAL") %>% 
  mutate(region = factor(region, levels = c("AI", "SBS", "EBS_SHELF", "EBS_SLOPE")))

write_csv(re_output, paste0(out_path, "/re_biomass_", YEAR, ".csv"))

# Calculate 95% CI by FMP (BS, AI), the level at which ABC/OFL are defined.
# Follow methods in tpl used to combine survey areas and calculate total
# biomass. TODO: ask PH if there is a reference for these formulas
fmp_biom <- area_sd %>% 
  left_join(re_area %>% select(species, year, region, re_est)) %>% 
  mutate(FMP = ifelse(region == "AI", "AI", "BS"),
         sd_numer = exp(2 * biomsd + biomsd_sd ^ 2) * (exp(biomsd_sd ^ 2) - 1),
         sd_denom = exp(biomsd + 0.5 * biomsd_sd ^ 2)) %>% 
  group_by(species, year, FMP) %>% 
  summarize(fmp_est = sum(re_est), # mean = sum of survey regions, on natural scale
            fmp_sd = sqrt(log( sum(sd_numer) / (sum(sd_denom) ^ 2) + 1)), # sd by fmp on log scale
            fmp_lci = exp(log(fmp_est) - 1.96 * fmp_sd), # lci/uci transformed to natural scale
            fmp_uci = exp(log(fmp_est) + 1.96 * fmp_sd))

write_csv(fmp_biom, paste0(out_path, "/re_biomass_fmp_", YEAR, ".csv"))

# ABC/OFL ----

# Tier 5: FOFL = M; FABC = 0.75 * M; OFL/ABC = FOFL/FABC * most recent biomass
# estimate

# SST and nonSST have different natural moralities. TODO: evaluate current M
# assumptions, create consistency with GOA Orox
SST_M <- 0.03
nonSST_M <- 0.09

specs <- re_total %>%
  group_by(species) %>% 
  filter(year == YEAR) %>%  # most recent years biomass
  mutate(M = ifelse(species == "SST", SST_M, nonSST_M),
         F_OFL = M,
         maxF_ABC = 0.75 * M, # currently maxF_ABC = F_ABC
         F_ABC = 0.75 * M,
         OFL = re_est * F_OFL,
         maxABC = re_est * maxF_ABC,
         ABC = re_est * F_ABC) %>% 
  # join in FMP-level biomass estimates. The fraction in each region determines
  # the AI and BS specific ABC and OFL
  left_join(fmp_biom %>% filter(year == YEAR)) %>% 
  mutate(fraction = fmp_est / re_est,
         FMP_ABC = fraction * ABC,
         FMP_OFL = fraction * OFL) %>% 
  select(species, FMP, M, Biomass = re_est, F_OFL, 
         maxF_ABC, F_ABC, OFL, maxABC, ABC, FMP_ABC, FMP_OFL) 

specs_clean <- specs %>% 
  select(-c(FMP, FMP_ABC, FMP_OFL)) %>% 
  pivot_longer(cols = -species) %>% 
  bind_rows(specs %>% 
              select(species, Area = FMP, ABC = FMP_ABC, OFL = FMP_OFL) %>% 
              pivot_longer(col = c("ABC", "OFL")) %>% 
              mutate(name = paste0(Area, " ", name)) %>%
              select(-Area)) %>% 
  select(Species = species, Variable = name, Value = value) %>% 
  arrange(Species)

write_csv(specs_clean, paste0(out_path, "/specs_table_", YEAR, ".csv"))

# Plot -----

# Region-specific biomass
ggplot() +
  geom_point(data = biom, aes(x = year, y = biomass)) +
  geom_errorbar(data = biom, aes(x = year, ymin = lci, ymax = uci)) +
  geom_line(data = re_area, aes(x = year, y = re_est), col = "darkgrey") +
  geom_ribbon(data = re_area, aes(x = year, ymin = re_lci, ymax = re_uci),
              fill = "grey80", alpha = 0.4) +
  facet_grid(region ~ species, scales = "free_y") +
  scale_y_continuous(labels = scales::comma) +
  labs(x = NULL, y = "Biomass (t)") +
  theme_bw()

ggsave(paste0(out_path, "/biomass_areaXspp_", YEAR, ".png"), 
       dpi=300, height=6, width=8, units="in")

# Biomass by FMP
ggplot(data = fmp_biom %>% filter(year >= 2000)) + #
  geom_line(aes(x = year, y = fmp_est), col = "black") +
  geom_ribbon(aes(x = year, ymin = fmp_lci, ymax = fmp_uci),
              fill = "grey80", alpha = 0.4) +
  # FMP ABC
  # geom_hline(data = specs %>% 
  #              select(species, FMP, ABC = FMP_ABC, OFL = FMP_OFL) %>% 
  #              pivot_longer(cols = c("ABC", "OFL"), names_to = "Area-specific values"), 
  #            aes(yintercept = value, col = `Area-specific values`, lty = `Area-specific values`)) + 
  facet_grid(species ~ FMP, scales = "free") +
  scale_y_continuous(labels = scales::comma) +
  labs(x = NULL, y = "Biomass (t)") +
  expand_limits(y = 0) +
  theme_bw()

ggsave(paste0(out_path, "/biomass_fmp_", YEAR, ".png"), 
       dpi=300, height=5, width=7, units="in")


# Total biomass
ggplot(data = re_total %>% filter(year >= 2000)) +
  geom_line(aes(x = year, y = re_est), col = "black") +
  geom_ribbon(aes(x = year, ymin = re_lci, ymax = re_uci),
              fill = "grey80", alpha = 0.4) +
  facet_wrap(~ species) +
  scale_y_continuous(labels = scales::comma) +
  labs(x = NULL, y = "Biomass (t)") +
  theme_bw()

ggsave(paste0(out_path, "/biomass_total_", YEAR, ".png"), 
       dpi=300, height=4, width=7, units="in")

# other ----

# other ways the complex has been split up...

# full_biom %>% 
#   mutate(full_complex = "All-OR",
#        complex_bySST = ifelse(species_code == 30020, "SST", "non-SST"),
#        complex_bySSTandDusky = ifelse(species_code == 30020, "SST",
#                                       ifelse(species_code == 30152, "Dusky", "non-SST-or-Dusky"))) %>% 
#   pivot_longer(cols = c("full_complex", "complex_bySST", "complex_bySSTandDusky"), 
#                names_to = "complex", values_to = "species") 