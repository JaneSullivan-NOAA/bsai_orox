# Survey and random effects model  %>% ass estimates
# Contact: jane.sullivan@noaa.gov
# Last updated: Nov 2020

# devtools::session_info()
# version  R version 4.0.2 (2020-06-22)
# os       Windows 10 x64              
# system   x86_64, mingw32             
# ui       RStudio 

# Set up ----

# Assessment year (most recent year with complete data set)
YEAR <- 2020

# date that catch was queried - important for long-term record keeping, tracking
# changes to db and tracing discrepancies in catch
access_date <- "10/13/2020" 

libs <- c("tidyverse", "R2admb", "viridis")
if(length(libs[which(libs %in% rownames(installed.packages()) == FALSE )]) > 0) {
  install.packages(libs[which(libs %in% rownames(installed.packages()) == FALSE)])}
lapply(libs, library, character.only = TRUE)

options(scipen = 999)

# User-defined fxns ----

source("R/admb_helper.R") # read in admb output
source("R/re_helper.R") # generalizes running random effects model

# Directory setup ----

root <- getwd()
dat_path <- paste0("data/", YEAR) # directory where source data is contained
out_path <- paste0("results/", YEAR) # directory for results/output
dir.create(out_path)
tpl_dir <- file.path("admb")

# If there's not an exe for the RE model, create one. If you update the tpl,
# you'll need to recompile it manually.

setwd(tpl_dir)

if(!file.exists("RE.exe")) {
  setup_admb()
  compile_admb("RE", verbose = TRUE)
}

if(!file.exists("re_test.exe")) {
  setup_admb()
  compile_admb("re_test", verbose = TRUE)
}

setwd(root)
files <- list.files(file.path(tpl_dir)) # searchable file list
tpl <- file.path(tpl_dir, files[which(grepl("RE.tpl", files))]) # find .tpl
exe <- file.path(tpl_dir, files[which(grepl("RE.exe", files))]) # find .exe
std <- paste0("RE.std")
# single survey version - this is the version on github...
stpl <- file.path(tpl_dir, files[which(grepl("RES.tpl", files))]) # find .tpl
sexe <- file.path(tpl_dir, files[which(grepl("RES.exe", files))]) # find .exe
# this is the single survey version cole emailed me 20210611
stpl <- file.path(tpl_dir, files[which(grepl("re_test.tpl", files))]) # find .tpl
sexe <- file.path(tpl_dir, files[which(grepl("re_test.exe", files))]) # find .exe

# Data ----

setwd(root)

# survey biomass: based on guidance from W Palsson 10/28/20 use EBS shelf
# (1982-present) and AI (1991-present) and the EBS slope which starts in 2002
full_biom <- read_csv(paste0(dat_path, "/bsai_orox_biomass_", YEAR, ".csv")) %>% 
  # year US trawl surveys were standardized (except for slope survey, which
  # began in 2002)
  filter(!(survey == "AI" & year < 1991))

# catch data for assessing overfishing status
catch <- read_csv(paste0(dat_path, "/bsai_orox_catch_2003_", YEAR, "_confidential.csv"))

# SST/non-SST biomass ----

glimpse(full_biom)

# FLAG: This assessment has historically summed variances across species to get
# the non-SST survey variance. PH and CT have informed me this is not correct,
# and we need to recalculate survey variance for all non-SST treated as a single
# stock.

# Stock structure for assessment: splits by SST and non-SST, also AI survey
# split into Southern Bering Sea and AI
data_sst_nonSST <- full_biom %>% 
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

# Run mod ----
out_sst_nonSST <- run_re_model(data = data_sst_nonSST, n_logsdlam = "single")

# Species glimpse ----

spp_eda <- full_biom %>% 
  # filter(common_name != "shortspine thornyhead") %>% 
  mutate(species = ifelse(species_code == 30020, "SST", "non-SST"),
         region = ifelse(survey == "AI" & area == "SBS", "SBS",
                         ifelse(survey == "AI" & grepl("AI", area), "AI",
                                survey)),
         fmp_subarea = ifelse(region == "AI", "AI", "EBS")) %>% 
  select(species, common_name, fmp_subarea, region, year, biomass) %>%
  mutate(species_name = case_when(common_name == "shortspine thornyhead" ~ "SST",
                                  common_name %in% c("dusky and dark rockfishes unid.",
                                                     "dusky rockfish") ~ "dusky",
                                  common_name %in% c("broadfin thornyhead", "longspine thornyhead") ~ "other thornyheads",
                                  common_name == "redbanded rockfish" ~ "redbanded",
                                  common_name == "sharpchin rockfish" ~ "sharpchin",
                                  common_name == "harlequin rockfish" ~ "harlequin",
                                  common_name == "yelloweye rockfish" ~ "yelloweye",
                                  common_name == "yellowmouth rockfish" ~ "yellowmouth",
                                  common_name == "rockfish unid." ~ "other"),
         species_name2 = ifelse(species_name %in% c("dusky", "SST", "other thornyheads"),
                                species_name, "other"))

spp_eda %>% 
  group_by(species, region, year) %>% 
  summarize(biomass = sum(biomass)) %>% 
  ungroup() %>% 
  complete(species, region, year, fill = list(biomass = 0)) %>% 
  ggplot(aes(x = year, y = biomass, fill = species)) +
  geom_bar(stat = "identity", alpha = 0.7, 
           # position = position_stack(reverse = TRUE),
           width = 1,
           colour = "grey") +
  scale_fill_grey(end = 0.7) +
  facet_wrap(~region, scales = "free_y") +
  theme_bw() +
  scale_y_continuous(labels = scales::comma) +
  # theme(legend.position = "bottom") +
  labs(x = NULL, y = "Biomass (t)", fill = "Complex")

ggsave(paste0(out_path, "/srvbiom_complexXregion_", YEAR, ".png"), 
       dpi=300, height=4, width=7, units="in")

nonsst_eda <- full_biom %>% 
  filter(common_name != "shortspine thornyhead") %>% 
  mutate(species = ifelse(species_code == 30020, "SST", "non-SST"),
         region = ifelse(survey == "AI" & area == "SBS", "SBS",
                         ifelse(survey == "AI" & grepl("AI", area), "AI",
                                survey))) %>% 
  select(species = common_name, region, year, biomass, var) %>%
  # combine dusky/dark and dusky
  mutate(species = ifelse(species == "dusky and dark rockfishes unid.",
                          "dusky rockfish", species)) %>% 
  # combine broadfin and longspine thornyhead
  mutate(species = ifelse(species %in% c("broadfin thornyhead", "longspine thornyhead"),
                          "other thornyhead", species)) %>% 
  group_by(species, region, year) %>% 
  summarize(biomass = sum(biomass),
            var = sum(var)) %>% 
  ungroup() %>% 
  complete(species, region, year, fill = list(biomass = 0, var = 0)) %>% 
  mutate(region = factor(region, levels = c("AI", "EBS_SHELF", "SBS", "EBS_SLOPE"), 
                         ordered = TRUE))

# library(viridis)
CUSTOMpalette <- c("#1a2ffa", "#0d177d", "#1a9ffa", "#fa751a", 
                   "#4b8e12", "#6fd21b", "#fae51a")#, "#c3b104", 
                   # "#f5df05", "#dcc805")
ggplot(data = nonsst_eda, aes(x = year, y = biomass, 
                              # fill = forcats::fct_rev(species))) +
                              fill = species)) +
  geom_bar(stat = "identity", alpha = 0.7, 
           position = position_stack(reverse = TRUE),
           size = 0.3,
           width = 1,
           colour = "black") +
  # scale_fill_viridis(discrete = TRUE, option = "D") +
  scale_fill_manual(values = CUSTOMpalette) +
  # scale_fill_grey(start = 1, end = 0) +
  facet_wrap(~region, scales = "free_y", ncol = 1) +
  # theme_minimal() +
  theme_bw() +
  scale_y_continuous(labels = scales::comma) +
  theme(legend.position = "bottom") +
  labs(x = NULL, y = "Biomass (t)", fill = "Non-SST species") +
  guides(fill = guide_legend(nrow = 3))

  
ggsave(paste0(out_path, "/srvbiom_nonSST_spp_", YEAR, ".png"), 
       dpi=400, height=8, width=6.5, units="in")

# nonsst_eda <- nonsst_eda %>% 
#   group_by(region, year) %>% 
#   mutate(total = sum(biomass),
#          prop = ifelse(total > 0, biomass / total, 0)) %>% 
#   ungroup()
# 
# nonsst2 <- nonsst_eda %>% 
#   filter(total > 0) %>% 
#   mutate(Year = factor(year),
#          Species = factor(species))

# ggplot(data = nonsst2 %>% 
#          filter(region %in% c("AI", "SBS")), aes(x = Year, y = prop, fill = Species)) +
#   geom_bar(stat = "identity", alpha = 0.6, width = 1, colour = "white") +
#   scale_fill_viridis(discrete = TRUE) +
#   facet_wrap(~region, scales = "free", ncol = 1) +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

# ggplot(data = nonsst2 %>% 
#          filter(region %in% c("AI", "SBS")),
#        aes(y = biomass, x = Year, fill = Species)) +
#   geom_bar(stat = "identity", position = position_stack(reverse = TRUE),
#            alpha = 0.6, width = 1, colour = "white") +
#   scale_fill_viridis(discrete = TRUE) +
#   facet_wrap(~region, scales = "free", ncol = 1) +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
#         legend.position = "top")

# Clean RE output ----

report <- c(paste0("Biomass estimates for BSAI Other Rockfish from AI/EBS trawl surveys and Random Effects Model", "\n",
                   "Contact: jane.sullivan@noaa.gov", "\n",
                   "Estimates generated ", Sys.Date()))

write.table(report, file = paste0(out_path, "/biomass_estimates_", YEAR, ".csv"), sep=",", quote = FALSE, row.names = FALSE, col.names = FALSE, eol = "\n")

# Survey biomass 95% confidence intervals (note that this is just replicating
# how RACE/GAP calculates them)
data_sst_nonSST <- data_sst_nonSST %>% 
  mutate(lci = biomass - 1.96 * sqrt(var),
         uci = biomass + 1.96 * sqrt(var),
         lci = ifelse(lci < 0, 0, lci),
         region = factor(region, levels = c("AI", "SBS", "EBS_SHELF", "EBS_SLOPE")))

f_srv_biomass <- data_sst_nonSST %>% 
  left_join(tibble(year = 1980:YEAR)) %>%
  complete(species, region, year, fill = list(biomass = NA, var = NA, cv = NA, lci = NA, uci = NA)) %>% 
  # Remove SST + EBS Shelf because they're never present there.
  filter(!c(species == "SST" & region == "EBS_SHELF"))

write.table(c(paste0("\n", "Survey biomass estimates (t) with estimated variance and 95% confidence intervals.
The reported CV was used as an input to random effects models. Biomass estimates of zero were treated as NAs in the model.
No SSTs occur on the EBS Shelf, so this species X region combination is not reported or estimated in subsequent RE models.")), 
            file = paste0(out_path, "/biomass_estimates_", YEAR, ".csv"), sep=",", quote = FALSE, row.names = FALSE, col.names = FALSE, eol = "\n", append = TRUE)
write.table(f_srv_biomass, file = paste0(out_path, "/biomass_estimates_", YEAR, ".csv"), sep=",", quote = FALSE, row.names = FALSE, col.names = TRUE, eol = "\n", append = TRUE)

# Alternative format for report:
alt <- f_srv_biomass %>% 
  select(species, region, year, biomass, cv) %>% 
  mutate(biomass = ifelse(is.na(biomass), "", 
                          prettyNum(formatC(biomass, format = "f", digits = 0), 
                                    big.mark = ",")),
         cv = ifelse(is.na(cv) | round(cv, 2) == 0, "", 
                     formatC(cv, format = "f", digits = 2)), 
         biomcv = ifelse((biomass == "" & cv == ""), "",
                         ifelse((biomass != "" & cv == ""),
                                biomass,
                                paste0(biomass, " (", cv, ")")))) %>% 
  select(-c(biomass, cv)) %>% 
  pivot_wider(id_cols = c(species, year), names_from = region, values_from = biomcv)

# write.table(c(paste0("\n", "ALTERNATIVE FORMAT: Survey biomass estimates (t) with CV used as inputs to random effects models. Biomass estimates of zero were treated as NAs in the model.
# No SSTs occur on the EBS Shelf, so this species X region combination is not reported or estimated in subsequent RE models.")), 
#             file = paste0(out_path, "/biomass_estimates_", YEAR, ".csv"), sep=",", quote = FALSE, row.names = FALSE, col.names = FALSE, eol = "\n", append = TRUE)
# write.table(alt, file = paste0(out_path, "/biomass_estimates_", YEAR, ".csv"), sep=",", quote = FALSE, row.names = FALSE, col.names = TRUE, eol = "\n", append = TRUE)

# write.csv(alt, paste0(out_path, "/survey_biomass_table.csv"))

# RE output
re_total <- out_sst_nonSST$re_output %>% filter(region == "TOTAL")
re_area <- out_sst_nonSST$re_output %>% filter(region != "TOTAL") %>% 
  mutate(region = factor(region, levels = c("AI", "SBS", "EBS_SHELF", "EBS_SLOPE")))

# results text
re_area %>% filter(year == YEAR)
re_total %>% filter(year == YEAR)

f_rebiom1 <- out_sst_nonSST$re_output %>% 
  arrange(species, region, year)

write.table(c(paste0("\n", "Biomass estimates (t) with 95% confidence intervals from the random effects model by species group,
region, and combined regions as region = TOTAL.")), 
            file = paste0(out_path, "/biomass_estimates_", YEAR, ".csv"), sep=",", quote = FALSE, row.names = FALSE, col.names = FALSE, eol = "\n", append = TRUE)
write.table(f_rebiom1, file = paste0(out_path, "/biomass_estimates_", YEAR, ".csv"), sep=",", quote = FALSE, row.names = FALSE, col.names = TRUE, eol = "\n", append = TRUE)
# write_csv(f_rebiom1, paste0(out_path, "/re_biomass_", YEAR, ".csv"))

# Calculate 95% CI by FMP (BS, AI), the level at which ABC/OFL are defined.
# Follow methods in tpl used to combine survey areas and calculate total
# biomass. TODO: ask PH if there is a reference for these formulas
fmp_biom <- out_sst_nonSST$area_sd %>% 
  left_join(re_area %>% select(species, year, region, re_est)) %>% 
  mutate(FMP = ifelse(region == "AI", "AI", "EBS"),
         sd_numer = exp(2 * biomsd + biomsd_sd ^ 2) * (exp(biomsd_sd ^ 2) - 1),
         sd_denom = exp(biomsd + 0.5 * biomsd_sd ^ 2)) %>% 
  group_by(species, year, FMP) %>% 
  summarize(fmp_est = sum(re_est), # mean = sum of survey regions, on natural scale
            fmp_sd = sqrt(log( sum(sd_numer) / (sum(sd_denom) ^ 2) + 1)), # sd by fmp on log scale
            fmp_lci = exp(log(fmp_est) - 1.96 * fmp_sd), # lci/uci transformed to natural scale
            fmp_uci = exp(log(fmp_est) + 1.96 * fmp_sd))

f_rebiom2 <- fmp_biom %>% 
  arrange(species, FMP, year)

write.table(c(paste0("\n", "Biomass estimates (t) with 95% confidence intervals from the random effects model by species group and FMP subarea.")), 
            file = paste0(out_path, "/biomass_estimates_", YEAR, ".csv"), sep=",", quote = FALSE, row.names = FALSE, col.names = FALSE, eol = "\n", append = TRUE)
write.table(f_rebiom2, file = paste0(out_path, "/biomass_estimates_", YEAR, ".csv"), sep=",", quote = FALSE, row.names = FALSE, col.names = TRUE, eol = "\n", append = TRUE)

# write_csv(f_rebiom2, paste0(out_path, "/re_biomass_fmp_", YEAR, ".csv"))

# Calculate 95% CI by species.
# Follow methods in tpl used to combine survey areas and calculate total
# biomass. TODO: ask PH if there is a reference for these formulas
spp_biom <- out_sst_nonSST$area_sd %>%
  left_join(re_area %>% select(species, year, region, re_est)) %>% 
  mutate(sd_numer = exp(2 * biomsd + biomsd_sd ^ 2) * (exp(biomsd_sd ^ 2) - 1),
         sd_denom = exp(biomsd + 0.5 * biomsd_sd ^ 2)) %>% 
  group_by(species, year) %>% 
  summarize(spp_est = sum(re_est), # mean = sum of survey regions, on natural scale
         spp_sd = sqrt(log( sum(sd_numer) / (sum(sd_denom) ^ 2) + 1)), # sd on log scale
         spp_lci = exp(log(spp_est) - 1.96 * spp_sd), # lci/uci transformed to natural scale
         spp_uci = exp(log(spp_est) + 1.96 * spp_sd)) %>% 
  ungroup()

# Same for all Other rockfish combined
all_orox_biom <- out_sst_nonSST$area_sd %>%
  left_join(re_area %>% select(species, year, region, re_est)) %>% 
  mutate(sd_numer = exp(2 * biomsd + biomsd_sd ^ 2) * (exp(biomsd_sd ^ 2) - 1),
         sd_denom = exp(biomsd + 0.5 * biomsd_sd ^ 2)) %>% 
  group_by(year) %>% 
  summarize(species = "Total Other Rockfish",
            spp_est = sum(re_est), # mean = sum of survey regions, on natural scale
            spp_sd = sqrt(log( sum(sd_numer) / (sum(sd_denom) ^ 2) + 1)), # sd on log scale
            spp_lci = exp(log(spp_est) - 1.96 * spp_sd), # lci/uci transformed to natural scale
            spp_uci = exp(log(spp_est) + 1.96 * spp_sd)) %>% 
  ungroup()

spp_biom <- bind_rows(spp_biom, all_orox_biom)
spp_biom %>% filter(year == max(year))

write.table(c(paste0("\n", "Biomass estimates (t) with 95% confidence intervals from the random effects
model by species group and all other rockfish combined.")), 
            file = paste0(out_path, "/biomass_estimates_", YEAR, ".csv"), sep=",", quote = FALSE, row.names = FALSE, col.names = FALSE, eol = "\n", append = TRUE)
write.table(spp_biom, file = paste0(out_path, "/biomass_estimates_", YEAR, ".csv"), sep=",", quote = FALSE, row.names = FALSE, col.names = TRUE, eol = "\n", append = TRUE)

# write_csv(spp_biom, paste0(out_path, "/re_biomass_bycomplex_", YEAR, ".csv"))

# ABC/OFL ----

# Tier 5: FOFL = M; FABC = 0.75 * M; OFL/ABC = FOFL/FABC * most recent biomass
# estimate

# SST and non-SST have different natural moralities. TODO: evaluate current M
# assumptions, create consistency with GOA Orox

natmat <- tibble(species = unique(re_total$species),
                 M = c(0.09, 0.03))

specs <- re_total %>%
  group_by(species) %>% 
  filter(year == YEAR) %>%  # most recent years biomass
  left_join(natmat) %>% 
  mutate(F_OFL = M,
         maxF_ABC = 0.75 * M, # maxF_ABC = F_ABC
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
  distinct() %>% 
  bind_rows(specs %>% 
              select(species, Area = FMP, ABC = FMP_ABC) %>% 
              pivot_longer(col = c("ABC")) %>% 
              mutate(name = paste0(Area, " ", name)) %>%
              select(-Area)) %>% 
  select(Species = species, Variable = name, Value = value) %>% 
  ungroup()

# Include summary for all OR - NOTE: Natural mortality and fishing mortality
# rates are specified separately for the SST and non-SST portions of the Other
# Rockfish complex.
specs_clean <- specs_clean %>% 
  bind_rows(specs_clean %>% 
  filter(Variable %in% c("Biomass", "OFL", "maxABC", "ABC", "AI ABC", "EBS ABC")) %>% 
  pivot_wider(id_cols = Species, names_from = "Variable", values_from = "Value", values_fill = 0) %>% 
  summarize_at(vars(-Species), list(sum)) %>% 
  mutate(Species = "All Other Rockfish",
         M = NA,
         F_OFL = NA,
         maxF_ABC = NA,
         F_ABC = NA) %>% 
  pivot_longer(cols = -Species, names_to = "Variable", values_to = "Value")) %>%
  mutate(Species = factor(Species, levels = c("SST", "non-SST", "All Other Rockfish"), ordered = TRUE),
         Variable = factor(Variable, levels = c("M", "Biomass", "F_OFL", "maxF_ABC",
                                                "F_ABC", "OFL", "maxABC", "ABC", "AI ABC", "EBS ABC"),
                           ordered = TRUE)) %>% 
  arrange(Species, Variable)

write_csv(specs_clean, paste0(out_path, "/specs_table_", YEAR, ".csv"))

specs_clean %>% filter(Variable == "Biomass")
re_total %>% filter(year == YEAR) # check, should match
specs_clean %>% filter(Variable == "OFL")
specs_clean %>% filter(Variable == "ABC")
specs_clean %>% filter(Variable == "AI ABC")
specs_clean %>% filter(Variable == "EBS ABC")
fmp_biom %>% 
  filter(year == YEAR) %>% 
  select(-c(fmp_sd, fmp_lci, fmp_uci), biom = fmp_est) %>% 
  group_by(species) %>% 
  mutate(biom_spp_sum = sum(biom)) %>% 
  group_by(FMP) %>% 
  mutate(biom_fmp_sum = sum(biom)) %>% 
  ungroup() %>% 
  mutate(prop_by_fmp = biom / biom_spp_sum)

# Theoretical FMP OFLs for total Other Rockfish
specs %>% 
  group_by(FMP) %>% 
  summarize(sum(FMP_OFL))
#Check that this adds up to total BSAI OFL
specs %>% 
  group_by(FMP) %>% 
  summarize(FMP_OFL = sum(FMP_OFL)) %>% 
  ungroup() %>% 
  summarize(OFL = sum(FMP_OFL))
# Biomass plots -----
  
# Region-specific biomass

# Flag for zero observations
biom <- data_sst_nonSST %>% mutate(zero_obs = ifelse(biomass == 0, TRUE, FALSE))

ggplot() +
  geom_point(data = biom, aes(x = year, y = biomass, col = zero_obs)) +
  geom_errorbar(data = biom, aes(x = year, ymin = lci, ymax = uci, col = zero_obs)) +
  geom_line(data = re_area %>% 
              filter(!(region == "EBS_SLOPE" & year < 2000)), 
            aes(x = year, y = re_est), col = "darkgrey") +
  geom_ribbon(data = re_area %>% 
                filter(!(region == "EBS_SLOPE" & year < 2000)), 
              aes(x = year, ymin = re_lci, ymax = re_uci),
              fill = "grey80", alpha = 0.4) +
  scale_colour_manual(values = c("black", "red")) +
  facet_grid(region ~ species, scales = "free_y") +
  # facet_grid(species ~ region, scales = "free_y") +
  scale_y_continuous(labels = scales::comma) +
  labs(x = NULL, y = "Biomass (t)") +
  theme_bw() +
  theme(legend.position = "none")

ggsave(paste0(out_path, "/biomass_areaXspp_", YEAR, ".png"), 
       dpi=300, height=6, width=8, units="in")

# ALT

# SST
tmp_biom <- biom %>% 
  mutate(region = factor(region, levels = c("EBS_SLOPE", "AI", "SBS"))) %>% 
  filter(species == "SST")

tmp_re <- re_area %>% 
  filter(species == "SST") %>% 
  filter(!(region == "EBS_SLOPE" & year < 2000) & 
           !(region %in% c("AI", "SBS") & year < 1991))
  # mutate(dat_flag = ifelse((region == "EBS_SLOPE" & year < 2002) &
  #                            (region == "AI" & year < 1991), 1, 0),
  #        re_lci = ifelse(dat_flag == 1, NA, re_lci),
  #        re_uci = ifelse(dat_flag == 1, NA, re_uci))
# filter(!(region == "EBS_SLOPE" & year < 2000))

# slope_fix <- data.frame(region = "EBS_SLOPE",
#                         year = min(full_biom$year):2002,
#                         yint = tmp_re %>% filter(dat_flag == 1) %>% distinct(re_est) %>% pull()) %>% 
#   mutate(region = factor(region, levels = c("EBS_SLOPE", "AI", "SBS")))

ggplot() +
  geom_point(data = tmp_biom, aes(x = year, y = biomass, col = zero_obs)) +
  geom_errorbar(data = tmp_biom, 
                aes(x = year, ymin = lci, ymax = uci, col = zero_obs),
                width = 0.5) +
  geom_line(data = tmp_re, # %>% filter(dat_flag == 0), 
            aes(x = year, y = re_est), col = "darkgrey") +
  # geom_line(data = slope_fix, aes(x = year, y = yint),
  #           col = "grey", lty = 2) +
  geom_ribbon(data = tmp_re, 
              aes(x = year, ymin = re_lci, ymax = re_uci),
              fill = "grey80", alpha = 0.4) +
  scale_colour_manual(values = c("black", "red")) +
  facet_wrap(~region, scales = "free_y", ncol = 1) +
  scale_y_continuous(labels = scales::comma) +
  labs(x = NULL, y = "Biomass (t)") +
  theme_bw(base_size = 14) +
  expand_limits(y = 0) +
  theme(legend.position = "none")

ggsave(paste0(out_path, "/biomass_SST_region_", YEAR, ".png"), 
       dpi=300, height=6, width=6, units="in")

# different SST plot for pres

ggplot() +
  geom_point(data = tmp_biom, 
             aes(x = year, y = biomass, col = region, fill = region)) +
  geom_errorbar(data = tmp_biom, 
                aes(x = year, ymin = lci, ymax = uci, 
                    col = region, fill = region),
                width = 0.5) +
  geom_line(data = tmp_re, 
            aes(x = year, y = re_est, col = region, fill = region)) +
  geom_ribbon(data = tmp_re, 
              aes(x = year, ymin = re_lci, ymax = re_uci, 
                  col = region, fill = region), 
              alpha = 0.3, col = "white") +
  scale_y_continuous(labels = scales::comma) +
  scale_fill_viridis(discrete = TRUE) +
  scale_colour_viridis(discrete = TRUE) +
  labs(x = NULL, y = NULL, col = NULL, fill = NULL,
       title = "SST biomass (t)") +
  theme_minimal(base_size = 14) +
  expand_limits(y = 0) +
  theme(legend.position = c(0.2, 0.8))

ggsave(paste0(out_path, "/pre_biomass_SST_region_", YEAR, ".png"), 
       dpi=300, height=4, width=4.5, units="in")


# non-SST
tmp_biom <- biom %>% 
  filter(species == "non-SST")

tmp_re <- re_area %>% 
  filter(species == "non-SST") %>% 
  filter(!(region == "EBS_SLOPE" & year < 2000) & 
           !(region %in% c("AI", "SBS") & year < 1991))
  # mutate(dat_flag = ifelse((region == "EBS_SLOPE" & year < 2002), 1, 0),
  #        re_lci = ifelse(dat_flag == 1, NA, re_lci),
  #        re_uci = ifelse(dat_flag == 1, NA, re_uci))
  # filter(!(region == "EBS_SLOPE" & year < 2000))

# slope_fix <- data.frame(region = "EBS_SLOPE",
#                         year = min(full_biom$year):2002,
#                         yint = tmp_re %>% filter(dat_flag == 1) %>% distinct(re_est) %>% pull()) %>% 
#   mutate(region = factor(region, levels = c("EBS_SLOPE", "AI", "SBS")))

ggplot() +
  geom_point(data = tmp_biom, aes(x = year, y = biomass, col = zero_obs)) +
  geom_errorbar(data = tmp_biom, 
                aes(x = year, ymin = lci, ymax = uci, col = zero_obs),
                width = 0.5) +
  geom_line(data = tmp_re, # %>% filter(dat_flag == 0),
            aes(x = year, y = re_est), 
            col = "darkgrey") +
  # geom_line(data = slope_fix, aes(x = year, y = yint),
  #           col = "darkgrey", lty = 2) +
  geom_ribbon(data = tmp_re, 
              aes(x = year, ymin = re_lci, ymax = re_uci),
              fill = "grey80", alpha = 0.4) +
  scale_colour_manual(values = c("black", "red")) +
  facet_wrap(~region, scales = "free_y", ncol = 2) +
  scale_y_continuous(labels = scales::comma) +
  labs(x = NULL, y = "Biomass (t)") +
  theme_bw(base_size = 14) +
  expand_limits(y = 0) +
  theme(legend.position = "none")

ggsave(paste0(out_path, "/biomass_nonSST_region_", YEAR, ".png"), 
       dpi=300, height=6, width=8, units="in")

# Biomass by FMP
ggplot(data = fmp_biom %>% 
         filter(year >= 2002) %>% 
         mutate(species = factor(species, levels = c("SST", "non-SST"), ordered = TRUE))) + #
  geom_line(aes(x = year, y = fmp_est), col = "black") +
  geom_ribbon(aes(x = year, ymin = fmp_lci, ymax = fmp_uci),
              fill = "grey80", alpha = 0.4) +
  facet_grid(species ~ FMP, scales = "free") +
  scale_y_continuous(labels = scales::comma) +
  labs(x = NULL, y = "Biomass (t)") +
  expand_limits(y = 0) +
  theme_bw()

ggsave(paste0(out_path, "/biomass_fmp_", YEAR, ".png"), 
       dpi=300, height=5, width=7, units="in")

# Total biomass
ggplot(data = re_total %>% 
         filter(year >= 2002) %>% 
         mutate(species = factor(species, levels = c("SST", "non-SST"), ordered = TRUE))) +
  geom_line(aes(x = year, y = re_est), col = "black") +
  geom_ribbon(aes(x = year, ymin = re_lci, ymax = re_uci),
              fill = "grey80", alpha = 0.4) +
  facet_wrap(~ species) +
  scale_y_continuous(labels = scales::comma) +
  labs(x = NULL, y = "Exploitable biomass (t)") +
  theme_bw(base_size = 14)

ggsave(paste0(out_path, "/biomass_total_", YEAR, ".png"), 
       dpi=300, height=3.5, width=5.5, units="in")

 Presentation total biomass ----
ggplot(data = re_total %>% 
         filter(year >= 2002) %>% 
         mutate(species = factor(species, levels = c("SST", "non-SST"), ordered = TRUE))) +
  geom_line(aes(x = year, y = re_est, col = species),
            size = 1) +
  geom_ribbon(aes(x = year, ymin = re_lci, ymax = re_uci, fill = species),
              alpha = 0.3, col = "white") +
  scale_color_manual(values = c("#6ba292", "#f0b74a")) +
  scale_fill_manual(values = c("#6ba292", "#f0b74a")) +
  scale_y_continuous(labels = scales::comma) +
  labs(x = NULL, y = NULL, title = "Exploitable biomass (t)",
       col = NULL, fill = NULL) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "top",
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

ggsave(paste0(out_path, "/pres_biom_", YEAR, ".png"), 
       dpi=300, height=3, width=4, units="in")

# Presentation catch ----
ggplot(data = catch %>%
         mutate(species = ifelse(species_name == "SST", "SST", "non-SST")) %>%
         group_by(year, species) %>%
         summarize(Catch = sum(tons)) %>% 
         mutate(species = factor(species, levels = c("SST", "non-SST"), ordered = TRUE))) +
  geom_area(aes(x = year, y = Catch, fill = species),
            alpha = 0.5 , size = 0.5, colour = "white") +
  scale_fill_manual(values = c("#6ba292", "#f0b74a")) +
  scale_y_continuous(labels = scales::comma) +
  labs(x = NULL, y = NULL, title = "Catch (t)", fill = NULL) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "top",
        plot.title = element_text(hjust = 0.5))

ggsave(paste0(out_path, "/pres_catch_", YEAR, ".png"), 
       dpi=300, height=3, width=4, units="in")

# Plot Catch/ABC/OFL  ----

# by FMP and species

specs %>%
  select(species, FMP, `ABC*` = FMP_ABC, `OFL*` = FMP_OFL) %>%
  left_join(catch %>%
  mutate(species = ifelse(species_name == "SST", "SST", "non-SST")) %>%
  group_by(year, FMP = fmp_subarea, species) %>%
  summarize(Catch = sum(tons))) %>% 
  pivot_longer(cols = c("ABC*", "OFL*", "Catch")) %>% 
  mutate(name = factor(name, levels = c("Catch", "ABC*", "OFL*"), ordered = TRUE),
         species = factor(species, levels = c("SST", "non-SST"), ordered = TRUE)) %>% #View()
  ggplot(aes(x = year, y = value, col = name, lty = name, size = name)) + #
  geom_line() +
  facet_grid(species ~ FMP, scales = "free") +
  scale_y_continuous(labels = scales::comma) +
  scale_colour_grey() +
  scale_size_manual(values = c(.5, 1, 1)) +
  labs(x = NULL, y = "Catch (t)", col = NULL, size = NULL, lty = NULL) +
  expand_limits(y = 0) +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "white"),
        panel.grid.minor = element_line(colour = "white"))

ggsave(paste0(out_path, "/catch_abc_ofl_fmp_", YEAR, ".png"), 
       dpi=300, height=5, width=7, units="in")

# Catch/ABC/OFL by species
specs %>%
  select(species, ABC, OFL) %>%
  left_join(catch %>%
              mutate(species = ifelse(species_name == "SST", "SST", "non-SST")) %>%
              group_by(year, species) %>%
              summarize(Catch = sum(tons))) %>% 
  pivot_longer(cols = c("ABC", "OFL", "Catch")) %>% 
  mutate(name = factor(name, levels = c("Catch", "ABC", "OFL"), ordered = TRUE)) %>% 
  ggplot(aes(x = year, y = value, col = name, lty = name, size = name)) + #
  geom_line() +
  facet_wrap(~ species, scales = "free", ncol = 1) +
  scale_y_continuous(labels = scales::comma) +
  scale_colour_grey() +
  scale_size_manual(values = c(.5, 1, 1)) +
  labs(x = NULL, y = "Catch (t)", col = NULL, size = NULL, lty = NULL) +
  expand_limits(y = 0) +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "white"),
        panel.grid.minor = element_line(colour = "white"))

ggsave(paste0(out_path, "/catch_abc_ofl_complex_", YEAR, ".png"), 
       dpi=300, height=5, width=4, units="in")

# Percent SST catch over time:

catch %>%
  mutate(species = ifelse(species_name == "SST", "SST", "non-SST")) %>%
  group_by(year, species) %>%
  summarize(Catch = sum(tons)) %>%
  group_by(year) %>% 
  mutate(Total_Catch = sum(Catch),
         percent = 100 * Catch / Total_Catch) %>% 
  filter(species == "SST" & year < 2020) %>% 
  ungroup() %>% 
  mutate(mean_percent = mean(percent),
         min_percent = min(percent),
         max_percent = max(percent))

re_total %>% 
  rename(re_biom = re_est) %>% 
  group_by(year) %>% 
  mutate(total_biom = sum(re_biom),
         percent = 100 * re_biom / total_biom) %>% 
  filter(species == "SST" & year < 2020 & year >= 2003) %>% 
  ungroup() %>% 
  mutate(mean_percent = mean(percent),
         min_percent = min(percent),
         max_percent = max(percent))

# Total Catch/ABC/OFL

specs %>%
  select(species, ABC, OFL) %>%
  ungroup() %>% 
  distinct() %>% 
  summarize(species = "All Other Rockfish",
            ABC = sum(ABC),
            OFL = sum(OFL)) %>% 
  left_join(catch %>%
              mutate(species = "All Other Rockfish") %>% 
              group_by(year, species) %>%
              summarize(Catch = sum(tons))) %>% 
  pivot_longer(cols = c("ABC", "OFL", "Catch")) %>% 
  mutate(name = factor(name, levels = c("Catch", "ABC", "OFL"), ordered = TRUE)) %>% 
  ggplot(aes(x = year, y = value, col = name, lty = name, size = name)) + #
  geom_line() +
  # facet_wrap(~ species, scales = "free", ncol = 1) +
  scale_y_continuous(labels = scales::comma) +
  scale_colour_grey() +
  scale_size_manual(values = c(.5, 1, 1)) +
  labs(x = NULL, y = NULL, title = "BSAI Other Rockfish",
       subtitle = "Comparison of catch to ABC/OFL (t)", 
       col = NULL, size = NULL, lty = NULL) +
  expand_limits(y = 0) +
  theme_bw(base_size = 14) +
  theme(panel.grid.major = element_line(colour = "white"),
        panel.grid.minor = element_line(colour = "white"),
        legend.position = c(0.8, 0.2))

ggsave(paste0(out_path, "/total_catch_abc_ofl_", YEAR, ".png"), 
       dpi=300, height=4, width=4.5, units="in")

# Catch to biomass ratio ----

refs <- data.frame(species = c("non-SST", "SST"),
                   yint = c(1, NA)) %>% 
  mutate(species = factor(species, 
                          levels = c("SST", "non-SST"), 
                          ordered = TRUE))
ratio <- specs %>%
  select(species, FMP, F_OFL) %>%
  left_join(catch %>%
              mutate(species = ifelse(species_name == "SST", "SST", "non-SST")) %>%
              group_by(year, FMP = fmp_subarea, species) %>%
              summarize(Catch = sum(tons))) %>% 
  left_join(fmp_biom %>% 
              select(species, year, FMP, biomass = fmp_est)) %>%
  mutate(ratio = Catch / biomass,
         species = factor(species, 
                          levels = c("SST", "non-SST"), 
                          ordered = TRUE)) 
ratio %>% 
  ggplot(aes(x = year, y = ratio, col = FMP, lty = FMP)) + #
  geom_line() +
  facet_wrap(~ species, scales = "free_y") +
  scale_colour_manual(values = c("grey", "black")) +
  scale_linetype_manual(values = c(1, 5)) +
  geom_hline(data = refs, aes(yintercept = yint),
             col = "red", lty = 3) +
  labs(x = NULL, y = "Catch/Biomass Ratio", col = NULL, size = NULL, lty = NULL) +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "white"),
        panel.grid.minor = element_line(colour = "white"))

ggsave(paste0(out_path, "/ratio_catch_biomass_", YEAR, ".png"), 
       dpi=300, height=3.75, width=6, units="in")

ratio %>% 
  group_by(FMP, species) %>% 
  summarize(mean_ratio = mean(ratio))

ratio %>% filter(year == YEAR & species == "non-SST" & FMP == "AI")

# Presentation catch/biomass----

# refs <- data.frame(FMP = c("AI", "EBS"),
#                    yint = c(1, NA)) %>% 
#   mutate(FMP = factor(FMP, 
#                       labels = c("Aleutian Islands", "Eastern Bering Sea"),
#                       levels = c("AI", "EBS"), 
#                       ordered = TRUE))
ratio %>% 
  mutate(FMP = factor(FMP, 
                      labels = c("Aleutian Islands", "Eastern Bering Sea"),
                      levels = c("AI", "EBS"), 
                      ordered = TRUE),
         species = factor(species, 
                          levels = c("SST", "non-SST"), 
                          ordered = TRUE)) %>% 
  # group_by(year, species) %>% 
  # summarize(ratio = sum(Catch) / sum(biomass)) %>% 
  ggplot(aes(x = year, y = ratio, col = species, lty = FMP)) +
  geom_line(size = 1) +
  labs(x = NULL, y = NULL, title = "Exploitation rate (catch/biomass)",
       # subtitle = "Note differences in y-axis scale",
       col = NULL, fill = NULL, lty = NULL) +
  geom_hline(data = refs, aes(yintercept = yint),
             col = "red", lty = 3) +
  # facet_wrap(~species, scales = "free_y") +
  # geom_hline(y = 1)
  scale_color_manual(values = c("#6ba292", "#f0b74a")) + #, guide = 'none') +
  theme_minimal(base_size = 13) +
  theme(#legend.position = "top",
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)#,
        # panel.grid.major = element_line(colour = "white"),
        # panel.grid.minor = element_line(colour = "white")
        )

ggsave(paste0(out_path, "/pres_ratio_", YEAR, ".png"), 
       dpi=300, height=4, width=7, units="in")


# SARA ----

# http://dev-afsc.nmfs.local/REFM/Stocks/Plan_Team/SARA/SARAdata.php

# Re-run as full complex (not split into SST/non-SST)
sara_dat <- full_biom %>% 
  mutate(species = "Other_Rockfish_SARA",
         region = survey) %>% 
  group_by(species, region, year) %>% 
  summarize(biomass = sum(biomass),
            var = sum(var)) %>% 
  ungroup() %>% 
  # Currently years with 0 biomass are removed from RE model
  mutate(cv = ifelse(biomass == 0, 0, sqrt(var) / biomass)) %>% 
  arrange(species, region, year) 

out_sara <- run_re_model(data = sara_dat)
sara_biomass <- out_sara$re_output %>% 
  filter(region == "TOTAL")
sara_biomass %>% filter(year == YEAR)
sara_biomass %>% pull(year)
sara_biomass %>% mutate(re_est = round(re_est, 0)) %>% pull(re_est)

catch %>%
  group_by(year) %>%
  summarize(catch = round(sum(tons), 0)) %>% 
  pull(catch)

# Survey biomass data
aisrv <- sara_dat %>% filter(region == "AI")
aisrv %>% pull(year)
round(aisrv %>% pull(biomass))

ebsshelfsrv <- sara_dat %>% filter(region == "EBS_SHELF")
ebsshelfsrv %>% pull(year)
round(ebsshelfsrv %>% pull(biomass))

ebsslopesrv <- sara_dat %>% filter(region == "EBS_SLOPE")
ebsslopesrv %>% pull(year)
round(ebsslopesrv %>% pull(biomass))

# logsdlam ----

# multi survey tpl with single logsdlam
msrv_slam <- run_re_model(data = data_sst_nonSST, n_logsdlam = "single")
names(msrv_slam)


# multi survey tpl with a logsdlam for each survey
msrv_mlam <- run_re_model(data = data_sst_nonSST, n_logsdlam = "multiple")

# single survey tpl
sout <- run_sre_model(data = data_sst_nonSST)

logsdlam <- msrv_slam$logsdlam %>% 
  left_join(msrv_slam$diagnostics) %>% 
  mutate(tpl = "multi-survey tpl",
         n_logsdlam = "single logsdlam for all surveys") %>% 
  bind_rows(msrv_mlam$logsdlam %>% 
    left_join(msrv_mlam$diagnostics) %>% 
    mutate(tpl = "multi-survey tpl",
           n_logsdlam = "logsdlam for each survey")) %>% 
  bind_rows(sout$logsdlam %>% 
    left_join(sout$diagnostics) %>% 
    mutate(tpl = "single survey tpl",
           n_logsdlam = "logsdlam for each survey"))

logsdlam <- logsdlam %>% 
  mutate(lci = logsdlam - 1.96 * logsdlam_se * logsdlam,
         uci = logsdlam + 1.96 * logsdlam_se * logsdlam,
         region = ifelse(!is.na(region), region,
                         ifelse(is.na(region) & !is.na(area), area,
                         "COMBINED"))) %>% 
  mutate(region = factor(region,
                         levels = c("AI", "SBS", "EBS_SLOPE", "EBS_SHELF", "COMBINED"),
                         ordered = TRUE)) %>% 
  select(-area)

logsdlam %>% 
  arrange(species, region, n_logsdlam) %>% 
  select(species, region, tpl, n_logsdlam, logsdlam, sd = logsdlam_se)

ggplot(logsdlam %>% filter(!c(species == "non-SST" & 
                            n_logsdlam == "logsdlam for each survey" &
                            region == "EBS_SLOPE")), 
       aes(x = region, y = logsdlam, 
           col = tpl,
           ymin = lci, ymax = uci)) +
  geom_point(position=position_dodge(width=0.2)) +
  geom_errorbar(width = 0.3, position=position_dodge(width=0.2)) +
  facet_grid(species ~ n_logsdlam) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  labs(x = NULL)
  
ggsave(paste0(out_path, "/logsdlam_by_tpl_bsai_orox.png"),
       dpi = 400, height = 4, width = 10, units = "in")

# other ----

# other ways the complex has been split up...

# full_biom %>% 
#   mutate(full_complex = "All-OR",
#        complex_bySST = ifelse(species_code == 30020, "SST", "non-SST"),
#        complex_bySSTandDusky = ifelse(species_code == 30020, "SST",
#                                       ifelse(species_code == 30152, "Dusky", "non-SST-or-Dusky"))) %>% 
#   pivot_longer(cols = c("full_complex", "complex_bySST", "complex_bySSTandDusky"), 
#                names_to = "complex", values_to = "species") 