# Survey and random effects model biomass estimates
# Contact: jane.sullivan@noaa.gov
# Last updated: Oct 2020

# devtools::session_info()
# version  R version 4.0.2 (2020-06-22)
# os       Windows 10 x64              
# system   x86_64, mingw32             
# ui       RStudio 

# Set up ----

# Assessment year (most recent year with complete data set)
YEAR <- 2019

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
files <- list.files(file.path(tpl_dir)) # searchable file list
tpl <- file.path(tpl_dir, files[which(grepl(".tpl", files))]) # find .tpl
exe <- file.path(tpl_dir, files[which(grepl(".exe", files))]) # find .exe
std <- paste0("RE.std")

# Executable ----

# If there's not an exe for the RE model, create one. If you update the tpl,
# you'll need to recompile it manually.

setwd(tpl_dir)

if(!file.exists("RE.exe")) {
  setup_admb()
  compile_admb("RE", verbose = TRUE)
}

# Data ----

setwd(root)

# survey biomass
full_biom <- read_csv(paste0(dat_path, "/bsai_orox_biomass_", YEAR, ".csv"))

# catch data for assessing overfishing status
catch <- read_csv(paste0(dat_path, "/bsai_orox_catch_2003_", YEAR+1, "_confidential.csv"))

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

out_sst_nonSST <- run_re_model(data = data_sst_nonSST)

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

library(viridis)
ggplot(data = nonsst_eda, aes(x = year, y = biomass, 
                              # fill = forcats::fct_rev(species))) +
                              fill = species)) +
  geom_bar(stat = "identity", alpha = 0.65, 
           position = position_stack(reverse = TRUE),
           # size = 0.1, 
           width = 1,
           colour = "grey") +
  scale_fill_viridis(discrete = TRUE, option = "A") +
  # scale_fill_grey(start = 1, end = 0) +
  facet_wrap(~region, scales = "free_y") +
  # theme_minimal() +
  theme_bw() +
  scale_y_continuous(labels = scales::comma) +
  theme(legend.position = "bottom") +
  labs(x = NULL, y = "Biomass (t)", fill = "Non-SST species")
  
ggsave(paste0(out_path, "/srvbiom_nonSST_spp_", YEAR, ".png"), 
       dpi=300, height=4, width=7, units="in")

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
re_area %>% filter(year == 2020)

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

specs_clean %>% filter(Variable == "OFL")
specs_clean %>% filter(Variable == "ABC")
specs_clean %>% filter(Variable == "AI ABC")
specs_clean %>% filter(Variable == "EBS ABC")

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
  mutate(dat_flag = ifelse((region == "EBS_SLOPE" & year < 2002), 1, 0),
         re_lci = ifelse(dat_flag == 1, NA, re_lci),
         re_uci = ifelse(dat_flag == 1, NA, re_uci))
# filter(!(region == "EBS_SLOPE" & year < 2000))

slope_fix <- data.frame(region = "EBS_SLOPE",
                        year = 1980:2002,
                        yint = tmp_re %>% filter(dat_flag == 1) %>% distinct(re_est) %>% pull()) %>% 
  mutate(region = factor(region, levels = c("EBS_SLOPE", "AI", "SBS")))


ggplot() +
  geom_point(data = tmp_biom, aes(x = year, y = biomass, col = zero_obs)) +
  geom_errorbar(data = tmp_biom, 
                aes(x = year, ymin = lci, ymax = uci, col = zero_obs),
                width = 0.5) +
  geom_line(data = tmp_re %>% filter(dat_flag == 0), 
            aes(x = year, y = re_est), col = "darkgrey") +
  geom_line(data = slope_fix, aes(x = year, y = yint),
            col = "grey", lty = 2) +
  geom_ribbon(data = tmp_re, 
              aes(x = year, ymin = re_lci, ymax = re_uci),
              fill = "grey80", alpha = 0.4) +
  scale_colour_manual(values = c("black", "red")) +
  facet_wrap(~region, scales = "free_y", ncol = 1) +
  scale_y_continuous(labels = scales::comma) +
  labs(x = NULL, y = "Biomass (t)") +
  theme_bw() +
  expand_limits(y = 0) +
  theme(legend.position = "none")

ggsave(paste0(out_path, "/biomass_SST_region_", YEAR, ".png"), 
       dpi=300, height=6, width=6, units="in")

# non-SST
tmp_biom <- biom %>% 
  filter(species == "non-SST")

tmp_re <- re_area %>% 
  filter(species == "non-SST") %>% 
  mutate(dat_flag = ifelse((region == "EBS_SLOPE" & year < 2002), 1, 0),
         re_lci = ifelse(dat_flag == 1, NA, re_lci),
         re_uci = ifelse(dat_flag == 1, NA, re_uci))
  # filter(!(region == "EBS_SLOPE" & year < 2000))

slope_fix <- data.frame(region = "EBS_SLOPE",
                        year = 1980:2002,
                        yint = tmp_re %>% filter(dat_flag == 1) %>% distinct(re_est) %>% pull()) %>% 
  mutate(region = factor(region, levels = c("EBS_SLOPE", "AI", "SBS")))

ggplot() +
  geom_point(data = tmp_biom, aes(x = year, y = biomass, col = zero_obs)) +
  geom_errorbar(data = tmp_biom, 
                aes(x = year, ymin = lci, ymax = uci, col = zero_obs),
                width = 0.5) +
  geom_line(data = tmp_re %>% filter(dat_flag == 0),
            aes(x = year, y = re_est), 
            col = "darkgrey") +
  geom_line(data = slope_fix, aes(x = year, y = yint),
            col = "darkgrey", lty = 2) +
  geom_ribbon(data = tmp_re, 
              aes(x = year, ymin = re_lci, ymax = re_uci),
              fill = "grey80", alpha = 0.4) +
  scale_colour_manual(values = c("black", "red")) +
  facet_wrap(~region, scales = "free_y", ncol = 1) +
  scale_y_continuous(labels = scales::comma) +
  labs(x = NULL, y = "Biomass (t)") +
  theme_bw() +
  expand_limits(y = 0) +
  theme(legend.position = "none")

ggsave(paste0(out_path, "/biomass_nonSST_region_", YEAR, ".png"), 
       dpi=300, height=7, width=6, units="in")

# Biomass by FMP
ggplot(data = fmp_biom %>% 
         filter(year >= 2000) %>% 
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
         species = factor(species, levels = c("SST", "non-SST"), ordered = TRUE)) %>% View()
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
  facet_wrap(~ species, scales = "free", ncol = 1) +
  scale_y_continuous(labels = scales::comma) +
  scale_colour_grey() +
  scale_size_manual(values = c(.5, 1, 1)) +
  labs(x = NULL, y = "Catch (t)", col = NULL, size = NULL, lty = NULL) +
  expand_limits(y = 0) +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "white"),
        panel.grid.minor = element_line(colour = "white"))

ggsave(paste0(out_path, "/total_catch_abc_ofl_", YEAR, ".png"), 
       dpi=300, height=4, width=5, units="in")

# Catch to biomass ratio ----

refs <- data.frame(species = c("non-SST", "SST"),
                   yint = c(1, NA)) %>% 
  mutate(species = factor(species, 
                          levels = c("SST", "non-SST"), 
                          ordered = TRUE))
specs %>%
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
                          ordered = TRUE)) %>% 
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

# other ----

# other ways the complex has been split up...

# full_biom %>% 
#   mutate(full_complex = "All-OR",
#        complex_bySST = ifelse(species_code == 30020, "SST", "non-SST"),
#        complex_bySSTandDusky = ifelse(species_code == 30020, "SST",
#                                       ifelse(species_code == 30152, "Dusky", "non-SST-or-Dusky"))) %>% 
#   pivot_longer(cols = c("full_complex", "complex_bySST", "complex_bySSTandDusky"), 
#                names_to = "complex", values_to = "species") 