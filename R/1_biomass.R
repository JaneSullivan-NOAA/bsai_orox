# Survey and random effects model biomass estimates
# Contact: jane.sullivan@noaa.gov
# Last updated: Sep 2020
# devtools::session_info()
# version  R version 4.0.2 (2020-06-22)
# os       Windows 10 x64              
# system   x86_64, mingw32             
# ui       RStudio 

# Set up ----

# Assessment year (most recent year with complete data set)
YEAR <- 2019

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

# catch data
catch <- read_csv(paste0(dat_path, "/bsai_orox_catch_2003_", YEAR+1, "_confidential.csv"))

# fishery catch is at the genus level for thornyheads, so use observer data to
# estimate proportion SST from other thornyheads
thorny <- read_csv(paste0(dat_path, "/obs_thornyheads_2003_", YEAR+1, ".csv"))

# ak land and nmfs area shapefiles (these are already cleaned up and saved as a
# dataframe)
load("data/map_ak_land.RData")
load("data/map_nmfs_areas.RData")

# SST/nonSST biomass ----

# Stock structure for assessment: splits by SST and non-SST, also AI survey
# split into Southern Bering Sea and AI
data_sst_nonsst <- full_biom %>% 
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

out_sst_nonsst <- run_re_model(data = data_sst_nonsst)

# Clean RE output ----

# Survey biomass 95% confidence intervals (note that this is just replicating
# how RACE/GAP calculates them)
data_sst_nonsst <- data_sst_nonsst %>% 
  mutate(lci = biomass - 1.96 * sqrt(var),
         uci = biomass + 1.96 * sqrt(var),
         lci = ifelse(lci < 0, 0, lci),
         region = factor(region, levels = c("AI", "SBS", "EBS_SHELF", "EBS_SLOPE")))

re_total <- out_sst_nonsst$re_output %>% filter(region == "TOTAL")
re_area <- out_sst_nonsst$re_output %>% filter(region != "TOTAL") %>% 
  mutate(region = factor(region, levels = c("AI", "SBS", "EBS_SHELF", "EBS_SLOPE")))

write_csv(out_sst_nonsst$re_output, paste0(out_path, "/re_biomass_", YEAR, ".csv"))

# Calculate 95% CI by FMP (BS, AI), the level at which ABC/OFL are defined.
# Follow methods in tpl used to combine survey areas and calculate total
# biomass. TODO: ask PH if there is a reference for these formulas
fmp_biom <- out_sst_nonsst$area_sd %>% 
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
              select(species, Area = FMP, ABC = FMP_ABC, OFL = FMP_OFL) %>% 
              pivot_longer(col = c("ABC", "OFL")) %>% 
              mutate(name = paste0(Area, " ", name)) %>%
              select(-Area)) %>% 
  select(Species = species, Variable = name, Value = value) %>% 
  arrange(Species) 

write_csv(specs_clean, paste0(out_path, "/specs_table_", YEAR, ".csv"))

# Catch ----

# Table 16.1 - percentage of catch of OR in AFSC research bottom trawl surveys
# and in observed fisheries hauls from a) 1991-2001, and b) 2004-2018.

# Table 16.2. -  Regulatory catch limits (OFL, ABC, and TAC) and total catch of
# Other Rockfish in the BSAI and split between BS and AI

catch_fmp <- catch %>% 
  group_by(year, fmp) %>% 
  summarize(catch = sum(tons))

catch_fmp_subarea <- catch %>% 
  group_by(year, fmp_subarea) %>% 
  summarize(catch = sum(tons)) %>% 
  arrange(fmp_subarea, year)

# Table 16.3 - Historical catch by foreign, joint venture, and domestic
# fisheries (STATIC)

# Tables 16.4 and 16.5. - Catch (t) of Other Rockfish species in the Aleutian
# Islands and Bering Sea 2004-2018. Species with catches < 1 ton of catch not
# shown. Also show % SST in thornyhead catch

# proportion of SSTs in thornyhead catch using observer data
prop_sst <- thorny %>% 
  group_by(year, fmp_subarea, thorny) %>% 
  summarize(w = sum(extrapolated_weight)) %>% 
  group_by(year, fmp_subarea) %>% 
  mutate(W = sum(w)) %>% 
  ungroup() %>% 
  mutate(prop_sst = w / W) %>% 
  filter(thorny == "shortspine thornyhead") %>% 
  select(year, fmp_subarea, prop_sst)

# create new species_names for SST and other thornyheads
catch <- catch %>% 
  filter(species_name == "thornyhead") %>%
  left_join(prop_sst) %>% 
  mutate(tons = tons * prop_sst,
         species_name = "SST") %>% 
  bind_rows(catch %>% 
              filter(species_name == "thornyhead") %>%
              left_join(prop_sst) %>% 
              mutate(tons = tons * (1 - prop_sst),
                     species_name = "other thornyheads")) %>% 
  select(- prop_sst) %>% 
  bind_rows(catch %>% 
              filter(species_name != "thornyhead")) 

# only report on the most common species
common_spp <- c("dusky", "SST", "other thornyheads",
                "harlequin", "sharpchin", "yelloweye", "redbanded",
                "redstripe", "black", "other")

# catch by species
catch_spp <- catch %>% 
  # new: define anything that's not a common species as "other" so that it's at
  # least reported somewhere
  mutate(species_name = ifelse(species_name %in% common_spp, species_name, "other")) %>% 
  group_by(year, fmp_subarea, species_name) %>% 
  summarize(catch = sum(tons)) %>% 
  ungroup() %>%
  # total annual catch
  bind_rows(catch_fmp_subarea %>% mutate(species_name = "total")) %>% 
  arrange(fmp_subarea, year, species_name)


catch_spp_wide <- catch_spp %>% 
  pivot_wider(id_cols = c("year","fmp_subarea"),
              names_from = "species_name", values_from = "catch", values_fill = 0) %>% 
  left_join(prop_sst) %>% 
  mutate(perc_sst = prop_sst * 100) %>% 
  select(-prop_sst)

# Final species table
f_catch_spp <- catch_spp_wide %>% 
  mutate(year = as.character(year)) %>% 
  # Average catch by species over all years
  bind_rows(catch_spp_wide %>% 
              group_by(fmp_subarea) %>% 
              summarize_at(vars(-year), list(mean)) %>% 
              mutate(year = "Average")) %>% 
  # Reformat columns
  select(Year = year, Area = fmp_subarea, `dusky rockfish` = dusky, SST, `other thornyheads`,
         `% SST in thornyhead catch` = perc_sst, `harlequin rockfish` = harlequin,
         `yelloweye rockfish` = yelloweye, `redbanded rockfish` = redbanded, 
         `redstripe rockfish` = redstripe, `black rockfish` = black, `other rockfish` = other,
         `Total (t)` = total) %>%   
  mutate_at(vars(-Year, -Area), list(~ formatC(., format = "f", digits = 1))) %>% 
  mutate(`% SST in thornyhead catch` = paste0(`% SST in thornyhead catch`, "%")) %>% 
  arrange(Area)

write_csv(f_catch_spp, paste0(out_path, "/catch_sppXarea_", YEAR, ".csv"))

# Catch by area, complex, and species

catch %>% 
  mutate(complex = ifelse(species_name == "SST", "SST", "nonSST")) %>% 
  group_by(year, fmp_subarea, complex) %>% 
  summarize(catch = sum(tons)) %>%
  ggplot(aes(x = year, y = catch, fill = complex)) +
  geom_area(alpha = 0.6 , size = 0.5, colour = "white") +
  scale_fill_grey() +
  facet_wrap(~ fmp_subarea) +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "grey95"),
        panel.grid.minor = element_line(colour = "grey95"),
        # axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5),
        axis.ticks.x = element_line(colour = "black")) +
  labs(x = NULL, y = "Catch (t)", fill = "Complex")

ggsave(paste0(out_path, "/catch_complexXarea_", YEAR, ".png"), 
       dpi=300, height=4, width=7, units="in")

tmp <- catch_spp %>% 
  filter(!species_name %in% c("total")) %>% 
  mutate(species_name2 = ifelse(species_name %in% c("dusky", "SST", "other thornyheads"),
                                species_name, "other")) 

tmp %>% 
  group_by(year, fmp_subarea, species_name2) %>% 
  summarize(catch = sum(catch)) %>% 
  ungroup() %>% 
  complete(year, species_name2, fmp_subarea, fill = list(catch = 0)) %>% 
  ggplot(aes(x = year, y = catch, fill = species_name2)) +
  geom_area(alpha = 0.6 , size = 0.5, colour = "white") +
  scale_fill_viridis(discrete = TRUE) +
  facet_wrap(~ fmp_subarea) +
  theme_bw() +
  labs(x = NULL, y = "Catch (t)", fill = "Species")

ggsave(paste0(out_path, "/catch_sppXarea_", YEAR, ".png"), 
       dpi=300, height=4, width=7, units="in")

# "Other" Other Rockfish

tmp %>% 
  filter(species_name2 == "other") %>% 
  group_by(year, fmp_subarea, species_name) %>% 
  summarize(catch = sum(catch)) %>% 
  ungroup() %>% 
  complete(year, species_name, fmp_subarea, fill = list(catch = 0)) %>% 
  ggplot(aes(x = year, y = catch, fill = species_name)) +
  geom_area(alpha = 0.6 , size = 0.5, colour = "white") +
  scale_fill_viridis(discrete = TRUE) +
  facet_wrap(~ fmp_subarea) +
  theme_bw() +
  labs(x = NULL, y = "Catch (t)", fill = "Species")

ggsave(paste0(out_path, "/catch_othersppXarea_", YEAR, ".png"), 
       dpi=300, height=4, width=7, units="in")

# Table 16.6. Percentage of the total catch (xx,xxx t) of Other Rockfish species
# from 2004-2018 by target fishery and gear type, over the entire Bering Sea and

catch %>% summarize(sum(tons)) # total catch since 2003

catch <- catch %>% 
  mutate(`Target fishery` = case_when(target == "Atka Mackerel" ~ "Atka mackerel",
                                      target %in% c("Greenland Turbot - BSAI", "Kamchatka Flounder - BSAI",
                                                    "Rock Sole - BSAI", "Flathead Sole", "Other Flatfish - BSAI",
                                                    "Yellowfin Sole - BSAI") ~ "Flatfish",
                                      target == "Pacific Cod" ~ "Pacific cod",
                                      target %in% c("Pollock - midwater", "Pollock - bottom") ~ "Pollock",
                                      target == "Rockfish" ~ "Rockfish",
                                      target == "Halibut" ~ "Halibut",
                                      target == "Sablefish" ~ "Sablefish",
                                      TRUE ~ "Other"),
         # in CAS, TRW is included in the NPT gear group
         `Gear` = case_when(gear %in% c("NPT", "TRW") ~ "Bottom trawl",
                            gear == "PTR" ~ "Pelagic trawl",
                            gear == "POT" ~ "Pot",
                            gear == "HAL" ~ "Longline",
                            gear == "JIG" ~ "Jig"))

# define fisheries by gear and target

complex_byfishery <- catch %>% 
  mutate(complex = ifelse(species_name == "SST", "SST", "nonSST"),
         fishery = paste(`Target fishery`, Gear, sep = " ")) 

complex_byfishery <- complex_byfishery %>% 
  group_by(fishery, complex) %>% 
  summarize(catch = sum(tons)) %>% 
  arrange(complex, -catch) %>%
  # group_by(complex) %>% 
  ungroup() %>% 
  top_n(n = 5, wt = catch) %>% 
  mutate(top = "top") %>% 
  right_join(complex_byfishery) %>% 
  mutate(plot_fishery = case_when(top == "top" ~ fishery,
                                  grepl("trawl", fishery) ~ "Other trawl",
                                  TRUE ~ "Other fixed gear")) %>% 
  group_by(year, nmfs_area, plot_fishery, complex) %>% 
  summarize(catch = sum(tons)) %>% 
  arrange(complex, -catch) 

complex_byfishery %>% 
  ungroup() %>% 
  complete(year, nmfs_area, plot_fishery, complex, fill = list(catch = 0)) %>% 
  filter(nmfs_area %in% c(517, 519, 521, 541, 542, 543)) %>% 
  mutate(nmfs_area2 = as.character(nmfs_area)) %>% 
  mutate(nmfs_area2 = factor(nmfs_area2,
                             levels = c("543", "542", "541",
                                        "517", "519", "521"),
                             labels = c("543 (WAI)", "542 (CAI)", "541 (EAI)",
                                        "517 (SBS)", "519 (SBS)",  "521 (EBS Slope)"),
                             ordered = TRUE)) %>% 
  ggplot(aes(x = year, y = catch, fill = plot_fishery)) +
  geom_area(alpha = 0.6 , size = 0.5, colour = "white") +
  scale_fill_viridis(discrete = TRUE) +
  # scale_fill_grey() +
  facet_grid(complex ~ nmfs_area2) +
  theme_minimal() +
  theme(panel.grid.major = element_line(colour = "grey95"),
        panel.grid.minor = element_line(colour = "grey95"),
        axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5),
        axis.ticks.x = element_line(colour = "black")) +
  labs(x = NULL, y = "Catch (t)", fill = "Fishery")
  
ggsave(paste0(out_path, "/catch_complexXfisheryXarea_", YEAR, ".png"), 
       dpi=300, height=5, width=12, units="in")

# since fishery dynamics seem to be consistent within regions, combine for ease
# of viewing
g <- complex_byfishery %>% 
  filter(nmfs_area %in% c(517, 519, 521, 541, 542, 543)) %>% 
  mutate(region = case_when(nmfs_area %in% c(543, 542, 541) ~ "AI",
                            nmfs_area %in% c(517, 519) ~ "SBS",
                            nmfs_area %in% c(521) ~ "EBS Slope")) %>%
  mutate(region2 = factor(region,
                          levels = c("AI", "SBS", "EBS Slope"),
                          ordered = TRUE)) %>% 
  group_by(year, complex, region2, plot_fishery) %>% 
  summarize(catch = sum(catch)) %>% 
  ungroup() %>% 
  complete(year, region2, plot_fishery, complex, fill = list(catch = 0)) %>% 
  ggplot(aes(x = year, y = catch, fill = plot_fishery)) +
  geom_area(alpha = 0.6 , size = 0.5, colour = "white") +
  scale_fill_viridis(discrete = TRUE) +
  # scale_fill_grey() +
  facet_grid(complex ~ region2) +
  theme_minimal() +
  theme(panel.grid.major = element_line(colour = "grey95"),
        panel.grid.minor = element_line(colour = "grey95"),
        axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5),
        axis.ticks.x = element_line(colour = "black")) +
  labs(x = NULL, y = "Catch (t)", fill = "Fishery")

# library(grid)
# gt <- ggplot_gtable(ggplot_build(g))
# library(gtable)
# gtable_show_layout(gt)
# gt$widths[5] = 1.67 * gt$widths[14]
# gt$widths[7] = 1.33 * gt$widths[14]
# gt$widths[9] = 1 * gt$widths[14]
# gt$widths[5] = 3 * gt$widths[14]
# gt$widths[7] = 2 * gt$widths[14]
# gt$widths[9] = 1 * gt$widths[14]
# 
# png(paste0(out_path, "/catch_complexXfisheryXarea_", YEAR, ".png"), 
#     height=375, width=1050)#res = 300
# grid.draw(gt)
# dev.off()
g
ggsave(paste0(out_path, "/catch_complexXfisheryXarea_", YEAR, "2.png"), 
       dpi=300, height=5, width=12, units="in")

catch_gear <- catch %>% 
  group_by(`Target fishery`, Gear) %>% 
  summarize(catch = sum(tons)) %>% 
  group_by(`Target fishery`, Gear) %>% 
  mutate(sum = sum(catch)) %>% 
  ungroup() %>% 
  mutate(total_sum = sum(catch),
         perc_catch = catch / total_sum * 100)

gear_wide <- catch_gear %>% 
  select(`Target fishery`, Gear, perc_catch) %>% 
  pivot_wider(id_cols = c(`Target fishery`),
              names_from = "Gear", values_from = "perc_catch", values_fill = 0) %>% 
  mutate(Total = `Bottom trawl` + `Pelagic trawl` + `Longline` + `Pot` + `Jig`)

# Final catch by gear/fishery target table
f_catch_gear <- gear_wide %>% 
  arrange(-Total) %>% 
  bind_rows(gear_wide %>% 
              summarize_at(vars(-`Target fishery`), list(sum)) %>% 
              mutate(`Target fishery` = "Total")) %>% 
  mutate_at(vars(-`Target fishery`), list(~ formatC(., format = "f", digits = 1, drop0trailing = TRUE)))

f_catch_gear
write_csv(f_catch_gear, paste0(out_path, "/catch_targetXgear_", YEAR, ".csv"))

# Catch by gear and area over time
catch %>% 
  mutate(Gear2 = Gear) %>% #ifelse(Gear %in% c("Bottom trawl", "Longline"), Gear, "Other")) %>% 
  group_by(year, Gear2, fmp_subarea) %>% 
  summarize(catch = sum(tons)) %>% 
  ungroup() %>% 
  complete(year, Gear2, fmp_subarea, fill = list(catch = 0)) %>% 
  ggplot(aes(x = year, y = catch, fill = Gear2)) +
  geom_area(alpha = 0.6 , size = 0.5, colour = "white") +
  scale_fill_viridis(discrete = TRUE) +
  facet_wrap(~ fmp_subarea) +
  theme_bw() +
  labs(x = NULL, y = "Catch (t)", fill = "Gear")

ggsave(paste0(out_path, "/catch_gearXarea_", YEAR, ".png"), 
       dpi=300, height=4, width=7, units="in")

# Catch by fishery target and area over time
catch %>% 
  mutate(target2 = ifelse(`Target fishery` %in% c("Pollock", "Halibut"), 
                          "Other", `Target fishery`)) %>% 
  group_by(year, target2, fmp_subarea) %>% 
  summarize(catch = sum(tons)) %>% 
  ungroup() %>% 
  complete(year, target2, fmp_subarea, fill = list(catch = 0)) %>% 
  ggplot(aes(x = year, y = catch, fill = target2)) +
  geom_area(alpha = 0.6 , size = 0.5, colour = "white") +
  scale_fill_viridis(discrete = TRUE) +
  facet_wrap(~ fmp_subarea) +
  theme_bw() +
  labs(x = NULL, y = "Catch (t)", fill = "Target\nfishery")

ggsave(paste0(out_path, "/catch_targetXarea_", YEAR, ".png"), 
       dpi=300, height=4, width=7, units="in")

# Table 16.7. - Sum of total catch (t) of EBS dusky rockfish from 2004-2018 by
# target fishery and gear type.

# Function to get catch by species - FMP subarea - Gear type - Target fishery -
# AND NMFS reporting area. Must provide properly filtered data.
detailed_catch_byarea <- function(data) {
  
  gear_area <- data %>% 
    group_by(Gear, `Target fishery`, nmfs_area) %>% 
    summarize(catch = sum(tons)) %>% 
    group_by(`Target fishery`, Gear) %>% 
    mutate(sum = sum(catch)) %>% 
    ungroup() %>% 
    mutate(total_sum = sum(catch),
           perc_tot = sum / total_sum * 100)
  
  # Go wide with NMFS reporting areas
  gear_area_w <- gear_area %>% 
    arrange(Gear, `Target fishery`, nmfs_area) %>% 
    select(`Target fishery`, Gear, nmfs_area, catch, total_sum, perc_tot) %>% 
    pivot_wider(id_cols = c(Gear, `Target fishery`),
                names_from = "nmfs_area", values_from = "catch", values_fill = 0) %>% 
    left_join(gear_area %>% distinct(Gear, `Target fishery`, `Total (t)` = sum, `% of total` = perc_tot))
  
  # Final catch by gear, fishery target, NMFS reporting area table
  f_gear_area <- gear_area_w %>% 
    arrange(-`Total (t)`) %>% 
    bind_rows(gear_area_w %>% 
                summarize_at(vars(-Gear, -`Target fishery`), list(sum)) %>% 
                mutate(Gear = "Total (t)")) %>% 
    mutate_at(vars(-Gear, -`Target fishery`), list(~ formatC(., format = "f", digits = 1, drop0trailing = TRUE)))
  
  return(f_gear_area)
}

# Duskies in the E Bering Sea
df <- catch %>% 
  filter(species_name %in% "dusky" & fmp_subarea == "BS")
dusky_bs <- detailed_catch_byarea(data = df)
dusky_bs
write_csv(dusky_bs, paste0(out_path, "/catch_dusky_BS_targetXgearXarea_2003_", YEAR, ".csv"))

# Duskies in the AI
df <- catch %>% 
  filter(species_name %in% "dusky" & fmp_subarea == "AI")
dusky_ai <- detailed_catch_byarea(data = df)
dusky_ai
write_csv(dusky_ai, paste0(out_path, "/catch_dusky_AI_targetXgearXarea_2003_", YEAR, ".csv"))

# SSTs - have to calculate proportions
sst_df <- catch %>% 
  filter(species_name %in% "thornyhead") %>%
  left_join(prop_sst) %>% 
  mutate(tons = tons * prop_sst,
         species_name = "SST")

# SST in BS
df <- sst_df %>% filter(fmp_subarea == "BS")
sst_bs <- detailed_catch_byarea(data = df)
sst_bs
write_csv(sst_bs, paste0(out_path, "/catch_SST_BS_targetXgearXarea_2003_", YEAR, ".csv"))

# SST in AI
df <- sst_df %>% filter(fmp_subarea == "AI")
sst_ai <- detailed_catch_byarea(data = df)
sst_ai
write_csv(sst_ai, paste0(out_path, "/catch_SST_AI_targetXgearXarea_2003_", YEAR, ".csv"))

# Map

# load("data/nmfs_areas.Rdata")

maps <- c("All Other Rockfish",
          "Non-SST or Thornyhead",
          "Thornyheads")

for(i in 1:length(maps)) {
  
  ii <- maps[i]
  
  # Subset to plot
  if(ii == "All Other Rockfish") {
    sub <- catch %>% filter(year == YEAR)
  } else if (ii == "Non-SST or Thornyhead") {
    sub <- catch %>% 
      filter(species_name != "thornyhead" & year == YEAR) 
  } else if (ii == "Thornyheads") {
    sub <- catch %>% 
      filter(species_name == "thornyhead" & year == YEAR) 
  }
    
  sub <- sub %>% 
    group_by(nmfs_area) %>% 
    summarize(catch = sum(tons)) %>% 
    rename(REP_AREA = nmfs_area)

  # join catch data to nmfs area shape files for mapping
  map_dat <- plyr::join(shp_nmfs, sub, by = c("REP_AREA"))
  
  ggplot() +
    #Land shapefile
    geom_polygon(data = alaska_land, aes(long, lat, group = group),
                 color = "grey60", fill="white") +
    #NMFS Reporting Areas shapefile (credit: Steve.Lewis@noaa.gov)
    geom_polygon(data = map_dat,
                 aes(long, lat, group = group, fill = catch),
                 color = "grey60") +
    coord_equal() +
    geom_text(data = nmfs_areas, aes(coords.x1, coords.x2, label = REP_AREA),
              check_overlap = TRUE, size=3, nudge_y=35000, fontface = "bold") +
    scale_fill_distiller(name="Catch (t)", palette = "YlGnBu", 
                         na.value = "white", guide = "colourbar", trans = "reverse", direction = -1) +# 
    xlim(c(-2540689, 500000)) +
    theme_void() +
    ggtitle(paste0(ii)) +
    theme(plot.title = element_text(hjust = 0.5))
  
  label <- ifelse(ii == "All Other Rockfish",
                  "allorox", ifelse(ii == "Non-SST or Thornyhead", 
                                    "nonthorny", "thorny"))
  
  ggsave(paste0(out_path, "/map_", label, "_", YEAR, ".png"), 
         dpi=300, height=5, width=7.5, units="in")
  
}

# By NMFS reporting area 
catch %>% 
  mutate(complex = ifelse(species_name == "SST", "SST", "nonSST")) %>% 
  group_by(year, nmfs_area, complex) %>% 
  mutate(catch = sum(tons)) %>% 
  distinct(year, complex, nmfs_area, catch) %>%  #, tot_catch)
  ungroup() %>% 
  complete(year, nmfs_area, complex, fill = list(catch = 0)) %>% 
  filter(nmfs_area %in% c(517, 519, 521, 541, 542, 543)) %>% 
  mutate(nmfs_area2 = as.character(nmfs_area)) %>% 
  mutate(nmfs_area2 = factor(nmfs_area2,
                             levels = c("543", "542", "541",
                                        "517", "519", "521"),
                             labels = c("543 (WAI)", "542 (CAI)", "541 (EAI)",
                                        "517 (SBS)", "519 (SBS)",  "521 (EBS Slope)"),
                             ordered = TRUE)) %>% 
  ggplot(aes(x = year, y = catch, fill = complex)) +
  geom_area(alpha = 0.6 , size = 0.5, colour = "white") +
  scale_fill_grey() +
  facet_wrap(~ nmfs_area2, ncol = 6) +
  theme_minimal() +
  theme(panel.grid.major = element_line(colour = "grey95"),
        panel.grid.minor = element_line(colour = "grey95"),
        axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5),
        axis.ticks.x = element_line(colour = "black")) +
  labs(x = NULL, y = "Catch (t)", fill = "Complex")

ggsave(paste0(out_path, "/catch_complexXnmfsarea_", YEAR, ".png"), 
       dpi=300, height=3, width=14, units="in")

# Plot -----
  
# Region-specific biomass

# Flag for zero observations
biom <- biom %>% mutate(zero_obs = ifelse(biomass == 0, TRUE, FALSE))

ggplot() +
  geom_point(data = biom, aes(x = year, y = biomass, col = zero_obs)) +
  geom_errorbar(data = biom, aes(x = year, ymin = lci, ymax = uci, col = zero_obs)) +
  geom_line(data = re_area, aes(x = year, y = re_est), col = "darkgrey") +
  geom_ribbon(data = re_area, aes(x = year, ymin = re_lci, ymax = re_uci),
              fill = "grey80", alpha = 0.4) +
  scale_colour_manual(values = c("black", "red")) +
  facet_grid(region ~ species, scales = "free_y") +
  scale_y_continuous(labels = scales::comma) +
  labs(x = NULL, y = "Biomass (t)") +
  theme_bw() +
  theme(legend.position = "none")

ggsave(paste0(out_path, "/biomass_areaXspp_", YEAR, ".png"), 
       dpi=300, height=6, width=8, units="in")

# Biomass by FMP
ggplot(data = fmp_biom %>% filter(year >= 2000)) + #
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
  select(species, FMP, ABC = FMP_ABC, OFL = FMP_OFL) %>%
  left_join(catch %>%
  mutate(species = ifelse(species_name == "SST", "SST", "nonSST")) %>%
  group_by(year, FMP = fmp_subarea, species) %>%
  summarize(Catch = sum(tons))) %>% 
  pivot_longer(cols = c("ABC", "OFL", "Catch")) %>% 
  mutate(name = factor(name, levels = c("Catch", "ABC", "OFL"), ordered = TRUE)) %>% 
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
              mutate(species = ifelse(species_name == "SST", "SST", "nonSST")) %>%
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
  mutate(species = ifelse(species_name == "SST", "SST", "nonSST")) %>%
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
  # group_by(year, species) %>%
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

# other ----

# other ways the complex has been split up...

# full_biom %>% 
#   mutate(full_complex = "All-OR",
#        complex_bySST = ifelse(species_code == 30020, "SST", "non-SST"),
#        complex_bySSTandDusky = ifelse(species_code == 30020, "SST",
#                                       ifelse(species_code == 30152, "Dusky", "non-SST-or-Dusky"))) %>% 
#   pivot_longer(cols = c("full_complex", "complex_bySST", "complex_bySSTandDusky"), 
#                names_to = "complex", values_to = "species") 