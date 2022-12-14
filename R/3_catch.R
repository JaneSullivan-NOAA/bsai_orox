# Catch tables
# Contact: jane.sullivan@noaa.gov
# Last updated: Sep 2022

# devtools::session_info()
# version  R version 4.0.2 (2020-06-22)
# os       Windows 10 x64              
# system   x86_64, mingw32             
# ui       RStudio 

# Set up ----

# Assessment year (most recent year with complete data set)
YEAR <- 2022

# date that catch was queried - important for long-term record keeping, tracking
# changes to db and tracing discrepancies in catch
access_date <- "10/03/2022" 

libs <- c("tidyverse", "R2admb", "viridis")
if(length(libs[which(libs %in% rownames(installed.packages()) == FALSE )]) > 0) {
  install.packages(libs[which(libs %in% rownames(installed.packages()) == FALSE)])}
lapply(libs, library, character.only = TRUE)

options(scipen = 999)

# Directory setup ----

root <- getwd()
dat_path <- paste0("data/", YEAR) # directory where source data is contained
out_path <- paste0("results/", YEAR) # directory for results/output
dir.create(out_path)

# Data ----

setwd(root)

# catch data
catch <- read_csv(paste0(dat_path, "/bsai_orox_catch_2003_", YEAR+1, "_confidential.csv"))

# proportion SST by FMP subarea
prop_sst <- read_csv(paste0(dat_path, "/prop_sst_catch_2003_", YEAR+1, ".csv"))

# ak land and nmfs area shapefiles (these are already cleaned up and saved as a
# dataframe)
load("data/map_ak_land.RData")
load("data/map_nmfs_areas.RData")

# directed atka catch
atka <- read_csv(paste0(dat_path, "/atka_ai_targeted_catch_2003_", YEAR, "_confidential.csv"))

# non-commercial catch - have to pull this from the AKFIN Answers dashboard
agency_catch <- read_csv(paste0(dat_path, "/Noncommercial_Fishery_Catch_AKFIN_DASH_20201019.csv"))

# Non-commercial catch ----

# the "Other Rockfish" species group is a little outdated.

# bsai_orox <- c(153, 154, 172, 148, 147, 157, 139, 158, 145, 176, 143, 142, 
#                 150, 156, 155, 175, 149, 159, 166, 146, 184, 137,
#                 138, 178, 182, 179)

agency_catch <- agency_catch %>% 
  select_all(~gsub("\\s+|\\.", "_", .)) %>% 
  select_all(tolower) %>% 
  rename(year = collection_year)
  
agency_catch %>% #distinct(species_code)
  # filter(species_code %in% bsai_orox) %>% 
  group_by(year, agency = collection_agency) %>% 
  dplyr::summarise(catch = sum(weight) / 1000) %>% 
  bind_rows(agency_catch %>% 
              group_by(year) %>% 
              dplyr::summarise(catch = sum(weight) / 1000) %>% 
              mutate(agency = "Total (t)")) %>% 
  ungroup() %>% 
  pivot_wider(id_cols = year, names_from = agency, values_from = catch)

# Catch by FMP ----

report <- c(paste0("Catch tables for BSAI Other Rockfish", "\n",
                   "Contact: jane.sullivan@noaa.gov", "\n",
                   "Catch assessed from AKFIN's COUNCIL.COMPREHENSIVE_BLEND_CA ", access_date))

write.table(report, file = paste0(out_path, "/catch_tables_", YEAR, ".csv"), sep=",", quote = FALSE, row.names = FALSE, col.names = FALSE, eol = "\n")

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

f_catcharea <- catch_fmp %>% 
  pivot_wider(id_cols = year, names_from = fmp, values_from = catch) %>% 
  left_join(catch_fmp_subarea %>% 
              pivot_wider(id_cols = year, names_from = fmp_subarea, values_from = catch)) %>% 
  rename(Year = year)

write.table(c(paste0("\n", "Total catch (t) by FMP and FMP subarea")), 
            file = paste0(out_path, "/catch_tables_", YEAR, ".csv"), sep=",", quote = FALSE, row.names = FALSE, col.names = FALSE, eol = "\n", append = TRUE)
write.table(f_catcharea, file = paste0(out_path, "/catch_tables_", YEAR, ".csv"), sep=",", quote = FALSE, row.names = FALSE, col.names = TRUE, eol = "\n", append = TRUE)

# Table 16.3 - Historical catch by foreign, joint venture, and domestic
# fisheries (STATIC)

# By species -----

# Tables 16.4 and 16.5. - Catch (t) of Other Rockfish species in the Aleutian
# Islands and Bering Sea 2004-2018. Species with catches < 1 ton of catch not
# shown. Also show % SST in thornyhead catch

# only report on the most common species
common_spp <- c("dusky", "SST", "other thornyheads",
                "harlequin", "sharpchin", "yelloweye", "redbanded",
                "redstripe", "black", "other")

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

write.table(c(paste0("\n", "Catch (t) of Other Rockfish species and percent shortspine thornyhead in the catch in the Aleutian
Islands and Bering Sea since 2003.")), 
            file = paste0(out_path, "/catch_tables_", YEAR, ".csv"), sep=",", quote = FALSE, row.names = FALSE, col.names = FALSE, eol = "\n", append = TRUE)
write.table(f_catch_spp, file = paste0(out_path, "/catch_tables_", YEAR, ".csv"), sep=",", quote = FALSE, row.names = FALSE, col.names = TRUE, eol = "\n", append = TRUE)

# Complex by NMFS area ----

catch %>% 
  mutate(complex = ifelse(species_name == "SST", "SST", "non-SST")) %>% 
  mutate(complex = factor(complex, 
                   levels = c("SST", "non-SST"), 
                   ordered = TRUE)) %>% 
  group_by(year, nmfs_area, complex) %>% 
  mutate(catch = sum(tons)) %>% 
  distinct(year, complex, nmfs_area, catch) %>% 
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
  # scale_fill_grey() +
  scale_fill_manual(values = c("#6ba292", "#f0b74a")) +
  theme_minimal(base_size = 14) +
  facet_wrap(~ nmfs_area2, ncol = 6) +
  theme_minimal() +
  theme(panel.grid.major = element_line(colour = "grey95"),
        panel.grid.minor = element_line(colour = "grey95"),
        axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5),
        axis.ticks.x = element_line(colour = "black")) +
  labs(x = NULL, y = "Catch (t)", fill = "Complex")

ggsave(paste0(out_path, "/catch_complexXnmfsarea_", YEAR, ".png"), 
       dpi=300, height=3, width=12, units="in")

# By area, complex, and species ----

catch %>% 
  mutate(complex = ifelse(species_name == "SST", "SST", "non-SST")) %>% 
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

catch %>% 
  mutate(complex = ifelse(species_name == "SST", "SST", "non-SST")) %>% 
  group_by(year, complex) %>% 
  summarize(catch = sum(tons)) %>% 
  pivot_wider(id_cols = c(year), 
              names_from = complex, values_from = catch, 
              values_fill = 0) %>% 
  mutate(prop_nonsst = `non-SST` / (`non-SST` + SST) * 100) %>% 
  ungroup() %>% 
  mutate(mean_prop = mean(prop_nonsst))

100-56.1

data.frame(complex = c("SST", "non-SST", "SST", "non-SST"),
           prop = c(43.9, 56.1, 95, 5),
           var = c("Catch", "Catch", "Biomass", "Biomass")) %>% 
  ggplot(aes(x = "", y = prop, fill = complex)) +
  geom_bar(stat = 'identity', width = 1, col = 'white', alpha = .7) +
  labs(fill = NULL) +
  scale_fill_manual(values = c("#f0b74a","#6ba292")) +
  coord_polar("y", start = 0) +
  theme_void(base_size = 18) +
  facet_wrap(~var)

ggsave(paste0(out_path, "/piechart_", YEAR, ".png"), 
       dpi=300, height=4, width=7, units="in")

# By area and major species
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

# By gear/target ----

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
  mutate(complex = ifelse(species_name == "SST", "SST", "non-SST"),
         fishery = paste(`Target fishery`, Gear, sep = " ")) 

# Top 5 fisheries by catch and nmfs area
complex_byfishery <- complex_byfishery %>% 
  group_by(fishery, complex) %>% 
  summarize(catch = sum(tons)) %>% 
  arrange(complex, -catch) %>%
  ungroup() %>% 
  top_n(n = 5, wt = catch) %>% 
  mutate(top = "top") %>% 
  right_join(complex_byfishery) %>% 
  mutate(plot_fishery = case_when(top == "top" ~ fishery,
                                  grepl("trawl", fishery) ~ "Other trawl",
                                  TRUE ~ "Other fixed gear"))

# Complex by fishery and FMP subarea
complex_byfishery_fmp <- complex_byfishery %>% 
  group_by(year, fmp_subarea, plot_fishery, complex) %>% 
  summarize(catch = sum(tons)) %>% 
  arrange(complex, year) %>% 
  ungroup() 

complex_byfishery_fmp %>% 
  complete(year, fmp_subarea, plot_fishery, complex, fill = list(catch = 0)) %>% 
  mutate(complex = factor(complex, levels = c("SST", "non-SST"), ordered = TRUE)) %>% 
  ggplot(aes(x = year, y = catch, fill = plot_fishery)) +
  geom_area(alpha = 0.6 , size = 0.5, colour = "white") +
  scale_fill_viridis(discrete = TRUE) +
  # facet_grid(fmp_subarea ~ complex) +
  facet_grid(complex ~ fmp_subarea) +
  theme_bw() +
  # theme_minimal() +
  theme(panel.grid.major = element_line(colour = "grey95"),
        panel.grid.minor = element_line(colour = "grey95"),
        axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5),
        axis.ticks.x = element_line(colour = "black")) +
  labs(x = NULL, y = "Catch (t)", fill = "Fishery")

ggsave(paste0(out_path, "/catch_complexXfisheryXfmp_", YEAR, ".png"), 
       dpi=300, height=6, width=7, units="in")

# Complex by fishery and NMFS area
complex_byfishery_nmfsarea <- complex_byfishery %>% 
  group_by(year, nmfs_area, plot_fishery, complex) %>% 
  summarize(catch = sum(tons)) %>% 
  arrange(complex, -catch) 

complex_byfishery_nmfsarea %>% 
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
  mutate(complex = factor(complex, levels = c("SST", "non-SST"), ordered = TRUE)) %>% 
  ggplot(aes(x = year, y = catch, fill = plot_fishery)) +
  geom_area(alpha = 0.6 , size = 0.5, colour = "white") +
  scale_fill_viridis(discrete = TRUE, option = "D") +
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
g <- complex_byfishery_nmfsarea %>% 
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
write.table(c(paste0("\n", "Percentage of the total catch of Other Rockfish species 
since 2003 by target fishery and gear type.")), 
            file = paste0(out_path, "/catch_tables_", YEAR, ".csv"), sep=",", quote = FALSE, row.names = FALSE, col.names = FALSE, eol = "\n", append = TRUE)
write.table(f_catch_gear, file = paste0(out_path, "/catch_tables_", YEAR, ".csv"), sep=",", quote = FALSE, row.names = FALSE, col.names = TRUE, eol = "\n", append = TRUE)

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

# Detailed species by NMFS area ----

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
  filter(species_name %in% "dusky" & fmp_subarea == "EBS")
dusky_bs <- detailed_catch_byarea(data = df)
dusky_bs

write.table(c(paste0("\n", "Sum of total catch (t) of EBS dusky rockfish since 2003 by target fishery and gear type.")), 
            file = paste0(out_path, "/catch_tables_", YEAR, ".csv"), sep=",", quote = FALSE, row.names = FALSE, col.names = FALSE, eol = "\n", append = TRUE)
write.table(dusky_bs, file = paste0(out_path, "/catch_tables_", YEAR, ".csv"), sep=",", quote = FALSE, row.names = FALSE, col.names = TRUE, eol = "\n", append = TRUE)

# Duskies in the AI
df <- catch %>% 
  filter(species_name %in% "dusky" & fmp_subarea == "AI")
dusky_ai <- detailed_catch_byarea(data = df)
dusky_ai

write.table(c(paste0("\n", "Sum of total catch (t) of AI dusky rockfish since 2003 by target fishery and gear type.")), 
            file = paste0(out_path, "/catch_tables_", YEAR, ".csv"), sep=",", quote = FALSE, row.names = FALSE, col.names = FALSE, eol = "\n", append = TRUE)
write.table(dusky_ai, file = paste0(out_path, "/catch_tables_", YEAR, ".csv"), sep=",", quote = FALSE, row.names = FALSE, col.names = TRUE, eol = "\n", append = TRUE)

# SST in BS
df <- catch %>% 
  filter(species_name %in% "SST" & fmp_subarea == "EBS")

sst_bs <- detailed_catch_byarea(data = df)
sst_bs

write.table(c(paste0("\n", "Sum of total catch (t) of EBS shortspine thornyhead since 2003 by target fishery and gear type.")), 
            file = paste0(out_path, "/catch_tables_", YEAR, ".csv"), sep=",", quote = FALSE, row.names = FALSE, col.names = FALSE, eol = "\n", append = TRUE)
write.table(sst_bs, file = paste0(out_path, "/catch_tables_", YEAR, ".csv"), sep=",", quote = FALSE, row.names = FALSE, col.names = TRUE, eol = "\n", append = TRUE)

# SST in AI
df <- catch %>% 
  filter(species_name %in% "SST" & fmp_subarea == "AI")
sst_ai <- detailed_catch_byarea(data = df)
sst_ai

write.table(c(paste0("\n", "Sum of total catch (t) of AI shortspine thornyhead since 2003 by target fishery and gear type.")), 
            file = paste0(out_path, "/catch_tables_", YEAR, ".csv"), sep=",", quote = FALSE, row.names = FALSE, col.names = FALSE, eol = "\n", append = TRUE)
write.table(sst_ai, file = paste0(out_path, "/catch_tables_", YEAR, ".csv"), sep=",", quote = FALSE, row.names = FALSE, col.names = TRUE, eol = "\n", append = TRUE)

# Discarded/retained ----

f_perc_discarded <- catch %>% 
  group_by(Area = fmp_subarea, Year = year, retained_or_discarded) %>% 
  summarize(catch = sum(tons)) %>% 
  mutate(retained_or_discarded = ifelse(retained_or_discarded == "R", 
                                        "Retained", "Discarded")) %>% 
  pivot_wider(id_cols = c(Area, Year), 
              names_from = retained_or_discarded, values_from = catch) %>% 
  mutate(`Total catch` = Retained + Discarded,
         `Percent discarded` = paste0(formatC(Discarded / `Total catch` * 100, 
                                              format = "f", digits = 1), "%")) %>% 
  mutate_at(vars(Discarded, Retained, `Total catch`), 
            list(~ formatC(., format = "f", digits = 0)))

f_perc_discarded

write.table(c(paste0("\n", "Retained and discarded catch of Other Rockfish species since 2003 in the Aleutian Islands and Bering Sea. ")), 
            file = paste0(out_path, "/catch_tables_", YEAR, ".csv"), sep=",", quote = FALSE, row.names = FALSE, col.names = FALSE, eol = "\n", append = TRUE)
write.table(f_perc_discarded, file = paste0(out_path, "/catch_tables_", YEAR, ".csv"), sep=",", quote = FALSE, row.names = FALSE, col.names = TRUE, eol = "\n", append = TRUE)

drates <- catch %>% 
  mutate(complex = ifelse(species_name == "SST", "SST", "non-SST"),
         retained_or_discarded = ifelse(retained_or_discarded == "D", "Discarded", "Retained")) %>% 
  group_by(complex, year) %>% 
  summarize(discard_rate = sum(tons[retained_or_discarded == "Discarded"]) / sum(tons)) %>% 
  group_by( complex) 

drates %>% 
  ggplot(aes(x = year, y = discard_rate, col = fmp_subarea)) +
  geom_line() +
  facet_wrap(~complex)

drates %>% 
  summarize(mean_discard_rate = mean(discard_rate),
         sd_discard_rate = sd(discard_rate)) %>% 
  print(n=Inf)

catch %>% 
  mutate(complex = ifelse(species_name == "SST", "SST", "non-SST"),
         retained_or_discarded = ifelse(retained_or_discarded == "D", "Discarded", "Retained")) %>% 
  group_by(fmp_subarea, complex, year, retained_or_discarded) %>% 
  summarize(catch = sum(tons)) %>% 
  distinct(fmp_subarea, complex, year, retained_or_discarded, catch) %>% 
  ungroup() %>% 
  complete(fmp_subarea, complex, year, retained_or_discarded, fill = list(catch = 0)) %>% 
  mutate(complex = factor(complex, levels = c("SST", "non-SST"), ordered = TRUE)) %>% 
  ggplot(aes(x = year, y = catch, fill = retained_or_discarded)) +
  geom_area(alpha = 0.6 , size = 0.5, colour = "white") +
  scale_fill_grey() +
  facet_grid(complex ~ fmp_subarea) +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "grey95"),
        panel.grid.minor = element_line(colour = "grey95"),
        axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5),
        axis.ticks.x = element_line(colour = "black")) +
  labs(x = NULL, y = "Catch (t)", fill = NULL)

ggsave(paste0(out_path, "/retained_or_discarded_complexXfmp_", YEAR, ".png"), 
       dpi=300, height=6, width=7, units="in")

# Atka in the AI ----

library(scales)
  
relate_atka <- atka %>% 
    mutate(complex = "Atka mackerel") %>% 
    bind_rows(catch %>% 
                filter(fmp_subarea == "AI") %>% # & target == "Atka Mackerel") %>% 
                mutate(complex = ifelse(species_name == "SST", "SST", "non-SST")) %>% 
                filter(complex == "non-SST")) %>%
                # filter(species_name == "dusky")) %>% 
    # filter(nmfs_area %in% c(541, 542, 543)) %>% 
    filter(nmfs_area %in% c(541)) %>%
    mutate(nmfs_area2 = as.character(nmfs_area)) %>% 
    mutate(nmfs_area2 = factor(nmfs_area2,
                               levels = c("541"),
                               labels = c("541 (EAI)"))) %>% 
  group_by(complex, year, nmfs_area2) %>% 
  summarize(catch = sum(tons)) %>% 
  complete(year, nmfs_area2, complex, fill = list(catch = 0)) %>% 
  group_by(complex, nmfs_area2) %>% 
  mutate(norm_catch = scales::rescale(catch, to = c(-1, 1))) 
  
relate_atka %>% 
  ggplot(aes(x = year, y = norm_catch, colour = complex, lty = complex)) +
  geom_line() +
  facet_wrap(~nmfs_area2) +
  theme_minimal() +
  labs(x = NULL, y = "Standardized catch by species group", colour = NULL, lty = NULL)
  

relate_atka %>% 
  ungroup() %>% 
  pivot_wider(id_cols = c("nmfs_area2", "year"), names_from = "complex", 
              values_from = "norm_catch") %>% 
  group_by(nmfs_area2) %>% 
  summarize(cor(`Atka mackerel`, `non-SST`))

# Export data set for JVB ----
catch %>% 
  select(year, date, nmfs_area, target_fishery = `Target fishery`)
  

# Map catch ----

maps <- c("All Other Rockfish")#,
          # "Non-SST",
          # "SST")

# Do you want all years or just current year? Toggle manually.
# map_yrs <- YEAR
map_yrs <- unique(catch$year)

for(i in 1:length(maps)) {
  
  ii <- maps[i]
  
  # Subset to plot
  if(ii == "All Other Rockfish") {
    sub <- catch %>% filter(year %in% map_yrs)
  } else if (ii == "Non-SST") {
    sub <- catch %>% 
      filter(species_name != "SST" & year %in% map_yrs) 
  } else if (ii == "SST") {
    sub <- catch %>% 
      filter(species_name == "SST" & year %in% map_yrs) 
  }
  
  sub <- sub %>% 
    group_by(nmfs_area) %>% 
    summarize(catch = sum(tons)) %>% 
    rename(REP_AREA = nmfs_area)
  
  # Labels
  label <- ifelse(ii == "All Other Rockfish",
                  "allorox", ifelse(ii == "SST", 
                                    "SST", "non-SST"))
  
  yr_labs <- ifelse(length(map_yrs) > 1, paste0(min(map_yrs), "_", max(map_yrs)),
                    YEAR)
  
  # join catch data to nmfs area shape files for mapping
  map_dat <- plyr::join(nmfs_areas, sub, by = c("REP_AREA"))
  
  ggplot() +
    #Land shapefile
    geom_polygon(data = alaska_land, aes(long, lat, group = group),
                 color = "grey60", fill="white") +
    #NMFS Reporting Areas shapefile (credit: Steve.Lewis@noaa.gov)
    geom_polygon(data = map_dat,
                 aes(long, lat, group = group, fill = catch),
                 color = "grey60", alpha = 0.8) +
    coord_equal() +
    geom_text(data = nmfs_areas, aes(coords.x1, coords.x2, label = REP_AREA),
              check_overlap = TRUE, size=3, nudge_y=35000, fontface = "bold") +
    scale_fill_distiller(name="Catch (t)", palette = "Greys", #palette = "YlGnBu", 
                         na.value = "white", guide = "colourbar", trans = "reverse", direction = -1) +# 
    xlim(c(-2540689, 500000)) +
    theme_void() +
    ggtitle(paste0(ii, "\n", yr_labs)) +
    theme(plot.title = element_text(hjust = 0.5))
  
  ggsave(paste0(out_path, "/map_", label, "_", yr_labs, ".png"), 
         dpi=300, height=5, width=7.5, units="in")
  
}

# detailed map ----
library(rgdal)
# sp_data <- readOGR(dsn = "data/shapefiles",layer = "P3_AK_All")
# sp_data <- readOGR(dsn = "data/shapefiles",layer = "NMFS_Zones_Clean")
sp_data <- readOGR(dsn = "data/shapefiles",layer = "gf95_nmfs_polygon")

dusky <- read_csv(paste0(raw_path, "/obs_dusky_2003_", YEAR, "_confidential.csv"))
names(dusky)

dusky_sm <- dusky %>% 
  filter(year %in% c(2015:2020)) %>% 
  # filter(nmfs_area %in% c(541, 542, 543, 610)) %>% 
  mutate(fmp = ifelse(nmfs_area > 600, 2, 4)) %>% 
  select(lon = londd_start, lat = latdd_start, haul_join, catch = extrapolated_weight, year, haul_join, nmfs_area, fmp) %>% 
  filter(!(is.na(lon) | is.na(lat)))

dusky_sm$lon = ifelse(dusky_sm$lon > 0, dusky_sm$lon - 360, dusky_sm$lon)

dusky_coords <- SpatialPointsDataFrame(coords = dusky_sm[, c(1, 2)],
                                       data = dusky_sm,
                                       proj4string = CRS("+proj=longlat"))
dusky_coords <- spTransform(dusky_coords, CRS(proj4string(sp_data)))

plot(sp_data, xlim = c(-2750000,-250000), ylim = c(0,2000000))
points(dusky_coords[ , c(1, 2)], pch=21, cex = dusky_coords$catch/2000, col = dusky_sm$fmp)

newproj <- "+proj=longlat +ellps=WGS84 +no_defs" # +no_defs" #CRS(proj4string(sp_data))
library(raster)
r1 <- marmap::as.raster(aleu)
r2 <- raster::projectRaster(r1, newproj)
aleu_proj <- as.bathy(r2)

title(main = paste0(unique(dusky_sm$year)))

# bathy -----

aleu <- getNOAA.bathy(170, -158, 50, 70, 
                      resolution = 10,
                      antimeridian = TRUE)

blues <- c("lightsteelblue4", "lightsteelblue3",
           "lightsteelblue2", "lightsteelblue1")

plot(aleu, image = TRUE, 
     # bpal = list(c(0,max(aleu),"lightgrey"),
     #             c(min(aleu),0,"#1f66e5","#dfe9fb")),
     # bpal = blues(100),
     bpal = list(c(0, max(aleu), "grey"),
                 c(min(aleu),0,blues)),
     land = T, lwd = 0.1, 
     col = c("lightgrey", "grey", "darkgrey", "grey", "lightgrey", "grey", "darkgrey", "grey", "lightgrey", "grey", "darkgrey", "grey", "lightgrey"),
     axes = FALSE)
# antimeridian.box(aleutians, 10)

dusky_sm <- dusky %>% 
  dplyr::filter(year %in% c(2015:2020)) %>% 
  # filter(nmfs_area %in% c(541, 542, 543, 610)) %>% 
  mutate(fmp = ifelse(nmfs_area == 610, 4, 6)) %>% 
  dplyr::select(lon = londd_start, lat = latdd_start, haul_join, catch = extrapolated_weight, year, haul_join, nmfs_area, fmp) %>% 
  filter(!(is.na(lon) | is.na(lat)))

dusky_sm$lon = ifelse(dusky_sm$lon < 0, dusky_sm$lon + 360, dusky_sm$lon)

points(dusky_sm[ , c(1, 2)], pch=21,
       cex = dusky_sm$catch/2500, 
       col = dusky_sm$fmp)

plot(sp_data)
nmfs <- fortify(sp_data)
summary(nmfs)
nmfs$long <- nmfs$long/1000
nmfs$lat <- nmfs$lat/1000
plot(aleutians, add = TRUE)
