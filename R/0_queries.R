# Queries for BSAI OROX
# Contact: jane.sullivan@noaa.gov
# Last updated: Feb 2022

# devtools::session_info()
# version  R version 4.1.2 (2021-11-01)
# os       Windows 10 x64 (build 19042)
# system   x86_64, mingw32
# ui       RStudio

# Note: also includes analysis using observer data to partition thornyhead catch
# to shortspine (SST) and other thornyhead

# Set up ----

# Assessment year (most recent year with complete data set)
YEAR <- 2022

libs <- c("tidyverse", "RODBC")
if(length(libs[which(libs %in% rownames(installed.packages()) == FALSE )]) > 0) {
  install.packages(libs[which(libs %in% rownames(installed.packages()) == FALSE)])}
lapply(libs, library, character.only = TRUE)

db <- read_csv("database.csv") # user-specific usernames and passwords, not tracked on github

database_afsc <- "afsc"
username_afsc <- db$username_afsc 
password_afsc <- db$password_afsc 
channel_afsc <- odbcConnect(database_afsc, uid = username_afsc, pwd = password_afsc, believeNRows=FALSE)

database_akfin <- "akfin" 
username_akfin <- db$username_akfin
password_akfin <- db$password_akfin
channel_akfin <- odbcConnect(database_akfin, uid = username_akfin, pwd = password_akfin, believeNRows=FALSE)

# Create a year subdirectory to store annual assessment data
dat_path <- paste0("data/", YEAR)
dir.create(dat_path)
raw_path <- paste0(dat_path, "/raw") # raw data
dir.create(raw_path) 

# Species lookup ----

query <- paste0("select   species_code, species_name, common_name
                 from     afsc.race_racespeciescodes")

spp <- sqlQuery(channel_akfin, query) %>% 
  rename_all(tolower) 

spp %>% filter(grepl("dusky", common_name))
spp %>% filter(grepl("thornyhead|rockfish", common_name))
spp %>% filter(grepl("Atka", common_name))

# Seven most common BSAI Orox
orox <- c("dusky rockfish", "shortspine thornyhead", "redstripe rockfish",
          "redbanded rockfish", "yelloweye rockfish", "harlequin rockfish",
          "sharpchin rockfish")

# Other Orox species that sometimes show up? catch tables
other_orox <- c("black rockfish", "darkblotched rockfish", "rosethorn rockfish",
                "silvergray rockfish", "rockfish unid.", "dusky and dark rockfishes unid.",
                "yellowmouth rockfish", "broadfin thornyhead", "longspine thornyhead")

codes <- spp %>% 
  filter(common_name %in% c(orox, other_orox)| grepl("thornyhead", common_name)) %>% 
  pull(species_code) 

codes_string <- toString(sprintf("'%s'", codes)) # allows you to pass vector into sql query

spp %>% 
  filter(species_code %in% codes) %>% 
  write_csv(paste0("data/bsai_orox_spp_lookup.csv"))

# Shortspine Thornyhead = SST = 30020 = 350
# Dusky = 30152 = 330

# AKR/ADFG (3 digits) uses different species codes than RACE (5 digits). You
# need the 3 digit codes to pull observer data

query <- "select   distinct species_code, common_name, scientific_name,
                   race_species_num as species_code_race
          from     norpac.atl_lov_species_code"

akr_spp <- sqlQuery(channel_afsc, query) %>%
  rename_all(tolower) %>%
  filter(species_code_race %in% codes)

# INPFC_AREA look up ----

# International North Pacific Fisheries Commission

query <- "select distinct   survey, inpfc_area as area, summary_area as area_code
          from              afsc.race_goastrataaigoa
          where             survey = 'AI'"

areas <- sqlQuery(channel_akfin, query) %>% rename_all(tolower) 

# AI Survey biomass ----

query <- "select   survey, year, summary_area as area_code, species_code, 
                   area_biomass as biomass, biomass_var as var
          from     afsc.race_biomassinpfcaigoa
          where    species_code in (%s) and
                   survey = 'AI' and
                   year >= 1991"

ai <- sqlQuery(channel_akfin, sprintf(query, codes_string)) %>% 
  write_csv(paste0(raw_path, "/ai_biomass_raw_", YEAR, ".csv"))

ai <- ai %>% 
  rename_all(tolower) %>% 
  left_join(areas) %>% 
  arrange(species_code, year) %>% 
  mutate(area = case_when(area == "Southern Bering Sea" ~ "SBS",
                          area == "Eastern Aleutians" ~ "Eastern AI",
                          area == "Western Aleutians" ~ "Western AI",
                          area == "Central Aleutians" ~ "Central AI"))

# EBS shelf ----

query <- "select   survey, year, species_code, 
                   biomass, varbio as var
          from     afsc.race_biomass_ebsshelf_standard
          where    species_code in (%s) and 
                   stratum in ('999') and 
                   year >= 1982"

ebs_shelf <- sqlQuery(channel_akfin, sprintf(query, codes_string)) %>% 
  write_csv(paste0(raw_path, "/ebs_shelf_biomass_raw_", YEAR, ".csv"))

ebs_shelf <- ebs_shelf %>% 
  rename_all(tolower) %>% 
  arrange(species_code, year) %>% 
  mutate(area = "EBS Shelf")

# EBS slope ----

query <- "select   year, species_code, 
                   stratum_biomass as biomass, bio_var as var
          from     afsc.race_biomass_ebsslope
          where    species_code in (%s) and 
                   stratum in ('999999')"

ebs_slope <- sqlQuery(channel_akfin, sprintf(query, codes_string)) %>% 
  write_csv(paste0(raw_path, "/ebs_slope_biomass_raw_", YEAR, ".csv"))

ebs_slope <- ebs_slope %>% 
  rename_all(tolower) %>% 
  arrange(species_code, year) %>% 
  mutate(survey = "EBS_SLOPE",
         area = "EBS Slope")

# Write biomass data ----

# Bind rows
biom <- ai %>% 
  bind_rows(ebs_shelf) %>% 
  bind_rows(ebs_slope) %>% 
  left_join(spp %>% select(species_code, common_name)) %>% # species names
  write_csv(paste0(dat_path, "/bsai_orox_biomass_", YEAR, ".csv"))  

# BSAI NMFS Area look up ----

# Look up table for NMFS areas in the BSAI FMP

query <- "select  distinct  fmp_area as fmp, fmp_subarea, 
                            reporting_area_code as nmfs_area
          from              council.comprehensive_blend_ca
          where             fmp_area = 'BSAI'"

bsai <- sqlQuery(channel_akfin, query) %>% rename_all(tolower)
bsai <- bsai %>% filter(!is.na(nmfs_area)) %>% write_csv("data/bsai_nmfs_area_lookup.csv")
bsai_lkup <- toString(sprintf("'%s'", pull(bsai, nmfs_area)))

# Fishery lengths ----

# Sampling started in 2002. Lengths only summarized for SSTs and Duskies

# Species are NMFS defn (not RACE);
# Shortspine Thornyhead = SST = 30020 = 350
# Dusky = 30152 = 330

query <- paste0("select   year, nmfs_area, gear, species as species_code, 
                          length, sex, frequency, sample_system
                 from     norpac.debriefed_length
                 where    species in ('350', '330') and
                          nmfs_area in (%s) and
                          year between 2002 and ", YEAR-1) 

fsh_len <- sqlQuery(channel_akfin, sprintf(query, bsai_lkup)) %>% 
  write_csv(paste0(raw_path, "/bsai_fishery_lengths_sstdusky_raw_", YEAR-1, ".csv"))

fsh_len <- fsh_len %>% 
  rename_all(tolower) %>% 
  mutate(species_code = ifelse(species_code == 350, 30020, 30152)) %>% 
  left_join(bsai) %>% 
  mutate(fmp_subarea = ifelse(fmp_subarea == "BS", "EBS", fmp_subarea)) %>% 
  arrange(species_code, year, length) 

# Survey lengths ----

# Only AI survey lengths reported in SAFE, probably due to sample size issues
# but not exactly sure why. Also only for SST (30020) and Dusky (30152)

# Lengths stored in mm, convert to cm

# AI
query <- "select   survey, year, stratum, length / 10 as length, species_code, 
                   total as frequency
          from     afsc.race_sizestratumaigoa
          where    species_code in ('30020','30152') and
                   survey = 'AI'"

ai_len <- sqlQuery(channel_akfin, query) %>% 
  write_csv(paste0(raw_path, "/ai_survey_lengths_sstdusky_raw_", YEAR, ".csv"))

ai_strata <- sqlQuery(channel_akfin,
                      "select distinct   survey, regulatory_area_name, inpfc_area as area,
                            stratum, min_depth, max_depth, area
          from              afsc.race_goastrataaigoa
          where             survey in ('AI')
          order by          regulatory_area_name asc, min_depth asc") %>% 
  rename_all(tolower) %>% 
  distinct(area, stratum, min_depth, max_depth, fmp = regulatory_area_name) %>% 
  mutate(description = paste0(area, ' (', min_depth, '-', max_depth, ' m)')) %>% 
  mutate(description = paste0(min_depth, '-', max_depth, ' m')) %>% 
  distinct(area, stratum, description)

ai_strata

ai_len <- ai_len %>% 
  rename_all(tolower) %>% 
  left_join(ai_strata) %>% 
  filter(area != 'Chirikof') %>% 
  mutate(area = ifelse(grepl('Aleutians', area), 'AI', 'SBS')) %>% 
  arrange(species_code, area, year, length) 

# EBS slope - there are no length data collected for duskies on the shelf (and
# no SST)
query <- "select   survey, year, stratum, length / 10 as length, species_code, 
                   total as frequency
          from     afsc.race_sizecomp_ebsslope
          where    species_code in ('30020','30152')"

ebs_slope_len <- sqlQuery(channel_akfin, query) %>% 
  write_csv(paste0(raw_path, "/ebs_slope_survey_lengths_sstdusky_raw_", YEAR, ".csv"))

ebs_slope_len <- ebs_slope_len %>% 
  rename_all(tolower) %>% 
  arrange(species_code, year, length)

# Catch ----

# 3 digit species codes - can't find a look up table for this one.
# northern rockfish (136) currently in species_group_name = "Other Rockfish"
# group but shouldn't be. removed from bsai_orox3 list

query <- "select   distinct agency_species_code, species_name, species_group_name
          from     council.comprehensive_blend_ca
          where    species_group_name = 'Other Rockfish'"
catch_spp <- sqlQuery(channel_akfin, query)

# dusky = (154, 172)
bsai_orox3 <- c(153, 154, 172, 148, 147, 157, 139, 158, 145, 176, 143, 142, 
                150, 156, 155, 175, 149, 159, 166, 146, 184, 137,
                138, 178, 182, 179)

# these are species classified in catch as "Other Rockfish" that occur primarily
# in GOA. Keep these in just in case any of these spp start to show up in the
# BSAI catch
goa_orox3 <- c(179, 182, 178, 138, 137, 184, 146, 149, 155, 156, 147, 148)

codes_string2 <- toString(sprintf("'%s'", c(bsai_orox3, goa_orox3))) # allows you to pass vector into sql query

query <- "select   *
          from     council.comprehensive_blend_ca
          where    agency_species_code in (%s) and
                   fmp_area = 'BSAI' and
                   year >= 2003"

catch <- sqlQuery(channel_akfin, sprintf(query, codes_string2)) %>% 
  rename_all(tolower) 

catch %>% write_csv(paste0(raw_path, "/bsai_orox_catch_raw_2003_", YEAR+1, "_confidential.csv"))

# Proportion SST ----

# Fishery catch is only reported as "thornyhead", so use observer data to
# estimate proportion SST from other thornyheads

thorny <- akr_spp %>% filter(grepl("THORNY", common_name))

thorny_string <- toString(sprintf("'%s'", thorny$species_code)) 

query <- "select    a.year, a.species, a.species_name, a.sample_number, a.sample_size,
                    a.sample_weight, a.extrapolated_weight, a.extrapolated_number, 
                    a.percent_retained, b.gear_type, b.latdd_start, b.londd_start,
                    b.nmfs_area, a.haul_join
                    
          from      obsint.debriefed_spcomp a
          
          join      obsint.debriefed_haul b on a.haul_join = b.haul_join
          
          where     a.species in (%s) and 
                    a.year >= 2003 and
                    b.nmfs_area between 500 and 543"

observed <- sqlQuery(channel_afsc, sprintf(query, thorny_string)) %>% 
  rename_all(tolower) 

observed %>% write_csv(paste0(raw_path, "/obs_thornyheads_2003_", YEAR+1, "_confidential.csv"))

# observed %>% filter(species == 350) %>% count(year, haul_join) %>% group_by(year) %>% summarize(tst = length(which(n >= 3)))
# observed %>% filter(species == 350) %>% count(haul_join) %>% filter(n > 10)
# observed %>% filter(species == 350) %>% count(year, gear_type, haul_join) %>% group_by(year, gear_type) %>% summarize(n_hauljoins_3plusrows = length(which(n >= 3))) #%>% View()
# observed %>% filter(haul_join == 24213002718000001024)

# Get catch by SST ----

catch_clean <- catch %>% 
  mutate(species = str_replace(species_name, "rockfish, ", ""),
         species = str_replace(species, fixed(" (red snapper)"), ""),
         species = str_replace(species, fixed(" (idiots)"), ""),
         # Re-level Bering Sea so it's consistent with SAFE terminology
         fmp_subarea = ifelse(fmp_subarea == "BS", "EBS", fmp_subarea)) %>% 
  select(year, date = catch_activity_date, fmp = fmp_area, fmp_subarea,
         nmfs_area = reporting_area_code, gear = agency_gear_code,
         retained_or_discarded, target = trip_target_name, 
         species_code = agency_species_code, species_group = species_group_name,
         species_name = species, tons = weight_posted)

# proportion of SSTs in thornyhead catch using observer data
prop_sst <- observed %>% 
  mutate(thorny = ifelse(species == 350, "shortspine thornyhead", "other thornyhead"),
         fmp_subarea = ifelse(between(nmfs_area, 500, 540), "EBS", "AI")) %>% 
  group_by(year, fmp_subarea, thorny) %>% 
  summarize(w = sum(extrapolated_weight)) %>% 
  group_by(year, fmp_subarea) %>% 
  mutate(W = sum(w)) %>% 
  ungroup() %>% 
  mutate(prop_sst = w / W) %>% 
  filter(thorny == "shortspine thornyhead") %>% 
  select(year, fmp_subarea, prop_sst)

# create new species_names for SST and other thornyheads
catch_clean <- catch_clean %>% 
  filter(species_name == "thornyhead") %>%
  left_join(prop_sst) %>% 
  mutate(tons = tons * prop_sst,
         species_name = "SST") %>% 
  bind_rows(catch_clean %>% 
              filter(species_name == "thornyhead") %>%
              left_join(prop_sst) %>% 
              mutate(tons = tons * (1 - prop_sst),
                     species_name = "other thornyheads")) %>% 
  select(- prop_sst) %>% 
  bind_rows(catch_clean %>% 
              filter(species_name != "thornyhead")) 

catch_clean %>% write_csv(paste0(dat_path, "/bsai_orox_catch_2003_", YEAR+1, "_confidential.csv"))

observed %>% 
  mutate(thorny = ifelse(species == 350, "SST", "other thornyhead"),
         fmp_subarea = ifelse(between(nmfs_area, 500, 540), "BS", "AI")) %>% 
  write_csv(paste0(dat_path, "/obs_thornyheads_2003_", YEAR+1, "_confidential.csv"))

prop_sst %>% write_csv(paste0(dat_path, "/prop_sst_catch_2003_", YEAR+1, ".csv"))

# Atka mackerel catch EDA ----

query <- "select   *
          from     council.comprehensive_blend_ca
          where    reporting_area_code in ('541', '542', '543') and 
                   agency_species_code in ('193')  and
                   year >= 2003"

atka <- sqlQuery(channel_akfin, query) %>% 
  rename_all(tolower) 

atka %>% write_csv(paste0(raw_path, "/atka_ai_catch_raw_2003_", YEAR, "_confidential.csv"))
names(atka)
unique(atka$sampling_strata_name)
atka %>%
  group_by(year, sampling_strata_name, reporting_area_code) %>% 
  summarise(catch = sum(weight_posted)) %>% 
  ungroup() %>% 
  tidyr::complete(year, nesting(sampling_strata_name,reporting_area_code ), fill = list(catch = 0)) %>% 
  ggplot(aes(x = year, y = catch, col = sampling_strata_name)) +
  geom_point(size = 2) +
  geom_line(size = 1) +
  facet_wrap(~reporting_area_code)

atka_clean <- atka %>% 
  mutate(species = str_replace(species_name, "greenling, ", "")) %>% 
  select(year, date = catch_activity_date, fmp = fmp_area, fmp_subarea,
         nmfs_area = reporting_area_code, gear = agency_gear_code,
         retained_or_discarded, target = trip_target_name, 
         species_code = agency_species_code, species_group = species_group_name,
         species_name = species, tons = weight_posted) %>% 
  # only include targeted catch
  filter(target == "Atka Mackerel")

atka_clean %>% write_csv(paste0(dat_path, "/atka_ai_targeted_catch_2003_", YEAR+1, "_confidential.csv"))

# Dusky bycatch EDA ----
query <- "select    a.year, a.species, a.species_name,
                    a.extrapolated_weight, b.gear_type, b.latdd_start, b.londd_start,
                    b.nmfs_area, a.haul_join
                    
          from      obsint.debriefed_spcomp a
          
          join      obsint.debriefed_haul b on a.haul_join = b.haul_join
          
          where     a.species in ('355', '330') and 
                    a.year >= 2003 and
                    b.nmfs_area in ('530','550','523','518','517',
                                    '513','509','516','512','508','514','524',
                                    '543','543', '542', '541', '610')"

dusky <- sqlQuery(channel_afsc, query) %>% 
  rename_all(tolower) 

dusky %>% write_csv(paste0(raw_path, "/obs_dusky_2003_", YEAR, "_confidential.csv"))

# LLS ----

lls <- sqlQuery(channel_akfin, query = ("
                select    *
                from      afsc.lls_area_rpn_all_strata
                where     species_code = '30020' and fmp_management_area = 'BSAI'
                order by  year asc
                ")) %>% 
  rename_all(tolower) 
unique(lls$geographic_area_name)

names(lls)

llssum <- lls %>% 
  filter(country == 'United States' & 
           !geographic_area_name %in% c('NW Aleutians slope', 'SW Aleutians slope')) %>% 
  mutate(strata = ifelse(geographic_area_name %in% c('NE Aleutians slope', 'SE Aleutians slope'),
                         'Eastern AI', 'EBS Slope')) %>% 
  group_by(strata, year) %>% 
  dplyr::summarise(cpue = sum(rpw, na.rm = TRUE),
                  cv = sqrt(sum(rpw_var, na.rm = TRUE)) / cpue)

llssum %>% 
  ggplot(aes(x = year, y = cpue)) +
  geom_line() +
  geom_point() +
  facet_wrap(~strata)

llssum %>% write_csv(paste0(dat_path, '/lls_rpw_sst.csv'))

lls %>% 
  filter(!is.na(rpn)) %>% 
  ggplot(aes(x = year, y = rpn, col = country)) +
  geom_point() +
  geom_line() +
  facet_wrap(~geographic_area_name, scales = 'free_y')

# LLS lengths ----

lls_len <- sqlQuery(channel_akfin, query = ("
                select    *
                from      afsc.lls_length_rpn_by_area_all_strata
                where     species_code = '30020' 
                order by  year asc
                ")) %>% 
  rename_all(tolower) 

lls_len_sum <- lls_len %>% 
  filter(grepl(c('Bering|Aleutians'), geographic_area_name)) %>% 
  group_by(year, council_sablefish_management_area, length) %>% 
  dplyr::summarise(rpw = sum(rpw, na.rm = TRUE)) %>% 
  write_csv(paste0(dat_path, '/lls_length_sst.csv'))

# Write lengths ----

comps <- fsh_len %>% 
  select(species_code, year, length, frequency, fmp_subarea) %>% 
  mutate(source = paste0(fmp_subarea, " fishery")) %>% 
  bind_rows(ai_len %>% 
              mutate(fmp_subarea = area,
                     source = "AI BTS") %>% 
              select(species_code, year, fmp_subarea, length, frequency, source)) %>% 
  bind_rows(ebs_slope_len %>% 
              select(species_code, year, length, frequency) %>% 
              mutate(fmp_subarea = "EBS",
                     source = "EBS slope BTS")) %>% 
  bind_rows(lls_len_sum %>% 
              ungroup() %>% 
              mutate(species_code = 30020,
                     fmp_subarea = ifelse(council_sablefish_management_area == 'Aleutians', 'AI', 'EBS'),
                     source = paste0(fmp_subarea, ' LLS')) %>% 
              select(species_code, year, fmp_subarea, length, frequency = rpw, source)) %>% 
  write_csv(paste0(dat_path, "/lengths_sstdusky_", YEAR, ".csv"))
