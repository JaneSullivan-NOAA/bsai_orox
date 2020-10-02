# Queries for BSAI OROX
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

libs <- c("tidyverse", "RODBC")
if(length(libs[which(libs %in% rownames(installed.packages()) == FALSE )]) > 0) {
  install.packages(libs[which(libs %in% rownames(installed.packages()) == FALSE)])}
lapply(libs, library, character.only = TRUE)

db <- read_csv("database.csv") # user-specific usernames and passwords, not tracked on github

database_afsc <- "afsc"
username_afsc <- db$username_afsc #"sullivanj"
password_afsc <- db$password_afsc #"ZwzkgB$7Es69LQ$3"
channel_afsc <- odbcConnect(database_afsc, uid = username_afsc, pwd = password_afsc, believeNRows=FALSE)

database_akfin <- "akfin" 
username_akfin <- db$username_akfin #"jsullivan"
password_akfin <- db$password_akfin #"sculja22"
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
                   survey = 'AI'"

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
                   stratum in ('999')"

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
                          year between 2002 and ", YEAR) 

fsh_len <- sqlQuery(channel_akfin, sprintf(query, bsai_lkup)) %>% 
  write_csv(paste0(raw_path, "/bsai_fishery_lengths_sstdusky_raw_", YEAR, ".csv"))

fsh_len <- fsh_len %>% 
  rename_all(tolower) %>% 
  mutate(species_code = ifelse(species_code == 350, 30020, 30152)) %>% 
  left_join(bsai) %>% 
  arrange(species_code, year, length) 

# Survey lengths ----

# Only AI survey lengths reported in SAFE, probably due to sample size issues
# but not exactly sure why. Also only for SST (30020) and Dusky (30152)

# Lengths stored in mm, convert to cm

query <- "select   survey, year, stratum, length / 10 as length, species_code, 
                   total as frequency
          from     afsc.race_sizestratumaigoa
          where    species_code in ('30020','30152') and
                   survey = 'AI'"

ai_len <- sqlQuery(channel_akfin, query) %>% 
  write_csv(paste0(raw_path, "/ai_survey_lengths_sstdusky_raw_", YEAR, ".csv"))

ai_len <- ai_len %>% 
  rename_all(tolower) %>% 
  arrange(species_code, year, length)

# Write lengths ----

comps <- fsh_len %>% 
  select(species_code, year, length, frequency) %>% 
  mutate(source = "BSAI fishery") %>% 
  bind_rows(ai_len %>% 
              select(species_code, year, length, frequency) %>% 
              mutate(source = "AI survey")) %>% 
  write_csv(paste0(dat_path, "/lengths_sstdusky_", YEAR, ".csv"))
 
# Catch ----

# 3 digit species codes - can't find a look up table for this one.
# northern rockfish (136) currently in species_group_name = "Other Rockfish"
# group but shouldn't be. removed from bsai_orox3 list

# query <- "select   distinct agency_species_code, species_name, species_group_name
#           from     council.comprehensive_blend_ca
#           where    species_group_name = 'Other Rockfish'"
# catch_spp <- sqlQuery(channel_akfin, query)

# dusky = (154, 172)
bsai_orox3 <- c(153, 154, 172, 148, 147, 157, 139, 158, 145, 176, 143, 142, 
                150, 156, 155, 175, 149, 159, 166, 146, 184, 137,
                138, 178, 182, 179)

# these are species classified in catch as "Other Rockfish" that occur primarily
# in GOA. Keep these in just in case any of these spp start to show up in the
# BSAI catch
goa_orox3 <- c(179, 182, 178, 138, 137, 184, 146, 149, 155, 156, 147, 148)

codes_string2 <- toString(sprintf("'%s'", c(bsai_orox3, goa_orox3))) # allows you to pass vector into sql query

# year >= 2019 and

query <- "select   *
          from     council.comprehensive_blend_ca
          where    agency_species_code in (%s) and
                   fmp_area = 'BSAI' and
                   year >= 2003"

catch <- sqlQuery(channel_akfin, sprintf(query, codes_string2)) %>% 
  rename_all(tolower) 

catch %>% write_csv(paste0(raw_path, "/bsai_orox_catch_raw_2003_", YEAR+1, "_confidential.csv"))

catch_clean <- catch %>% 
  mutate(species = str_replace(species_name, "rockfish, ", ""),
         species = str_replace(species, fixed(" (red snapper)"), ""),
         species = str_replace(species, fixed(" (idiots)"), "")) %>% 
  select(year, date = catch_activity_date, fmp = fmp_area, fmp_subarea,
         nmfs_area = reporting_area_code, gear = agency_gear_code,
         retained_or_discarded, target = trip_target_name, 
         species_code = agency_species_code, species_group = species_group_name,
         species_name = species, tons = weight_posted)

catch_clean %>% write_csv(paste0(dat_path, "/bsai_orox_catch_2003_", YEAR+1, "_confidential.csv"))

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

observed %>% write_csv(paste0(raw_path, "/obs_thornyheads_2003_", YEAR+1, ".csv"))

observed %>% 
  mutate(thorny = ifelse(species == 350, "shortspine thornyhead", "other thornyhead"),
         fmp_subarea = ifelse(between(nmfs_area, 500, 540), "BS", "AI")) %>% 
  write_csv(paste0(dat_path, "/obs_thornyheads_2003_", YEAR+1, ".csv"))

# observed %>% filter(species == 350) %>% count(year, haul_join) %>% group_by(year) %>% summarize(tst = length(which(n >= 3)))
# observed %>% filter(species == 350) %>% count(haul_join) %>% filter(n > 10)
# observed %>% filter(species == 350) %>% count(year, gear_type, haul_join) %>% group_by(year, gear_type) %>% summarize(n_hauljoins_3plusrows = length(which(n >= 3))) #%>% View()
# observed %>% filter(haul_join == 24213002718000001024)

# # EBS biomass EDA ----
# 
# # There are multiple EBS Biomass options in AFKIN.AFSC.
# # RACE_BIOMASS_EBSS_PLUSNW_GRPD and RACE_BIOMASS_EBSS_STANDARD_GRP tables return
# # empty dataframes for our spp. I ultimately went with
# # RACE_BIOMASS_EBSSHELF_STANDARD (#3 here) because it includes 1984 and 1985
# # data and otherwise matches RACE_BIOMASS_EBSSHELF_PLUSNW. It also matches
# # Ingrid Spies' 2018 All_survey_data_RE.xlsx sheet titled "EBS"
# 
# # EBS shelf 1: RACE_BIOMASS_EBSSHELF
# query <- "select   survey, year, species_code, 
#                    stratum_biomass as biomass, bio_var as var
#           from     afsc.race_biomass_ebsshelf
#           where    species_code in (%s) and 
#                    stratum in ('999')"
# 
# ebs_shelf <- sqlQuery(channel_akfin, sprintf(query, codes_string)) %>% 
#   rename_all(tolower) %>% 
#   arrange(species_code, year) %>% 
#   mutate(table_name = "RACE_BIOMASS_EBSSHELF")
# 
# ebs_shelf <- ebs_shelf %>% group_by(table_name, year) %>% summarize(biomass = sum(biomass), var = sum(var)) %>% print(n = Inf)
# 
# # EBS shelf 2: RACE_BIOMASS_EBSSHELF_PLUSNW
# query <- "select   survey, year, species_code, 
#                    biomass, varbio as var
#           from     afsc.race_biomass_ebsshelf_plusnw
#           where    species_code in (%s) and 
#                    stratum in ('999')"
# 
# ebs_shelf2 <- sqlQuery(channel_akfin, sprintf(query, codes_string)) %>% 
#   rename_all(tolower) %>% 
#   arrange(species_code, year) %>% 
#   mutate(table_name = "RACE_BIOMASS_EBSSHELF_PLUSNW")
# 
# ebs_shelf2 %>% dim()
# ebs_shelf2 <- ebs_shelf2 %>% group_by(table_name, year) %>% summarize(biomass = sum(biomass), var = sum(var)) %>% print(n = Inf)
# 
# # EBS shelf 3: RACE_BIOMASS_EBSSHELF_STANDARD
# query <- "select   survey, year, species_code, 
#                    biomass, varbio as var
#           from     afsc.race_biomass_ebsshelf_standard
#           where    species_code in (%s) and 
#                    stratum in ('999')"
# 
# ebs_shelf3 <- sqlQuery(channel_akfin, sprintf(query, codes_string)) %>% 
#   rename_all(tolower) %>% 
#   arrange(species_code, year) %>% 
#   mutate(table_name = "RACE_BIOMASS_EBSSHELF_STANDARD")
# 
# ebs_shelf3 %>% dim()
# ebs_shelf3 <- ebs_shelf3 %>% group_by(table_name, year) %>% summarize(biomass = sum(biomass), var = sum(var)) %>% print(n = Inf)
# 
# ebs_shelf %>% 
#   bind_rows(ebs_shelf2) %>% 
#   bind_rows(ebs_shelf3) %>% 
#   ggplot(aes(x = year, y = biomass, fill = table_name)) +
#   geom_bar(stat = "identity", position = "dodge", width = 0.8, col = "black") 
# 
# ebs_shelf3 %>% filter(year %in% c(1984, 1985))
# 
