# Queries for BSAI OROX
# Contact: jane.sullivan@noaa.gov
# Last updated: Sep 2020

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

spp %>% filter(grepl("thornyhead|rockfish", common_name))

# Seven most common BSAI Orox
orox <- c("dusky rockfish", "shortspine thornyhead", "redstripe rockfish",
          "redbanded rockfish", "yelloweye rockfish", "harlequin rockfish",
          "sharpchin rockfish")

# Other Orox species that sometimes show up? catch tables
other_orox <- c("black rockfish", "darkblotched rockfish", "rosethorn rockfish",
                "silvergray rockfish", "rockfish unid.", "dusky and dark rockfishes unid.")

codes <- spp %>% 
  filter(common_name %in% c(orox, other_orox)| grepl("thornyhead", common_name)) %>% 
  pull(species_code) 

codes_string <- toString(sprintf("'%s'", codes)) # allows you to pass vector into sql query

spp %>% filter(species_code %in% codes) %>% write_csv(paste0("data/bsai_orox_spp_lookup.csv"))

# Shortspine Thornyhead = SST = 30020 = 350
# Dusky = 30152 = 330

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
