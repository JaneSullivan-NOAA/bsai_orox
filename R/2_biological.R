# Biological data (length comps) for SST and Duskies
# Contact: jane.sullivan@noaa.gov
# Last updated: Oct 2020

# devtools::session_info()
# version  R version 4.0.2 (2020-06-22)
# os       Windows 10 x64              
# system   x86_64, mingw32             
# ui       RStudio 

# Set up ----

# Assessment year (most recent year with complete data set)
YEAR <- 2022
dat_path <- paste0("data/", YEAR) # directory where source data is contained
out_path <- paste0("results/", YEAR) # directory for results/output
dir.create(out_path)

libs <- c("tidyverse", "ggridges")
if(length(libs[which(libs %in% rownames(installed.packages()) == FALSE )]) > 0) {
  install.packages(libs[which(libs %in% rownames(installed.packages()) == FALSE)])}
lapply(libs, library, character.only = TRUE)

# Data ----

lengths <- read_csv(paste0(dat_path, "/lengths_sstdusky_", YEAR, ".csv"))
spp <- read_csv(paste0("data/bsai_orox_spp_lookup.csv"))

# Get length comps ----

# Just for SSTs combine BSAI fisheries
lengths <- lengths %>% 
  filter(species_code == 30020 & source == "fishery") %>% 
  mutate(fmp_subarea = "BSAI") %>% 
  bind_rows(lengths %>% 
              filter(!(species_code == 30020 & source == "fishery"))) #%>% 
  # distinct(species_code, fmp_subarea, source)

comps <- lengths %>% 
  group_by(source, fmp_subarea, species_code, year, length) %>% 
  summarize(n = sum(frequency)) %>% 
  mutate(prop = n / sum(n)) %>% 
  left_join(spp %>% select(species_code, common_name)) %>% 
  group_by(source, fmp_subarea, species_code, year) %>% 
  mutate(N = paste0("N = ", sum(n)))

# Comp test: should all be 1!
comps %>% group_by(source, fmp_subarea, species_code, year) %>% summarize(tst = sum(prop)) %>% pull(tst) 
comps %>% filter(source == "EBS fishery") %>% 
  group_by(species_code, year) %>% summarize(tst = sum(prop)) %>% pull(tst) 

comps <- comps %>% 
  mutate(source = paste0(fmp_subarea, " ", source))

# Plot dusky ----

ggplot(data = comps %>% 
         filter(year >= 2002 & species_code == 30152 & between(length, 20, 50) & 
                  source %in% c("AI fishery", "AI survey")) %>% 
         mutate(source = factor(source, levels = c("AI survey", "AI fishery"))), 
       aes(x = length, y = factor(year), height = prop, fill = source)) +
  geom_density_ridges(stat = "identity", col = "lightgrey", alpha = 0.5, 
                      panel_scaling = TRUE, size = 0.5) +
  # geom_hline(yintercept = factor(2002:YEAR), col = "lightgrey") +
  scale_fill_manual(values = c("grey30", "#00BFC4")) +
  labs(x = "Length (cm)", y = NULL, fill = NULL) + #, title = "Dusky rockfish") +
  theme_light() +
  theme(legend.position = "top")

ggsave(paste0(out_path, "/lencomps_dusky_", YEAR, ".png"), 
       dpi=300, height=7, width=5, units="in")

# Plot SST ----

ggplot(data = comps %>% 
         filter(year >= 2002 & species_code == 30020 & between(length, 10, 70)) %>% #View() 
         filter(source %in% c("AI survey", "BSAI fishery")),
         # filter(source %in% c("AI fishery", "EBS fishery", "AI survey")), #, "EBS slope survey"
       aes(x = length, y = factor(year), height = prop, fill = source)) +
  geom_density_ridges(stat = "identity", col = "white", alpha = 0.4, #alpha = 0.8
                      panel_scaling = TRUE, size = 0.5) +
  scale_fill_manual(values = c("grey30", "#00BFC4", "#daa520", "#da2055")) +
  # geom_hline(yintercept = factor(2002:YEAR), col = "lightgrey") +
  labs(x = "Length (cm)", y = NULL, fill = NULL) + #, title = "Shortspine thornyhead") +
  theme_light() +
  # facet_wrap(~ fmp_subarea) +
  theme(legend.position = "top",
        strip.background = element_blank())

ggsave(paste0(out_path, "/lencomps_sst_", YEAR, ".png"), 
       dpi=300, height=7, width=5, units="in")
