# rema

# present models 20 and 22

# set up ----

# assessment year
YEAR <- 2022

libs <- c('readr', 'dplyr', 'tidyr', 'ggplot2', 'cowplot')
if(length(libs[which(libs %in% rownames(installed.packages()) == FALSE )]) > 0) {install.packages(libs[which(libs %in% rownames(installed.packages()) == FALSE)])}
lapply(libs, library, character.only = TRUE)

# install.packages("devtools")
# devtools::install_github("afsc-assessments/rema", dependencies = TRUE, build_vignettes = TRUE)
library(rema)

# folder set up
dat_path <- paste0("data/", YEAR); dir.create(dat_path)
out_path <- paste0("results/", YEAR); dir.create(out_path)

ggplot2::theme_set(cowplot::theme_cowplot(font_size = 13) +
                     cowplot::background_grid() +
                     cowplot::panel_border())

# new data ----

biomass_dat <- read_csv(paste0(dat_path, "/bsai_orox_biomass_", YEAR, ".csv")) %>% #distinct(area)
  mutate(area = ifelse(area %in% c('Eastern AI', 'Central AI', 'Western AI'), 'AI', area),
         species = ifelse(species_code == 30020, 'SST', 'non-SST')) %>% 
  group_by(species, year, area) %>% 
  summarize(biomass = sum(biomass),
            var = sum(var)) %>% 
  ungroup() %>% 
  mutate(cv = sqrt(var)/biomass)

sst_biomass_dat <- biomass_dat %>% 
  filter(species == 'SST') %>% 
  select(strata = area, year, biomass, cv)

nonsst_biomass_dat <- biomass_dat %>% 
  filter(species == 'non-SST') %>% 
  select(strata = area, year, biomass, cv) %>% 
  mutate(cv = ifelse(biomass == 0 & is.nan(cv), NA, cv))

nonsst_biomass_dat %>% print(n = Inf)
nonsst_biomass_dat %>% filter(biomass == 0)
nonsst_biomass_dat %>% filter(strata == 'EBS Shelf') %>% nrow()

cpue_dat <- read_csv((paste0(dat_path, '/lls_rpw_sst.csv'))) %>% 
  filter(strata == 'EBS Slope')

# Projection year ----
YEAR <- 2023

# M 20 TMB SST ----

input <- prepare_rema_input(model_name = 'Model 20',
                            biomass_dat = sst_biomass_dat,
                            end_year = YEAR,
                            PE_options = list(pointer_PE_biomass = c(1, 1, 1)))
m20_sst <- fit_rema(input)

input <- prepare_rema_input(model_name = 'Model 20',
                            biomass_dat = nonsst_biomass_dat,
                            end_year = YEAR,
                            PE_options = list(pointer_PE_biomass = c(1, 1, 1, 1)),
                            zeros = list(assumption = 'NA'))
m20_nonsst <- fit_rema(input)

# M22 add LLS in EBS slope for SST ----

input <- prepare_rema_input(model_name = 'Model 22',
                            biomass_dat = sst_biomass_dat,
                            sum_cpue_index = TRUE,
                            cpue_dat = cpue_dat,
                            end_year = YEAR,
                            multi_survey = TRUE,
                            PE_options = list(pointer_PE_biomass = c(1, 1, 1)),
                            q_options = list(pointer_biomass_cpue_strata = c(NA, 1, NA)))

m22_sst <- fit_rema(input)
m22_sst_out <- tidy_rema(m22_sst)
m22_sst_plots <- plot_rema(tidy_rema = m22_sst_out)
m22_sst_plots$cpue_by_strata

# compare M20 and M22 ----

# sst
compare_sst <- compare_rema_models(list(m20_sst,
                                        m22_sst))

p1 <- ggplot(data = compare_sst$output$biomass_by_strata,
             aes(x = year, y = pred,
                 col = model_name)) +
  geom_ribbon(aes(ymin = pred_lci, ymax = pred_uci,
                  fill = model_name), col = NA,
              alpha = 0.25) +
  geom_line(aes(lty = model_name), size = 0.7) +
  facet_wrap(~strata, nrow = NULL, scales = 'free_y') +
  geom_point(aes(x = year, y = obs), col = 'black') +
  geom_errorbar(aes(x = year, ymin = obs_lci, ymax = obs_uci), col = 'black') +
  scale_y_continuous(labels = scales::comma, expand = c(0.01, 0), limits = c(0, NA)) +
  labs(x = NULL, y = NULL, 
       title = 'Shortspine thornyhead (SST)',
       subtitle = 'Model fits to bottom trawl survey biomass (t) by region',
       fill = NULL, colour = NULL, lty = NULL, shape = NULL) +
  ggplot2::scale_fill_viridis_d(direction = 1) +
  ggplot2::scale_colour_viridis_d(direction = 1)
p1
p2 <- compare_sst$plots$total_predicted_biomass + 
  geom_line(aes(lty = model_name), size = 0.7) +
  theme(legend.position = 'none') +
  labs(lty = NULL, subtitle = 'Total predicted biomass (t)', y = NULL)
p2

compare_sst <- compare_rema_models(list(m22_sst), cpue_ylab = NULL)
p3 <- compare_sst$plots$cpue_by_strata +
  labs(subtitle = 'Model 22 fit to the NMFS longline survey RPWs',
       fill = NULL, colour = NULL, lty = NULL, shape = NULL) +
  ggplot2::scale_fill_viridis_d(direction = -1) +
  ggplot2::scale_colour_viridis_d(direction = -1) +
  theme(legend.position = 'none')

cowplot::plot_grid(p1, p3, p2, ncol = 1, 
                   rel_heights = c(0.4, 0.3, 0.3))

ggsave(paste0(out_path, '/sst_m20_m22.png'), units = 'in', bg = 'white',
       height = 8.5, width = 8, dpi = 300)

# M22 nonSST (same as M20) ----
input <- prepare_rema_input(model_name = 'Model 22',
                            biomass_dat = nonsst_biomass_dat,
                            end_year = YEAR,
                            zeros = list(assumptions = 'NA'),
                            PE_options = list(pointer_PE_biomass = c(1, 1, 1, 1)))
m22_nonsst <- fit_rema(input)
m22_nonsst_out <- tidy_rema(m22_nonsst)

compare_nonsst <- compare_rema_models(list(m20_nonsst, m22_nonsst), 
                                      biomass_ylab = 'Biomass (t)')

p1 <- ggplot(data = compare_nonsst$output$biomass_by_strata,
             aes(x = year, y = pred,
                 col = model_name)) +
  geom_ribbon(aes(ymin = pred_lci, ymax = pred_uci,
                  fill = model_name), col = NA,
              alpha = 0.25) +
  geom_line(aes(lty = model_name), size = 0.7) +
  facet_wrap(~strata, nrow = NULL, scales = 'free_y') +
  geom_point(aes(x = year, y = obs), col = 'black') +
  geom_errorbar(aes(x = year, ymin = obs_lci, ymax = obs_uci), col = 'black') +
  scale_y_continuous(labels = scales::comma, expand = c(0.01, 0), limits = c(0, NA)) +
  labs(x = NULL, y = NULL, 
       title = 'Other non-SST rockfish',
       subtitle = 'Model fits to bottom trawl survey biomass (t) by region',
       fill = NULL, colour = NULL, lty = NULL, shape = NULL) +
  ggplot2::scale_fill_viridis_d(direction = 1) +
  ggplot2::scale_colour_viridis_d(direction = 1)
p1
p2 <- compare_nonsst$plots$total_predicted_biomass + 
  geom_line(aes(lty = model_name), size = 0.7) +
  # scale_x_continuous(limits = c(1991, 2020)) +
  scale_y_continuous(labels = scales::comma, expand = c(0.01, 0)) +
  theme(legend.position = 'none') +
  labs(lty = NULL, subtitle = 'Total predicted biomass (t)', y = NULL)
p2

cowplot::plot_grid(p1, p2, ncol = 1, rel_heights = c(0.7, 0.3))

ggsave(paste0(out_path, '/nonsst_m20_m22.png'), units = 'in', bg = 'white',
       height = 8, width = 7, dpi = 300)

# write output ----

m20_tot <- tidy_rema(m20_sst)$total_predicted_biomass %>% 
  mutate(species_group = 'SST') %>% 
  bind_rows(tidy_rema(m20_nonsst)$total_predicted_biomass %>% 
              mutate(species_group = 'non-SST')) %>% 
  select(species_group, model_name, year, pred, pred_lci, pred_uci) 

m22_tot <- tidy_rema(m22_sst)$total_predicted_biomass %>% 
  mutate(species_group = 'SST') %>% 
  bind_rows(tidy_rema(m22_nonsst)$total_predicted_biomass %>% 
              mutate(species_group = 'non-SST')) %>% 
  select(species_group, model_name, year, pred, pred_lci, pred_uci) 

tot <- m20_tot %>% 
  bind_rows(m22_tot) 

tot %>%
  pivot_wider(id_cols = c(species_group, year), names_from = model_name, values_from = pred) %>% 
  bind_rows(tot %>% 
              group_by(model_name, year) %>% 
              summarise(biomass = sum(pred)) %>% 
              pivot_wider(id_cols = year, names_from = model_name, values_from = biomass) %>% 
              arrange(year) %>% 
              mutate(species_group = 'Total')) %>% 
  write_csv(paste0(out_path, '/tot_biomass.csv'))

# total abc/ofl ---

natmat <- data.frame(species_group = c('SST', 'non-SST'),
                     M = c(0.03, 0.09))
spp_abc <- tot %>% 
  filter(year == 2023) %>% 
  dplyr::select(species_group, model_name, year, biomass = pred) %>% 
  left_join(natmat) %>% 
  mutate(OFL = biomass * M,
         maxABC = biomass * (0.75 * M)) 

spp_abc %>% write_csv(paste0(out_path, '/abc_ofl_by_sppgroup.csv'))

sumtable <- spp_abc %>% 
  group_by(model_name, year) %>% 
  dplyr::summarise(biomass = sum(biomass),
            OFL = sum(OFL),
            maxABC = sum(maxABC)) 

statquo_biomass <- sumtable %>% filter(model_name == 'Model 20') %>% pull(biomass)
statquo_OFL <- sumtable %>% filter(model_name == 'Model 20') %>% pull(OFL)
statquo_maxABC <- sumtable %>% filter(model_name == 'Model 20') %>% pull(maxABC)

sumtable %>%
  mutate(biomass_percent_change = (biomass - statquo_biomass)/statquo_biomass * 100,
         OFL_maxABC_percent_change = (OFL - statquo_OFL)/statquo_OFL * 100) %>%
  write_csv(paste0(out_path, '/abc_ofl_summary.csv'))

# apportionment ----

m20_appo <- tidy_rema(m20_sst)$biomass_by_strata %>% 
  mutate(species_group = 'SST') %>% 
  bind_rows(tidy_rema(m20_nonsst)$biomass_by_strata %>% 
              mutate(species_group = 'non-SST')) %>% 
  select(species_group, model_name, strata, year, biomass = pred) 

m22_appo <- tidy_rema(m22_sst)$biomass_by_strata %>% 
  mutate(species_group = 'SST') %>% 
  bind_rows(tidy_rema(m22_nonsst)$biomass_by_strata %>% 
              mutate(species_group = 'non-SST')) %>% 
  select(species_group, model_name, strata, year, biomass = pred) 

appo <- m20_appo %>%
  bind_rows(m22_appo) %>% 
  filter(year == YEAR) %>% 
  mutate(fmp = ifelse(strata == 'AI', 'AI', 'EBS')) 

apposum <- appo  %>% 
  group_by(species_group, model_name, fmp) %>% 
  summarise(biomass = sum(biomass)) %>% 
  left_join(natmat) %>% 
  mutate(species_fmp_ABC = biomass * (0.75 * M),
         species_fmp_OFL = biomass * M) %>% 
  group_by(model_name, species_group) %>% 
  mutate(total_biomass = sum(biomass),
         prop_biomass = biomass / total_biomass) 

apposum %>% 
  write_csv(paste0(out_path, '/abc_apportionment.csv'))  

# parameter estimates ----

m20_par <- tidy_rema(m20_sst)$parameter_estimates %>% 
  mutate(species_group = 'SST') %>% 
  bind_rows(tidy_rema(m20_nonsst)$parameter_estimates %>% 
              mutate(species_group = 'non-SST')) %>% 
  select(species_group, model_name, parameter, estimate, std_err, lci, uci) 

m22_par <- tidy_rema(m22_sst)$parameter_estimates %>% 
  mutate(species_group = 'SST') %>% 
  bind_rows(tidy_rema(m22_nonsst)$parameter_estimates %>% 
              mutate(species_group = 'non-SST')) %>% 
  select(species_group, model_name, parameter, estimate, std_err, lci, uci) 

pars <- m20_par %>% 
  bind_rows(m22_par) %>% 
  arrange(species_group) %>% 
  select(species_group, model_name, parameter, estimate, std_err, lci, uci) %>% 
  write_csv(paste0(out_path, '/sst_nonsst_parameters.csv'))

pars

# LLS data summary tbl ----

read_csv((paste0(dat_path, '/lls_rpw_sst.csv'))) %>% 
  mutate(cpue2 = ifelse(is.na(cpue), NA,
                        paste0(prettyNum(cpue, big.mark = ',', digits = 1, trim = TRUE),
                               " (", prettyNum(round(cv, 2), digits = 2), ")"))) %>%
  pivot_wider(id_cols = year, names_from = strata, values_from = cpue2) %>% 
  arrange(year) %>% 
  write_csv(paste0(out_path, '/rpws.csv'))

# BTS data summary tbl ----

data.frame(year = 1982:YEAR) %>% 
  left_join(nonsst_biomass_dat %>% 
              mutate(biomass2 = ifelse(is.na(biomass), NA,
                                       paste0(prettyNum(biomass, big.mark = ',', digits = 1, trim = TRUE),
                                              " (", prettyNum(round(cv, 2), digits = 2), ")"))) %>%
              pivot_wider(id_cols = year, names_from = strata, values_from = biomass2) %>% 
              arrange(year)) %>%
  write_csv(paste0(out_path, '/nonsst_bts_biomass.csv'))

data.frame(year = 1982:YEAR) %>% 
  left_join(sst_biomass_dat %>% 
              mutate(biomass2 = ifelse(is.na(biomass), NA,
                                       paste0(prettyNum(biomass, big.mark = ',', digits = 1, trim = TRUE),
                                              " (", prettyNum(round(cv, 2), digits = 2), ")"))) %>%
              pivot_wider(id_cols = year, names_from = strata, values_from = biomass2) %>% 
              arrange(year)) %>% 
  write_csv(paste0(out_path, '/sst_bts_biomass.csv'))

# Presentation total biomass ----
ggplot(data = tot %>% 
         filter(year >= 2002 & model_name == 'Model 22') %>% 
         mutate(species_group = factor(species_group, levels = c("SST", "non-SST"), ordered = TRUE))) +
  geom_line(aes(x = year, y = pred, col = species_group),
            size = 1) +
  geom_ribbon(aes(x = year, ymin = pred_lci, ymax = pred_uci, fill = species_group),
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

catch <- read_csv(paste0(dat_path, "/bsai_orox_catch_2003_", YEAR, "_confidential.csv"))

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
apposum %>% 
  ungroup() %>%
  filter(model_name == 'Model 22') %>% 
  dplyr::select(species = species_group, FMP = fmp, `ABC*` = species_fmp_ABC,
         `OFL*` = species_fmp_OFL) %>%
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
spp_abc %>%
  filter(model_name == 'Model 22') %>% 
  dplyr::select(species = species_group, ABC = maxABC, OFL) %>%
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
  filter(species == "SST" & year < YEAR) %>% 
  ungroup() %>% 
  mutate(mean_percent = mean(percent),
         min_percent = min(percent),
         max_percent = max(percent))

tot %>% 
  group_by(year) %>% 
  mutate(total_biom = sum(pred),
         percent = 100 * pred / total_biom) %>% 
  filter(species_group == "SST" & year < YEAR & year >= 2003) %>% 
  ungroup() %>% 
  mutate(mean_percent = mean(percent),
         min_percent = min(percent),
         max_percent = max(percent))

# Total Catch/ABC/OFL
spp_abc %>%
  filter(model_name == 'Model 22') %>% 
  dplyr::select(species = species_group, ABC = maxABC, OFL) %>%
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
ratio <- apposum %>% 
  ungroup() %>%
  filter(model_name == 'Model 22') %>% 
  dplyr::select(species = species_group, FMP = fmp, F_OFL = M, ) %>% 
  left_join(catch %>%
              mutate(species = ifelse(species_name == "SST", "SST", "non-SST")) %>%
              group_by(year, FMP = fmp_subarea, species) %>%
              summarize(Catch = sum(tons))) %>% 
  left_join(m22_appo %>% 
              mutate(FMP = ifelse(strata == 'AI', 'AI', 'EBS')) %>% 
              dplyr::select(species = species_group, year, FMP, biomass) %>%
              group_by(species, year, FMP) %>% 
              dplyr::summarise(biomass = sum(biomass))) %>% 
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
  labs(x = NULL, y = "Catch/biomass", col = NULL, size = NULL, lty = NULL) +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "white"),
        panel.grid.minor = element_line(colour = "white"))

ggsave(paste0(out_path, "/ratio_catch_biomass_", YEAR, ".png"), 
       dpi=300, height=3.75, width=6, units="in")

ratio %>% 
  group_by(FMP, species) %>% 
  summarize(mean_ratio = mean(ratio))

ratio %>% filter(year == YEAR & species == "non-SST" & FMP == "AI")

ratio %>% 
  group_by(species, year) %>% 
  summarise(ratio = (sum(Catch) / sum(biomass))*100) %>% 
  # filter(year >= YEAR-3) 
  print(n=Inf)

# Presentation catch/biomass----

ratio %>% 
  mutate(FMP = factor(FMP, 
                      labels = c("Aleutian Islands", "Eastern Bering Sea"),
                      levels = c("AI", "EBS"), 
                      ordered = TRUE),
         species = factor(species, 
                          levels = c("SST", "non-SST"), 
                          ordered = TRUE)) %>% 
  ggplot(aes(x = year, y = ratio, col = species, lty = FMP)) +
  geom_line(size = 1) +
  labs(x = NULL, y = NULL, title = "Exploitation rate (catch/biomass)",
       # subtitle = "Note differences in y-axis scale",
       col = NULL, fill = NULL, lty = NULL) +
  geom_hline(data = refs, aes(yintercept = yint),
             col = "red", lty = 3) +
  scale_color_manual(values = c("#6ba292", "#f0b74a")) + #, guide = 'none') +
  theme_minimal(base_size = 13) +
  theme(#legend.position = "top",
        legend.position = c(.8, .8),
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5))

ggsave(paste0(out_path, "/pres_ratio_", YEAR, ".png"), 
       dpi=300, height=4, width=7, units="in")

# detailed non-SST ----

CUSTOMpalette <- c("#1a2ffa", "#0d177d", "#1a9ffa", "#fa751a", 
                   "#4b8e12", "#6fd21b", "#fae51a")#, "#c3b104", 
# "#f5df05", "#dcc805")

nonsst_eda <- read_csv(paste0(dat_path, "/bsai_orox_biomass_", YEAR-1, ".csv")) %>% #distinct(area)
  mutate(area = ifelse(area %in% c('Eastern AI', 'Central AI', 'Western AI'), 'AI', area)) %>% 
  filter(species_code != 30020) %>% 
  # combine dusky/dark and dusky
  mutate(common_name = ifelse(common_name == "dusky and dark rockfishes unid.",
                          "dusky and dark rockfish", common_name)) %>% 
  mutate(common_name = ifelse(common_name == "dusky rockfish",
                          "dusky and dark rockfish", common_name)) %>% 
  # combine broadfin and longspine thornyhead
  mutate(common_name = ifelse(common_name %in% c("broadfin thornyhead", "longspine thornyhead"),
                          "other thornyhead", common_name)) %>% 
  group_by(common_name, area, year) %>% 
  summarize(biomass = sum(biomass),
            var = sum(var)) %>% 
  ungroup() %>% 
  complete(common_name, area, year, fill = list(biomass = 0, var = 0)) %>% 
  mutate(area = factor(area, levels = c("AI", "EBS Shelf", "SBS", "EBS Slope"), 
                         ordered = TRUE))
         
ggplot(data = nonsst_eda, aes(x = year, y = biomass, 
                              fill = common_name)) +
  geom_bar(stat = "identity", alpha = 0.7, 
           position = position_stack(reverse = TRUE),
           size = 0.3,
           width = 1,
           colour = "black") +
  # scale_fill_viridis(discrete = TRUE, option = "D") +
  scale_fill_manual(values = CUSTOMpalette) +
  # scale_fill_grey(start = 1, end = 0) +
  facet_wrap(~area, scales = "free_y", ncol = 1) +
  # theme_minimal() +
  theme_bw() +
  scale_y_continuous(labels = scales::comma) +
  theme(legend.position = "bottom") +
  labs(x = NULL, y = "Biomass (t)", fill = "Non-SST species") +
  guides(fill = guide_legend(nrow = 3))

ggsave(paste0(out_path, "/srvbiom_nonSST_spp_", YEAR, ".png"), 
       dpi=400, height=8, width=6.5, units="in")

# SARA ----

# base weighted M on ratio of biomass estimates between SST/non-SST in terminal
# year (same as projection year)
sara_pred <- m22_tot %>% 
  filter(year >= 1991) %>% 
  pivot_wider(id_cols = year, names_from = species_group, values_from = pred) %>% 
  mutate(total = SST + `non-SST`,
         p_SST = SST / total,
         SST_M = 0.03,
         nonSST_M = 0.09,
         weighted_M = (SST_M * p_SST) + (nonSST_M * (1-p_SST))) 

write_csv(sara_pred, paste0(out_path, '/sara_biomass_pred.csv'))

weighted_M <- sara_pred %>% 
  filter(year == YEAR) %>% 
  pull(weighted_M)
  
weighted_M

# previously the SARA biomass estimates combine biomass, re-run model. this is how it was done previously, sticking
# with status quo

# based on combined 
sara_biomass_dat <- biomass_dat %>% 
  mutate(area = ifelse(area == 'SBS', 'AI', area)) %>% 
  group_by(year, area) %>% 
  dplyr::summarise(biomass = sum(biomass),
                   var = sum(var)) %>% 
  mutate(cv = sqrt(var) / biomass) %>% 
  select(-var)

# sara_biomass_dat <- sara_biomass_dat %>% 
#   filter(year >= 1991) %>% 
#   pivot_wider(id_cols = year, names_from = area, values_from = biomass,
#               values_fill = NA)
# 
# cpue_dat %>% 
#   mutate(area = 'LLS') %>% 
#   select(-strata, -cv, biomass = cpue)

sara_biomass_dat <- sara_biomass_dat %>% 
  filter(year >= 1991) %>% 
  pivot_wider(id_cols = area, names_from = year, values_from = biomass,
              values_fill = NA)

write_csv(sara_biomass_dat, paste0(out_path, '/sara_survey_biomass_data.csv'))

sara_exploitation <- catch %>% 
  group_by(year) %>% 
  dplyr::summarise(catch = sum(tons)) %>% 
  left_join(sara_pred %>% 
              select(year, biomass = total)) %>% 
  mutate(exploitation = catch / biomass)
  
write_csv(sara_exploitation, paste0(out_path, '/sara_exploitation_rates.csv'))


