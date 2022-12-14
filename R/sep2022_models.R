# Bridge from ADMB to TMB

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
  mutate(cv = sqrt(var)/biomass) %>% 
  filter(year >= 1982)

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

# Change after September -----

YEAR <- 2020

# M 20 and 20.a sst ----

admb_sst <- read_admb_re(filename = 'results/2020/RE_SST/RE_SST_rwout.rep',
                        model_name = 'Model 20',
                        biomass_strata_names = c('AI', 'EBS Slope', 'SBS'))

input <- prepare_rema_input(model_name = 'Model 20.a',
                            admb_re = admb_sst,
                            # shared PE
                            PE_options = list(pointer_PE_biomass = c(1, 1, 1)))
mod_sst <- fit_rema(input)
out_sst <- tidy_rema(mod_sst)

compare_sst <- compare_rema_models(list(mod_sst), admb_re = admb_sst)

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
  # scale_x_continuous(limits = c(1991, 2020)) +
  # scale_y_continuous(limits = c(0, 36000), labels = scales::comma, expand = c(0.01, 0)) +
  theme(legend.position = 'none') +
  labs(lty = NULL, subtitle = 'Total predicted biomass (t)', y = NULL)
# p2

cowplot::plot_grid(p1, p2, ncol = 1, rel_heights = c(0.55, 0.45))

ggsave(paste0(out_path, '/sst_m20_m20a.png'), units = 'in', bg = 'white',
       height = 6, width = 7.5, dpi = 300)

compare_sst$output$biomass_by_strata %>% 
  pivot_wider(id_cols = c('strata', 'year'), names_from = model_name, values_from = pred) %>% 
  mutate(percent_difference = abs((`Model 20.a` - `Model 20`) / `Model 20` * 100)) %>% 
  pull(percent_difference) %>% 
  summary()

compare_sst$output$total_predicted_biomass %>%
  pivot_wider(id_cols = c('year'), names_from = model_name, values_from = pred) %>% 
  mutate(percent_difference = abs((`Model 20.a` - `Model 20`) / `Model 20` * 100)) %>% 
  pull(percent_difference) %>% 
  summary()

# M 20 and 20.a nonsst ----

admb_nonsst <- read_admb_re(filename = 'results/2020/RE_non-sst/RE_non-SST_rwout.rep',
                        model_name = 'Model 20',
                        biomass_strata_names = c('AI', 'EBS Shelf', 'EBS Slope', 'SBS'))

input <- prepare_rema_input(model_name = 'Model 20.a',
                            admb_re = admb_nonsst,
                            zeros = list(assumptions = 'NA'),
                            PE_options = list(pointer_PE_biomass = c(1, 1, 1, 1)))
mod_nonsst <- fit_rema(input)
out_nonsst <- tidy_rema(mod_nonsst)

compare_nonsst <- compare_rema_models(list(mod_nonsst), 
                                      admb_re = admb_nonsst,
                                      biomass_ylab = 'Biomass (t)')

p1 <- ggplot(data = compare_nonsst$output$biomass_by_strata %>% 
         filter(year >= 1990),
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

p2 <- compare_nonsst$plots$total_predicted_biomass + 
  geom_line(aes(lty = model_name), size = 0.7) +
  scale_x_continuous(limits = c(1991, 2020)) +
  scale_y_continuous(limits = c(0, 36000), labels = scales::comma, expand = c(0.01, 0)) +
  theme(legend.position = 'none') +
  labs(lty = NULL, subtitle = 'Total predicted biomass (t)', y = NULL)
# p2

cowplot::plot_grid(p1, p2, ncol = 1, rel_heights = c(0.7, 0.3))

ggsave(paste0(out_path, '/nonsst_m20_m20a.png'), units = 'in', bg = 'white',
       height = 8, width = 7, dpi = 300)

compare_nonsst$output$biomass_by_strata %>% 
  pivot_wider(id_cols = c('strata', 'year'), names_from = model_name, values_from = pred) %>% 
  mutate(percent_difference = abs((`Model 20.a` - `Model 20`) / `Model 20` * 100)) %>% 
  pull(percent_difference) %>% 
  summary()

compare_nonsst$output$total_predicted_biomass %>%
  pivot_wider(id_cols = c('year'), names_from = model_name, values_from = pred) %>% 
  mutate(percent_difference = ((`Model 20.a` - `Model 20`) / `Model 20` * 100)) %>% 
  pull(percent_difference) %>% 
  summary()

# SST on EBS Slope ----

compare_nonsst$output$total_predicted_biomass %>% 
  bind_rows(compare_sst$output$total_predicted_biomass) %>% 
  filter(model_name == 'Model 20') %>% 
  group_by(year) %>% 
  summarize(total_biom = sum(pred)) %>% 
  left_join(compare_sst$output$biomass_by_strata %>% 
              filter(model_name == 'Model 20' & strata == 'EBS Slope') %>% 
              select(year, sst_ebsslope = pred)) %>% 
  mutate(p = sst_ebsslope / total_biom) %>% print(n=Inf)
  
# M 20.a TMB + update data ----

input <- prepare_rema_input(model_name = 'Model 20.a',
                            biomass_dat = sst_biomass_dat,
                            end_year = YEAR,
                            PE_options = list(pointer_PE_biomass = c(1, 1, 1)))
m20a_sst <- fit_rema(input)

input <- prepare_rema_input(model_name = 'Model 20.a',
                            biomass_dat = nonsst_biomass_dat %>% 
                              # remove the updated data based on recommendation from PH
                              filter(!c(strata == 'EBS Shelf' & year == 2021)),
                            end_year = YEAR,
                            PE_options = list(pointer_PE_biomass = c(1, 1, 1, 1)),
                            zeros = list(assumption = 'NA'))
m20a_nonsst <- fit_rema(input)

# M22 (22.1) add LLS in EBS slope ----

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

# compare M20b and M22 ----

# sst
compare_sst <- compare_rema_models(list(m20a_sst,
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
# p2

compare_sst <- compare_rema_models(list(m22_sst), cpue_ylab = NULL)
p3 <- compare_sst$plots$cpue_by_strata +
  labs(subtitle = 'Model 22 fit to the NMFS longline survey RPWs',
       fill = NULL, colour = NULL, lty = NULL, shape = NULL) +
  ggplot2::scale_fill_viridis_d(direction = -1) +
  ggplot2::scale_colour_viridis_d(direction = -1) +
  theme(legend.position = 'none')

cowplot::plot_grid(p1, p3, p2, ncol = 1, 
                   rel_heights = c(0.4, 0.3, 0.3))

ggsave(paste0(out_path, '/sst_m20a_m22.png'), units = 'in', bg = 'white',
       height = 8.5, width = 8, dpi = 300)

# M22 nonSST (same as M20a) ----
input <- prepare_rema_input(model_name = 'Model 22',
                            biomass_dat = nonsst_biomass_dat %>% 
                              filter(!c(strata == 'EBS Shelf' & year == 2021)),
                            end_year = YEAR,
                            zeros = list(assumptions = 'NA'),
                            PE_options = list(pointer_PE_biomass = c(1, 1, 1, 1)))
m22_nonsst <- fit_rema(input)
m22_nonsst_out <- tidy_rema(m22_nonsst)

# write output ----

m20_tot <- admb_sst$admb_re_results$total_predicted_biomass %>% 
  mutate(species_group = 'SST') %>% 
  bind_rows(admb_nonsst$admb_re_results$total_predicted_biomass %>% 
              mutate(species_group = 'non-SST')) %>% 
  select(species_group, model_name, year, pred, pred_lci, pred_uci) 
  
m20a_tot <- tidy_rema(m20a_sst)$total_predicted_biomass %>% 
  mutate(species_group = 'SST') %>% 
  bind_rows(tidy_rema(m20a_nonsst)$total_predicted_biomass %>% 
              mutate(species_group = 'non-SST')) %>% 
  select(species_group, model_name, year, pred, pred_lci, pred_uci) 

m22_tot <- tidy_rema(m22_sst)$total_predicted_biomass %>% 
  mutate(species_group = 'SST') %>% 
  bind_rows(tidy_rema(m22_nonsst)$total_predicted_biomass %>% 
              mutate(species_group = 'non-SST')) %>% 
  select(species_group, model_name, year, pred, pred_lci, pred_uci) 

tot <- m20_tot %>% 
  bind_rows(m20a_tot) %>% 
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
sumtable <- tot %>% 
  filter(year == 2020) %>% 
  select(species_group, model_name, year, biomass = pred) %>% 
  left_join(natmat) %>% 
  mutate(OFL = biomass * M,
         maxABC = biomass * (0.75 * M)) %>% 
  group_by(model_name, year) %>% 
  summarise(biomass = sum(biomass),
            OFL = sum(OFL),
            maxABC = sum(maxABC)) %>% 
  # rounding error in 2020 biomass, just fix value to the value in the table for
  # clarity
  mutate(biomass = ifelse(model_name == 'Model 20', 53248, biomass))

statquo_biomass <- sumtable %>% filter(model_name == 'Model 20') %>% pull(biomass)
statquo_OFL <- sumtable %>% filter(model_name == 'Model 20') %>% pull(OFL)
statquo_maxABC <- sumtable %>% filter(model_name == 'Model 20') %>% pull(maxABC)

sumtable %>%
  mutate(biomass_percent_change = (biomass - statquo_biomass)/statquo_biomass * 100,
         OFL_maxABC_percent_change = (OFL - statquo_OFL)/statquo_OFL * 100) %>%
  write_csv(paste0(out_path, '/abc_ofl_summary.csv'))

# apportionment ----

m20_appo <- admb_sst$admb_re_results$biomass_by_strata %>% 
  mutate(species_group = 'SST') %>% 
  bind_rows(admb_nonsst$admb_re_results$biomass_by_strata %>% 
              mutate(species_group = 'non-SST')) %>% 
  select(species_group, model_name, strata, year, biomass = pred) 

m20a_appo <- tidy_rema(m20a_sst)$biomass_by_strata %>% 
  mutate(species_group = 'SST') %>% 
  bind_rows(tidy_rema(m20a_nonsst)$biomass_by_strata %>% 
              mutate(species_group = 'non-SST')) %>% 
  select(species_group, model_name, strata, year, biomass = pred) 

m22_appo <- tidy_rema(m22_sst)$biomass_by_strata %>% 
  mutate(species_group = 'SST') %>% 
  bind_rows(tidy_rema(m22_nonsst)$biomass_by_strata %>% 
              mutate(species_group = 'non-SST')) %>% 
  select(species_group, model_name, strata, year, biomass = pred) 

appo <- m20_appo %>%
  bind_rows(m20a_appo) %>% 
  bind_rows(m22_appo) %>% 
  filter(year == YEAR) %>% 
  mutate(fmp = ifelse(strata == 'AI', 'AI', 'EBS')) 
  
appo  %>% 
  group_by(species_group, model_name, fmp) %>% 
    summarise(biomass = sum(biomass)) %>% 
    # group_by(species_group, model_name) %>% 
    # mutate(fmp_biomass = sum(biomass)) %>% 
  left_join(natmat) %>% 
  mutate(species_fmp_ABC = biomass * (0.75 * M)) %>% 
  group_by(model_name, fmp) %>% 
  summarise(fmp_ABC = sum(species_fmp_ABC)) %>% 
  group_by(model_name) %>% 
  mutate(ABC = sum(fmp_ABC),
         perc_ABC = fmp_ABC / ABC * 100) %>% 
  pivot_wider(id_cols = c(model_name), names_from = fmp, values_from = perc_ABC) %>% 
  write_csv(paste0(out_path, '/abc_apportionment.csv'))

# parameter estimates ----

m20a_par <- tidy_rema(m20a_sst)$parameter_estimates %>% 
  mutate(species_group = 'SST') %>% 
  bind_rows(tidy_rema(m20a_nonsst)$parameter_estimates %>% 
              mutate(species_group = 'non-SST')) %>% 
  select(species_group, model_name, parameter, estimate, std_err, lci, uci) 

m22_par <- tidy_rema(m22_sst)$parameter_estimates %>% 
  mutate(species_group = 'SST') %>% 
  bind_rows(tidy_rema(m22_nonsst)$parameter_estimates %>% 
              mutate(species_group = 'non-SST')) %>% 
  select(species_group, model_name, parameter, estimate, std_err, lci, uci) 

# Model 20 parameters 

# from .std files. note i needed to rerun the admb model because i had
# accidentally deleted the std
# sst logSdLam -2.1095e+00 2.3650e-01
# nonsst logSdLam -3.5468e-01 1.8466e-01

m20_par <- data.frame(model_name = 'Model 20',
           parameter = c('process error'),
           est = c(-2.1095e+00),
           se = c(2.3650e-01),
           species_group = 'SST') %>%
  bind_rows(data.frame(model_name = 'Model 20',
                       parameter = c('process error'),
                       est = c(-3.5468e-01),
                       se = c(1.8466e-01),
                       species_group = 'non-SST')) %>% 
  mutate(estimate = exp(est),
         std_err = exp(est) * se,
         lci = exp(est - qnorm(1 - 0.05/2) * se),
         uci = exp(est + qnorm(1 - 0.05/2) * se)) %>%
  select(-est, - se) 

pars <- m20_par %>% 
  bind_rows(m20a_par) %>% 
  bind_rows(m22_par) %>% 
  arrange(species_group) %>% 
  select(species_group, model_name, parameter, estimate, std_err, lci, uci) %>% 
  write_csv(paste0(out_path, '/sst_nonsst_parameters.csv'))

pars

# data summary tbl ----
read_csv((paste0(dat_path, '/lls_rpw_sst.csv'))) %>% 
  mutate(cpue2 = ifelse(is.na(cpue), NA,
                      paste0(prettyNum(cpue, big.mark = ',', digits = 1, trim = TRUE),
                             " (", prettyNum(round(cv, 2), digits = 2), ")"))) %>%
  pivot_wider(id_cols = year, names_from = strata, values_from = cpue2) %>% 
  arrange(year) %>% 
  write_csv(paste0(out_path, '/rpws.csv'))

# correlation between biomass and RPW ----

df <- expand.grid(year = min(cpue_dat$year):2021) %>% 
  left_join(cpue_dat %>% 
              select(-cv, -strata)) %>% 
  left_join(sst_biomass_dat %>% 
              filter(strata == 'EBS Slope') %>% 
              select(-cv, -strata))
cor.test(df$biomass, df$cpue)

df %>% 
  summarize(cor(biomass, cpue))

# STOP ----

# everything below was exploratory and not in the report.

# M 22.2a extra biomass obs error ----

input <- prepare_rema_input(model_name = 'Model 22.2.a',
                            biomass_dat = sst_biomass_dat,
                            sum_cpue_index = TRUE,
                            cpue_dat = cpue_dat,
                            end_year = YEAR,
                            multi_survey = TRUE,
                            PE_options = list(pointer_PE_biomass = c(1, 1, 1)),
                            q_options = list(pointer_biomass_cpue_strata = c(NA, 1, NA)),
                            extra_biomass_cv = list(assumption = 'extra_cv'))
m22.2a_sst <- fit_rema(input)

input <- prepare_rema_input(model_name = 'Model 22.2.a',
                            biomass_dat = nonsst_biomass_dat,
                            end_year = YEAR,
                            PE_options = list(pointer_PE_biomass = c(1, 1, 1, 1)),
                            zeros = list(assumption = 'NA'),
                            extra_biomass_cv = list(assumption = 'extra_cv'))
m22.2a_nonsst <- fit_rema(input)

# compare 1 ----

compare_sst <- compare_rema_models(list(m20b_sst, m22a_sst))
compare_sst$plots$total_predicted_biomass
compare_sst$plots$biomass_by_strata +
  facet_wrap(~strata, scales = 'free_y')

compare_sst$output$parameter_estimates # mod 22a NOT converged for SST (CI hitting bounds)

compare_nonsst <- compare_rema_models(list(m20b_nonsst, m22a_nonsst))
compare_nonsst$aic
compare_nonsst$plots$total_predicted_biomass
compare_nonsst$plots$biomass_by_strata +
  facet_wrap(~strata, scales = 'free_y')
compare_nonsst$output$parameter_estimates # mod 22a NOT converged for non-SST (est/CI hitting bounds)

# M 22b extra cpue obs error ----

input <- prepare_rema_input(model_name = 'Model 22.2.b',
                            biomass_dat = sst_biomass_dat,
                            sum_cpue_index = TRUE,
                            cpue_dat = cpue_dat,
                            end_year = YEAR,
                            multi_survey = TRUE,
                            PE_options = list(pointer_PE_biomass = c(1, 1, 1)),
                            q_options = list(pointer_biomass_cpue_strata = c(NA, 1, NA)),
                            extra_cpue_cv = list(assumption = 'extra_cv'))
m22.2b_sst <- fit_rema(input)

# M 22c extra biomass + cpue obs error ----

input <- prepare_rema_input(model_name = 'Model 22.2.c',
                            biomass_dat = sst_biomass_dat,
                            sum_cpue_index = TRUE,
                            cpue_dat = cpue_dat,
                            end_year = YEAR,
                            multi_survey = TRUE,
                            PE_options = list(pointer_PE_biomass = c(1, 1, 1)),
                            q_options = list(pointer_biomass_cpue_strata = c(NA, 1, NA)),
                            extra_biomass_cv = list(assumption = 'extra_cv'),
                            extra_cpue_cv = list(assumption = 'extra_cv'))
m22.2c_sst <- fit_rema(input)

compare_sst <- compare_rema_models(list(m22a_sst,
                                        m22.2a_sst, m22.2b_sst,
                                        m22.2c_sst))
compare_sst$aic
compare_sst$plots$total_predicted_biomass
compare_sst$plots$biomass_by_strata
compare_sst$output$parameter_estimates # mod 22.2.a and 22.2.c NOT converged for SST (CI hitting bounds)

compare_sst <- compare_rema_models(list(m22a_sst,m22.2b_sst))
compare_sst$aic
compare_sst$plots$total_predicted_biomass
compare_sst$plots$biomass_by_strata
compare_sst$output$parameter_estimates 

# SST results summary 

# Alternative models that estimated additional observation error (i.e., using
# eqn 7 in the main document) were explored during model development. However,
# these models either exhibited convergence issues (e.g., hitting parameter
# bounds) or were not statistically supported via AIC. We therefore do not
# present results for these more complex models.

# compare M20b and M22a ----
compare_sst <- compare_rema_models(list(m20b_sst, m22a_sst))
compare_sst$aic 
