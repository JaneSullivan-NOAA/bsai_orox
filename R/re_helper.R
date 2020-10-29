
# Random effects model helper functions


# Generalized function to run random effects model. Using input data with
# defined "species" and "region" values, this fxn creates a new admb
# subdirectory for each "species", write an admb dat file, runs the model, and
# summarizes the results and convergence criteria.

run_re_model <- function(data) {
  
  setwd(root)
  
  # expand to include all years, species, regions for input to RE model
  input_biom <- data %>% 
    expand(year = min(data$year):YEAR, species, region) %>% 
    left_join(data)
  
  input_biom <- input_biom %>% 
    pivot_wider(id_cols = c("species", "year"), names_from = "region", 
                values_from = "biomass", values_fill = 0) %>% 
    mutate(value = "mu") %>% 
    bind_rows(input_biom %>% 
                pivot_wider(id_cols = c("species", "year"), names_from = "region", 
                            values_from = "cv", values_fill = 0) %>% 
                mutate(value = "cv")) %>% 
    arrange(species, value, year) %>% 
    replace(is.na(.), 0)
  
  # Create a look up table of species x area combinations to include in the RE
  # model runs (remove areas where species are not sampled, e.g. SST in EBS_SHELF)
  model_lkup <- input_biom %>% 
    filter(value == "mu") %>% 
    pivot_longer(cols = unique(data$region), names_to = "area", values_to = "biomass") %>% 
    group_by(species, area) %>% 
    summarize(tot_biom = sum(biomass)) %>% 
    ungroup() %>% 
    filter(tot_biom > 0) %>% 
    select(species, area)
  
  # Run model ----
  
  spp <- unique(model_lkup$species)
  yrs <- unique(input_biom$year)
  process_error <- 1 # TODO: talk with PH about this assumption and alternative options
  
  diagnostics <- list() # for model diagnostics (convergence, mgc)
  re_output <- list() # random effect biomass estimate output
  area_sd <- list() # sd of the estimated biomass sd generated, property of sdreport obj
  
  for(i in 1:length(spp)) {
    
    name <- paste0("RE_", spp[i]) # will be tpl name
    areas <- model_lkup %>% filter(species == spp[i]) %>% pull(area) # areas in model
    error <- rep(process_error, length(areas))
    
    # Create new directory and copy over tpl and exe files. Rename them to
    # specific species/complex
    setwd(root)
    spp_out <- file.path(paste0(out_path, "/RE_", spp[i]))
    dir.create(spp_out)
    file.copy(from = tpl, to = spp_out, overwrite = TRUE)
    file.copy(from = exe, to = spp_out, overwrite = TRUE)
    setwd(spp_out)
    file.rename("RE.tpl", paste0(name, ".tpl"))
    file.rename("RE.exe", paste0(name, ".exe"))
    
    # Biomass estimates
    srv_est <- input_biom %>% 
      filter(species == spp[i] & value == "mu") %>% 
      select(-c(species, year, value)) %>% 
      select(areas) %>% 
      as.matrix()
    srv_est <- unname(srv_est) # gets rid of column names
    srv_est[srv_est == 0] <- "-9" # ADMB flag for NA or zero data
    
    # Biomass variability
    srv_cv <- input_biom %>% 
      filter(species == spp[i] & value == "cv") %>% 
      select(-c(species, year, value)) %>% 
      select(areas) %>% 
      as.matrix()
    srv_cv <- unname(srv_cv) # gets rid of column names
    srv_cv[srv_cv == 0] <- "-9" # ADMB flag for NA or zero data
    
    # Write dat file
    cat("# Model start and end years","\n", min(yrs), max(yrs), "\n",
        "# Number of survey indices fit (i.e., regions/depth strata)","\n", length(areas), "\n",
        "# Number or process error parameters","\n",1,"\n",
        "# Process error index","\n", error, "\n",
        "# Number of surveys","\n", length(yrs), "\n",
        "# Survey years","\n", yrs, "\n",
        "#", areas, "\n", 
        "# Survey biomass","\n", sep = " ", file = paste0(name, ".dat"))
    
    write.table(srv_est, file = paste0(name, ".dat"), sep = " ", append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE)
    write.table(paste0("# Survey biomass SE"), file = paste0(name, ".dat"), sep = " ", append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE)
    write.table(srv_cv, file = paste0(name, ".dat"), sep = " ", append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE)
    
    # remove pre-existing std file if it exists and run admb
    if (file.exists(paste0(name, ".std"))) {file.remove(paste0(name, ".std"), showWarnings = FALSE)}
    run_admb(name, verbose = TRUE)
    file.rename("rwout.rep", paste0(name, "_rwout.rep"))
    
    # full admb output
    df <- read_admb(tplfile = name, repfile = paste0(name, "_rwout"))
    
    # dataframe with species complex and year
    id <- tibble(species = spp[i], 
                 year = df$yrs)
    
    # total biomass + 95% ci
    re_output[[i]] <- id %>% 
      mutate(region = "TOTAL",
             re_est = df$biom_TOT,
             re_lci = df$biom_TOT_LCI,
             re_uci = df$biom_TOT_UCI) %>% 
      # area-specific biomass estimates + 95% ci
      bind_rows(id %>% 
                  bind_cols(setNames(as_tibble(df$biomA), areas)) %>% # colnames = areas where survey occurred
                  pivot_longer(cols = areas, names_to = "region", values_to = "re_est") %>% 
                  left_join(id %>% 
                              bind_cols(setNames(as_tibble(df$LCI), areas)) %>%  
                              pivot_longer(cols = areas, names_to = "region", values_to = "re_lci")) %>% 
                  left_join(id %>% 
                              bind_cols(setNames(as_tibble(df$UCI), areas)) %>%  
                              pivot_longer(cols = areas, names_to = "region", values_to = "re_uci")))
    
    # Area-specific sd of biomass sd - used to calculate BS or AI wide RE estimates
    area_sd[[i]] <- id %>% 
      bind_cols(setNames(as_tibble(df$biomsd), areas)) %>% 
      pivot_longer(cols = areas, names_to = "region", values_to = "biomsd") %>% 
      left_join(id %>% 
                  bind_cols(setNames(as_tibble(df$biomsd.sd), areas)) %>% 
                  pivot_longer(cols = areas, names_to = "region", values_to = "biomsd_sd"))
    
    # diagnostics - just mgc and test for positive definite hessian
    diagnostics[[i]] <- tibble(species = spp[i],
                               logDetHess = df$fit$logDetHess,
                               mgc = df$fit$maxgrad) %>% 
      mutate(converged = mgc < 0.001 & logDetHess > 0)
  }
  
  re_output <- do.call(rbind, re_output)
  area_sd <- do.call(rbind, area_sd)
  diagnostics <- do.call(rbind, diagnostics)
  
  print(diagnostics)
  
  setwd(root)
  
  out <- list(re_output = re_output, area_sd = area_sd, diagnostics = diagnostics)
  return(out)
}
