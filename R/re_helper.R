
# Random effects model helper functions


# Generalized function to run random effects model. Using input data with
# defined "species" and "region" values, this fxn creates a new admb
# subdirectory for each "species", write an admb dat file, runs the model, and
# summarizes the results and convergence criteria.

run_re_model <- function(data, 
                         # one overall process error or multiple?
                         n_logsdlam = c("single", "multiple")) {
  
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
  logsdlam <- list() # process error estimates
  
  for(i in 1:length(spp)) {
    
    name <- paste0("RE_", spp[i]) # will be tpl name
    areas <- model_lkup %>% filter(species == spp[i]) %>% pull(area) # areas in model
    if(n_logsdlam == "single") {
      error <- 1
      pe <- rep(1, length(areas)) } else {
        error <- 1:length(areas)
        pe <- length(areas)
      }
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
        "# Number or process error parameters","\n", pe,"\n",
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
    
    # process error estimates
    fit <- read_fit(tplfile = name)
    
    if (n_logsdlam == "single") {
      logsdlam[[i]] <- id %>% 
        distinct(species) %>% 
        mutate(logsdlam = fit$est[which(fit$names %in% "logSdLam")],
               logsdlam_se = fit$std[which(fit$names %in% "logSdLam")])
    } else {
    logsdlam[[i]] <- id %>% 
      left_join(model_lkup) %>% 
      distinct(species, area) %>% 
      mutate(logsdlam = fit$est[which(fit$names %in% "logSdLam")],
             logsdlam_se = fit$std[which(fit$names %in% "logSdLam")])
    }
  }
  
  re_output <- do.call(rbind, re_output)
  area_sd <- do.call(rbind, area_sd)
  diagnostics <- do.call(rbind, diagnostics)
  logsdlam <- do.call(rbind, logsdlam)
  
  print(diagnostics)
  
  setwd(root)
  
  out <- list(re_output = re_output, area_sd = area_sd, diagnostics = diagnostics,
              logsdlam = logsdlam)
  return(out)
}

# single species version ----

run_sre_model <- function(data,
                          # use all years (same as multi-survey vs)?
                          all_years = TRUE) {
    
  dir.create(paste0(out_path, "/single_survey_RE"))
  
  setwd(root)
  # data = data_sst_nonSST; all_years = TRUE
  input_biom <- data 
  syr <- min(input_biom$year)
  spp <- unique(input_biom$species)

  diagnostics <- list() # for model diagnostics (convergence, mgc)
  re_output <- list() # random effect biomass estimate output
  logsdlam <- list() # process error
  
  # i=1;j=1
  for(i in 1:length(spp)) {
    
    # SST is only in 3 surveys, non-SST in 4
    regions <- input_biom %>% filter(species == spp[i]) %>% distinct(region) %>% pull()
    
    for(j in 1:length(regions)) {
      
      name <- paste0("RE_", spp[i], "_", regions[j]) # will be tpl name
      
      # Create new directory and copy over tpl and exe files. Rename them to
      # specific species/complex
      setwd(root)
      spp_out <- file.path(paste0(out_path, "/single_survey_RE/RE_", spp[i], "_", regions[j]))
      dir.create(spp_out)
      
      file.copy(from = stpl, to = spp_out, overwrite = TRUE)
      file.copy(from = sexe, to = spp_out, overwrite = TRUE)
      setwd(spp_out)
      # file.rename("RES.tpl", paste0(name, ".tpl")) # file on github
      # file.rename("RES.exe", paste0(name, ".exe"))
      file.rename("re_test.tpl", paste0(name, ".tpl")) # file cole emailed me 20210611
      file.rename("re_test.exe", paste0(name, ".exe"))
      
      # prep for dat file
      sppdat <- input_biom %>% 
        filter(species == spp[i] & region == regions[j]) %>% 
        # zeros are removed, treated as NAs
        filter(biomass > 0) %>% 
        arrange(year)
      
      # fill in missing years with NAs
      if(isTRUE(all_years)) {
        sppdat <- tibble(year = syr:YEAR) %>% 
          left_join(sppdat) 
      }
      
      nobs <- length(unique(sppdat$year))
      yrs <- sppdat %>% pull(year)
      srv_est <- sppdat %>% pull(biomass)
      srv_est[is.na(srv_est)] <- "-9" # ADMB flag for NA or zero data
      srv_cv <- sppdat %>% pull(cv)
      srv_cv[is.na(srv_cv)] <- "-9" # ADMB flag for NA or zero data
      
      # Write dat file
      cat("# Single survey RE for", spp[i], regions[j], "\n",
          "# Model start and end years","\n", min(yrs), max(yrs), "\n",
          "# nobs","\n", nobs, "\n",
          "# Years","\n", yrs,"\n",
          "# Biomass estimates","\n", srv_est, "\n",
          "# Biomass CV","\n", srv_cv, "\n", sep = " ", file = paste0(name, ".dat"))
      
      # remove pre-existing std file if it exists and run admb
      if (file.exists(paste0(name, ".std"))) {file.remove(paste0(name, ".std"), showWarnings = FALSE)}
      run_admb(name, verbose = TRUE)
      file.rename("rwout.rep", paste0(name, "_rwout.rep"))
      
      # full admb output
      df <- read_admb(tplfile = name, repfile = paste0(name, "_rwout"))
      
      # dataframe with species complex and year
      id <- tibble(species = spp[i], 
                   region = regions[j],
                   year = df$yrs)

      res <- id %>% 
        mutate(re_est = df$biomA,
               re_lci = df$LCI,
               re_uci = df$UCI)
      
      # process error estimates
      fit <- read_fit(tplfile = name)
      param <- id %>% 
        distinct(species, region) %>% 
        mutate(logsdlam = fit$est[which(fit$names %in% "logSdLam")],
               logsdlam_se = fit$std[which(fit$names %in% "logSdLam")])
      
      
      # diagnostics - just mgc and test for positive definite hessian
      diag <- tibble(species = spp[i],
                     region = regions[j],
                     logDetHess = df$fit$logDetHess,
                     mgc = df$fit$maxgrad) %>% 
        mutate(converged = mgc < 0.001 & logDetHess > 0)
      
      if(j == 1) {
        res_f <- res
        param_f <- param  
        diag_f <- diag
      } else {
        res_f <- bind_rows(res, res_f)
        param_f <- bind_rows(param, param_f)
        diag_f <- bind_rows(diag, diag_f)
      } 
    }
    re_output[[i]] <- res_f
    diagnostics[[i]] <- diag_f
    logsdlam[[i]] <- param_f
    rm(res_f, param_f, diag_f)
    
  }
  
  re_output <- do.call(rbind, re_output)
  logsdlam <- do.call(rbind, logsdlam)
  diagnostics <- do.call(rbind, diagnostics)
  
  setwd(root)
  
  out <- list(re_output = re_output, logsdlam = logsdlam, diagnostics = diagnostics)
  return(out)
}
