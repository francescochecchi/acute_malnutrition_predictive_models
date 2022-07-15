#..........................................................................................
### +++++++++++++++ SMALL-AREA PREDICTION OF ACUTE MALNUTRITION BURDEN ++++++++++++++++ ###
#..........................................................................................

#..........................................................................................
## ----- R CODE TO EVALUATE REGRESSION MODELS TO PREDICT ACUTE MALNUTRITION BURDEN ----- ##
#..........................................................................................

                                       # Written by Francesco Checchi, LSHTM (November 2020)
                                       # francesco.checchi@lshtm.ac.uk 



#.........................................................................................
### Reading in required files
#.........................................................................................
    
  #...................................
  ## Model data (from earlier script)
  ch_obs <- read.csv(paste(country, "_ch_obs_model.csv", sep = ""), sep = "," )


#.........................................................................................
### Prepare and explore the dependent variable
#.........................................................................................

  #...................................
  ## Plot prevalences if binary
  if (y_hat %in% c("samz", "gamz", "samm", "gamm")) {
    x1 <- aggregate(ch_obs[, y_hat], by = ch_obs[, c(part_unit, "tm")], FUN = mean, na.rm = TRUE)
    colnames(x1)[1] <- paste(y_hat, "_prevalence", sep = "")
    f_hist(paste(y_hat, "_prevalence", sep = ""), x1, c(NA, NA))}

  #...................................
  ## Plot means if continuous
  if (! y_hat %in% c("samz", "gamz", "samm", "gamm")) { f_hist(y_hat, ch_obs, c(NA, NA)) }

  #...................................
  ## Partition data randomly into "training" sample and "hold-out" sample for validation
    
    # If the partitioning method is to split the data by period...
    if (part_method == "period") {
      # figure out the time point at which the desired fraction of the data is reached
      x1 <- ch_obs
      x1[, "unit"] <- 1
      x1 <- as.data.frame(table(x1[, c("tm", "unit")]) )
      x1[, "cum"] <- cumsum(x1$Freq) / sum(x1$Freq)
      x2 <- x1[findInterval(1 - part_fraction, x1[, "cum"]) , "tm"]
      x2 <- x1[which(x1[, "tm"] == x2), "tm"]
      
      # label observations in origin dataset as holdout or training
      ch_obs[, "holdout"] <- FALSE
      ch_obs[, "holdout"] <- ifelse(ch_obs[, "tm"] > as.integer(as.character(x2)), TRUE, FALSE)
          
    }
  
    # If the partitioning method is to select a sample of the data across its entire time span...
    if (part_method == "sample") {
      # select holdout units at random without replacement
      holdout_units <- sample(unique(ch_obs[, part_unit]), 
        round(part_fraction * length(unique(ch_obs[, part_unit])), digits = 0), replace = FALSE )
        
      # label observations in origin dataset as holdout or training
      ch_obs[, "holdout"] <- FALSE
      ch_obs[ch_obs[, part_unit] %in% holdout_units, "holdout"] <- TRUE
        
    }

   # If no partitioning is to occur, and instead cross-validation is to be used...
    if (part_method == "cv") {

      # label observations in origin dataset as holdout or training
      ch_obs[, "holdout"] <- FALSE
        
    }
  
  
  #...................................
  ## Compare the average of the malnutrition index within the holdout and training datasets
  if (part_method %in% c("period", "sample") ) {  
    x1 <- aggregate(ch_obs[, y_hat], by = list(ch_obs[, "holdout"]), FUN = mean, na.rm = TRUE)
    colnames(x1) <- c("holdout", y_hat)
    print(x1)
  }

   
  
#.........................................................................................
### Exploring and preparing predictor variables
#.........................................................................................
      
  #...................................
  ## Identify "parent" variables (i.e. from which lags were derived)
    # Identify variables that are not parent
    x1 <- na.omit(grep("_lag", colnames(ch_obs)))
    vars_in <- colnames(ch_obs)[-x1]
    
    # Also eliminate variables that are not potential predictor variables...
    vars_in <- vars_in[! vars_in %in% c("survey_id", "stratum", "admin1")]

    #...any that are not included in the list of variable parameters...
    vars_in <- vars_in[vars_in %in% var_pars$variable]
    
    # ...and any that are being forced out of the model
    vars_in <- vars_in[! vars_in %in% subset(var_pars, force=="out")[, "variable"] ]
 
  #...................................
  ## Explore distributions of predictors

    # For each parent variable...
    for (i in vars_in) { 
      # full range  
      f_hist(i, ch_obs, c(NA, NA)) 
      
      # zooming in
      f_hist(i, ch_obs, c(NA, quantile(ch_obs[, i], 0.90, na.rm = TRUE))) 
      
    }

        
  #...................................
  ## Categorise predictors if needed (including lags)
    
    # For each variable in the var_pars data frame...
    for (i in 1:nrow(var_pars) ) {
      
      # if the variable is among the parent variables and should be categorised...  
      if (var_pars[i, "variable"] %in% vars_in & is.na(var_pars[i, "cat_only"]) == FALSE ) {
       
        # work out the indices of the variable and any lags
        ind <- grep(var_pars[i, "variable"], colnames(ch_obs) )
        names <- colnames(ch_obs)[ind]
          
        # if the variable should be categorised according to cut-offs...
        if (var_pars[i, "cat_method"] == "cut") {
          
          # work out cut-offs and corresponding labels
          x1 <- as.numeric(unlist(strsplit(var_pars[i, "cat_cutoffs"], split=",")))
          x2 <- unlist(strsplit(var_pars[i, "cat_labels"], split=", "))
          
          # for the variable and each of its lags, if any...
          for (j in names) {
            # create new categorical variables
            ch_obs[, paste(j , "_cat", sep="")] <- cut(ch_obs[, j], breaks = x1, labels = x2, include.lowest = TRUE)
          } 
          print(paste("categorisation of variable ", var_pars[i, "variable"], " :", sep=""))
          print(table(ch_obs[, paste(var_pars[i, "variable"] , "_cat", sep="")] ))
            
        }
        
        # if the variable should be categorised according to values...
        if (var_pars[i, "cat_method"] == "values") {
          
          # work out categories for values
          x1 <- cbind(unlist(strsplit(var_pars[i, "cat_values"], split=",")),
            unlist(strsplit(var_pars[i, "cat_labels"], split=", ")) )
          
          # for the variable and each of its lags, if any...
          for (j in names) {
            
            # assign corresponding column names
            colnames(x1) <- c(j, paste(j , "_cat", sep="") )
            
            # create new categorical variable in training dataset
            ch_obs <- merge(ch_obs, x1, by = paste(colnames(x1)[1]), all.x = TRUE)
          }
          print(paste("categorisation of variable ", var_pars[i, "variable"], " :", sep=""))
          print(table(ch_obs[, paste(var_pars[i, "variable"] , "_cat", sep="")] ))
          
        }

        # add categorical variable to list of parent variables to be considered in model
        vars_in <- c(vars_in, paste(names , "_cat", sep=""))
        
        # if only the categorical variables are to be retained...
        if ( var_pars[i, "cat_only"] == "Y"  ) {
          #...remove the continuous ones and from list of variables to be considered in model
          vars_in <- vars_in[! vars_in %in% names ]
          
        }
        
      }
      
    }


  #...................................
  ## Convert categorical variables to factors
  for (i in colnames(ch_obs)) {
    if (typeof(ch_obs[, i]) == "character") {ch_obs[, i] <- as.factor(ch_obs[, i])}
  }  
    
  #...................................
  ## Change reference categories, if needed
    # For each variable in the var_pars data frame...
    for (i in 1:nrow(var_pars) ) {
      # if the variable is among the parent variables and should be categorised...  
      if (var_pars[i, "variable"] %in% vars_in & is.na(var_pars[i, "cat_only"]) == FALSE ) {
        
        # work out the indices of the categorised variable and any lags
        ind_1 <- grep(var_pars[i, "variable"], colnames(ch_obs))
        ind_2 <- grep("_cat", colnames(ch_obs))
        ind <- intersect(ind_1, ind_2)
        names <- colnames(ch_obs)[ind]
          
        # if the reference category should be specified...
        if (is.na(var_pars[i, "cat_ref"]) == FALSE) {
          # for the variable and each of its lags, if any...
          for (j in names) {
            # specify reference category in the training dataset
            ch_obs[, j] <- relevel(ch_obs[, j], ref = var_pars[i, "cat_ref"])
          } 
        }
      }
  }

    
  #...................................
  ## Update the list of predictors to include in modelling, to make sure lags are still included
  x1 <- c()
  for (i in vars_in) {x1 <- c(x1, grep(i, colnames(ch_obs))) }
  vars_in <- colnames(ch_obs)[x1]
  vars_in <- unique(vars_in)

  
  #...................................
  ## Restrict lags as desired
    
    # Eliminate all lags below the desired time horizon for prediction
    if (horizon > 0) {

      # eliminate lag 0
      x1 <- grep(paste("_lag", horizon, sep = ""), vars_in, value = TRUE) 
      x1 <- gsub(paste("_lag", horizon, sep = ""), "", x1)
      vars_in <- vars_in[! vars_in %in% x1]

      # eliminate lags > 0 and less than the horizon
      for (i in 0:(horizon - 1)) { 
        x1 <- grep(paste("_lag", i, sep = ""), vars_in)
        if (length(x1) > 0) {vars_in <- vars_in[-x1] }
      }
    }

    # Remove lags outside a specified range, if desired
    x1 <- subset(var_pars, ! is.na(restrict_lags) & ! force %in% c("out"))[, c("variable", "restrict_lags")]
    x2 <- rep(TRUE, times = length(vars_in))
    for (i in x1$variable) {
      x2[grep(i, vars_in)] <- FALSE
      x3 <- x1[x1$variable == i, "restrict_lags"]
      x3 <- suppressWarnings(as.numeric(unlist(strsplit(x3, ","))) )
      if (is.na(x3[1]) ) {x3[1] <- 0}
      if (is.na(x3[2]) ) {x3[2] <- 10}
      if (0 %in% x3) {x2[grepl(i, vars_in) & ! grepl("_lag", vars_in)] <- TRUE}
      if (x3[2] != 0) {
        for (j in x3[1]:x3[2]) {
          x2[grepl(i, vars_in) & grepl(paste("_lag", j, sep=""), vars_in)] <- TRUE
        }
      }
    }
    
    vars_in <- vars_in[x2]
      
    # Retain any specific lags if desired
    x1 <- subset(var_pars, ! is.na(force_lag) & ! force %in% c("out"))[, c("variable", "force_lag")]
    x2 <- rep(TRUE, times = length(vars_in))
    for (i in x1$variable) {
      x3 <- x1[x1$variable == i, "force_lag"]
      if(grepl(",", x3)) {x4 <- as.integer(unlist(strsplit(x3, ",")))}
      if(! grepl(",", x3)) {x4 <- as.integer(x3)}

      x2[grep(i, vars_in)] <- FALSE
      if (0 %in% x4) {x2[grepl(i, vars_in) & ! grepl("_lag", vars_in)] <- TRUE}
      for (j in 1:12) {
        if (j %in% x4) {x2[grepl(i, vars_in) & grepl(paste("_lag", j, sep=""), vars_in)] <- TRUE}
      }
    }
    
    vars_in <- vars_in[x2]

  
#.........................................................................................
### Univariate analysis
#.........................................................................................

  #...................................
  ## Split data into holdout and training sets
  
    # Hold-out set
    ch_obs_h <- subset(ch_obs, holdout == TRUE)
    
    # Training set (everything else)
    ch_obs_t <- subset(ch_obs, holdout == FALSE)
  
  
  #...................................
  ## Observe collinearity among pairs of predictors
    
    # Select all continuous variables and continuous versions of those that have been categorised
    x1 <- unique(c(vars_in[- which(sapply(ch_obs_t[, vars_in], is.factor))], gsub("_cat", "", vars_in)) )
        
    # Collinearity heatmap
    plot <- ggcorr(ch_obs_t[, x1], method = c("pairwise.complete.obs", "pearson"),
      low = "steelblue", mid = "grey90", high = "darkred", geom = "tile", nbreaks = 5, min_size = 0, max_size = 6, 
      label = TRUE, label_size = 3, label_round = 2, size = 3, hjust = 0.75, layout.exp = 1, legend.position = "off")
    print(plot)
    ggsave(paste(country, "_collinearity_heatmap.png", sep=""), plot, height = 45, width = 45, units = "cm", dpi = "print")    
        
                    
  #...................................
  ## Univariate analysis
    
    # Fit univariate models for each predictor
    stats <- c("aic_glm", "mse_glm", "p_value_glm", "mse_cv_glm", "prop_folds_ok_glm", "mse_holdout_glm",
      "aic_glmm", "mse_glmm", "p_value_glmm", "mse_cv_glmm", "prop_folds_ok_glmm", "mse_holdout_glmm")
    out <- data.frame(matrix(NA, nrow = length(vars_in), ncol = 1 + length(stats))) 
    colnames(out) <- c("predictor", stats)
    out[, "predictor"] <- vars_in
    
    for (i in 1:nrow(out)) {
      print (paste("now fitting univariate models for predictor ", out[i, "predictor"], sep="") )
      out[i, stats] <- 
        f_val(y_hat, out[i, "predictor"], c(), ch_obs_t, NA, re, fit_family, FALSE, part_unit, k_folds, FALSE)
    }

    # Select best-fitting lags and categorical vs continuous versions, for the same base variable
      # first need to identify base variables
      out[, "base_var"] <- NA
      for (i in var_pars$variable) { out[grep(i, out$predictor), "base_var"] <- i }
    
      # then identify version with lowest deviance p-test
      x1 <- tapply(out$p_value_glm, out$base_var, which.min)
      out[, "keep"] <- FALSE
      for (i in unique(out$base_var)) {
        x2 <- rep(FALSE, nrow(out[out$base_var == i, ]) )
        x2[x1[names(x1) == i]] <- TRUE
        out[out$base_var == i, "keep"] <- x2
      }
      
      # keep all lags in if desired for any predictor
      x1 <- subset(var_pars, multiple_lags == "Y")$variable
      if (length(x1) > 0) {
        for (i in x1) {
          x2 <- out[out$base_var == i, "predictor"]
          if (length(grep("_cat", out[out$predictor %in% x2 & out$keep == TRUE, "predictor"]) ) > 0) {
            out[out$predictor %in% x2 & grepl("_cat", out$predictor) & out$keep == FALSE, "keep"] <- TRUE
          }
          if (length(grep("_cat", out[out$predictor %in% x2 & out$keep == TRUE, "predictor"]) ) == 0) {
            out[out$predictor %in% x2 & ! grepl("_cat", out$predictor) & out$keep == FALSE, "keep"] <- TRUE
          }
          
        }
      }
      
      vars_in <- out[out$keep == TRUE, "predictor"]
      
    # Screen out variables that don't pass the univariate screening p-value threshold test
    vars_in <- vars_in[vars_in %in% out[out$p_value_glm < f_univar, "predictor"]]

    # Save output
    write.csv(out, paste(country, "_est_", y_hat, "_univar_models.csv", sep = ""), row.names = FALSE)
      
    
#.........................................................................................
### Fit, evaluate and select among candidate models
#.........................................................................................
    
  #...................................
  ## Specify data frame of all possible candidate models
    
    # Create all possible models
    mods <- as.data.frame(matrix(c(0,1), nrow = 2, ncol = length(vars_in)))
    colnames(mods) <- vars_in
    mods <- expand.grid(mods[, vars_in])
      # remove first row (null model)
      mods <- mods[-1, ]
    
      # remove models where variables to be forced in are left out (i.e. == 0)
      if (nrow(subset(var_pars, force=="in")) > 0) {
        x1 <- grep(paste(subset(var_pars, force=="in")[, "variable"], collapse="|"), colnames(mods), value=FALSE)
        x2 <- apply(mods, 1, function(x, x1) {if (all(x[x1]==1)) return(TRUE) else return(FALSE)}, x1 ) 
        mods <- mods[which(x2 == TRUE), ]
      }
      
    # Control statement  
    print(paste("the number of possible models at this stage is: ", nrow(mods), sep=""))

  #...................................  
  ## Fit all possible fixed effects models and retain the most promising few
    ### NOTE: COMPUTATIONALLY INTENSIVE IF LARGE NUMBER OF CANDIDATE MODELS ARE FIT
    
    # Fit all models
    out <- as.data.frame(matrix(NA, ncol = ncol(mods) + length(stats), nrow = nrow(mods)))
    colnames(out) <- c(colnames(mods), stats)
    for (i in 1:nrow(mods)) {
      print(paste("now fitting model ", i, " of ", nrow(mods), ":", sep="") )
      x1 <- names(mods)[mods[i, ] == 1]
      out[i, colnames(mods)] <- mods[i, ]
      out[i, stats] <- f_val(y_hat, x1, c(), ch_obs_t, ch_obs_h, re, fit_family, FALSE, part_unit, k_folds, FALSE)
    }  

    # Save output
    write.csv(out, paste(country, "_est_", y_hat, "_all_models.csv", sep = ""), row.names = FALSE)
    
    # Keep the best models, based on their holdout mse score being in the bottom n% (f_multivar) of all models
    x1 <- quantile(out$mse_holdout_glm, f_multivar)
    out <- out[order(out[, "mse_holdout_glm"]), ]
    out_best <- out[out$mse_holdout_glm <= x1, ]
    View(out_best)
    
   
  #...................................  
  ## Evaluate best candidate models for overfitting and predictive power on cross-validation and holdout dataset
    ### NOTE: COMPUTATIONALLY INTENSIVE IF LARGE NUMBER OF CANDIDATE MODELS ARE PUT THROUGH CV, OR IF LOOCV IS USED
    
    # Prepare output
    out_best_stats <- as.data.frame(matrix(NA, ncol = ncol(mods) + length(stats), nrow = nrow(out_best)))
    colnames(out_best_stats) <- c(colnames(mods), stats)
    
    # For each of the models...
    for (i in 1:nrow(out_best)) {
      print("################################################")
      print(paste("now fitting model ", i, " of ", nrow(out_best), ":", sep="") )
      
      # specify model formula 
      x1 <- names(mods)[out_best[i, colnames(mods)] == 1]
      out_best_stats[i, colnames(mods)] <- out_best[i, colnames(mods)]
      
      # fit model with cross-validation and predicting on holdout data, and compute fit statistics
      out_best_stats[i, stats] <- f_val(y_hat, x1, c(), ch_obs_t, ch_obs_h, 
        re, fit_family, TRUE, part_unit, k_folds, TRUE)
      
      # compute penalty in mse score on cross-validation
      out_best_stats[i, "mse_penalty_cv_glm"] <- out_best_stats[i, "mse_cv_glm"] - 
        out_best_stats[i, "mse_glm"]
      out_best_stats[i, "mse_penalty_cv_glmm"] <- out_best_stats[i, "mse_cv_glmm"] - 
        out_best_stats[i, "mse_glmm"]
      
      # compute penalty in mse score on holdout
      out_best_stats[i, "mse_penalty_holdout_glm"] <- out_best_stats[i, "mse_holdout_glm"] - 
        out_best_stats[i, "mse_glm"]
      out_best_stats[i, "mse_penalty_holdout_glmm"] <- out_best_stats[i, "mse_holdout_glmm"] - 
        out_best_stats[i, "mse_glmm"]
      
    }  
    
    View(out_best_stats)
    
    # Save output
    write.csv(out_best_stats, paste(country, "_est_", y_hat, "_best_models.csv", sep = ""), row.names = FALSE)

  #...................................  
  ## Evaluate best fixed-effects model
        
    # Fit and plot best fixed-effects model
      # identify best model based on mse score on CV
      x1 <- which.min(out_best_stats$mse_cv_glm)
        # or specify x1 = row number of desired model
      
      # select predictors
      vars_in <- names(mods)[out_best_stats[x1, colnames(mods)] == 1]      
      
      # fixed effect GLM with cross-validation and on holdout sample
      scores_glm <- out_best_stats[x1, c("aic_glm", "mse_glm", "p_value_glm", "mse_cv_glm", "prop_folds_ok_glm", 
        "mse_holdout_glm", "aic_glmm", "mse_glmm", "p_value_glmm", "mse_cv_glmm", "prop_folds_ok_glmm",
        "mse_holdout_glmm")]
      form <- as.formula( paste(y_hat, "~", paste(vars_in, collapse = "+"), sep = "")  ) 
      fit_glm <- glm(form, family = fit_family, data = ch_obs_t, weights = wt )

      # observations and predictions on 
        # training sample, 
        x1 <- data.frame(ch_obs_t[, part_unit], ch_obs_t[, "tm"], ch_obs_t[, y_hat], 
          predict(fit_glm, newdata = ch_obs_t, type = "response", allow.new.levels = TRUE))
        colnames(x1) <- c(part_unit, "tm", "observed", "predicted")
        x1 <- aggregate(x1[, c("observed", "predicted")], by = x1[, c(part_unit, "tm")], FUN = mean, na.rm=TRUE)
        colnames(x1) <- c(part_unit, "tm", "observed", "predicted")
        
        # cross-validation,
        x2 <- f_cv(fit_glm, ch_obs_t, part_unit, k_folds, FALSE)
        colnames(x2) <- c(part_unit, "tm", "observed", "predicted")
 
        # and holdout sample
        x3 <- data.frame(ch_obs_h[, part_unit], ch_obs_h[, "tm"], ch_obs_h[, y_hat], 
          predict(fit_glm, newdata = ch_obs_h, type = "response", allow.new.levels = TRUE))
        colnames(x3) <- c(part_unit, "tm", "observed", "predicted")
        x3 <- aggregate(x3[, c("observed", "predicted")], by = x3[, c(part_unit, "tm")], FUN = mean, na.rm=TRUE)
        colnames(x3) <- c(part_unit, "tm", "observed", "predicted")
      
      # combined plot    
      f_plot_prev(country, y_hat, x1, x2, x3, "glm")
            

#.........................................................................................
### Exploring interactions and random effects; finalising choice of model
#.........................................................................................
            
  #...................................
  ## Explore plausible interactions and update fixed-effects fit accordingly

    # Identify any predictors involved in plausible two-way interactions
    int_vars <- subset(var_pars, is.na(interactions) == FALSE )[, c("variable", "interactions")]

    # Check whether these predictors are in the dataset and identify the right variant of the predictors
    if (nrow(int_vars) > 0) {int_vars[, "actual_variable"] <- NA}
    for (i in int_vars[, "variable"]) {
      
      # find all the variants of this variable in the dataset
      x1 <- grep(i, colnames(ch_obs), value = TRUE )
      
      # if the variable is not in the dataset, break
      if (length(x1) == 0) next;
      
      # if there is only one variant, use that
      if (length(x1) == 1) {int_vars[int_vars$variable == i, "actual_variable"] <- x1}
      
      # if there are several variants...
        # if one is in the model formula, choose that, but continuous version
        if (length(intersect(x1, all.vars(formula(fit_glm)) )) > 0 ) {
          int_vars[int_vars$variable == i, "actual_variable"] <- x1[x1 %in% all.vars(formula(fit_glm))]
          int_vars[int_vars$variable == i, "actual_variable"] <- 
            gsub("_cat", "", int_vars[int_vars$variable == i, "actual_variable"])
        }
        
        # if none is in the model formula...
        if (length(intersect(x1, all.vars(formula(fit_glm)) )) == 0 ) {
          # if a specific lag has been specified, choose that
          if (! is.na(var_pars[var_pars$variable == i, "force_lag"]) &  var_pars[var_pars$variable == i, "force_lag"] != 0) {
            x2 <- grep(paste("lag", var_pars[var_pars$variable == i, "force_lag"], sep =""), x1, value = TRUE)
            x2 <- gsub("_cat", "", x2)
            int_vars[int_vars$variable == i, "actual_variable"] <- unique(x2)
          }
          if (! is.na(var_pars[var_pars$variable == i, "force_lag"]) &  var_pars[var_pars$variable == i, "force_lag"] == 0) {
            x2 <- x1[! grepl("lag", x1)]
            x2 <- gsub("_cat", "", x2)
            int_vars[int_vars$variable == i, "actual_variable"] <- unique(x2)
          }
          # otherwise choose the variant without lag
          if (is.na(var_pars[var_pars$variable == i, "force_lag"]) ) {
            x2 <- x1[! grepl("lag", x1)]
            x2 <- gsub("_cat", "", x2)
            int_vars[int_vars$variable == i, "actual_variable"] <- unique(x2)
          }
        }
    }
    
    # Test each interaction in turn and record fit statistics
    if (nrow(int_vars) > 0) {
     
      # prepare
      x1 <- unique(int_vars$interactions)
      interaction_terms <- c()
      
      # for each potential interaction term...
      for (i in x1) {
        # create interaction term
        x2 <- paste(as.vector(subset(int_vars, interactions == i)$actual_variable ), collapse=":")
        print(paste("now fitting model with interaction term  ", x2, sep = "") )
  
        # evaluate model
        scores_int <- f_val(y_hat, vars_in, x2, ch_obs_t, ch_obs_h, FALSE, fit_family, TRUE, part_unit, k_folds, TRUE)

        # write the model formula
        form <- as.formula(paste(paste(terms(fit_glm))[2], " ~ ", paste(terms(fit_glm))[3], " + ", x2, sep = ""))

        # fit GLM
        fit_int <- update(fit_glm, formula = form, data = ch_obs_t)
          # print model summary
          print(tidy(fit_int, exponentiate = ifelse(fit_family == "gaussian", FALSE, TRUE) ))
          print(glance(fit_int))
          
        # keep model with interaction term, if desired and if it improves fit on both CV and holdout
        if (interactions_in == TRUE & scores_int["prop_folds_ok_glm"] >= 0.90 & scores_glm["prop_folds_ok_glm"] >= 0.90 &
            scores_int["mse_cv_glm"] < scores_glm["mse_cv_glm"] &
            scores_int["mse_holdout_glm"] < scores_glm["mse_holdout_glm"]) { 
          fit_glm <- fit_int
          interaction_terms <- c(interaction_terms, x2)
          }
      }
    }
    
    # Best model after testing for interactions
    tidy(fit_glm, exponentiate = ifelse(fit_family == "gaussian", FALSE, TRUE))
    glance(fit_glm)

        
  #...................................  
  ## Evaluate mixed model option, if desired
    ### NOTE: COMPUTATIONALLY INTENSIVE, ESPECIALLY IF THE DATASET IS LARGE
  
  if (force_glm != "glm") {    
    # Fit mixed model
    scores_glmm <- f_val(y_hat, vars_in, interaction_terms, ch_obs_t, ch_obs_h, TRUE, fit_family, TRUE, part_unit, 10, TRUE)      
    form <- as.formula( paste(y_hat, "~", paste(c(vars_in, interaction_terms), collapse = "+"), "+ (1|", part_unit, ")",  sep = "")  )
    fit_glmm <- glmer(form, family = gsub("quasi", "", fit_family), data = ch_obs_t) # note no weights (fitting problems)
    summary(fit_glmm)
  } 
    
    # Select between fixed-effects only and mixed model and save fit statistics
      # if the choice should be based on relative mse score...
      if (force_glm == "either") {
        fit_best <- fit_glm
        fit_best_stats <- scores_glm
          
        if (scores_glmm["prop_folds_ok_glmm"] >= 0.90) {
          if (scores_glmm["mse_cv_glmm"] < scores_glmm["mse_cv_glm"]) 
          {fit_best <- fit_glmm; fit_best_stats <- scores_glmm}
        }
      } 
      
      # if the GLM (fixed-effects only) model should be adopted...
      if (force_glm == "glm") {
        fit_best <- fit_glm; 
        fit_best_stats <- scores_glm
      } 
      
      # if the GLMM (mixed-effects) model should be adopted...
      if (force_glm == "glmm") {
        fit_best <- fit_glmm; 
        fit_best_stats <- scores_glmm
      } 
      
  #...................................  
  ## Calculate robust standard errors (only if a fixed-effect model is selected)
  
  if (! part_unit %in% all.vars(formula(fit_best)) ) {      
    # Generate variance-covariance matrix by specifying the desired (unique) cluster variable(s)
    ch_obs_t[, cluster_vars] <- paste(ch_obs_t[, "survey_id"], "_", ch_obs_t[, cluster_vars], sep = "")
    ch_obs_h[, cluster_vars] <- paste(ch_obs_h[, "survey_id"], "_", ch_obs_h[, cluster_vars], sep = "")
    vcov_cl <- cluster.vcov(fit_best, ch_obs_t[, cluster_vars])   
    
    # Compute model output with robust standard errors
    out <- f_rob(fit_best, vcov_cl, "response")
    print(out)
  }

  #...................................  
  ## Calculate additional model performance metrics on cross-validation and holdout
    ### NOTE: COMPUTATIONALLY INTENSIVE IF RANDOM EFFECTS MODEL IS SELECTED
    
    # Thresholds to be used for classification accuracy
    thresholds <- as.numeric(unlist(strsplit(thresholds, ",")))
      
    # Do cross-validation and compute metrics by partition unit
    cv_metrics <- f_cv_metrics(fit_best, ch_obs_t, part_unit, NA, thresholds, TRUE)
    
    # Predict on holdout and compute metrics by partition unit
    holdout_metrics <- f_holdout_metrics(fit_best, ch_obs_h, part_unit, thresholds, TRUE)
    
    # Plot performance accuracy
      # observations and predictions on 
        # training sample, 
        ys <- ch_obs_t[, y_hat]
        x1 <- data.frame(ch_obs_t[, part_unit], ch_obs_t[, "tm"], ys, 
          predict(fit_best, newdata = ch_obs_t, type = "response", allow.new.levels = TRUE))
        colnames(x1) <- c(part_unit, "tm", "observed", "predicted")
        x1 <- aggregate(x1[, c("observed", "predicted")], by = x1[, c(part_unit, "tm")], FUN = mean, na.rm=TRUE)
        colnames(x1) <- c(part_unit, "tm", "observed", "predicted")
        
        # cross-validation,
        x2 <- f_cv(fit_best, ch_obs_t, part_unit, k_folds, FALSE)
        colnames(x2) <- c(part_unit, "tm", "observed", "predicted")
 
        # and holdout sample
        ys <- ch_obs_h[, y_hat]        
        x3 <- data.frame(ch_obs_h[, part_unit], ch_obs_h[, "tm"], ys, 
          predict(fit_best, newdata = ch_obs_h, type = "response", allow.new.levels = TRUE))
        colnames(x3) <- c(part_unit, "tm", "observed", "predicted")
        x3 <- aggregate(x3[, c("observed", "predicted")], by = x3[, c(part_unit, "tm")], FUN = mean, na.rm=TRUE)
        colnames(x3) <- c(part_unit, "tm", "observed", "predicted")
      
      # combined plot    
      f_plot_prev(country, y_hat, x1, x2, x3, "glm")
            
      
  #...................................  
  ## Write final model to file
        
    # Write output to file
    x1 <- paste(country, "_est_", y_hat, "_final_model", ".csv", sep ="")
    if (part_unit %in% all.vars(formula(fit_best)) ) {
      write.table(cbind(broom.mixed::tidy(fit_best, exponentiate = TRUE), exp(confint(fit_best))),
        x1, sep = ",", row.names = FALSE)
      write.table("--------------------", x1, sep = ",", col.names = FALSE, row.names = FALSE, append = TRUE)
      write.table("Fit statistics:", x1, sep = ",", col.names = FALSE, row.names = FALSE, append = TRUE)
      write.table(rbind(c("aic", "p_value", "mse", "mse_cv", "mse_holdout") , 
        fit_best_stats[c("aic_glmm", "p_value_glmm", "mse_glmm", "mse_cv_glmm", "mse_holdout_glmm")]), 
        x1, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
    }
    if (! part_unit %in% all.vars(formula(fit_best)) ) {
      write.table(out, x1, sep = ",", row.names = FALSE) 
      write.table("--------------------", x1, sep = ",", col.names = FALSE, row.names = FALSE, append = TRUE)
      write.table("Fit statistics:", x1, sep = ",", col.names = FALSE, row.names = FALSE, append = TRUE)
      write.table(rbind(c("aic", "p_value", "mse", "mse_cv", "mse_holdout") , 
        fit_best_stats[c("aic_glm", "p_value_glm", "mse_glm", "mse_cv_glm", "mse_holdout_glm")]), 
        x1, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
    }
    write.table("--------------------", x1, sep = ",", col.names = FALSE, row.names = FALSE, append = TRUE)
    write.table("Additional performance metrics...", x1, sep = ",", col.names = FALSE, row.names = FALSE, append = TRUE)
    write.table("  on cross-validation:", x1, sep = ",", col.names = FALSE, row.names = FALSE, append = TRUE)
    write.table(rbind(names(cv_metrics), cv_metrics), x1, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
    write.table("  on holdout dataset:", x1, sep = ",", col.names = FALSE, row.names = FALSE, append = TRUE)
    write.table(rbind(names(holdout_metrics), holdout_metrics), x1, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)

  # Save fit object
  saveRDS(fit_best, paste(country, "_", y_hat, "_final_model", ".rds", sep="")) 


#.........................................................................................
### ENDS
#.........................................................................................
