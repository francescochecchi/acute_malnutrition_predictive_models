#..........................................................................................
### +++++++++++++ SMALL-AREA ESTIMATION OF CRISIS-ATTRIBUTABLE MORTALITY ++++++++++++++ ###
#..........................................................................................

#..........................................................................................
## --- R CODE TO EVALUATE RANDOM FOREST METHODS TO PREDICT ACUTE MALNUTRITION BURDEN --- ##
#..........................................................................................

                                       # Written by Francesco Checchi, LSHTM (May 2021)
                                       # francesco.checchi@lshtm.ac.uk 



#.........................................................................................
### Reading in required files
#.........................................................................................
    
  #...................................
  ## Model data (from earlier script)
  ch_obs <- read.csv(paste(country, "_ch_obs_model.csv", sep = ""), sep = "," )


#.........................................................................................
### Preparing the data for random forest
#.........................................................................................

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
    # NOTE: "cv" option should not be selected alongside random forest

 
  #...................................
  ## Identify "parent" variables (i.e. from which lags were derived)
    # Identify variables that are not parent
    x1 <- na.omit(grep("_lag", colnames(ch_obs)))
    vars_in <- colnames(ch_obs)[-x1]
    
    # Also eliminate variables that are not potential predictor variables...
    vars_in <- vars_in[! vars_in %in% c("survey_id", "stratum", "cluster", "admin1")]

    #...any that are not included in the list of variable parameters...
    vars_in <- vars_in[vars_in %in% var_pars$variable]
    
    # ...and any that are being forced out of the model
    vars_in <- vars_in[! vars_in %in% subset(var_pars, force=="out")[, "variable"] ]
  
  #...................................
  ## Update the list of predictors to include in random forest, to make sure lags are still included
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

    # Add observation weights and cluster IDs to variables to be retained   
    vars_in <- c(vars_in, "wt", "cluster")
    
  #...................................  
  ## Prepare data
    
    # Factorise outcome if binary (needed for classification)  
    if (length(unique(na.omit(ch_obs[, y_hat]))) == 2) {
      ch_obs[, y_hat] <- as.factor(ch_obs[, y_hat]) 
    }    

    # Split data into holdout and training sets
      # hold-out set
      ch_obs_h <- subset(ch_obs, holdout == TRUE)
      
      # training set (everything else)
      ch_obs_t <- subset(ch_obs, holdout == FALSE)
    
    # Impute data, if desired; otherwise just select complete cases
    df_t <- f_rf_prep(y_hat, vars_in, ch_obs_t, FALSE, part_unit)
    df_h <- f_rf_prep(y_hat, vars_in, ch_obs_h, FALSE, part_unit)
    
    # Aggregate data by partition unit
      # training data
     if (is.factor(df_t[, y_hat])) {df_t[, y_hat] <- 
        as.numeric(levels(df_t[, y_hat]))[df_t[, y_hat]]}
      x1 <- aggregate(df_t[, y_hat], by = df_t[, c(part_unit, "tm")], FUN = mean, na.rm = TRUE)
      colnames(x1) <- c(part_unit, "tm", y_hat)
      x2 <- unique(df_t[, c(part_unit, "tm")])
      x3 <- c()
      for (i in 1:nrow(x2)) {x3 <- rbind(x3, df_t[df_t[, part_unit] == x2[i, 1] & df_t[, "tm"] == x2[i, 2], ][1, c(part_unit, "tm", vars_in)]) }
      df_t <- merge(x1, x3, by = c(part_unit, "tm"))
      
      # holdout data
     if (is.factor(df_h[, y_hat])) {df_h[, y_hat] <- 
        as.numeric(levels(df_h[, y_hat]))[df_h[, y_hat]]}
      x1 <- aggregate(df_h[, y_hat], by = df_h[, c(part_unit, "tm")], FUN = mean, na.rm = TRUE)
      colnames(x1) <- c(part_unit, "tm", y_hat)
      x2 <- unique(df_h[, c(part_unit, "tm")])
      x3 <- c()
      for (i in 1:nrow(x2)) {x3 <- rbind(x3, df_h[df_h[, part_unit] == x2[i, 1] & df_h[, "tm"] == x2[i, 2], ][1, c(part_unit, "tm", vars_in)]) }
      df_h <- merge(x1, x3, by = c(part_unit, "tm"))

#.........................................................................................
### Growing random forest and evaluating its performance
#.........................................................................................
        
  #...................................  
  ## Grow random forest       
  print("now growing forest...")
  
    # # Classification for binary outcomes
    # if (length(unique(df_t[, y_hat])) == 2) {
    #   
    #   # # # grow random forest
    #   # fit <- randomForest(as.formula(paste(y_hat, "~", ".", sep = "")),
    #   #   data = df_t[complete.cases(df_t[, c(y_hat, vars_in)]), c(y_hat, vars_in)],
    #   #   sampsize = c(1000, 1000), strata = factor(df_t[, y_hat]),
    #   #   na.action = na.fail, ntree = 100, mtry = 3, localImp = TRUE)
    #   # print(fit)
    #   
    #   # grow and tune random forest
    #   fit <- tuneRF(df_t[complete.cases(df_t[, c(y_hat, vars_in)]), vars_in], df_t[complete.cases(df_t[, c(y_hat, vars_in)]), y_hat],
    #     mtryStart = 3, ntreeTry = 50, stepFactor = 1.5, improve = 0.05, trace = TRUE, plot = TRUE, doBest = TRUE,
    #     na.action = na.fail, sampsize = c(1000, 1000), strata = factor(df_t[, y_hat]), localImp = TRUE )
    #   print(fit)
    # }
    
    # # Regression for continuous outcomes
    # if (length(unique(df_t[, y_hat])) > 2) {
      
      # # grow random forest
      # fit <- randomForest(as.formula(paste(y_hat, "~", ".", sep = "")), data = df_t[, c(y_hat, vars_in)],
      #   na.action = na.fail, ntree = 100, mtry = 4, maxnodes = 5, localImp = TRUE)
      # print(fit)
      
      # grow and tune random forest
      # fit <- tuneRF(df_t[complete.cases(df_t[, c(y_hat, vars_in)]), vars_in], df_t[complete.cases(df_t[, c(y_hat, vars_in)]), y_hat],
      #   mtryStart = 3, ntreeTry = 50, stepFactor = 1.5, improve = 0.05, trace = TRUE, plot = TRUE, doBest = TRUE,
      #   na.action = na.fail, localImp = TRUE )
      
    # }
  
  form <- as.formula(paste(y_hat, " ~ ", paste(vars_in[! vars_in %in% c("wt", part_unit, "cluster")], collapse = " + "), sep = ""))
  fit <- ranger(formula = form, data = df_t, case.weights = "wt", keep.inbag = TRUE, importance = "permutation", num.trees = 1000)
  print(fit)
      
  
  #...........................
  ## Generate fit performance and variable influence statistics
  
    # Predictions on training data
    pred <- predict(fit, data = df_t)
    
      # aggregate predictions
      x1 <- as.data.frame(cbind(df_t[, c(y_hat, part_unit, "tm")], pred$predictions))
      colnames(x1) <- c("observed", part_unit, "tm", "predicted")
      if (is.factor(x1[, "observed"])) {x1[, "observed"] <- 
        as.numeric(levels(x1[, "observed"]))[x1[, "observed"]]}
      if (is.factor(x1[, "predicted"])) {x1[, "predicted"] <- 
        as.numeric(levels(x1[, "predicted"]))[x1[, "predicted"]]}
      pred_t <- aggregate(x1[, c("observed", "predicted")], by = x1[, c(part_unit, "tm")], FUN = mean, na.rm = TRUE)
      colnames(pred_t) <- c(part_unit, "tm", "observed", "predicted")
 
    # Predictions on cross-validation
    pred_cv <- f_cv_rf(fit, df_t, part_unit, k_folds)
    
    # Predictions on holdout data
    pred <- predict(fit, data = df_h )
    
      # aggregate predictions
      x1 <- as.data.frame(cbind(df_h[, c(y_hat, part_unit, "tm")], pred$predictions))
      colnames(x1) <- c("observed", part_unit, "tm", "predicted")
      if (is.factor(x1[, "observed"])) {x1[, "observed"] <- 
        as.numeric(levels(x1[, "observed"]))[x1[, "observed"]]}
      if (is.factor(x1[, "predicted"])) {x1[, "predicted"] <- 
        as.numeric(levels(x1[, "predicted"]))[x1[, "predicted"]]}
      pred_h <- aggregate(x1[, c("observed", "predicted")], by = x1[, c(part_unit, "tm")], FUN = mean, na.rm = TRUE)
      colnames(pred_h) <- c(part_unit, "tm", "observed", "predicted")

    # Plot predictions
    f_plot_prev("som", y_hat, pred_t, pred_cv, pred_h, "rf")
  
    # Compute MSEs
    mses <- c(
      mean((pred_t$predicted - pred_t$observed )^2, na.rm = TRUE),
      mean((pred_cv$predicted - pred_cv$observed )^2, na.rm = TRUE),
      mean((pred_h$predicted - pred_h$observed )^2, na.rm = TRUE)
    )
    names(mses) <- c("mse_training", "mse_cv", "mse_holdout")
    
    # Generate performance metrics on cross-validation and on holdout dataset
    thresholds <- as.numeric(unlist(strsplit(thresholds, ",")))
    cv_metrics <- f_cv_metrics_rf(fit, df_t, ch_obs_t, part_unit, k_folds, thresholds, TRUE)
    holdout_metrics <- f_holdout_metrics_rf(fit, df_h, ch_obs_h, part_unit, thresholds, TRUE)
    
    # Generate variable importance metrics
    print("now evaluating variable importance...")
    x1 <- importance_pvalues(fit, formula = formula(fit), data = df_t, method = "altmann")
    out <- data.frame("predictor" = row.names(x1), x1)
    
  #...........................
  ## Write output to file
  x1 <- paste(country, "_est_", y_hat, "_final_rf", ".csv", sep ="")
  write.table(capture.output(fit, type = "output"), x1, sep = ",", col.names = FALSE, row.names = FALSE)
  write.table("--------------------", x1, sep = ",", col.names = FALSE, row.names = FALSE, append = TRUE)
  write.table("Mean square errors:", x1, sep = ",", col.names = FALSE, row.names = FALSE, append = TRUE)
  write.table(mses, x1, sep = ",", col.names = FALSE, append = TRUE)
  write.table("--------------------", x1, sep = ",", col.names = FALSE, row.names = FALSE, append = TRUE)  
  write.table("Additional performance metrics...", x1, sep = ",", col.names = FALSE, row.names = FALSE, append = TRUE)
  write.table("  on cross-validation:", x1, sep = ",", col.names = FALSE, row.names = FALSE, append = TRUE)
  write.table(rbind(names(cv_metrics), cv_metrics), x1, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
  write.table("  on holdout dataset:", x1, sep = ",", col.names = FALSE, row.names = FALSE, append = TRUE)
  write.table(rbind(names(holdout_metrics), holdout_metrics), x1, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
  write.table("--------------------", x1, sep = ",", col.names = FALSE, row.names = FALSE, append = TRUE)
  write.table("Variable importance:", x1, sep = ",", col.names = FALSE, row.names = FALSE, append = TRUE)
  write.table(out, x1, sep = ",", row.names = FALSE, append = TRUE)

  # Save fit object
  saveRDS(fit, paste(country, "_", y_hat, "_final_rf", ".rds", sep="")) 

    
#.........................................................................................
### ENDS
#.........................................................................................

  
   
  