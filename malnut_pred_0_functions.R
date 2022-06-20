#..........................................................................................
### +++++++++++++++ SMALL-AREA PREDICTION OF ACUTE MALNUTRITION BURDEN ++++++++++++++++ ###
#..........................................................................................

#..........................................................................................
## --------------------------------- FUNCTIONS FOR ANALYSIS ---------------------------- ##
#..........................................................................................

                                          # Written by Francesco Checchi, LSHTM (Nov 2020)
                                          # francesco.checchi@lshtm.ac.uk 



#.........................................................................................
### Function that calculates the proportion of a given month's days that is covered by a survey's recall period
#.........................................................................................

f_calc_days <- function(f_surveys_cov, f_df) {
  # select survey
  s <- subset(f_df, survey_id == f_surveys_cov["survey_id"])
  tm_now <- as.integer(f_surveys_cov["tm"])
  c1 <- as.integer(s["tm_recall_start"]) - tm_now
  c2 <- as.integer(s["tm_recall_end"])- tm_now
  x1 <- 0
  
  # calculate proportion of the month's days that are covered by the survey's recall period
  if (c1 > 0) { x1 <- 0.0 }
  if (c1 == 0) { x1 <- (as.integer(s["days_in_month_start"]) - as.integer(s["day_start"]) ) / 
    as.integer(s["days_in_month_start"]) }
  if (c1 < 0 & c2 > 0) { x1 <- 1.0 }
  if (c2 == 0) { x1 <- as.integer(s["day_end"]) / as.integer(s["days_in_month_end"]) }
  if (c2 < 0) { x1 <- 0.0 }  
    
  return(x1)
  }
  

#.........................................................................................
### Function that calculates and formats median, range and number of observations for any quantity, overall and by year
#.........................................................................................

f_calc_svy <- function(f_quantity, f_df, f_digits, f_years) {
    
    # output of n elements, where n = 1 (total) + number of years
    x1 <- rep(NA, times=length(f_years + 1))
  
    # overall (all years)
    x1[1] <- paste( 
      round(median(f_df[, paste(f_quantity)], na.rm=TRUE), digits=f_digits), " (", 
      round(min(f_df[, paste(f_quantity)], na.rm=TRUE), digits=f_digits), " to ",
      round(max(f_df[, paste(f_quantity)], na.rm=TRUE), digits=f_digits), ", ", 
      length(na.omit(f_df[, paste(f_quantity)]) ), ")", 
      sep="" )
    
    x1_pos <- 1
    # for each year...
    for (i in f_years) {
      x1_pos <- x1_pos + 1
      x1[x1_pos] <- paste(
        round(median(subset(f_df, year_survey == i)[, paste(f_quantity)], na.rm=TRUE), digits=f_digits), " (", 
        round(min(subset(f_df, year_survey == i)[, paste(f_quantity)], na.rm=TRUE), digits=f_digits), " to ",
        round(max(subset(f_df, year_survey == i)[, paste(f_quantity)], na.rm=TRUE), digits=f_digits), ", ", 
        length(na.omit(subset(f_df, year_survey == i)[, paste(f_quantity)]) ), ")", 
        sep="" )
      }
    return(x1)
} 


#.........................................................................................
### Function to perform K-fold cross-validation, compute predictive scores and plot results
#.........................................................................................

f_cv <- function(f_fit, f_data, f_part_unit, f_k_folds, f_scores) {
  
  # Remove missing observations from data
  f_data <- f_data[complete.cases(f_data[, all.vars(formula(f_fit)) ]), ]
  
  # Determine number of folds if f_k_folds = NA (i.e. LOOCV case)
  if (is.na(f_k_folds) ) {f_k_folds <- nrow(unique(f_data[, c(f_part_unit, "tm")])) }
    
  # Shuffle dataset
    # attribute a random value to each unique partition unit
    x1 <- unique(f_data[, c(f_part_unit, "tm")])
    x2 <- data.frame(x1, sample(c(1:nrow(x1)), nrow(x1), replace = FALSE) )
    colnames(x2) <- c(f_part_unit, "tm", "rand_rank")
    f_data <- merge(f_data, x2, by = c(f_part_unit, "tm") )
    
    # sort dataset on the basis of that random value
    f_data <- f_data[order(f_data[, "rand_rank"]), ]
    
  # Split data into K folds
    # split unique partition units into the desired number of folds
    x1 <- split(x2[, "rand_rank"], sort(x2[, "rand_rank"]%%f_k_folds))
    
    # reshape as dataframe with fold ID as factor and partition unit rank
    x2 <- data.frame(fold = rep(names(x1), sapply(x1, length)), rand_rank = unlist(x1))
    
    # merge with dataset
    f_data <- merge(f_data, x2, by="rand_rank")
    
    # now split dataset by fold as the factor
    folds <- split(f_data, f_data$fold)

  # Fit model on all the unfolded sets and track predictive scores of model fit on each fold  
    # vector to hold whether the fold fitted successfully
    folds_ok <- rep(1, length(folds))
    
    # vector to hold the partition units
    part_units <- c()
    
    # vector to hold the time units
    tms <- c()
    
    # vector to hold the observations
    obs <- c()
    
    # vector to hold the predictions
    pred <- c()
    
    # database to hold the components to calculate MSE
    scores <- data.frame()
    
  for (i in 1:length(folds) ) {	
    # control statement
    print(paste("now working on fold  ", i, " of  ", length(folds), sep=""))
    
    # fit on all data but the fold
    cv_fit <- update(f_fit, formula = formula(f_fit),  family = family(f_fit)[[1]], 
      data = do.call(rbind, folds[-i]))
    
    # calculate components of predictive scores of model when predicting fold data
    lambdas <- try(predict(cv_fit, newdata = folds[[i]], type = "response", allow.new.levels = TRUE) )
      # if fold prediction doesn't work (usually because of new levels), record this and skip to next fold
      if (class(lambdas) == "try-error") {folds_ok[i] <- 0; next}
    ys <- folds[[i]][complete.cases(folds[[i]][, all.vars(formula(f_fit))]), all.vars(formula(f_fit))[1] ]
    x3 <- cbind(ys, lambdas)
    
    # add fold results to output vectors
    part_units <- c(part_units, folds[[i]][, f_part_unit])
    tms <- c(tms, folds[[i]][, "tm"])
    obs <- c(obs, ys)
    pred <- c(pred, lambdas)
    scores <- rbind(scores , x3 )

  }

  # Aggregate data and predictions by partition units
  x1 <- as.data.frame(cbind(part_units, tms, scores))
  colnames(x1) <- c(f_part_unit, "tm", "ys", "lambdas")
  x1 <- aggregate(x1[, c("ys", "lambdas")], by = x1[, c(f_part_unit, "tm")], FUN = mean, na.rm = TRUE)
  colnames(x1) <- c(f_part_unit, "tm", "ys", "lambdas")
  
  # Return mean MSE across all folds, if desired
  if (f_scores == TRUE) {
    out <- c(mean((x1$lambdas - x1$ys )^2, na.rm = TRUE), sum(folds_ok) / length(folds) )
    if (f_part_unit %in% all.vars(formula(f_fit))) {names(out) <- c("mse_cv_glmm", "prop_folds_ok_glm")}
    if (! f_part_unit %in% all.vars(formula(f_fit))) {names(out) <- c("mse_cv_glm", "prop_folds_ok_glmm")}    
    return(out)
  }
    
  # Alternatively, return the cross-validation predictions and observations for each partition unit
  if (f_scores == FALSE) { return(x1) }
  
}
  

#.........................................................................................
### Function to generate additional metrics of fit performance, based on cross-validation results
#.........................................................................................
 
f_cv_metrics <- function(f_fit, f_data, f_part_unit, f_k_folds, f_thresholds, f_overall) {
  
  # Remove missing observations from data
  f_data <- f_data[complete.cases(f_data[, all.vars(formula(f_fit)) ]), ]
  
  # Determine number of folds if f_k_folds = NA (i.e. LOOCV case)
  if (is.na(f_k_folds) ) {f_k_folds <- nrow(unique(f_data[, c(f_part_unit, "tm")])) }
    
  # Shuffle dataset
    # attribute a random value to each unique partition unit
    x1 <- unique(f_data[, c(f_part_unit, "tm")])
    x2 <- data.frame(x1, sample(c(1:nrow(x1)), nrow(x1), replace = FALSE) )
    colnames(x2) <- c(f_part_unit, "tm", "rand_rank")
    f_data <- merge(f_data, x2, by = c(f_part_unit, "tm") )
    
    # sort dataset on the basis of that random value
    f_data <- f_data[order(f_data[, "rand_rank"]), ]
    
  # Split data into K folds
    # split unique partition units into the desired number of folds
    x1 <- split(x2[, "rand_rank"], sort(x2[, "rand_rank"]%%f_k_folds))
    
    # reshape as dataframe with fold ID as factor and partition unit rank
    x2 <- data.frame(fold = rep(names(x1), sapply(x1, length)), rand_rank = unlist(x1))
    
    # merge with dataset
    f_data <- merge(f_data, x2, by="rand_rank")
    
    # now split dataset by fold as the factor
    folds <- split(f_data, f_data$fold)
        
  # Fit model on all the unfolded sets and compute point estimates and confidence intervals for the predictions  
    
    # vector to hold whether the fold fitted successfully
    folds_ok <- rep(1, length(folds))

    # dataframe to hold results
    pred <- c()
    f_y_hat <- all.vars(formula(f_fit))[1]

  for (i in 1:length(folds) ) {	
    # control statement
    print(paste("now working on fold  ", i, " of  ", length(folds), sep=""))
    
    # fit on all data but the fold
    cv_fit <- update(f_fit, formula = formula(f_fit),  family = family(f_fit)[[1]], 
      data = do.call(rbind, folds[-i]))
    
    # partition units
    part_units <- as.character(folds[[i]][complete.cases(folds[[i]][, all.vars(formula(f_fit))]), f_part_unit ])
    tms <- folds[[i]][complete.cases(folds[[i]][, all.vars(formula(f_fit))]), "tm" ]
    
    # calculate point estimates, standard errors and bootstrap 80% and 95% confidence intervals of predictions
    lambdas <- try(predict(cv_fit, newdata = folds[[i]], type = "link", allow.new.levels = TRUE) )
      # if fold prediction doesn't work (usually because of new levels), record this and skip to next fold
      if (class(lambdas) == "try-error") {folds_ok[i] <- 0; next}
    
    ses <- unlist(predict(cv_fit, newdata = folds[[i]], se.fit = TRUE, allow.new.levels = TRUE)["se.fit"])

    x1 <- t(apply(cbind(lambdas, ses), 1, FUN = function(x) {rnorm(1000, mean = x[1], sd = x[2]) } ) ) 
    if (as.character(family(f_fit))[1] != "gaussian") {x1 <- inv.logit(x1)}
    x1 <- t(apply(x1, 1, sort) )
    x2 <- t(apply(x1, 1, quantile, c(0.025, 0.10, 0.90, 0.975)))
    x2 <- data.frame(part_units, tms, predict(cv_fit, newdata = folds[[i]], type = "response", 
      allow.new.levels = TRUE), x2)
    
    # add fold results to output vectors
    pred <- rbind(pred, x2)

  }
 
  # Aggregate predictions by partition unit
  colnames(pred) <- c(f_part_unit, "tm", "pred", "pred_l95", "pred_l80", "pred_u80", "pred_u95")
  pred <- aggregate(pred[, c("pred", "pred_l95", "pred_l80", "pred_u80", "pred_u95")], 
    by = pred[, c(f_part_unit, "tm")], FUN = mean, na.rm = TRUE)
  colnames(pred)[c(1,2)] <- c(f_part_unit, "tm")
    
  # Aggregate observed data by partition unit and compute cluster-adjusted confidence intervals
    # unique instances of partition units
    obs <- unique(pred[, c(f_part_unit, "tm")])
    
    # set up output
    obs[, c("obs", "obs_l95", "obs_l80", "obs_u80", "obs_u95")] <- NA
    
  for (i in 1:nrow(obs)) {
    # capture dataset
    df <- f_data[f_data[, f_part_unit] == obs[i, 1] & f_data[, "tm"] == obs[i, 2], ]

    # specify survey design
      
      # if all cluster values == 99999, assume SRS or exhaustive (latter to be treated as an SRS)
      if (length(unique(df[, "cluster"])) == 1) {
        survey_design <-  suppressWarnings(svydesign(id = ~0, data = df) )
      }

      # otherwise, assume cluster design
      if (length(unique(df[, "cluster"])) > 1) {
        survey_design <- suppressWarnings(svydesign(id = ~cluster, data = subset(df, ! is.na(cluster) ) ) )
      }

    # estimate point estimates and 95%CIs for continuous variables 
    if (f_y_hat %in% c("wfhz", "muac", "mfaz")) {
 
      # fit survey GLM
      fit <- svyglm(as.formula(paste(f_y_hat, "~", "NULL", sep = " ")), survey_design, family = "gaussian")
      
      # compute point estimate and CIs      
      obs[i, "obs"] <- summary(fit)$coefficients[[1]]
      obs[i, "obs_l95"] <- summary(fit)$coefficients[[1]] - 1.96 * summary(fit)$coefficients[[2]]
      obs[i, "obs_u95"] <- summary(fit)$coefficients[[1]] + 1.96 * summary(fit)$coefficients[[2]]
      obs[i, "obs_l80"] <- summary(fit)$coefficients[[1]] - 1.28 * summary(fit)$coefficients[[2]]
      obs[i, "obs_u80"] <- summary(fit)$coefficients[[1]] + 1.28 * summary(fit)$coefficients[[2]]
    }
    
    # estimate point estimates and 95%CIs for binary variables
    if (f_y_hat %in% c("gamz", "samz", "gamm", "samm")) {
 
      # fit survey GLM
      fit <- svyglm(as.formula(paste(f_y_hat, "~", "NULL", sep = " ")), survey_design, family = "binomial")
      
      # compute point estimate and CIs      
      obs[i, "obs"] <- inv.logit(summary(fit)$coefficients[[1]] )
      obs[i, "obs_l95"] <- inv.logit(summary(fit)$coefficients[[1]] - 1.96 * summary(fit)$coefficients[[2]] )
      obs[i, "obs_u95"] <- inv.logit(summary(fit)$coefficients[[1]] + 1.96 * summary(fit)$coefficients[[2]] )
      obs[i, "obs_l80"] <- inv.logit(summary(fit)$coefficients[[1]] - 1.28 * summary(fit)$coefficients[[2]] )
      obs[i, "obs_u80"] <- inv.logit(summary(fit)$coefficients[[1]] + 1.28 * summary(fit)$coefficients[[2]] )
    }
        
  }

  # Merge predictions and observations
  out <- merge(obs, pred, by = c(f_part_unit, "tm"))
    
  # Calculate performance metrics
  out[, "rel_bias"] <- (out[, "pred"] - out[, "obs"]) / out[, "obs"]
  out[, "rel_precision95"] <- (abs(out[, "pred_u95"] - out[, "pred_l95"] ) / 2) / abs(out[, "pred"])
  out[, "coverage95"] <- data.table::between(out[, "pred"], out[, "obs_l95"], out[, "obs_u95"])
  out[, "coverage80"] <- data.table::between(out[, "pred"], out[, "obs_l80"], out[, "obs_u80"])
  out[, "denom"] <- ifelse(is.na(out[, "pred"]), FALSE, TRUE)
  
  for (i in f_thresholds) {
    out[, paste("sensitivity_", i, sep = "")] <- NA
    out[, paste("denom_sens_", i, sep = "")] <- NA
    if (f_y_hat %in% c("gamz", "samz", "gamm", "samm")) {    
      out[, paste("sensitivity_", i, sep = "")] <- ifelse(out[, "obs"] >= i & out[, "pred"] >= i, TRUE, 
        out[, paste("sensitivity_", i, sep = "")] )
      out[, paste("sensitivity_", i, sep = "")] <- ifelse(out[, "obs"] >= i & out[, "pred"] < i, FALSE, 
        out[, paste("sensitivity_", i, sep = "")] )
      out[, paste("denom_sens_", i, sep = "")] <- ifelse(out[, "obs"] >= i, TRUE, 
        out[, paste("denom_sens_", i, sep = "")] )
      
      out[, paste("specificity_", i, sep = "")] <- NA
      out[, paste("specificity_", i, sep = "")] <- ifelse(out[, "obs"] < i & out[, "pred"] < i, TRUE, 
        out[, paste("specificity_", i, sep = "")] )
      out[, paste("specificity_", i, sep = "")] <- ifelse(out[, "obs"] < i & out[, "pred"] >= i, FALSE, 
        out[, paste("specificity_", i, sep = "")] )
      out[, paste("denom_spec_", i, sep = "")] <- NA
      out[, paste("denom_spec_", i, sep = "")] <- ifelse(out[, "obs"] < i, TRUE, 
        out[, paste("denom_spec_", i, sep = "")] )
    }
  }
    

  # Return performance metrics
    # by partition unit if desired...
    if (f_overall == FALSE) {return(out)}
  
    # ...or overall otherwise
    if (f_overall == TRUE) {
      x2 <- c(
        round(100 * colMeans(out[, c("rel_bias", "rel_precision95", "coverage95", "coverage80", 
          grep("sensitivity", colnames(out), value = TRUE), grep("specificity", colnames(out), value = TRUE))], 
          na.rm = TRUE), 1),
        as.integer(colSums(out[, grep("denom", colnames(out), value = TRUE)], na.rm = TRUE)),
        round(100 * sum(folds_ok) / length(folds), 1)
      )
      names(x2) <- c(
        "rel_bias", "rel_precision95", "coverage95", "coverage80", 
        grep("sensitivity", colnames(out), value = TRUE), grep("specificity", colnames(out), value = TRUE),
        grep("denom", colnames(out), value = TRUE),
        "prop_folds_ok"
      )
      return(x2)
    }
}


#.........................................................................................
### Function to generate additional metrics of random forest performance, based on cross-validation results
#.........................................................................................
 
f_cv_metrics_rf <- function(f_fit, f_data_agg, f_data, f_part_unit, f_k_folds, f_thresholds, f_overall) {
  
  # Remove missing observations from data (not aggregated)
  f_data_raw <- f_data[complete.cases(f_data[, all.vars(formula(f_fit)) ]), ]
  
  # Choose between aggregated and not aggregated dataset for prediction
  if (is.null(f_data_agg)) {f_data <- f_data_raw}
  if (! is.null(f_data_agg)) {f_data <- f_data_agg}
  
  # Determine number of folds if f_k_folds = NA (i.e. LOOCV case)
  if (is.na(f_k_folds) ) {f_k_folds <- nrow(unique(f_data[, c(f_part_unit, "tm")])) }
    
  # Shuffle dataset
    # attribute a random value to each unique partition unit
    x1 <- unique(f_data[, c(f_part_unit, "tm")])
    x2 <- data.frame(x1, sample(c(1:nrow(x1)), nrow(x1), replace = FALSE) )
    colnames(x2) <- c(f_part_unit, "tm", "rand_rank")
    f_data <- merge(f_data, x2, by = c(f_part_unit, "tm") )
    
    # sort dataset on the basis of that random value
    f_data <- f_data[order(f_data[, "rand_rank"]), ]
    
  # Split data into K folds
    # split unique partition units into the desired number of folds
    x1 <- split(x2[, "rand_rank"], sort(x2[, "rand_rank"]%%f_k_folds))
    
    # reshape as dataframe with fold ID as factor and partition unit rank
    x2 <- data.frame(fold = rep(names(x1), sapply(x1, length)), rand_rank = unlist(x1))
    
    # merge with dataset
    f_data <- merge(f_data, x2, by="rand_rank")
    
    # now split dataset by fold as the factor
    folds <- split(f_data, f_data$fold)
        
  # Fit model on all the unfolded sets and compute point estimates and confidence intervals for the predictions  
    
    # vector to hold whether the fold fitted successfully
    folds_ok <- rep(1, length(folds))

    # dataframe to hold results
    pred <- c()
    f_y_hat <- all.vars(formula(f_fit))[1]

  for (i in 1:length(folds) ) {	
    # control statement
    print(paste("now working on fold  ", i, " of  ", length(folds), sep=""))
    
    # fit on all data but the fold
    #cv_fit <- update(f_fit, formula = formula(f_fit), data = do.call(rbind, folds[-i]))
    cv_fit <- ranger(formula = formula(f_fit), data = do.call(rbind, folds[-i]), case.weights = "wt", 
      keep.inbag = TRUE)
    
    # partition units
    part_units <- as.character(folds[[i]][complete.cases(folds[[i]][, all.vars(formula(f_fit))]), f_part_unit ])
    tms <- folds[[i]][complete.cases(folds[[i]][, all.vars(formula(f_fit))]), "tm" ]
    
    # calculate point estimates, standard errors and jackknife 80% and 95% confidence intervals of predictions
    pred_i <- try(predict(cv_fit, data = folds[[i]], type = "se", se.method = "jack") )
      # if fold prediction doesn't work (usually because of new levels), record this and skip to next fold
      if (class(pred_i) == "try-error") {folds_ok[i] <- 0; next}
    x2 <- c(
      part_units,
      tms,
      pred_i$predictions, 
      pred_i$predictions - pred_i$se * 1.96,
      pred_i$predictions - pred_i$se * 1.28,
      pred_i$predictions + pred_i$se * 1.28,
      pred_i$predictions + pred_i$se * 1.96
    )
    # add fold results to output vectors
    pred <- rbind(pred, unlist(x2) )

  }
 
  # Aggregate predictions by partition unit
  pred <- as.data.frame(pred)
  colnames(pred) <- c(f_part_unit, "tm", "pred", "pred_l95", "pred_l80", "pred_u80", "pred_u95")
  pred[, 2:ncol(pred)] <- lapply(pred[, 2:ncol(pred)], as.numeric)
  pred <- aggregate(pred[, c("pred", "pred_l95", "pred_l80", "pred_u80", "pred_u95")], 
    by = pred[, c(f_part_unit, "tm")], FUN = mean, na.rm = TRUE)
  colnames(pred)[c(1,2)] <- c(f_part_unit, "tm")
    
  # Aggregate observed data by partition unit and compute cluster-adjusted confidence intervals
    # unique instances of partition units
    obs <- unique(pred[, c(f_part_unit, "tm")])
    
    # set up output
    obs[, c("obs", "obs_l95", "obs_l80", "obs_u80", "obs_u95")] <- NA
    
  for (i in 1:nrow(obs)) {
    # capture dataset
    df <- f_data_raw[f_data_raw[, f_part_unit] == obs[i, 1] & f_data_raw[, "tm"] == obs[i, 2], ]

    # specify survey design
      
      # if all cluster values == 99999, assume SRS or exhaustive (latter to be treated as an SRS)
      if (length(unique(df[, "cluster"])) == 1) {
        survey_design <-  suppressWarnings(svydesign(id = ~0, data = df) )
      }

      # otherwise, assume cluster design
      if (length(unique(df[, "cluster"])) > 1) {
        survey_design <- suppressWarnings(svydesign(id = ~cluster, data = subset(df, ! is.na(cluster) ) ) )
      }

    # estimate point estimates and 95%CIs for continuous variables 
    if (f_y_hat %in% c("wfhz", "muac", "mfaz")) {
 
      # fit survey GLM
      fit <- svyglm(as.formula(paste(f_y_hat, "~", "NULL", sep = " ")), survey_design, family = "gaussian")
      
      # compute point estimate and CIs      
      obs[i, "obs"] <- summary(fit)$coefficients[[1]]
      obs[i, "obs_l95"] <- summary(fit)$coefficients[[1]] - 1.96 * summary(fit)$coefficients[[2]]
      obs[i, "obs_u95"] <- summary(fit)$coefficients[[1]] + 1.96 * summary(fit)$coefficients[[2]]
      obs[i, "obs_l80"] <- summary(fit)$coefficients[[1]] - 1.28 * summary(fit)$coefficients[[2]]
      obs[i, "obs_u80"] <- summary(fit)$coefficients[[1]] + 1.28 * summary(fit)$coefficients[[2]]
    }
    
    # estimate point estimates and 95%CIs for binary variables
    if (f_y_hat %in% c("gamz", "samz", "gamm", "samm")) {
 
      # fit survey GLM
      fit <- svyglm(as.formula(paste(f_y_hat, "~", "NULL", sep = " ")), survey_design, family = "binomial")
      
      # compute point estimate and CIs      
      obs[i, "obs"] <- inv.logit(summary(fit)$coefficients[[1]] )
      obs[i, "obs_l95"] <- inv.logit(summary(fit)$coefficients[[1]] - 1.96 * summary(fit)$coefficients[[2]] )
      obs[i, "obs_u95"] <- inv.logit(summary(fit)$coefficients[[1]] + 1.96 * summary(fit)$coefficients[[2]] )
      obs[i, "obs_l80"] <- inv.logit(summary(fit)$coefficients[[1]] - 1.28 * summary(fit)$coefficients[[2]] )
      obs[i, "obs_u80"] <- inv.logit(summary(fit)$coefficients[[1]] + 1.28 * summary(fit)$coefficients[[2]] )
    }
        
  }

  # Merge predictions and observations
  out <- merge(obs, pred, by = c(f_part_unit, "tm"))
    
  # Calculate performance metrics
  out[, "rel_bias"] <- (out[, "pred"] - out[, "obs"]) / out[, "obs"]
  out[, "rel_precision95"] <- (abs(out[, "pred_u95"] - out[, "pred_l95"] ) / 2) / abs(out[, "pred"])
  out[, "coverage95"] <- data.table::between(out[, "pred"], out[, "obs_l95"], out[, "obs_u95"])
  out[, "coverage80"] <- data.table::between(out[, "pred"], out[, "obs_l80"], out[, "obs_u80"])
  out[, "denom"] <- ifelse(is.na(out[, "pred"]), FALSE, TRUE)
  
  for (i in f_thresholds) {
    out[, paste("sensitivity_", i, sep = "")] <- NA
    out[, paste("denom_sens_", i, sep = "")] <- NA
    if (f_y_hat %in% c("gamz", "samz", "gamm", "samm")) {    
      out[, paste("sensitivity_", i, sep = "")] <- ifelse(out[, "obs"] >= i & out[, "pred"] >= i, TRUE, 
        out[, paste("sensitivity_", i, sep = "")] )
      out[, paste("sensitivity_", i, sep = "")] <- ifelse(out[, "obs"] >= i & out[, "pred"] < i, FALSE, 
        out[, paste("sensitivity_", i, sep = "")] )
      out[, paste("denom_sens_", i, sep = "")] <- ifelse(out[, "obs"] >= i, TRUE, 
        out[, paste("denom_sens_", i, sep = "")] )
      
      out[, paste("specificity_", i, sep = "")] <- NA
      out[, paste("specificity_", i, sep = "")] <- ifelse(out[, "obs"] < i & out[, "pred"] < i, TRUE, 
        out[, paste("specificity_", i, sep = "")] )
      out[, paste("specificity_", i, sep = "")] <- ifelse(out[, "obs"] < i & out[, "pred"] >= i, FALSE, 
        out[, paste("specificity_", i, sep = "")] )
      out[, paste("denom_spec_", i, sep = "")] <- NA
      out[, paste("denom_spec_", i, sep = "")] <- ifelse(out[, "obs"] < i, TRUE, 
        out[, paste("denom_spec_", i, sep = "")] )
    }
  }
    
  # Return performance metrics
    # by partition unit if desired...
    if (f_overall == FALSE) {return(out)}
  
    # ...or overall otherwise
    if (f_overall == TRUE) {
      x2 <- c(
        round(100 * colMeans(out[, c("rel_bias", "rel_precision95", "coverage95", "coverage80", 
          grep("sensitivity", colnames(out), value = TRUE), grep("specificity", colnames(out), value = TRUE))], 
          na.rm = TRUE), 1),
        as.integer(colSums(out[, grep("denom", colnames(out), value = TRUE)], na.rm = TRUE)),
        round(100 * sum(folds_ok) / length(folds), 1)
      )
      names(x2) <- c(
        "rel_bias", "rel_precision95", "coverage95", "coverage80", 
        grep("sensitivity", colnames(out), value = TRUE), grep("specificity", colnames(out), value = TRUE),
        grep("denom", colnames(out), value = TRUE),
        "prop_folds_ok"
      )
      return(x2)
    }
}


#.........................................................................................
### Function to perform cross-validation of a random forest model and return predictions
#.........................................................................................
 
f_cv_rf <- function(f_fit, f_data, f_part_unit, f_k_folds) {
  
  # Determine number of folds if f_k_folds = NA (i.e. LOOCV case)
  if (is.na(f_k_folds) ) {f_k_folds <- nrow(unique(f_data[, c(f_part_unit, "tm")])) }
    
  # Shuffle dataset
    # attribute a random value to each unique partition unit
    x1 <- unique(f_data[, c(f_part_unit, "tm")])
    x2 <- data.frame(x1, sample(c(1:nrow(x1)), nrow(x1), replace = FALSE) )
    colnames(x2) <- c(f_part_unit, "tm", "rand_rank")
    f_data <- merge(f_data, x2, by = c(f_part_unit, "tm") )
    
    # sort dataset on the basis of that random value
    f_data <- f_data[order(f_data[, "rand_rank"]), ]
    
  # Split data into K folds
    # split unique partition units into the desired number of folds
    x1 <- split(x2[, "rand_rank"], sort(x2[, "rand_rank"]%%f_k_folds))
    
    # reshape as dataframe with fold ID as factor and partition unit rank
    x2 <- data.frame(fold = rep(names(x1), sapply(x1, length)), rand_rank = unlist(x1))
    
    # merge with dataset
    f_data <- merge(f_data, x2, by="rand_rank")
    
    # now split dataset by fold as the factor
    folds <- split(f_data, f_data$fold)
        
  # Fit model on all the unfolded sets and compute point estimates for the predictions  
    
    # vector to hold whether the fold fitted successfully
    folds_ok <- rep(1, length(folds))

    # dataframe to hold results
    pred <- c()
    f_y_hat <- all.vars(formula(f_fit))[1]

  for (i in 1:length(folds) ) {	
    # control statement
    print(paste("now working on fold  ", i, " of  ", length(folds), sep=""))
    
    # fit on all data but the fold
    #cv_fit <- update(f_fit, formula = formula(f_fit), data = do.call(rbind, folds[-i]))
    cv_fit <- ranger(formula = formula(f_fit), data = do.call(rbind, folds[-i]), case.weights = "wt", 
      keep.inbag = TRUE)
    
    # partition units
    part_units <- as.character(folds[[i]][complete.cases(folds[[i]][, all.vars(formula(f_fit))]), f_part_unit ])
    tms <- folds[[i]][complete.cases(folds[[i]][, all.vars(formula(f_fit))]), "tm" ]
    
    # observations
    ys <- folds[[i]][complete.cases(folds[[i]][, all.vars(formula(f_fit))]), f_y_hat ]
    
    # calculate point estimates of predictions
    pred_i <- try(predict(cv_fit, data = folds[[i]]) )
      # if fold prediction doesn't work (usually because of new levels), record this and skip to next fold
      if (class(pred_i) == "try-error") {folds_ok[i] <- 0; next}
    x2 <- c(
      part_units,
      tms,
      ys,
      pred_i$predictions
    )
    # add fold results to output vectors
    pred <- rbind(pred, unlist(x2) )

  }
 
  # Aggregate predictions by partition unit
  pred <- as.data.frame(pred)
  colnames(pred) <- c(f_part_unit, "tm", "observed", "predicted")
  pred[, 2:ncol(pred)] <- lapply(pred[, 2:ncol(pred)], as.numeric)
  pred <- aggregate(pred[, c("observed", "predicted")], 
    by = pred[, c(f_part_unit, "tm")], FUN = mean, na.rm = TRUE)
  colnames(pred) <- c(f_part_unit, "tm", "observed", "predicted")
    
  # Return output
  return(pred)
}



#.........................................................................................
### Function to plot histograms of variables
#.........................................................................................  

f_hist <- function(f_var, f_data, f_lims) {
    
  plot <- ggplot(f_data)
      
    # if the variable has >= 20 unique values...
      if (length(unique(na.omit(f_data[, f_var]))) >= 20) {
        plot <- plot + geom_histogram(aes(x = as.numeric(f_data[, f_var]) ), 
          color="seagreen", fill="seagreen3", alpha = 0.5 ) +
          theme_bw() + xlab(f_var) + scale_x_continuous(expand = c(0, 0), limits = f_lims )
      }
 
    # otherwise...
      if (length(unique(na.omit(f_data[, f_var]))) < 20) {
        plot <- plot + geom_histogram(aes(x = as.factor(f_data[, f_var]) ), stat="count", 
          color="seagreen", fill="seagreen3", alpha = 0.5) +
          theme_bw() + xlab(f_var)
      }
        
    print(plot)
  }


#.........................................................................................
### Function to generate additional metrics of fit performance, based on predicting on holdout data
#.........................................................................................
 
f_holdout_metrics <- function(f_fit, f_data_holdout, f_part_unit, f_thresholds, f_overall) {
  
  # Remove missing observations from data
  f_data_holdout <- f_data_holdout[complete.cases(f_data_holdout[, all.vars(formula(f_fit)) ]), ]
  
  # Predict on holdout data
  lambdas <- try(predict(f_fit, newdata = f_data_holdout, type = "link", allow.new.levels = TRUE) )

  # Calculate 80% and 95% confidence intervals of predictions
  if (f_part_unit %in% all.vars(formula(f_fit)) ) {
    lambdas <- try(predict(f_fit, newdata = f_data_holdout, type = "link", re.form = NA, allow.new.levels = TRUE) )
    mm <- model.matrix(terms(f_fit), f_data_holdout)
    pvar1 <- diag(mm %*% tcrossprod(vcov(f_fit), mm))
    tvar1 <- pvar1 + as.numeric(VarCorr(f_fit)[f_part_unit])
    x2 <- data.frame(
      lambdas - 1.96 * sqrt(tvar1), 
      lambdas - 1.28 * sqrt(tvar1), 
      lambdas + 1.28 * sqrt(tvar1), 
      lambdas + 1.96 * sqrt(tvar1)
    )
    if (as.character(family(f_fit))[1] != "gaussian") {x2 <- inv.logit(x2)}
  }
  
  if (! f_part_unit %in% all.vars(formula(f_fit)) ) {
    ses <- unlist(predict(f_fit, newdata = f_data_holdout, se.fit = TRUE, allow.new.levels = TRUE)["se.fit"])
    x1 <- t(apply(cbind(lambdas, ses), 1, FUN = function(x) {rnorm(1000, mean = x[1], sd = x[2]) } ) ) 
    if (as.character(family(f_fit))[1] != "gaussian") {x1 <- inv.logit(x1)}
    x1 <- t(apply(x1, 1, sort) )
    x2 <- t(apply(x1, 1, quantile, c(0.025, 0.10, 0.90, 0.975)))
  }
  
  # Aggregate predictions by partition unit
  pred <- as.data.frame(cbind(f_data_holdout[, c(f_part_unit, "tm")], lambdas, x2) )
  colnames(pred) <- c(f_part_unit, "tm", "pred", "pred_l95", "pred_l80", "pred_u80", "pred_u95")
  pred <- aggregate(pred[, c("pred", "pred_l95", "pred_l80", "pred_u80", "pred_u95")], 
    by = pred[, c(f_part_unit, "tm")], FUN = mean, na.rm = TRUE)
  colnames(pred)[c(1,2)] <- c(f_part_unit, "tm")
  if (as.character(family(f_fit))[1] != "gaussian") {pred[, "pred"] <- inv.logit(pred[, "pred"])}
  
  # Aggregate observed data by partition unit and compute its cluster-adjusted confidence intervals
    # unique instances of partition units
    obs <- unique(f_data_holdout[, c(f_part_unit, "tm")])
    
    # set up output
    obs[, c("obs", "obs_l95", "obs_l80", "obs_u80", "obs_u95")] <- NA
    f_y_hat <- all.vars(formula(f_fit))[1]
    
  for (i in 1:nrow(obs)) {
    # capture dataset
    df <- f_data_holdout[f_data_holdout[, f_part_unit] == obs[i, 1] & f_data_holdout[, "tm"] == obs[i, 2], ]
  
    # specify survey design
      
      # if all cluster values == 99999, assume SRS or exhaustive (latter to be treated as an SRS)
      if (length(unique(df[, "cluster"])) == 1) {
        survey_design <-  suppressWarnings(svydesign(id = ~0, data = df) )
      }

      # otherwise, assume cluster design
      if (length(unique(df[, "cluster"])) > 1) {
      survey_design <- suppressWarnings(
        svydesign(id = ~cluster, data = subset(df, ! is.na(cluster) ) ) )
      }

    # estimate point estimates and 95%CIs for continuous variables 
    if (f_y_hat %in% c("wfhz", "muac", "mfaz")) {
 
      # fit survey GLM
      fit <- svyglm(as.formula(paste(f_y_hat, "~", "NULL", sep = " ")), survey_design, family = "gaussian")
      
      # compute point estimate and CIs      
      obs[i, "obs"] <- summary(fit)$coefficients[[1]]
      obs[i, "obs_l95"] <- summary(fit)$coefficients[[1]] - 1.96 * summary(fit)$coefficients[[2]]
      obs[i, "obs_u95"] <- summary(fit)$coefficients[[1]] + 1.96 * summary(fit)$coefficients[[2]]
      obs[i, "obs_l80"] <- summary(fit)$coefficients[[1]] - 1.28 * summary(fit)$coefficients[[2]]
      obs[i, "obs_u80"] <- summary(fit)$coefficients[[1]] + 1.28 * summary(fit)$coefficients[[2]]
    }
    
    # estimate point estimates and 95%CIs for binary variables
    if (f_y_hat %in% c("gamz", "samz", "gamm", "samm")) {
 
      # fit survey GLM
      fit <- svyglm(as.formula(paste(f_y_hat, "~", "NULL", sep = " ")), survey_design, family = "binomial")
      
      # compute point estimate and CIs      
      obs[i, "obs"] <- inv.logit(summary(fit)$coefficients[[1]] )
      obs[i, "obs_l95"] <- inv.logit(summary(fit)$coefficients[[1]] - 1.96 * summary(fit)$coefficients[[2]] )
      obs[i, "obs_u95"] <- inv.logit(summary(fit)$coefficients[[1]] + 1.96 * summary(fit)$coefficients[[2]] )
      obs[i, "obs_l80"] <- inv.logit(summary(fit)$coefficients[[1]] - 1.28 * summary(fit)$coefficients[[2]] )
      obs[i, "obs_u80"] <- inv.logit(summary(fit)$coefficients[[1]] + 1.28 * summary(fit)$coefficients[[2]] )
    }
        
  }

  # Merge predictions and observations
  out <- merge(obs, pred, by = c(f_part_unit, "tm"))
  
  # Compute performance metrics
  out[, "rel_bias"] <- (out[, "pred"] - out[, "obs"]) / out[, "obs"]
  out[, "rel_precision95"] <- (abs(out[, "pred_u95"] - out[, "pred_l95"] ) / 2) / abs(out[, "pred"])
  out[, "coverage95"] <- data.table::between(out[, "pred"], out[, "obs_l95"], out[, "obs_u95"])
  out[, "coverage80"] <- data.table::between(out[, "pred"], out[, "obs_l80"], out[, "obs_u80"])
  out[, "denom"] <- ifelse(is.na(out[, "pred"]), FALSE, TRUE)
  
  for (i in f_thresholds) {
    out[, paste("sensitivity_", i, sep = "")] <- NA
    out[, paste("denom_sens_", i, sep = "")] <- NA
    if (f_y_hat %in% c("gamz", "samz", "gamm", "samm")) {    
      out[, paste("sensitivity_", i, sep = "")] <- ifelse(out[, "obs"] >= i & out[, "pred"] >= i, TRUE, 
        out[, paste("sensitivity_", i, sep = "")] )
      out[, paste("sensitivity_", i, sep = "")] <- ifelse(out[, "obs"] >= i & out[, "pred"] < i, FALSE, 
        out[, paste("sensitivity_", i, sep = "")] )
      out[, paste("denom_sens_", i, sep = "")] <- ifelse(out[, "obs"] >= i, TRUE, 
        out[, paste("denom_sens_", i, sep = "")] )
      
      out[, paste("specificity_", i, sep = "")] <- NA
      out[, paste("specificity_", i, sep = "")] <- ifelse(out[, "obs"] < i & out[, "pred"] < i, TRUE, 
        out[, paste("specificity_", i, sep = "")] )
      out[, paste("specificity_", i, sep = "")] <- ifelse(out[, "obs"] < i & out[, "pred"] >= i, FALSE, 
        out[, paste("specificity_", i, sep = "")] )
      out[, paste("denom_spec_", i, sep = "")] <- NA
      out[, paste("denom_spec_", i, sep = "")] <- ifelse(out[, "obs"] < i, TRUE, 
        out[, paste("denom_spec_", i, sep = "")] )
    }
  }
    
  # Return performance metrics
    # by partition unit if desired...
    if (f_overall == FALSE) {return(out)}
  
    # ...or overall otherwise
    if (f_overall == TRUE) {
      x2 <- c(
        round(100 * colMeans(out[, c("rel_bias", "rel_precision95", "coverage95", "coverage80", 
          grep("sensitivity", colnames(out), value = TRUE), grep("specificity", colnames(out), value = TRUE))], 
          na.rm = TRUE), 1),
        as.integer(colSums(out[, grep("denom", colnames(out), value = TRUE)], na.rm = TRUE))
      )
      names(x2) <- c(
        "rel_bias", "rel_precision95", "coverage95", "coverage80", 
        grep("sensitivity", colnames(out), value = TRUE), grep("specificity", colnames(out), value = TRUE),
        grep("denom", colnames(out), value = TRUE) )
      return(x2)
    }
}


#.........................................................................................
### Function to generate additional metrics of random forest performance, based on predicting on holdout data
#.........................................................................................
 
f_holdout_metrics_rf <- function(f_fit, f_data_holdout_agg, f_data_holdout, f_part_unit, f_thresholds, f_overall) {
  
  # Remove missing observations from data
  f_data_holdout <- f_data_holdout[complete.cases(f_data_holdout[, all.vars(formula(f_fit)) ]), ]
  
  # Predict on holdout data
  if (! is.null(f_data_holdout_agg)) {
    x1 <- try(predict(fit, data = f_data_holdout_agg, type = "se" ) )
  }
  if (is.null(f_data_holdout_agg)) {
    x1 <- try(predict(fit, data = f_data_holdout, type = "se" ) )
  }

  # Calculate 80% and 95% confidence intervals of predictions
  if (! is.null(f_data_holdout_agg))  {
    pred <- data.frame(f_data_holdout_agg[, c(part_unit, "tm")], 
      x1$predictions, 
      x1$predictions - x1$se * 1.96,
      x1$predictions - x1$se * 1.28,
      x1$predictions + x1$se * 1.28,
      x1$predictions + x1$se * 1.96
    )
    colnames(pred) <- c(f_part_unit, "tm", "pred", "pred_l95", "pred_l80", "pred_u80", "pred_u95")
  }
  
  if (is.null(f_data_holdout_agg)) {
    pred <- data.frame(f_data_holdout[, c(part_unit, "tm")], 
      x1$predictions, 
      x1$predictions - x1$se * 1.96,
      x1$predictions - x1$se * 1.28,
      x1$predictions + x1$se * 1.28,
      x1$predictions + x1$se * 1.96
    )
    colnames(pred) <- c(f_part_unit, "tm", "pred", "pred_l95", "pred_l80", "pred_u80", "pred_u95")
  }
  
  # Aggregate predictions by partition unit
  pred <- aggregate(pred[, c("pred", "pred_l95", "pred_l80", "pred_u80", "pred_u95")], 
    by = pred[, c(f_part_unit, "tm")], FUN = mean, na.rm = TRUE)
  colnames(pred)[c(1,2)] <- c(f_part_unit, "tm")
  
  # Aggregate observed data by partition unit and compute its cluster-adjusted confidence intervals
    # unique instances of partition units
    obs <- unique(f_data_holdout[, c(f_part_unit, "tm")])
    
    # set up output
    obs[, c("obs", "obs_l95", "obs_l80", "obs_u80", "obs_u95")] <- NA
    f_y_hat <- all.vars(formula(f_fit))[1]
    
  for (i in 1:nrow(obs)) {
    # capture dataset
    df <- f_data_holdout[f_data_holdout[, f_part_unit] == obs[i, 1] & f_data_holdout[, "tm"] == obs[i, 2], ]
  
    # specify survey design
      
      # if all cluster values == 99999, assume SRS or exhaustive (latter to be treated as an SRS)
      if (length(unique(df[, "cluster"])) == 1) {
        survey_design <-  suppressWarnings(svydesign(id = ~0, data = df) )
      }

      # otherwise, assume cluster design
      if (length(unique(df[, "cluster"])) > 1) {
      survey_design <- suppressWarnings(
        svydesign(id = ~cluster, data = subset(df, ! is.na(cluster) ) ) )
      }

    # estimate point estimates and 95%CIs for continuous variables 
    if (f_y_hat %in% c("wfhz", "muac", "mfaz")) {
 
      # fit survey GLM
      fit <- svyglm(as.formula(paste(f_y_hat, "~", "NULL", sep = " ")), survey_design, family = "gaussian")
      
      # compute point estimate and CIs      
      obs[i, "obs"] <- summary(fit)$coefficients[[1]]
      obs[i, "obs_l95"] <- summary(fit)$coefficients[[1]] - 1.96 * summary(fit)$coefficients[[2]]
      obs[i, "obs_u95"] <- summary(fit)$coefficients[[1]] + 1.96 * summary(fit)$coefficients[[2]]
      obs[i, "obs_l80"] <- summary(fit)$coefficients[[1]] - 1.28 * summary(fit)$coefficients[[2]]
      obs[i, "obs_u80"] <- summary(fit)$coefficients[[1]] + 1.28 * summary(fit)$coefficients[[2]]
    }
    
    # estimate point estimates and 95%CIs for binary variables
    if (f_y_hat %in% c("gamz", "samz", "gamm", "samm")) {
 
      # fit survey GLM
      fit <- svyglm(as.formula(paste(f_y_hat, "~", "NULL", sep = " ")), survey_design, family = "binomial")
      
      # compute point estimate and CIs      
      obs[i, "obs"] <- inv.logit(summary(fit)$coefficients[[1]] )
      obs[i, "obs_l95"] <- inv.logit(summary(fit)$coefficients[[1]] - 1.96 * summary(fit)$coefficients[[2]] )
      obs[i, "obs_u95"] <- inv.logit(summary(fit)$coefficients[[1]] + 1.96 * summary(fit)$coefficients[[2]] )
      obs[i, "obs_l80"] <- inv.logit(summary(fit)$coefficients[[1]] - 1.28 * summary(fit)$coefficients[[2]] )
      obs[i, "obs_u80"] <- inv.logit(summary(fit)$coefficients[[1]] + 1.28 * summary(fit)$coefficients[[2]] )
    }
        
  }

  # Merge predictions and observations
  out <- merge(obs, pred, by = c(f_part_unit, "tm"))
  
  # Compute performance metrics
  out[, "rel_bias"] <- (out[, "pred"] - out[, "obs"]) / out[, "obs"]
  out[, "rel_precision95"] <- (abs(out[, "pred_u95"] - out[, "pred_l95"] ) / 2) / abs(out[, "pred"])
  out[, "coverage95"] <- data.table::between(out[, "pred"], out[, "obs_l95"], out[, "obs_u95"])
  out[, "coverage80"] <- data.table::between(out[, "pred"], out[, "obs_l80"], out[, "obs_u80"])
  out[, "denom"] <- ifelse(is.na(out[, "pred"]), FALSE, TRUE)
  
  for (i in f_thresholds) {
    out[, paste("sensitivity_", i, sep = "")] <- NA
    out[, paste("denom_sens_", i, sep = "")] <- NA
    if (f_y_hat %in% c("gamz", "samz", "gamm", "samm")) {    
      out[, paste("sensitivity_", i, sep = "")] <- ifelse(out[, "obs"] >= i & out[, "pred"] >= i, TRUE, 
        out[, paste("sensitivity_", i, sep = "")] )
      out[, paste("sensitivity_", i, sep = "")] <- ifelse(out[, "obs"] >= i & out[, "pred"] < i, FALSE, 
        out[, paste("sensitivity_", i, sep = "")] )
      out[, paste("denom_sens_", i, sep = "")] <- ifelse(out[, "obs"] >= i, TRUE, 
        out[, paste("denom_sens_", i, sep = "")] )
      
      out[, paste("specificity_", i, sep = "")] <- NA
      out[, paste("specificity_", i, sep = "")] <- ifelse(out[, "obs"] < i & out[, "pred"] < i, TRUE, 
        out[, paste("specificity_", i, sep = "")] )
      out[, paste("specificity_", i, sep = "")] <- ifelse(out[, "obs"] < i & out[, "pred"] >= i, FALSE, 
        out[, paste("specificity_", i, sep = "")] )
      out[, paste("denom_spec_", i, sep = "")] <- NA
      out[, paste("denom_spec_", i, sep = "")] <- ifelse(out[, "obs"] < i, TRUE, 
        out[, paste("denom_spec_", i, sep = "")] )
    }
  }
    
  # Return performance metrics
    # by partition unit if desired...
    if (f_overall == FALSE) {return(out)}
  
    # ...or overall otherwise
    if (f_overall == TRUE) {
      x2 <- c(
        round(100 * colMeans(out[, c("rel_bias", "rel_precision95", "coverage95", "coverage80", 
          grep("sensitivity", colnames(out), value = TRUE), grep("specificity", colnames(out), value = TRUE))], 
          na.rm = TRUE), 1),
        as.integer(colSums(out[, grep("denom", colnames(out), value = TRUE)], na.rm = TRUE))
      )
      names(x2) <- c(
        "rel_bias", "rel_precision95", "coverage95", "coverage80", 
        grep("sensitivity", colnames(out), value = TRUE), grep("specificity", colnames(out), value = TRUE),
        grep("denom", colnames(out), value = TRUE) )
      return(x2)
    }
}




#.........................................................................................
### Function to standardise livelihood type nomenclature
#.........................................................................................   

f_liv <- function(f_ts, f_livelihood_substrings) {
  # Agriculturalists
  if (length(grep(paste(f_livelihood_substrings$agriculturalists, collapse="|"), f_ts )) > 0 ) 
    return( paste(names(livelihood_substrings)[2]) )
  
  # Pastoralists (slightly different to avoid confusion with agropastoralists)
  if (length(grep(paste(f_livelihood_substrings$pastoralists, collapse="|"), f_ts )) > 0  & 
    length(grep(paste(f_livelihood_substrings$agropastoralists, collapse="|"), f_ts )) == 0
    ) 
    return( paste(names(livelihood_substrings)[3]) )
  
  # Agropastoralists       
  if (length(grep(paste(f_livelihood_substrings$agropastoralists, collapse="|"), f_ts )) > 0 ) 
    return( paste(names(livelihood_substrings)[4]) )
  
  # Riverine        
  if (length(grep(paste(f_livelihood_substrings$riverine, collapse="|"), f_ts )) > 0 ) 
    return( paste(names(livelihood_substrings)[5]) )
  
  # Fishing       
  if (length(grep(paste(f_livelihood_substrings$fishing, collapse="|"), f_ts )) > 0 ) 
    return( paste(names(livelihood_substrings)[6]) )
  
  # Urban        
  if (length(grep(paste(f_livelihood_substrings$urban, collapse="|"), f_ts )) > 0 ) 
    return( paste(names(livelihood_substrings)[7]) )
  
  # Displaced        
  if (length(grep(paste(f_livelihood_substrings$displaced, collapse="|"), f_ts )) > 0 ) 
    return( paste(names(livelihood_substrings)[8]) )
  
  # Refugee       
  if (length(grep(paste(f_livelihood_substrings$refugee, collapse="|"), f_ts )) > 0 ) 
    return( paste(names(livelihood_substrings)[9]) )
  
  # Otherwise return NA
  if (length(grep(paste(unlist(f_livelihood_substrings), collapse="|"), f_ts )) == 0 ) 
    return( NA )
    
}      


#.........................................................................................
### Function to graph prevalence predictions versus observations (on training data, on CV, on holdout)
  # each of the input data frames should contained an 'observed' and 'predicted' column
  # f_model_type = "glm" for GLM, "rf" for random forest1
#.........................................................................................

f_plot_prev <- function(f_country, f_y_hat, f_training, f_cv_out, f_holdout, f_model_type) {
  
  # Prepare for plotting
  if (! is.data.frame(f_cv_out) ) {
    dfs <- list(f_training, f_holdout);
    lim_max <- 0.01 + max(
      f_training[, "observed"], f_training[, "predicted"], 
      f_holdout[, "observed"], f_holdout[, "predicted"],
      na.rm = TRUE)
    lim_min <- min(
      f_training[, "observed"], f_training[, "predicted"], 
      f_holdout[, "observed"], f_holdout[, "predicted"],
      na.rm = TRUE) - 0.01   

  }
  if (is.data.frame(f_cv_out) ) {
    dfs <- list(f_training, f_cv_out, f_holdout);
    lim_max <- 0.01 + max(
      f_training[, "observed"], f_training[, "predicted"], 
      f_cv_out[, "observed"], f_cv_out[, "predicted"],
      f_holdout[, "observed"], f_holdout[, "predicted"],
      na.rm = TRUE)
    lim_min <- min(
      f_training[, "observed"], f_training[, "predicted"], 
      f_cv_out[, "observed"], f_cv_out[, "predicted"],
      f_holdout[, "observed"], f_holdout[, "predicted"],
      na.rm = TRUE) - 0.01   
  }
  
  if (f_y_hat %in% c("samz", "samm")) {
    errors <- c(0.01, 0.02, 0.03) # note these will be absolute percentage errors
    thresholds <- c(0.02, 0.05)
    breaks_step <- 0.01
    error_labs <- paste("\u00B1", errors[1:3] * 100, "%", sep = "")
  }
  if (f_y_hat %in% c("gamz", "gamm")) {
    errors <- c(0.02, 0.05, 0.10) # note these will be absolute percentage errors
    thresholds <- c(0.15, 0.20)
    breaks_step <- 0.05
    error_labs <- paste("\u00B1", errors[1:3] * 100, "%", sep = "")
  }
  if (! f_y_hat %in% c("samz", "gamz", "samm", "gamm") ) {
    if (f_y_hat == "muac") {
      thresholds <- c(13.5, 13.0) # not very useful, in mm
      errors <- c(0.10, 0.25, 0.50)
      breaks_step <- round(abs(lim_max - lim_min) / 10, digits = 1)
    }  
    if (f_y_hat != "muac") {
      thresholds <- c(-0.5, -1.0) # not very useful, in z-scores
      errors <- c(0.05, 0.10, 0.20)
      breaks_step <- round(abs(lim_max - lim_min) / 10, digits = 1)
    }    
    error_labs <- paste("\u00B1", errors[1:3], sep = "")
  }
  
#  max_error <- max(errors)

  if (f_model_type == "glm") {colours <- c(palette_cb[7], palette_cb[4])}
  if (f_model_type == "rf") {colours <- c(palette_cb[6], palette_cb[8])}
  
  # Produce plots
  for (i in 1:length(dfs) ) {
    
    # prepare data
    df <- dfs[[i]]
    df <- na.omit(df)
    df[, "y_perfect"] <- df[, "observed"]
    df <- rbind(df, c(NA, NA, 0, NA, 0), c(NA, NA, lim_min, NA, lim_min) , c(NA, NA, lim_max, NA, lim_max))
    for (j in 1:length(errors[3:1]) ) {df[, paste("y_min", j, sep = "")] <- df[, "y_perfect"] - errors[j]}
    for (j in 1:length(errors[1:3]) ) {df[, paste("y_max", j, sep = "")] <- df[, "y_perfect"] + errors[j]}
    # text_angle <- 45 * abs(lim_max - lim_min) / (abs(lim_max - lim_min) + max(errors))
    text_angle <- 45

    # plot
    plot <- ggplot(df) +
      geom_point(aes(x = observed, y = predicted), size = 4, colour = colours[1], 
        fill = colours[1], alpha = 0.5) + 
      theme_bw() +
      geom_line (aes(y = y_perfect, x = observed), colour = colours[2], size = 1) +
      geom_ribbon(aes(ymin = y_min3, ymax = y_max3, x = observed), fill = colours[2], alpha = 0.1) +
      geom_ribbon(aes(ymin = y_min2, ymax = y_max2, x = observed), fill = colours[2], alpha = 0.2) +
      geom_ribbon(aes(ymin = y_min1, ymax = y_max1, x = observed), fill = colours[2], alpha = 0.3) +
      theme(axis.title = element_text(colour="grey20", size = 9)) +
      # annotate("text", x = lim_max - breaks_step, y = lim_max - breaks_step + errors[1] / 2, 
      annotate("text", x = lim_max - breaks_step - errors[1] / 2, y = lim_max - breaks_step, 
        label = error_labs[1], colour = colours[2], size = 3, angle = text_angle) +
      # annotate("text", x = lim_max - breaks_step, y = lim_max - breaks_step + (errors[2] + errors[1]) / 2, 
      annotate("text", x = lim_max - breaks_step - (errors[1] + errors[2]) / 2, y = lim_max - breaks_step, 
        label = error_labs[2], colour = colours[2], size = 3, angle = text_angle) +
      # annotate("text", x = lim_max - breaks_step, y = lim_max - breaks_step + (errors[3] + errors[2]) / 2,
      annotate("text", x = lim_max - breaks_step - (errors[3] + errors[2]) / 2, y = lim_max - breaks_step,
        label = error_labs[3], colour = colours[2], size = 3, angle = text_angle) +
      geom_hline(yintercept = thresholds, linetype = c("dotted", "dashed", "twodash")[1:length(thresholds)], 
        alpha = 0.5, size = 1, colour = colours[1] ) +
      geom_vline(xintercept = thresholds, linetype = c("dotted", "dashed", "twodash")[1:length(thresholds)], 
        alpha = 0.5, size = 1, colour = colours[1] )

    if (f_y_hat %in% c("samz", "gamz", "samm", "gamm") ) { 
      plot <- plot +
        scale_x_continuous("observed",
          breaks = seq(lim_min, lim_max, by = breaks_step),
          expand = c(0, 0), labels = scales::percent_format(accuracy = 1) ) +
        scale_y_continuous("predicted",
          # breaks = seq(lim_min, lim_max + max_error, by = breaks_step),
          breaks = seq(lim_min, lim_max, by = breaks_step),
          expand = c(0, 0), labels = scales::percent_format(accuracy = 1) ) +
        #coord_cartesian(ylim = c(lim_min, lim_max + max_error), xlim = c(lim_min, lim_max))
        coord_cartesian(ylim = c(lim_min, lim_max), xlim = c(lim_min, lim_max))
    }    
    
    if (! f_y_hat %in% c("samz", "gamz", "samm", "gamm") ) { 
      plot <- plot +
        scale_x_continuous("observed",
          breaks = seq(round(lim_min, 1), round(lim_max, 1), by = breaks_step),
          expand = c(0, 0) ) +
        scale_y_continuous("predicted",
          # breaks = seq(round(lim_min, 1), round(lim_max + max_error, 1), by = breaks_step),
          breaks = seq(round(lim_min, 1), round(lim_max, 1), by = breaks_step),
          expand = c(0, 0) ) +
        # coord_cartesian(ylim = c(lim_min, lim_max + max_error), xlim = c(lim_min, lim_max))         
        coord_cartesian(ylim = c(lim_min, lim_max), xlim = c(lim_min, lim_max))
    }
    
    print(plot)
    assign(paste("plot", i, sep = ""), plot)
    
  }

  # Produce combined plot
  if (exists("plot3")) {
    plot <- ggarrange(plot1, 
      plot2 + theme(axis.title.y = element_blank() ), 
      plot3 + theme(axis.title.y = element_blank() ), 
      ncol = 3, align = "v",
      labels = c("training data", "cross-validation", "holdout data"), vjust = 2.5, hjust = -1,
      font.label = list(size = 10, colour = "grey20", face = "plain")
      )
    print(plot)
    ggsave(paste(f_country, "_", f_y_hat, "_est_performance_", f_model_type, ".png", sep=""),
      height = 10, width = 30, units = "cm", dpi = "print")
  }

  if (! exists("plot3")) {
    plot <- ggarrange(plot1, 
      plot2 + theme(axis.title.y = element_blank() ), 
      ncol = 2, align = "v",
      labels = c("training data", "holdout data"), vjust = 2.5, hjust = -1,
      font.label = list(size = 10, colour = "grey20", face = "plain")
      )
    print(plot)
    ggsave(paste(f_country, "_", f_y_hat, "_est_performance_", f_model_type, ".png", sep=""),
      height = 10, width = 20, units = "cm", dpi = "print")
  }
  
}


#.........................................................................................
### Function to predict on new data, given a robust variance-covariance matrix (for fixed-effects models only)
  # (based on https://stackoverflow.com/questions/3790116/using-clustered-covariance-matrix-in-predict-lm?noredirect=1&lq=1 )
#.........................................................................................

f_predict <- function(f_fit, f_vcov_cl, f_newdata, f_se_fit) {
    
  # if new data are missing, revert to predicting on training dataset
  if (missing(f_newdata)) { f_newdata <- f_fit$model }
    
  # identify terms of the model and remove the response term
  fit_terms <- delete.response( terms(f_fit) )
    
  # construct model matrix
  m_mat<- model.matrix(fit_terms, model.frame(fit_terms, f_newdata, na.action = "na.pass"))
    
  # access model coefficients
  m_coef <- f_fit$coef
    
  # generate predictions on the desired scale, as well as standard errors for the predictions
  if (family(fit_glm)[2] == "log") {fit <- as.vector(m_mat %*% f_fit$coef)}
  if (family(fit_glm)[2] == "linear") {fit <- inv.logit(as.vector(m_mat %*% f_fit$coef))}
  se_fit <- sqrt(diag(m_mat %*% f_vcov_cl %*% t(m_mat)))
    
  # return predictions or standard errors, as desired      
  if (f_se_fit == TRUE) {return(se_fit)} else {return(fit)}
}



#.........................................................................................
### Function to prepare data for random forest (only complete cases or imputation, as desired)
#.........................................................................................

f_rf_prep <- function(f_y_hat, f_vars, f_data, f_impute, f_part_unit) {

  # Select complete cases or perform imputation 
  if (f_impute == FALSE) {
    # select only complete cases
    out <- f_data[complete.cases(f_data[, c(f_part_unit, "tm", f_y_hat, f_vars)]), 
      c(f_part_unit, "tm", f_y_hat, f_vars)]
  }
  
  if (f_impute == TRUE) {
    # select cases for which the outcome is complete
    x1 <- f_data[complete.cases(f_data[, f_y_hat]), c(f_y_hat, f_vars)]
    
    # perform imputation by proximity method for missing independent variable values
    out <- rfImpute(as.formula(paste(f_y_hat, "~", ".", sep = "")), x1, iter = 50)
    out <- cbind(f_data[, c(f_part_unit, "tm")], out)
  }

  # Output prepared dataframe            
  return(out)
}



#.........................................................................................
### Function to estimate robust standard errors for model coefficients, given a robust variance-covariance matrix
#.........................................................................................
  
f_rob <- function(f_fit, f_vcov_cl, f_scale) {
  
  # Calculate robust standard errors
  std.err <- sqrt(diag(f_vcov_cl))
    
  # Output on the link scale
  if (f_scale == "link") {
    r_est <- data.frame("Relative risk" = round(coef(f_fit), 3),
      "95%CI - lower" = round(coef(f_fit) - 1.96 * std.err, 3),
      "95%CI - upper" = round(coef(f_fit) + 1.96 * std.err, 3),
      "Pr(>|z|)" = round(2 * pnorm(abs(coef(f_fit)/std.err), lower.tail = FALSE), 3)
    )
  }
  
  # Output on the response scale
  if (f_scale == "response" & as.character(family(f_fit))[1] != "gaussian") {
    r_est <- data.frame("Relative risk" = round(exp(coef(f_fit)), 3),
      "95%CI - lower" = round(exp(coef(f_fit) - 1.96 * std.err), 3),
      "95%CI - upper" = round(exp(coef(f_fit) + 1.96 * std.err), 3),
      "Pr(>|z|)" = round(2 * pnorm(abs(coef(f_fit)/std.err), lower.tail = FALSE), 3)
    )
  }

  if (f_scale == "response" & as.character(family(f_fit))[1] == "gaussian") {
    r_est <- data.frame("Relative risk" = round(coef(f_fit), 3),
      "95%CI - lower" = round(coef(f_fit) - 1.96 * std.err, 3),
      "95%CI - upper" = round(coef(f_fit) + 1.96 * std.err, 3),
      "Pr(>|z|)" = round(2 * pnorm(abs(coef(f_fit)/std.err), lower.tail = FALSE), 3)
    )
  }
  
    
  # Return output
  r_est <- cbind(names(coef(f_fit)), r_est)
  colnames(r_est) <- c("coefficient", "relative risk", "95%CI - lower", "95%CI - upper", "p-value")
  rownames(r_est) <- c()
  return(r_est)
}


#.........................................................................................
### Function to execute sections of code 
  # (based on https://stackoverflow.com/questions/26245554/execute-a-set-of-lines-from-another-r-file )
#.........................................................................................

f_source_part <- function(f_script, f_start_tag, f_end_tag) {

  # Identify lines with start and end tags
  st <- grep(f_start_tag, f_script)
  en <- grep(f_end_tag, f_script)
  
  # Set up a connection
  tc <- textConnection(f_script[(st + 1):(en - 1)])
  
  # Run the script
  source(tc)
  
  # Close the connection
  close(tc)
}



#.........................................................................................
### Function to fit any model formula (with/out random effects), generate model fit statistics and plot CV 
#.........................................................................................

f_val <- function(f_y_hat, f_preds, f_interaction_terms, f_data, f_data_holdout, f_re, f_family, f_cv, 
  f_part_unit, f_k_folds, f_show_output) {

  # Identify variables needed within the dataset
   # variables that are part of interaction terms
    int_terms <- c()
    if (length(f_interaction_terms) > 0) {
      int_terms <- unlist(strsplit(f_interaction_terms, ":"))
    }
  
  # Select non-missing data
  f_data <- f_data[complete.cases(f_data[, c(f_y_hat, f_preds, int_terms)]), ]
  if (is.data.frame(f_data_holdout) ) {f_data_holdout <- 
    f_data_holdout[complete.cases(f_data_holdout[, c(f_y_hat, f_preds, int_terms)]), ] }

  # Prepare output
  out <- rep(NA, times = 12)
  
  # Fixed effects only model (always)
    # write the model formula
   form <- as.formula( paste(f_y_hat, "~", paste(c(f_preds, f_interaction_terms), collapse = "+"), sep = "")  )
      
    # fit GLM
    fit <- glm(form, data = f_data, family = f_family, weights = wt )
      # print model summary
      if (f_show_output == TRUE) {
        print(tidy(fit, exponentiate = ifelse(f_family == "gaussian", FALSE, TRUE) ));
        print(glance(fit))
      }
      
      # record AIC
      out[1] <- fit$aic
      
      # generate elements to calculate predictive scores
      lambdas <- predict(fit, type = "response")
      ys <- f_data[, f_y_hat]

      # aggregate data and predictions
      x1 <- as.data.frame(cbind(f_data[, c(f_part_unit, "tm")], ys, lambdas))
      colnames(x1) <- c(f_part_unit, "tm", "ys", "lambdas")
      x1 <- aggregate(x1[, c("ys", "lambdas")], by = x1[, c(f_part_unit, "tm")], FUN = mean, na.rm = TRUE)
      colnames(x1) <- c(f_part_unit, "tm", "ys", "lambdas")

      # calculate mean square error (MSE)
      out[2] <- mean((x1$lambdas - x1$ys )^2, na.rm = TRUE)
      
      # calculate model Chi-square or F-test p-value
      if (f_family == "gaussian") {out[3] <- anova(fit, test = "F")[2, "Pr(>F)"] }
      if (! f_family == "gaussian") {out[3] <- anova(fit, test = "Chisq")[2, "Pr(>Chi)"] }
      
    # if desired, do cross-validation and extract MSE, plus proportion of folds that fitted OK
    if (f_cv == TRUE) { out[c(4, 5)] <- f_cv(fit, f_data, f_part_unit, f_k_folds, TRUE) }

    # if desired, predict on holdout data and compute MSE
    if (is.data.frame(f_data_holdout) ) {
      
      # generate elements to calculate predictive scores
      lambdas <- predict(fit, newdata = f_data_holdout, type = "response")
      ys <- f_data_holdout[, f_y_hat]

      # aggregate data and predictions
      x1 <- as.data.frame(cbind(f_data_holdout[, c(f_part_unit, "tm")], ys, lambdas))
      colnames(x1) <- c(f_part_unit, "tm", "ys", "lambdas")
      x1 <- aggregate(x1[, c("ys", "lambdas")], by = x1[, c(f_part_unit, "tm")], FUN = mean, na.rm = TRUE)
      colnames(x1) <- c(f_part_unit, "tm", "ys", "lambdas")

      # calculate mean square error
      out[6] <- mean((x1$lambdas - x1$ys )^2, na.rm = TRUE)
    }        
      
  # Random effects model, if desired: fitted without weights, as the latter cause singular fits / fitting problems
  if (f_re == TRUE) {
    # write the model formula
      form <- as.formula( paste(f_y_hat, "~", paste(c(f_preds, f_interaction_terms), collapse = "+"), 
        "+ (1|", f_part_unit, ")", sep = "")  )
    
    # fit GLM with random effect
    # fit <- glmer(form, data = f_data, family = gsub("quasi", "", f_family), weights = wt )
    fit <- glmer(form, data = f_data, family = gsub("quasi", "", f_family) )
      
      # print model summary
      if (f_show_output == TRUE) {
        print(summary(fit))
        print(broom.mixed::tidy(fit, exponentiate = ifelse(f_family == "gaussian", FALSE, TRUE) ))
        print(broom.mixed::glance(fit))
      }
    
      # record AIC
      out[7] <- AIC(fit)
      
      # generate elements to calculate predictive scores        
      lambdas <- predict(fit, type = "response")
      ys <- f_data[, f_y_hat]
      
      # aggregate data and predictions
      x1 <- as.data.frame(cbind(f_data[, c(f_part_unit, "tm")], ys, lambdas))
      colnames(x1) <- c(f_part_unit, "tm", "ys", "lambdas")
      x1 <- aggregate(x1[, c("ys", "lambdas")], by = x1[, c(f_part_unit, "tm")], FUN = mean, na.rm = TRUE)
      colnames(x1) <- c(f_part_unit, "tm", "ys", "lambdas")

      # calculate mean square error
      out[8] <- mean((x1$lambdas - x1$ys )^2, na.rm = TRUE)
      
      # calculate model p-value / not straightforward for mixed models
      out[9] <- NA
    
    # if desired, do cross-validation and extract MSE, plus proportion of folds that fitted OK
    if (f_cv == TRUE) { out[c(10,11)] <- f_cv(fit, f_data, f_part_unit, f_k_folds, TRUE) }
      
    # if desired, predict on holdout data and compute MSE 
    if (is.data.frame(f_data_holdout) ) {
      
      # generate elements to calculate predictive scores
      lambdas <- predict(fit, newdata = f_data_holdout, type = "response", allow.new.levels = TRUE)
      ys <- f_data_holdout[, f_y_hat]

      # aggregate data and predictions
      x1 <- as.data.frame(cbind(f_data_holdout[, c(f_part_unit, "tm")], ys, lambdas))
      colnames(x1) <- c(f_part_unit, "tm", "ys", "lambdas")
      x1 <- aggregate(x1[, c("ys", "lambdas")], by = x1[, c(f_part_unit, "tm")], FUN = mean, na.rm = TRUE)
      colnames(x1) <- c(f_part_unit, "tm", "ys", "lambdas")

      # calculate mean square error
      out[12] <- mean((x1$lambdas - x1$ys )^2, na.rm = TRUE)
    }    
      
  }
      
  # Return results
  names(out) <- c("aic_glm", "mse_glm", "p_value_glm", "mse_cv_glm", "prop_folds_ok_glm", 
    "mse_holdout_glm", "aic_glmm", "mse_glmm", "p_value_glmm", "mse_cv_glmm", "prop_folds_ok_glmm",
    "mse_holdout_glmm")
  return(out)
}    
  
    
  
#.........................................................................................
### ENDS
#.........................................................................................
