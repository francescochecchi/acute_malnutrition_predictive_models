#..........................................................................................
### +++++++++++++++ SMALL-AREA PREDICTION OF ACUTE MALNUTRITION BURDEN ++++++++++++++++ ###
#..........................................................................................

#..........................................................................................
## --------------- R CODE TO RE-ANALYSE AND PREPARE PREDICTOR VARIABLES ---------------- ##
#..........................................................................................

                                          # Written by Francesco Checchi, LSHTM (Nov 2020)
                                          # francesco.checchi@lshtm.ac.uk 


#.........................................................................................                            
### Reading files from previous scripts
#.........................................................................................

  #...................................   
  ## Child observations
  ch_obs <- read.csv(paste(country, "_ch_obs.csv", sep = ""), sep = ",")  


#.........................................................................................                            
### Merging population onto time series
#.........................................................................................
    
  #...................................   
  ## Merge population data with time series
  
    # identify relevant columns in population dataframe
    x1 <- colnames(pop)[colnames(pop) %in% c("stratum", "m", "y", "pop_average", "pop_average_u5")]
    
    # merge
    ts <- merge(ts, pop[, x1], by = c("stratum", "y", "m"), all.x = TRUE, sort = TRUE)

    
#.........................................................................................                            
### Merging the other predictors onto the time series
#.........................................................................................
    
  #...................................   
  ## Aggregate any predictors by variables of interest prior to merging
    # (so as to be left with unique values for each geo-time unit)
    
  for (i in 1:nrow(predictors) )  {
    # If the predictor needs to be aggregated...
    if (is.na( predictors[i, "aggregate_operation"]) == FALSE ) {
    
      # get predictor dataset that is being reshaped
      x1 <- get(predictors[i, "worksheet"] )
  
      # variables to aggregate over
      aggregate_over <- c(predictors[i, "aggregate_over"] )
        # exclude variable(s) to aggregate over and missing observations
        x1 <- x1[, ! colnames(x1) %in% aggregate_over]
        x1 <- na.omit(x1)
        
      # operation to perform (e.g. sum, mean, etc.)
      op <- as.name(predictors[i, "aggregate_operation"] )
      
      # time_unit of predictor
      timeunit <- predictors[i, "time_unit"]
      
      # geo_unit of predictor
      geounit <- predictors[i, "geo_unit"]
        # work out generic geo_unit of predictor if needed
        if (geounit %in% c(admin2_name, admin1_name) ) {
          x3 <- ifelse(geounit == admin2_name, "stratum", "admin1")
          geounit <- x3
        }   
      
      # work out which variables the dataset should be aggregated by, based on geo_unit and time_unit
        # if geo_unit = admin1...
        if (geounit == "admin1") {
          if (timeunit == "m") {by_vars <- list(admin1 = x1$admin1, y = x1$y, m = x1$m) }
          if (timeunit == "y") {by_vars <- list(admin1 = x1$admin1, y = x1$y) }
          if (timeunit == "static") {by_vars <- list(admin1 = x1$admin1) }
          }
        # if geo_unit = stratum and admin1 is not part of the predictor dataset...
        if (geounit == "stratum" & ! "admin1" %in% colnames(x1) ) {
          if (timeunit == "m") {by_vars <- list(stratum = x1$stratum, y = x1$y, m = x1$m) }
          if (timeunit == "y") {by_vars <- list(stratum = x1$stratum, y = x1$y) }
          if (timeunit == "static") {by_vars <- list(stratum = x1$stratum) }
          }
      # if geo_unit = stratum and admin1 is also part of the predictor dataset...
      if (geounit == "stratum" & "admin1" %in% colnames(x1) ) {
        if (timeunit =="m") {by_vars <- list(admin1 = x1$admin1, stratum = x1$stratum, y = x1$y, m = x1$m) }
        if (timeunit =="y") {by_vars <- list(admin1 = x1$admin1, stratum = x1$stratum, y = x1$y) }
        if (timeunit =="static") {by_vars <- list(admin1 = x1$admin1, stratum = x1$stratum) }
      }      
        # if there is any other variable to aggregate by...
        if (is.na( predictors[i, "aggregate_by"]) == FALSE ) { 
          by_vars <- c(by_vars, as.list(x1[, c(unlist(predictors[i, "aggregate_by"]))] ) ) 
          }
      
      # aggregate:
      x1 <- aggregate(x1[ , ! colnames(x1) %in% names(by_vars)], by = by_vars, FUN = op )
      
      # replace the predictor dataset with the aggregated one
      assign(predictors[i, "worksheet"], x1 )
    
    }
  }


  #...................................       
  ## Reshape any predictors prior to merging
  
    for (i in 1:nrow(predictors) )  {
    
      # If the predictor needs to be aggregated...
      if (is.na( predictors[i, "reshape_needed"]) == FALSE ) {
        
        # get predictor dataset that is being reshaped
        x1 <- get(predictors[i, "worksheet"] )
        x2 <- predictors[i, "worksheet"]
        
        # set unique time unit to reshape on
          
          # what is the time_unit of the predictor?
          timeunit <- predictors[predictors$worksheet == x2 , "time_unit"]
          
          # if time_unit is year, use that; if it is month, need to merge with tm to have unique time unit
          if (timeunit == "y") { time_reshape <- "y" }
          if (timeunit == "m") { x1 <- merge(x1, t_units, by = c("y", "m"), y.all = FALSE, sort = TRUE) ;
                                 time_reshape <- "tm" ;
                                 x1 <- x1[ , ! colnames(x1) %in% c("y", "m")] 
                               }
          
      # Set unique geo_unit to reshape on
        
        # what is the geo_unit of the predictor?
        geo_reshape <- predictors[predictors$worksheet == x2 , "geo_unit"]
        
        # work out generic geo_unit of predictor if needed
        if (geo_reshape %in% c(admin2_name, admin1_name) ) {
          x3 <- ifelse(geo_reshape == admin2_name, "stratum", "admin1")
          geo_reshape <- x3
        }
        
      # reshape wide:
      x1 <- reshape(x1, idvar = c(geo_reshape, time_reshape), 
        timevar = paste(subset(predictors, worksheet == x2 )[, "reshape_by"]), direction = "wide", sep = "_")
      
      # merge time units back in if the time unit is months (as both m and y would have been replaced by tm only)
      if (time_reshape == "tm") { 
        x1 <- merge(x1, t_units , by="tm", all.y=FALSE, sort=TRUE) ;
        x1 <- subset(x1, select = -tm) # otherwise tm will be merged twice later on 
      }
      
      # replace the predictor dataset with the reshaped one
      assign(predictors[i, "worksheet"], x1 )
        
    }
  }
 
       
  #...................................    
  ## Merge predictors onto the time series 
    
    # Create a vector in which to store names of predictor datasets that cannot be merged on first pass
    later <- vector()
    
    # First pass at merging:
    for (i in 1:nrow(predictors) )  {
      # get predictor dataset that is being merged
      x1 <- get(predictors[i, "worksheet"] )
      
      # work out which time units (y, m) and geo units (admin1, stratum, other) are in the dataset
      by_vars <- colnames(x1)[colnames(x1) %in% c("y", "m", "admin1", "stratum", predictors[i, "geo_unit"])]

      # if the predictor's time and geo_units are all in the master dataset (x), then go ahead with the merge...
      if( all(by_vars %in% colnames(ts)) ) {
        ts <- merge(ts, x1, by = by_vars, all.x=TRUE, sort=TRUE) 
        print(paste("after merging predictor ", paste(predictors[i, "worksheet"]), 
          ", the number of rows is ", paste(nrow(ts)), sep="") )
                
      }
      
      # if the predictor's geo_unit is not (yet) in the master dataset, then store predictor name for second pass...
        # by which time the geo_unit should be in the master dataset
      if( ! all(by_vars %in% colnames(ts)) ) {
        later <- c(later, paste(predictors[i, "worksheet"]))
      }
      
    # Close first merging pass loop      
    }    

          
    # Second pass at merging with any predictors that couldn't yet be merged:
      # (for each predictor that couldn't be merged the first time...)
    if (length(later) > 0 ) {
    
      for (i in 1:length(later) )  {
        # get predictor dataset that is being merged
        x1 <- get(later[i] )
        x2 <- which(predictors[, "worksheet"] == later[i])
        
        # work out which time units (y, m) and geo units (admin1, stratum, other) are in the dataset
        by_vars <- colnames(x1)[colnames(x1) %in% c("y", "m", "admin1", "stratum", predictors[x2, "geo_unit"])]
        
        # if the predictor's time and geo_units are all in the master dataset (x), then go ahead with the merge...
        if( all(by_vars %in% colnames(ts)) ) {
          ts <- merge(ts, x1, by = by_vars , all.x=TRUE, sort=TRUE) 
          print(paste("after merging predictor ", predictors[i, "worksheet"], 
            ", the number of rows is ", nrow(ts), sep="") )
          
          }
        
        else print( paste("cannot merge predictor ", x2, " even on second pass", sep="") )
        
      # Close second merging pass loops
      }   
    }
    

#.........................................................................................                            
### Preparing predictor data for analysis, depending on their specific features / requirements
#.........................................................................................
    
  #...................................
  ## Standardise livelihood nomenclature, if there is at least one livelihood variable
    
    # Find the index/indices of the livelihood variable(s), if there are any (if not, NA)
    indices <- grep(paste(livelihood_substrings$livelihood, collapse="|"), colnames(ts))
    
    # If any livelihood variable is included...
    if (length(indices) > 0) {
      
      # for each variable...
      for (i in indices) {
        
        # if the variable contains digits, then it should not be changed (likely to be livelihood ID)
        if (length(grep(paste(1:10, collapse="|"), ts[, i] ))  > 0 ) {next}
        
        # otherwise, continue...
        ts[, i] <- sapply(ts[, i], f_liv, livelihood_substrings)
        
      }
    }
    
    
  #...................................
  ## Set missing values to zero/missing/informative missing, for predictors where this is the assumed interpretation of missingness
    # but only within data availability periods
  
  for (i in 1:nrow(predictors) )  {
    # If the predictor needs missing values set to zero or 'informative missing'...
    if (predictors[i, "missing_meaning"] %in% c("zero", "informative_missing") ) {
  
      # which predictor dataset is being worked on?
      x1 <- predictors[i, "worksheet"]   
    
      # identify columns belonging to the predictor
      x_vars <- subset(dictionary, worksheet == x1 & used_in_analysis == "Y" )$variable
        # exclude geo_unit and time_unit columns
        x_vars <- x_vars[! x_vars %in% c(admin1_name, admin2_name, "y", "m")]
        # also exclude columns that may have already been aggregated out or otherwise removed
        x_vars <- x_vars[x_vars %in% colnames(ts)]
        
      # identify start and end points for the predictor's data availability (if no dates are given, assume start/end of analysis period)
      m_start <- ifelse(is.na(predictors[i, "date_start"]) == TRUE, m_analysis_start, month(predictors[i, "date_start"]) )
      m_end <- ifelse(is.na(predictors[i, "date_end"]) == TRUE, m_analysis_end, month(predictors[i, "date_end"]) )
      y_start <- ifelse(is.na(predictors[i, "date_start"]) == TRUE, y_analysis_start, year(predictors[i, "date_start"]) )
      y_end <- ifelse(is.na(predictors[i, "date_end"]) == TRUE, y_analysis_end, year(predictors[i, "date_end"]) )
        
        # convert to corresponding tm values 
        x2 <- as.Date.character(paste(y_start, m_start, "01", sep="-")) - 
          as.Date.character(paste(t_units[t_units$tm == min(t_units$tm), "y"], t_units[t_units$tm == min(t_units$tm), "m"], "01", sep="-"))
        if (x2 <= 0) {tm_start <- min(t_units$tm)}
        if (x2 > 0) {tm_start <- t_units[t_units$y == y_start & t_units$m == m_start, "tm"] }
        x2 <- as.Date.character(paste(y_end, m_end, "01", sep="-")) - 
          as.Date.character(paste(t_units[t_units$tm == max(t_units$tm), "y"], t_units[t_units$tm == max(t_units$tm), "m"], "01", sep="-"))
        if (x2 <= 0) {tm_end <- t_units[t_units$y == y_end & t_units$m == m_end, "tm"] }
        if (x2 > 0) {tm_end <- max(t_units$tm) }

      # set NA values to desired value within period of data availability; leave NA values as NA outside this period
        # if want 0's...
        if (predictors[i, "missing_meaning"] == "zero") {
          ts[ts$tm >= tm_start & ts$tm <= tm_end, x_vars] <- na.replace(ts[ts$tm >= tm_start & ts$tm <= tm_end, x_vars], 0)
        }
        
        # if want informative missing...
        if (predictors[i, "missing_meaning"] == "informative_missing") {
          ts[ts$tm >= tm_start & ts$tm <= tm_end, x_vars] <- na.replace(ts[ts$tm >= tm_start & ts$tm <= tm_end, x_vars], "missing/none")
        }
        
    } 
  }

   
  #...................................
  ## Visualise and describe predictor completeness; remove predictors or years that are below completeness cutoff
    
    # For each variable of each predictor, set value = 0 if NA and 1 otherwise
      
      # identify columns that this operation should not apply to
      cols_include <- colnames(ts)[! colnames(ts) %in% c("admin1", "stratum", "y", "m", "tm", colnames(strata), "pop_average", "pop_average_u5", predictors$geo_unit)]
      cols_exclude <- colnames(ts)[! colnames(ts) %in% cols_include]
      
      # generate completeness database
      complete <- apply(ts[, cols_include], c(1,2), function(x) { if (is.na(x) ) return(0) else return(1) })
      complete <- cbind(ts[, cols_exclude], complete)
    
    # Restrict analysis of completeness to the actual period of analysis + maximum lag
      
      # remove observations for the first few months
      complete <- subset(complete, tm > (burn_in_period * 12 - max(predictors$lags) ) )
    
    # Compute completeness by stratum
    complete_by_stratum <- aggregate(complete[, cols_include], by=list(admin1 <- complete$admin1, stratum <- complete$stratum), FUN=mean)
    colnames(complete_by_stratum)[1:2] <- c("admin1", "stratum")

      # plot completeness by stratum
      x1 <- melt(complete_by_stratum, id.vars = c("admin1", "stratum"), measure.vars = cols_include )

      plot <- ggplot(x1, aes(x = stratum, y = variable) )
      plot <- plot + geom_tile(aes(fill=value), colour = "grey80", show.legend = TRUE) + 
        scale_x_discrete(paste(admin2_name, ", ", admin1_name, sep=""), expand=c(0,0)) + 
        scale_fill_gradientn(colours = c("grey90", "red", "green"), 
                             values = c(0, 0.0000001, 1 )) +
        facet_grid(~admin1, space="free", scales="free", switch="x") +
        theme_bw() +
        labs(fill = "completeness") +
        theme(axis.text.x = element_text(angle = 90, vjust=0.2)) +
        theme(axis.title = element_text(colour="grey20")) +
        theme(strip.placement = "outside",
              strip.background = element_rect(fill=NA, colour="grey50"),
              panel.spacing=unit(0,"cm"), strip.text.x = element_text(angle = 90, colour="grey20"))
      
      plot
      ggsave(paste(country, "_prd_completeness_stratum.png", sep=""), height = 30, width = 45, units = "cm", dpi = "print")
       
      # apply tolerance cut-offs for missingness
      complete_by_stratum_pass <- apply(complete_by_stratum[, cols_include], c(1,2),
        function(f_x, f_cutoff) { if (f_x >= f_cutoff) return(1) else return(0) } , tol_completeness_geo )

      # identify which variables pass the cut-off
      complete_by_stratum_pass <- colMeans(complete_by_stratum_pass)
      complete_by_stratum_pass <- complete_by_stratum_pass[which(complete_by_stratum_pass >= tol_completeness_geo)]

    # Compute completeness by month
    complete_by_month <- aggregate(complete[, cols_include], by=list(tm <- complete$tm), FUN=mean)
    colnames(complete_by_month)[1] <- "tm"
        
      # plot completeness by month
      x1 <- melt(complete_by_month, id.vars = "tm", measure.vars = cols_include )
      x1 <- merge(x1, t_units, by = "tm", x.all = TRUE, sort = TRUE)
      x1[, "date"] <- lubridate::ymd(paste(x1[, "y"], x1[, "m"], 15, sep = "-"))

      plot <- ggplot(x1, aes(x = date, y = variable) )
      plot <- plot + geom_tile(aes(fill=value), colour = "grey80", show.legend = TRUE) + 
        scale_x_date("year" ) +
        scale_fill_gradientn(colours = c("grey90", "red", "green"), values = c(0, 0.0000001, 1 )) +
        facet_grid(~y, space="free", scales="free", switch="x") +
        theme_bw() +
        theme(axis.title = element_text(colour="grey20")) +
        labs(fill = "completeness") +
        theme(strip.placement = "outside",
          strip.background = element_rect(fill=NA, colour="grey50"),
          panel.spacing=unit(0,"cm"), strip.text.y = element_text(angle = 0), strip.text.x = element_text(vjust=0)) +
        theme(axis.text.x = element_blank(), axis.ticks.x = element_blank() )
      
      plot      
      ggsave(paste(country, "_prd_completeness_time.png", sep=""), height = 30, width = 45, units = "cm", dpi = "print")

      # apply tolerance cut-offs for missingness

        # re-compute completeness by month with further restriction to actual period of estimation
        complete_by_month <- aggregate(complete[, cols_include], by=list(tm <- complete$tm), FUN=mean)
        colnames(complete_by_month)[1] <- "tm"

        # apply function
        complete_by_month_pass <- apply(complete_by_month[, cols_include], c(1,2),
          function(f_x, f_cutoff) { if (f_x >= f_cutoff) return(1) else return(0) } , tol_completeness_tm )

        # identify which variables pass the cut-off (at least x1 complete for at least complete months)
        complete_by_month_pass <- colMeans(complete_by_month_pass)
        complete_by_month_pass <- complete_by_month_pass[which(complete_by_month_pass >= tol_completeness_tm)]

    # Eliminate variables that don't meet completeness cutoffs for either the stratum or the time dimension
    cols_keep <- intersect(names(complete_by_month_pass), names(complete_by_stratum_pass))
    ts <- ts[, c(cols_exclude, cols_keep)]


  #...................................  
  ## Perform any manual imputations as per manual_imputations table

  # If there are any manual imputations to be done...  
  if (nrow(manual_imputations) > 0) {
    
    # for each manual imputation to do...
    for (i in 1:nrow(manual_imputations)) {
    
      # identify predictor...
      pred_imp <- manual_imputations[i, "worksheet"]
        # if the predictor has already been eliminated from analysis, skip to next loop
        if (! pred_imp %in% c(predictors$worksheet)) next; 
  
      #...and corresponding variables:
      vars_imp <- c(subset(dictionary, worksheet == paste(pred_imp) & used_in_analysis == "Y")[, "variable"] )
        # but exclude any geo_unit or time_unit variables:
        vars_imp <- unlist(vars_imp)
        vars_imp <- vars_imp[! vars_imp %in% c("admin1", "stratum", "y", "m", "tm", admin1_name, admin2_name, 
          colnames(strata), colnames(pop), predictors$geo_unit)]
        
        # also exclude variables that do not pass the completeness tests
        vars_imp <- vars_imp[vars_imp %in% cols_keep]
        
        # if all of the variables have already been eliminated from analysis, skip to next loop
        if (length(vars_imp) == 0) next;    
      
      # identify strata to be attributed an imputed value:
        # if the admin1 column is filled in, take all the strata in those admin1...
        if (is.na(manual_imputations[i, "admin1"]) == FALSE) {
          admin1_imp <- unlist(strsplit(paste(manual_imputations[i, "admin1"]), ", "))
          admin1_all <- unique(strata[, c("admin1", "stratum")])
          strata_imp <- admin1_all[admin1_all[, "admin1"] %in% admin1_imp, ][, "stratum"]
          }
        
        #if the stratum column is filled in, only consider those strata after all...
        if (is.na(manual_imputations[i, "stratum"]) == FALSE) {
          strata_imp <- unlist(strsplit(paste(manual_imputations[i, "stratum"]), ", ")) 
          }
          
      # identify start and end points for the imputation (if no dates are given, assume start/end of analysis period)
      m_start <- ifelse(is.na(manual_imputations[i, "date_start"]) == TRUE, m_analysis_start, month(manual_imputations[i, "date_start"]) )
      m_end <- ifelse(is.na(manual_imputations[i, "date_end"]) == TRUE, m_analysis_end, month(manual_imputations[i, "date_end"]) )
      y_start <- ifelse(is.na(manual_imputations[i, "date_start"]) == TRUE, y_analysis_start, year(manual_imputations[i, "date_start"]) )
      y_end <- ifelse(is.na(manual_imputations[i, "date_end"]) == TRUE, y_analysis_end, year(manual_imputations[i, "date_end"]) )
      
        # convert to corresponding tm values 
        x2 <- as.Date.character(paste(y_start, m_start, "01", sep="-")) - 
          as.Date.character(paste(t_units[t_units$tm == min(t_units$tm), "y"], t_units[t_units$tm == min(t_units$tm), "m"], "01", sep="-"))
        if (x2 <= 0) {tm_start <- min(t_units$tm)}
        if (x2 > 0) {tm_start <- t_units[t_units$y == y_start & t_units$m == m_start, "tm"] }
        x2 <- as.Date.character(paste(y_end, m_end, "01", sep="-")) - 
          as.Date.character(paste(t_units[t_units$tm == max(t_units$tm), "y"], t_units[t_units$tm == max(t_units$tm), "m"], "01", sep="-"))
        if (x2 <= 0) {tm_end <- t_units[t_units$y == y_end & t_units$m == m_end, "tm"] }
        if (x2 > 0) {tm_end <- max(t_units$tm) }
      
      # work out length of imputation period, constrained by limits of total data availability/analysis period
      period_all <- unique(ts[,"tm"])
      period_imp <- ( min(tm_end, max(period_all)) - max(tm_start, min(period_all)) + 1)
      
      # if the actual imputation value(s) are specified...
      if (is.na(manual_imputations[i, "imputation_value"]) == FALSE) { 
        
        # attribute imputed values
        ts <- ts[order(ts[, "admin1"], ts[, "stratum"], ts[, "tm"]), ]
        ts1 <- subset(ts[ts[, "stratum"] %in% strata_imp, ], tm >= tm_start & tm <= tm_end)
        ts2 <- subset(ts[ts[, "stratum"] %in% strata_imp, ], tm < tm_start | tm > tm_end)
        ts3 <- ts[! ts[, "stratum"] %in% strata_imp, ]
        ts1[, c(vars_imp)] <- manual_imputations[i, "imputation_value"]
        ts <- rbind(ts1, ts2, ts3)
      }
  
      # if the imputation value(s) are to be weighted means of values from other geo_units...
      if (is.na(manual_imputations[i, "ref_geo_unit1"]) == FALSE) {
        # identify level of geo_units to take imputation values from
        level_ref <- ifelse(manual_imputations[i, "ref_geo_unit1"] %in% c(unique(strata$stratum)) , "stratum", "admin1")
        
        # identify geo_units in different reference categories, if the reference geo_unit is
          # NOT specified as "all others" (i.e. any not yet mentioned)
        if (paste(manual_imputations[i, "ref_geo_unit1"]) != "all others") { 
           refunits1 <- unlist(strsplit(paste(manual_imputations[i, "ref_geo_unit1"]), ", ")) }
        if (paste(manual_imputations[i, "ref_geo_unit2"]) != "all others") { 
          refunits2 <- unlist(strsplit(paste(manual_imputations[i, "ref_geo_unit2"]), ", ")) }
        if (paste(manual_imputations[i, "ref_geo_unit3"]) != "all others") { 
          refunits3 <- unlist(strsplit(paste(manual_imputations[i, "ref_geo_unit3"]), ", ")) }
        
        # resolve instances in which the reference geo_unit is specified as "all others" (i.e. any not yet mentioned)
        if (paste(manual_imputations[i, "ref_geo_unit1"]) == "all others") { 
          all <- unique(strata[, level_ref])
          refunits1 <- all[! all %in% manual_imputations[i, level_ref]]
          refunits1 <- as.character(refunits1)
          }
        if (paste(manual_imputations[i, "ref_geo_unit2"]) == "all others") { 
          all <- unique(strata[, level_ref])
          refunits2 <- all[! all %in% c(manual_imputations[i, level_ref], refunits1)]
          refunits2 <- as.character(refunits2)
          }
        if (paste(manual_imputations[i, "ref_geo_unit3"]) == "all others") { 
          all <- unique(strata[, level_ref])
          refunits3 <- all[! all %in% c(manual_imputations[i, level_ref], refunits1, refunits2)]
          refunits3 <- as.character(refunits3)
          }
        
        # for each variable to be imputed...
        for (j in vars_imp) {
          # extract imputation values and arrange them as columns of strata, with n rows = imputation period
          x1 <- ts[ts[, paste(level_ref)] %in% refunits1, ]
            x1 <- x1[order(x1[, "admin1"], x1[, "stratum"], x1[, "tm"]), ]
            x1 <- matrix(subset(x1, tm >= tm_start & tm <= tm_end)[, j], nrow = period_imp )
          x2 <- ts[ts[, paste(level_ref)] %in% refunits2, ]
            x2 <- x2[order(x2[, "admin1"], x2[, "stratum"], x2[, "tm"]), ]
            x2 <- matrix(subset(x2, tm >= tm_start & tm <= tm_end)[, j], nrow = period_imp )
          x3 <- ts[ts[, paste(level_ref)] %in% refunits3, ]
            x3 <- x3[order(x3[, "admin1"], x3[, "stratum"], x3[, "tm"]), ]
            x3 <- matrix(subset(x3, tm >= tm_start & tm <= tm_end)[, j], nrow = period_imp )
          
          # take arithmetic means of each set of imputation values
          x1 <- rowMeans(x1,  na.rm=TRUE) 
          x2 <- rowMeans(x2,  na.rm=TRUE) 
          x3 <- rowMeans(x3,  na.rm=TRUE) 
          
          # imputation weights, structured them same way as the imputation values
          w1 <- rep(manual_imputations[i, "ref_geo_unit1_weight"], times=length(x1) )
          w2 <- rep(manual_imputations[i, "ref_geo_unit2_weight"], times=length(x2) )
          w3 <- rep(manual_imputations[i, "ref_geo_unit3_weight"], times=length(x3) )
          
          # compute weighted mean(s) of values
          tmp <- cbind(x1, x2, x3, w1, w2, w3)
          imp_values <- apply(tmp, 1, function(f_x) {weighted.mean(f_x[1:3], f_x[4:6], na.rm=TRUE)} )
          imp_values <- na.replace(imp_values, NA)
            # for each of the strata being imputed, recycle the same vector of imputed values
            imp_values <- rep(imp_values, times = length(strata_imp) )
         
          # attribute imputed values
          ts <- ts[order(ts[, "admin1"], ts[, "stratum"], ts[, "tm"]), ]
          ts1 <- subset(ts[ts[, "stratum"] %in% strata_imp, ], tm >= tm_start & tm <=tm_end)
          ts2 <- subset(ts[ts[, "stratum"] %in% strata_imp, ], tm < tm_start | tm > tm_end)
          ts3 <- ts[! ts[, "stratum"] %in% strata_imp, ]
          ts1[, j] <- as.numeric(imp_values)
          ts <- rbind(ts1, ts2, ts3)
        
        }
      }
      
    }    
  }

  
  #...................................        
  ## Perform automatic imputations to deal with the remaining missingness
    
    # Identify variables to be imputed automatically
    x_vars <- c()
    
    for (i in 1:nrow(predictors) )  {
    
      # if the predictor needs to be interpolated...
      if (predictors[i, "imputation"] == "Y" ) {
        
        # identify predictor...
        x1 <- predictors[i, "worksheet"]
          # if the predictor has already been eliminated from analysis, skip to next loop
          if (! x1 %in% c(predictors$worksheet)) {next;} 
          
        #...and corresponding variables:
        x2 <- subset(dictionary, worksheet == x1 & used_in_analysis == "Y")[, "variable"]
          # but exclude any geo_unit or time_unit variables:
          x2 <- x2[! x2 %in% c("admin1", "stratum", "y", "m", "tm", colnames(strata), colnames(pop), predictors$geo_unit)]
        
        x_vars <- c(x_vars, grep(paste(x2, collapse="|"), colnames(ts), value=TRUE) )
      }
    }  

    # Perform automatic imputation
    
      # restrict dataset to columns without missingness and timeframe
      cols_include <- colnames(ts)[! colnames(ts) %in% c("admin1", "stratum", "y", "m", "tm", colnames(strata), "pop_average", "pop_average_u5", predictors$geo_unit)]
      x1 <- ts[ts$tm %in% c(tm_analysis_start:tm_analysis_end), cols_include]
      x1 <- ts[ts$tm %in% c(tm_analysis_start:tm_analysis_end), ]
      
      # declare methods for imputation (if "", no imputation to be done)
      mice_methods <- rep("", ncol(x1))
      mice_methods[colnames(x1) %in% x_vars] <- mice_method

      # perform imputation (will be done for all missing values across dataset)
      x3 <- mice(x1, method = mice_method, m = mice_it)
      x3 <- complete(x3)    

      # merge back columns to be imputed into main database
      ts <- merge(ts[, ! colnames(ts) %in% x_vars], x3[, c("stratum", "tm", x_vars)],
        by = c("stratum", "tm"), all.x = TRUE)
      
      # check that imputation has worked
      summary(subset(ts, tm %in% c(tm_analysis_start:tm_analysis_end) ))

      
  #...................................  
  ## Interpolate indicators as needed    
  for (i in 1:nrow(predictors) )  {
    # If the predictor needs to be interpolated...
    if (predictors[i, "interpolation"] == "Y" ) {
      
      # identify predictor...
      x1 <- predictors[i, "worksheet"]
        # if the predictor has already been eliminated from analysis, skip to next loop
        if (! x1 %in% c(predictors$worksheet)) {next;} 
        
      #...and corresponding variables:
      x2 <- subset(dictionary, worksheet == x1 & used_in_analysis == "Y")[, "variable"]
        # but exclude any geo_unit or time_unit variables:
        x2 <- x2[! x2 %in% c("admin1", "stratum", "y", "m", "tm", colnames(strata), colnames(pop), predictors$geo_unit)]
      x_vars <- grep(paste(c(x2), collapse="|"), colnames(ts), value=TRUE)
        # if all of the variables have already been eliminated from analysis, skip to next loop
        if (length(x_vars) == 0) {next;}         
      
      # interpolate each variable for each geo_unit, create a new interpolated variable:
      # for each variable to be interpolated...
      for (j in 1:length(x_vars)) {
        # data frame to hold output of interpolation
        ipol <- data.frame("stratum" <- c(), "tm"<- c() )
          # column to hold interpolated variable
          ipol[, paste(x_vars[j], "_ipol", sep="")] <- c() 
          
        # for each stratum...        
        for (k in 1:length(unique(strata[, "stratum"])) ) {
          # sort ts just to be sure
          ts <- ts[order(ts[, "stratum"], ts[, "tm"]), ]
          
          # select time series for stratum k...
          x2 <- subset(ts, stratum == unique(strata[, "stratum"])[k] )[, c("stratum", "tm", x_vars[j])]
          
          # restrict to period
          x2 <- subset(x2, tm >= tm_analysis_start & tm <= tm_analysis_end)
          
          # interpolate the non-missing time series values
            # if the variable is numeric...
            if (typeof(x2[, x_vars[j] ]) != "character") {x3 <- na.approx(x2[, x_vars[j] ], rule = 2) }
            # if the variable is character...
            if (typeof(x2[, x_vars[j] ]) == "character") {
              # default output is the non-interpolated series (which could be all NA)
              x3 <- x2[, x_vars[j] ]
              # if there is at least one non-missing value...
              x4 <- na.omit(x2[, x_vars[j] ])
              if (length(x4) > 0) {
                # start interpolated time series with the first non-missing value
                x3[1] <- x4[1]
                # interpolate the rest by filling missing values with previous non-missing value
                for (k in 1:length(x3) ) {
                  if (is.na(x3[k])) {x3[k] <- x3[k-1]}
                }
              }
            
            }
          
          # merge output with output data frame
          ipol <- rbind(ipol, cbind(x2[, c("stratum", "tm")], x3) )
          
          }
        
        # merge output with main time series
        colnames(ipol) <- c("stratum", "tm", paste(x_vars[j], "_ipol", sep=""))
        ts <- merge(ts, ipol, by=c("stratum", "tm"), all.x=TRUE, sort=TRUE)

        # plot interpolated and original variables, by stratum (if variable is numeric)
        if (typeof(x2[, x_vars[j] ]) != "character") {  
          x1 <- ts
          
          # code if instead wish to plot by admin1
          # x1 <- aggregate(ts[, c(paste(x_vars[j]), paste(x_vars[j], "_ipol", sep=""))],
          #                    by = as.list(ts[, c("admin1", "y", "m", "tm")]), FUN=mean)
          
          # identify the y axis title
          axis_title <- subset(dictionary, variable == x_vars[j])$axis_title
        
          # create a date variable for the x axis
          x1[, "date"] <- dmy(paste("1", x1[, "m"], x1[, "y"], sep="/"))
          
          # create breaks for years
          year_breaks <- subset(x1, m==1)[, "date"]
          
          # draw plot
          plot <- ggplot(x1, aes(x = date) ) +
            geom_point(aes(y = x1[, x_vars[j] ] ), colour = "grey80") +
            geom_line(aes(y = x1[, paste(x_vars[j], "_ipol", sep="")] ), colour = "indianred3") +
            scale_y_continuous(axis_title, labels = function(x) format(x, big.mark = ",", scientific = FALSE)) +
            theme_bw() + theme(plot.margin = unit(c(0.5, 0.5, 1, 0.5), "cm") ) +
            labs(x = "\nmonth", y = paste(axis_title) ) +
            geom_vline(xintercept = year_breaks, color="grey50") +
            facet_wrap(~stratum, ncol=10, scales = "free_y" ) +
            theme(axis.text.x=element_text(angle = -90, vjust=0.5)) +
            scale_x_date("\nmonth - year", 
              limits = c(dmy(paste("1", m_analysis_start, y_analysis_start, sep="/")), 
                dmy(paste("1", m_analysis_end, y_analysis_end, sep="/"))) ,
              expand = c(0,0) , minor_breaks=NULL, date_breaks="3 months", date_labels = "%b-%Y") +
            theme(plot.title = element_text(color="grey30"), legend.title = element_text(color="grey30"),
              axis.title.x = element_text(color="grey30"), axis.title.y = element_text(color="grey30") )
          
          # call plot
          print(plot)
        }
        
        # remove original variable from the main time series
        ts <- ts[, colnames(ts) != x_vars[j] ] 
          
        }
    }
  }
      
        
  #...................................  
  ## Smooth indicators as needed    
  for (i in 1:nrow(predictors) )  {
    
    # If the predictor needs to be smoothed...
    if (predictors[i, "smoothing"] == "Y" ) {
      
      # identify predictor...
      x1 <- predictors[i, "worksheet"]
        # if the predictor has already been eliminated from analysis, skip to next loop
        if (! x1 %in% c(predictors$worksheet)) {next;} 
        
      #...and corresponding variables:
      x2 <- subset(dictionary, worksheet == x1 & used_in_analysis == "Y")[, "variable"]
        # but exclude any geo_unit or time_unit variables:
        x2 <- x2[! x2 %in% c("admin1", "stratum", "y", "m", "tm", colnames(strata), colnames(pop), predictors$geo_unit)]
      x_vars <- grep(paste(x2, collapse="|"), colnames(ts), value=TRUE)
        # if all of the variables have already been eliminated from analysis, skip to next loop
        if (length(x_vars) == 0) {next;}         
      
      # identify start and end points for the predictor's data availability (if no dates are given, assume start/end of analysis period)
      m_start <- ifelse(is.na(predictors[i, "date_start"]) == TRUE, m_analysis_start, month(predictors[i, "date_start"]) )
      m_end <- ifelse(is.na(predictors[i, "date_end"]) == TRUE, m_analysis_end, month(predictors[i, "date_end"]) )
      y_start <- ifelse(is.na(predictors[i, "date_start"]) == TRUE, y_analysis_start, year(predictors[i, "date_start"]) )
      y_end <- ifelse(is.na(predictors[i, "date_end"]) == TRUE, y_analysis_end, year(predictors[i, "date_end"]) )
      
        # convert to corresponding tm values 
        x2 <- as.Date.character(paste(y_start, m_start, "01", sep="-")) - 
          as.Date.character(paste(t_units[t_units$tm == min(t_units$tm), "y"], t_units[t_units$tm == min(t_units$tm), "m"], "01", sep="-"))
        if (x2 <= 0) {tm_start <- min(t_units$tm)}
        if (x2 > 0) {tm_start <- t_units[t_units$y == y_start & t_units$m == m_start, "tm"] }
        x2 <- as.Date.character(paste(y_end, m_end, "01", sep="-")) - 
          as.Date.character(paste(t_units[t_units$tm == max(t_units$tm), "y"], t_units[t_units$tm == max(t_units$tm), "m"], "01", sep="-"))
        if (x2 <= 0) {tm_end <- t_units[t_units$y == y_end & t_units$m == m_end, "tm"] }
        if (x2 > 0) {tm_end <- max(t_units$tm) }
              
      # smooth each variable for each geo_unit, create a new smoothed variable:
      # for each variable to be smoothed...
      for (j in 1:length(x_vars)) {
        # data frame to hold output of smoothing
        smooth <- data.frame("stratum" <- c(), "tm"<- c() )
          # column to hold smoothed variable
          smooth[, paste(x_vars[j], "_smooth", sep="")] <- c() 
          
        # for each stratum...        
        for (k in 1:length(unique(strata[, "stratum"])) ) {
          # sort ts just to be sure
          ts <- ts[order(ts[, "stratum"], ts[, "tm"]), ]
          # select time series for stratum k...
          x2 <- subset(ts, stratum == paste(unique(strata[, "stratum"])[k]) )[, c("stratum", "tm", x_vars[j])]
          # restrict to period of data availability
          x2 <- subset(x2, tm >= tm_start & tm <= tm_end)
          # smooth the non-missing time series values and use the smoothing function to predict across all tm points, 
            # including ones with missing values; only do smoothing if there are at least three values, otherwise all NA
          if (length(na.omit(x2[, x_vars[j]])) > 2) {
            x3 <- predict(smooth.spline(as.matrix(na.omit(x2[, c("tm", x_vars[j])]) ), spar = value_spar), x2[, "tm"] )
          }
          if (length(na.omit(x2[, x_vars[j]])) < 3) {
            x3 <- rep(NA, times = length(x2[, "tm"]) )
          }
          
          # merge output with output data frame
          smooth <- rbind(smooth, cbind(x2[, c("stratum", "tm")], x3$y) )
          
          }
        
        # merge output with main time series
        colnames(smooth) <- c("stratum", "tm", paste(x_vars[j], "_smooth", sep=""))
        ts <- merge(ts, smooth, by=c("stratum", "tm"), all.x=TRUE, sort=TRUE)

        # plot smoothed and unsmoothed variables, by stratum
        
          x1 <- ts
          
          # code if instead wish to plot by admin1
          # x1 <- aggregate(ts[, c(paste(x_vars[j]), paste(x_vars[j], "_smooth", sep=""))],
          #                    by = as.list(ts[, c("admin1", "y", "m", "tm")]), FUN=mean)
          
          # identify the y axis title
          axis_title <- x_vars[j]
        
          # create a date variable for the x axis
          x1[, "date"] <- dmy(paste("1", x1[, "m"], x1[, "y"], sep="/"))
          
          # create breaks for years
          year_breaks <- subset(x1, m==1)[, "date"]
          
          # draw plot
          plot <- ggplot(x1, aes(x = date) ) +
            geom_point(aes(y = x1[, paste(x_vars[j])] ), colour = "grey80") +
            geom_line(aes(y = x1[, paste(x_vars[j], "_smooth", sep="")] ), colour = "indianred3") +
            scale_y_continuous(axis_title, labels = function(x) format(x, big.mark = ",", scientific = FALSE)) +
            theme_bw() + theme(plot.margin = unit(c(0.5, 0.5, 1, 0.5), "cm") ) +
            labs(x = "\nmonth", y = paste(axis_title) ) +
            geom_vline(xintercept = year_breaks, color="grey50") +
            facet_wrap(~stratum, ncol=10, scales = "free_y" ) +
            theme(axis.text.x = element_text(angle = -90, vjust=0.5)) +
            scale_x_date("\nmonth - year", 
              limits= c(dmy(paste("1", m_start, y_start, sep="/")), dmy(paste("1", m_end, y_end, sep="/"))) ,
              expand=c(0,0) , minor_breaks=NULL, date_breaks="3 months", date_labels = "%b-%Y") +
            theme(plot.title = element_text(color="grey30"), legend.title = element_text(color="grey30"),
              axis.title.x = element_text(color="grey30"), axis.title.y = element_text(color="grey30") )
          
          # call plot
          print(plot)
      
        # remove unsmoothed variable from the main time series
        ts <- ts[, colnames(ts) != x_vars[j] ] 
          
        }
    }
  }
      
    
  #...................................
  ## Calculate per-capita rates for predictors that can be expressed this way    
  for (i in 1:nrow(predictors) )  {
    # Identify predictor...
    x1 <- predictors[i, "worksheet"]
    
    # If the predictor can be expressed as a rate...
    if (predictors[i, "analyse_as_rate"] == "Y" ) {
      # if the predictor has already been eliminated from analysis, skip to next loop
      if (! x1 %in% c(predictors$worksheet)) {next;}
      
      # identify columns belonging to the predictor
      x_vars <- subset(dictionary, worksheet == x1 & used_in_analysis == "Y" )$variable
      # exclude geo_unit and time_unit columns
      x_vars <- x_vars[! x_vars %in% c(admin1_name, admin2_name, "y", "m")]
      # include any smoothed variables
      x2 <- na.omit(grep(paste(c(x_vars), collapse="|"), colnames(ts)) )
      x_vars <- colnames(ts)[x2]
      # final check that variables are in fact in the time series
      x_vars <- x_vars[x_vars %in% colnames(ts)]
        # if all of the variables have already been eliminated from analysis, skip to next loop
        if (length(x_vars) == 0) {next;}  
        
      # for each variable that can be converted into a rate...
      for (j in 1:length(x_vars) ) {
        # generate by dividing by population and multiplying by scaling unit
        ts[, paste(x_vars[j], "_rate", sep="")] <- ts[, x_vars[j] ] * predictors[i, "rate_scaling_unit"] / ts[, "pop_average"]
        
        # if the predictor should ONLY be expressed as a rate, remove absolute value variable
        if (predictors[i, "analyse_as_is"] == "N" )
        { ts <- ts[, colnames(ts) != x_vars[j]] }
        
      }
        
    } 
  }
  
    
#.........................................................................................      
### Create secondary variables from rolling means, relative within-stratum levels and to explore lag effects
#.........................................................................................
    
  #...................................    
  ## Calculate running means for indicators for which this is warranted
  for (i in 1:nrow(predictors) )  {
    # identify predictor...
    x1 <- predictors[i, "worksheet"]
      
    # if the predictor needs to be recalculated as a running mean...
    if (predictors[i, "running_mean"] > 1 ) {
      # if the predictor has already been eliminated from analysis, skip to next loop
      if (! x1 %in% c(predictors$worksheet)) {next;}
      
      # identify columns belonging to the predictor
      x_vars <- subset(dictionary, worksheet == x1 & used_in_analysis == "Y" )$variable
      # exclude geo_unit and time_unit columns
      x_vars <- x_vars[! x_vars %in% c(admin1_name, admin2_name, "y", "m")]
      # include any smoothed or rate variables
      x2 <- na.omit(grep(paste(c(x_vars), collapse="|"), colnames(ts)) )
      x_vars <- colnames(ts)[x2]
      # final check that variables are in fact in the time series
      x_vars <- x_vars[x_vars %in% colnames(ts)]
      # if all of the variables have already been eliminated from analysis, skip to next loop
      if (length(x_vars) == 0) {next;}  
      
      # set number of months for mean calculation
      k <- predictors[i, "running_mean"]

      # sort time series
      ts <- ts[order(ts[, "stratum"], ts[, "tm"]), ]
              
      # for each variable that needs to be converted into a running mean...
      for (j in 1:length(x_vars) ) {
        # output vector
        out <-c()
       
        # for each stratum...
        for (l in 1:length(sort(strata[, "stratum"])) ) {
          # select time series for stratum l...
          x1 <- subset(ts, stratum == sort(strata[, "stratum"])[l] )[, x_vars[j] ]
          # calculate rolling mean (leave a burn-in period as NA at the start, until k values are reached)
          out <- c(out, x1[1:(k - 1)], rollmean(x1, k, align="right") )
        }
     
      # add new variable to time series
      ts[, paste(x_vars[j] ,"_rollmean", sep="")] <- out
    
      # remove old variable from the time series
      ts <- ts[, colnames(ts) != x_vars[j] ] 
          
      } 
    }
  }  
    
  #...................................
  ## Create relative values for indicators (compared to value at start of time series) for which this is warranted
  
  for (i in 1:nrow(predictors) )  {
    # Identify predictor...
    x1 <- predictors[i, "worksheet"]
      
    # If a relative value of the predictor needs to be generated...
    if (predictors[i, "relative"] == "Y" ) {
      # if the predictor has already been eliminated from analysis, skip to next loop
        if (! x1 %in% c(predictors$worksheet)) {next;}
        
        # identify columns belonging to the predictor
        x_vars <- subset(dictionary, worksheet == x1 & used_in_analysis == "Y" )$variable
        # exclude geo_unit and time_unit columns
        x_vars <- x_vars[! x_vars %in% c(admin1_name, admin2_name, "y", "m")]
        # include any smoothed, running mean or rate variables
        x2 <- na.omit(grep(paste(c(x_vars), collapse="|"), colnames(ts)) )
        x_vars <- colnames(ts)[x2]
        # final check that variables are in fact in the time series
        x_vars <- x_vars[x_vars %in% colnames(ts)]
        # if all of the variables have already been eliminated from analysis, skip to next loop
        if (length(x_vars) == 0) {next;}  
        
        # sort time series
        ts <- ts[order(ts[, "stratum"], ts[, "tm"]), ]
                
        # for each variable that needs to be converted into a relative value...
        for (j in 1:length(x_vars) ) {
          # output vector
          out <-c()
          
          # for each stratum...
          for (k in 1:length(sort(strata[, "stratum"])) ) {
            # select time series for stratum k...
            x1 <- subset(ts, stratum == sort(strata[, "stratum"])[k] )[, x_vars[j] ]
            # value at start of time series (or first non-missing value)
            x2 <- na.omit(x1)[1]
            # calculate ratio of values:reference and add it to the output vector
            out <- c(out, x1 / x2 )
          }
        # add new variable to time series
        ts[, paste(x_vars[j] ,"_rel", sep="")] <- out
       
        # plot relative variable, by admin1
          x1 <- aggregate(ts[, paste(x_vars[j], "_rel", sep="")], by = as.list(ts[, c("admin1", "y", "m", "tm")]), FUN=mean)
          colnames(x1) <- c("admin1", "y", "m", "tm", paste(x_vars[j] ,"_rel", sep=""))
          
          # create a date variable for the x axis
          x1[, "date"] <- dmy(paste("1", x1[, "m"], x1[, "y"], sep="/"))
          
          # create breaks for years
          year_breaks <- subset(x1, m==1)[, "date"]
          
          # identify the y axis title
          axis_title <- subset(dictionary, 
            worksheet == predictors[i, "worksheet"] & used_in_analysis == "Y" & ! variable %in% c(admin1_name, admin2_name, "y", "m"))$axis_title
          
          # draw plot
          plot <- ggplot(x1, aes(x = date) ) +
            geom_line(aes(y = x1[, paste(x_vars[j], "_rel", sep="")] ), colour = "indianred3") +
            theme_bw() + theme(plot.margin = unit(c(0.5, 0.5, 1, 0.5), "cm") ) +
            labs(x = "\nmonth", y = paste(x_vars[j], "_rel", sep="") ) +
            scale_y_continuous(axis_title, labels = function(x) format(x, big.mark = ",", scientific = FALSE)) +
            geom_vline(xintercept = year_breaks, color="grey50") +
            facet_wrap(~admin1, ncol=3, scales = "free_y") +
            theme(axis.text.x=element_text(angle = -90, vjust=0.5)) +
            scale_x_date("\nmonth, year", 
              limits= c(dmy(paste("1", m_analysis_start, y_analysis_start + pars[which(pars$parameter == "burn.in.period"), "value"], sep="/")),
              dmy(paste("1", m_analysis_end, y_analysis_end, sep="/"))) ,
              expand=c(0,0) , minor_breaks=NULL, date_breaks="3 months", date_labels = "%b-%Y") +
            theme(plot.title = element_text(color="grey30"), legend.title = element_text(color="grey30"),
              axis_title.x = element_text(color="grey30"), axis.title.y = element_text(color="grey30") )
          
          # call plot
          print(plot)

      
        # remove old variable from the time series
        ts <- ts[, colnames(ts) != x_vars[j]] 
          
      } 
    }
  }  
    

  #...................................
  ## Create lags for indicators for which this is warranted
    
    for (i in 1:nrow(predictors) )  {
      # identify predictor...
      x1 <- predictors[i, "worksheet"]
      
      # if the predictor should be analysed as lagged...
      if (predictors[i, "lags"] > 0 ) {
        # if the predictor has already been eliminated from analysis, skip to next loop
        if (predictors[i, "used_in_analysis"] == "N") {next;}
        
        # identify columns belonging to the predictor
        x_vars <- subset(dictionary, worksheet == x1 & used_in_analysis == "Y" )$variable
        # exclude geo_unit and time_unit columns
        x_vars <- x_vars[! x_vars %in% c(admin1_name, admin2_name, "y", "m")]
        # include any smoothed, relative, rate or rolling mean variables
        x_vars <- grep(paste(c(x_vars), collapse="|"), colnames(ts), value=TRUE) 
        # final check that the selected variables are actually in the time series
        x_vars <- x_vars[x_vars %in% colnames(ts)]
          # if all of the variables have already been eliminated from analysis, skip to next loop
          if (length(x_vars) == 0) {next;}  
        
        # sort time series
        ts <- ts[order(ts[, "stratum"], ts[, "tm"]), ]
        
        # create lags
          # name lagged columns
          x_vars_lags <- c()
          for (j in 1:length(x_vars)) {
            x_vars_lags <- c(x_vars_lags, paste(x_vars[j], "_lag", c(1:as.numeric(predictors[i, "lags"])), sep="") )
            }
          
          # use data.table shift function
          ts <- data.table(ts)
          ts[, (x_vars_lags) :=  shift(.SD, 1:as.numeric(predictors[i, "lags"]), fill=NA ), by=stratum, .SDcols=x_vars]
          ts <- as.data.frame(ts)
          
      }
    }  


#.........................................................................................      
### Merge with child observations to create modelling dataset; add a few additional predictors
#.........................................................................................
      
  #...................................               
  ## Add any additional geographic levels from the analysis strata dataset (could be used as predictors or random effects)
  ts <- merge(ts, strata, by=c( colnames(strata)[colnames(strata) %in% colnames(ts)] ), sort = TRUE )  
    

  #...................................               
  ## Add stratum, year and month to child observations
  ch_obs <- merge(ch_obs, surveys[, c("survey_id", "stratum", "year_survey", 
    "month_end")], by = "survey_id", all.x = TRUE)
  colnames(ch_obs)[colnames(ch_obs) %in% c("year_survey", "month_end")] <- c("y", "m")
    
  #...................................               
  ## Merge predictors with child observations
  ch_obs <- merge(ch_obs, ts, by = c("stratum", "y", "m"), all.x = TRUE)  
  
  #...................................               
  ## Add any additional variables pertaining to the survey or its sampling frame, that could be used as predictors   
    
    # Livelihood of surveyed population: if non-missing, override livelihood variable (otherwise assume
      # the general population was surveyed)
  
      # find the index/indices of the original livelihood variable(s), if there are any (if not, NA)
      indices <- grep(paste(livelihood_substrings$livelihood, collapse="|"), colnames(ch_obs) )
      
      # if there is a livelihood variable...
      if (length(indices) > 0 ) {
        # if more than two variables, exclude the one that just contains the ID of the livelihood type/zone
        for (i in 1:length(indices))
            { if (length(grep(paste(1:10, collapse="|"), paste(ch_obs[, indices[i] ]) ))  > 0 ) {indices[i] <- NA} }
        indices <- na.omit(indices)  
        
        # merge in livelihood survey population variable from survey meta-data
        ch_obs <- merge(ch_obs, surveys[, c("survey_id", "livelihood_survey_pop")], 
          by = "survey_id", all.x = TRUE)
        
        # standardise livelihood terms of livelihood survey population variable
        ch_obs[, "livelihood_survey_pop"] <- sapply(ch_obs[, "livelihood_survey_pop"], 
          f_liv, livelihood_substrings)
        
        # override original variable
        ch_obs[, indices] <- ifelse( is.na(ch_obs[, "livelihood_survey_pop"]) == TRUE, ch_obs[, indices], 
          ch_obs[, "livelihood_survey_pop"] )
        
        # now can delete the variable from surveys
        ch_obs <- subset(ch_obs, select = -livelihood_survey_pop)
        
        # if livelihood is displaced, set proportion of IDPs at 100%
        if (length(grep("prop_idp", colnames(ch_obs))) > 0 ) {
          x1 <- grep("prop_idp", colnames(ch_obs))
          ch_obs[ch_obs[, "lhz"] == "displaced", x1] <- 1
        }
      }  
  
      
  #...................................
  ## Add model weights
  ch_obs <- merge(ch_obs, surveys[, c("survey_id", "quality_score", "sampling_coverage")], 
    by = "survey_id", all.x = TRUE)
  ch_obs[, "wt"] <- ch_obs[, "quality_score"] * ch_obs[, "sampling_coverage"]
      
        
  #...................................   
  ## Write data
  write.csv(ch_obs, paste(country, "_ch_obs_model.csv", sep=""), row.names = FALSE)  
        

#.........................................................................................
### ENDS
#.........................................................................................
  
