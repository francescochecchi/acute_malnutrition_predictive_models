#..........................................................................................
### +++++++++++++++ SMALL-AREA PREDICTION OF ACUTE MALNUTRITION BURDEN ++++++++++++++++ ###
#..........................................................................................

#..........................................................................................
## ------- CODE TO RE-ANALYSE AND DESCRIBE TRENDS IN SMART ANTHROPOMETRIC SURVEYS ------ ##
#..........................................................................................

                        # Written by Severine Frison and Francesco Checchi, LSHTM (Nov 2020)
                        # francesco.checchi@lshtm.ac.uk 



#.........................................................................................                            
### Cleaning and appending surveys with datasets to the overall child dataset
#.........................................................................................
    
  #...................................    
  ## Set working directory for this part of the code
  setwd(dir_surveys)

  #...................................  
  ## Create (empty) dataset that will contain all child observations from all surveys with datasets
  cols_ch_obs <- c("survey_id", "cluster", "age_months", "gender", "weight", "height", "oedema", "muac")

    # initialise dataset of all survey observations
    ch_obs <- data.frame(matrix(nrow = 0, ncol = length(cols_ch_obs) ))
    colnames(ch_obs) <- cols_ch_obs
  
    
for (i in 1:nrow(surveys) ) {
  
  #...................................  
  ## Read survey id and decide whether it's analysable
  svy_id <- surveys[i, "survey_id"]
    
    # Control message showing progress of loop
    print(paste ("now cleaning data for survey...", i, " of ", nrow(surveys), "...ID ", svy_id, sep=""))

    # If the survey is excluded, skip to next survey
    if (surveys[i, "exclude"] == "Y") { next }


  #...................................
  ## Read survey dataset
  df <- suppressWarnings(read_excel(paste(svy_id, ".xlsx", sep = ""),
    range = cell_cols(1:17), col_types = c("date", "numeric", "numeric", "numeric", "numeric", 
      "text", "date", "numeric", "numeric", "numeric", "text",
      "numeric", "numeric", "numeric", "numeric", "numeric", "numeric")
  ))
    # remove tibble
    df <- as.data.frame(df)
  

  #...................................
  ## Clean survey dataset
    # Rename columns
    colnames(df) <- tolower(colnames(df))
    colnames(df)[colnames(df) == "sex"] <- "gender"
    colnames(df)[colnames(df) == "edema"] <- "oedema"
    colnames(df)[colnames(df) == "months"] <- "age_months"
      
    # Restrict columns
    df <- df[, c("cluster", "gender", "age_months", "weight", "height", "oedema", "muac")]
    
    # Add survey ID
    df[, "survey_id"] <- svy_id
    
    # Correct gender and remove implausible age entries
    df[, "gender"] <- ifelse(df[, "gender"] %in% c(2, "F", "f"), "f", df[, "gender"] )
    df[, "gender"] <- ifelse(df[, "gender"] %in% c(1, "M", "m"), "m", df[, "gender"] )
    df[, "gender"] <- ifelse(! df[, "gender"] %in% c("f", "m"), NA, df[, "gender"] )
    
    df[, "age_months"] <- ifelse(! df[, "age_months"] %in% 6:59, NA, df[, "age_months"] )
    
    # Correct / remove implausible anthropometric readings
    df[, "oedema"] <- ifelse(df[, "oedema"] %in% c(2, "N", "No", "no", "nn", "n", "0=No"), "n", df[, "oedema"] )    
    df[, "oedema"] <- ifelse(df[, "oedema"] %in% c(1, "Y", "Yes", "yes", "yy", "y", "1=Yes"), "y", df[, "oedema"] )    
    df[, "oedema"] <- ifelse(! df[, "oedema"] %in% c("y", "n"), NA, df[, "oedema"] )
    
    df[, "muac"] <- suppressWarnings(as.numeric(df[, "muac"]) )
    df[, "muac"] <- ifelse(! df[, "muac"] %in% 85:200, NA, df[, "muac"] )
    df[, "muac"] <- df[, "muac"] / 10 # anthro needs it in cm

    df[, "height"] <- suppressWarnings(as.numeric(df[, "height"]) )
    df[, "weight"] <- suppressWarnings(as.numeric(df[, "weight"]) )
    
    # Resolve cluster variable: if missing across the entire database, assume exhaustive or random sampling
    df[, "cluster"] <- suppressWarnings(as.integer(df[, "cluster"]) )
      # if missing across the entire database, assume exhaustive or random sampling
      if (length(unique(na.omit(df[, "cluster"]))) <= 1) { df[, "cluster"] <- paste(svy_id, 99999, sep = "")}
      # otherwise, create a unique cluster ID across all surveys
      if (length(unique(na.omit(df[, "cluster"]))) > 1) { df[, "cluster"] <- paste(svy_id, df[, "cluster"], sep = "")}
      

  #...................................
  ## Append to child observations
  ch_obs <- rbind(ch_obs, df[, cols_ch_obs])  
    
}    


#.........................................................................................                            
### Generating anthropometric indices and classifications
#.........................................................................................
     
  #...................................
  ## Generate anthropometric indices
    # Generate indices using anthro package
    indices <- with(ch_obs, 
      anthro_zscores( sex = gender, age = age_months, is_age_in_month = TRUE,
        weight = weight, lenhei = height, oedema = oedema, armc = muac
      )
    )

    # Attribute to database, including flags  
    ch_obs[, c("wfhz", "hfaz", "wfaz", "mfaz")]  <- indices[, c("zwfl", "zlen", "zwei", "zac")]
    ch_obs[, paste("flag_", c("wfhz", "hfaz", "wfaz", "mfaz"), sep = "")]  <- 
      indices[, c("fwfl", "flen", "fwei", "fac")]
    
    # Rectify flags to also include any instance of Z scores < -5
    for (j in c("wfhz", "hfaz", "wfaz", "mfaz") ) {
      ch_obs[, paste("flag_", j, sep = "")] <- ifelse(ch_obs[, j] < (-5), 1, ch_obs[, paste("flag_", j, sep = "")] )
    }
    
    # Identify all observations with at least 1 flag
    ch_obs[, "flag"] <- rowSums(ch_obs[ grep("flag", colnames(ch_obs))], na.rm = TRUE)
    ch_obs[, "flag"] <- ifelse(ch_obs[, "flag"] >= 1, "y", "n")

  
  #...................................
  ## Create GAM and SAM classifications
    
    # GAMZ and SAMZ (i.e. according to WFH and/or oedema)
    ch_obs[, "gamz"] <- NA
    ch_obs[, "gamz"] <- ifelse(ch_obs[, "wfhz"] < (-2) | ch_obs[, "oedema"] == "y", 1, ch_obs[, "gamz"])
    ch_obs[, "gamz"] <- ifelse(ch_obs[, "wfhz"] >= (-2) & ch_obs[, "oedema"] == "n", 0, ch_obs[, "gamz"])
  
    ch_obs[, "samz"] <- NA
    ch_obs[, "samz"] <- ifelse(ch_obs[, "wfhz"] < (-3) | ch_obs[, "oedema"] == "y", 1, ch_obs[, "samz"])
    ch_obs[, "samz"] <- ifelse(ch_obs[, "wfhz"] >= (-3) & ch_obs[, "oedema"] == "n", 0, ch_obs[, "samz"])
  
    # GAMM and SAMM (i.e. according to MUAC and/or oedema)
    ch_obs[, "gamm"] <- NA
    ch_obs[, "gamm"] <- ifelse(ch_obs[, "muac"] < 12.5 | ch_obs[, "oedema"] == "y", 1, ch_obs[, "gamm"])
    ch_obs[, "gamm"] <- ifelse(ch_obs[, "muac"] >= 12.5 & ch_obs[, "oedema"] == "n", 0, ch_obs[, "gamm"])
  
    ch_obs[, "samm"] <- NA
    ch_obs[, "samm"] <- ifelse(ch_obs[, "muac"] < 11.5 | ch_obs[, "oedema"] == "y", 1, ch_obs[, "samm"])
    ch_obs[, "samm"] <- ifelse(ch_obs[, "muac"] >= 11.5 & ch_obs[, "oedema"] == "n", 0, ch_obs[, "samm"])
  

#.........................................................................................                            
### Re-analysing surveys to generate point estimates and confidence intervals
#.........................................................................................
      
  #...................................  
  ## Assign extra columns of meta-data dataframe that will store re-analysis results for each survey
  out_col_names <- c(
    # survey design based on re-analysis
    "survey_design",
    # GAM prevalence based on whz + oedema, according to anthro (WHO flags)
    "gamz_est", "gamz_lci", "gamz_uci", 
    # SAM prevalence based on whz + oedema, according to anthro (WHO flags)
    "samz_est", "samz_lci", "samz_uci",
    # GAM prevalence based on MUAC + oedema
    "gamm_est", "gamm_lci", "gamm_uci", 
    # SAM prevalence based on MUAC + oedema
    "samm_est", "samm_lci", "samm_uci",
    # mean WHZ, according to anthro (WHO flags)
    "wfhz_est", "wfhz_lci", "wfhz_uci",
    # mean MUAC
    "muac_est", "muac_lci", "muac_uci",
    # mean MUAC-for-age Z score
    "mfaz_est", "mfaz_lci", "mfaz_uci",
    # proportion of flags
    "prop_flag",
    # sample size
    "sample_size"
    )
  
  surveys[, out_col_names] <- NA
  
  
for (i in 1:nrow(surveys) ) {
  
  #...................................  
  ## Read survey id and decide whether it's analysable
  svy_id <- surveys[i, "survey_id"]
    
    # Control message showing progress of loop
    print(paste ("now re-analysing survey...", i, " of ", nrow(surveys), "...ID ", svy_id, sep=""))

    # If the survey is excluded, skip to next survey
    if (surveys[i, "exclude"] == "Y") { next }


  #...................................
  ## Capture and prepare survey dataset
    
    # Capture dataset
    df <- subset(ch_obs, survey_id == svy_id)
  
    # Specify survey design (always exclude observations with any flags)
      
      # if all cluster values == 99999, assume SRS or exhaustive (latter to be treated as an SRS)
      if (length(unique(df[, "cluster"])) == 1) {
        surveys[i, "survey_design"] <- "SRS or exhaustive";
        survey_design <-  suppressWarnings(svydesign(id = ~0, data = subset(df, flag == "n")) )
      }

      # otherwise, assume cluster design
      if (length(unique(df[, "cluster"])) > 1) {
      surveys[i, "survey_design"] <- "multi-stage cluster";
      survey_design <- suppressWarnings(
        svydesign(id = ~cluster, data = subset(df, ! is.na(cluster) & flag == "n") ) )
      }

    # Estimate point estimates and 95%CIs for continuous variables 
    for (j in c("wfhz", "muac", "mfaz")) {
 
      # fit survey GLM
      fit <- svyglm(as.formula(paste(j, "~", "NULL", sep = " ")), survey_design, family = "gaussian")
      
      # compute point estimate and 95%CIs      
      surveys[i, paste(j, "_est", sep = "")] <- summary(fit)$coefficients[[1]]
      surveys[i, paste(j, "_lci", sep = "")] <- summary(fit)$coefficients[[1]] - 1.96 * summary(fit)$coefficients[[2]]
      surveys[i, paste(j, "_uci", sep = "")] <- summary(fit)$coefficients[[1]] + 1.96 * summary(fit)$coefficients[[2]]
    }
    
    # Estimate point estimates and 95%CIs for binary variables (always exclude observations with any flags)
    for (j in c("gamz", "samz", "gamm", "samm")) {

      # fit survey GLM
      fit <- svyglm(as.formula(paste(j, "~", "NULL", sep = " ")), survey_design, family = "binomial")
      
      # compute point estimate and 95%CIs      
      surveys[i, paste(j, "_est", sep = "")] <- exp(summary(fit)$coefficients[[1]])
      surveys[i, paste(j, "_lci", sep = "")] <- 
        exp(summary(fit)$coefficients[[1]] - 1.96 * summary(fit)$coefficients[[2]] )
      surveys[i, paste(j, "_uci", sep = "")] <- 
        exp(summary(fit)$coefficients[[1]] + 1.96 * summary(fit)$coefficients[[2]] )
    }
    
    # Compute the proportion of flags
    surveys[i, "prop_flag"] <- table(df[, "flag"])["y"] / nrow(df)
    if (length(grep("y", df[, "flag"])) == 0) {surveys[i, "prop_flag"] <- 0}
    
    # Compute the survey's sample size
    surveys[i, "sample_size"] <- nrow(df)
    
}    
 

#.........................................................................................
### Saving output of reanalysis
#.........................................................................................  
  
  #...................................    
  ## Write child observations database
    # Set working directory
    setwd(dirname(current_path ))
    
    # Remove flags
    ch_obs <- subset(ch_obs, flag == "n")
    ch_obs <- ch_obs[, - grep("flag", colnames(ch_obs))]
    
    # Write
    write.csv(ch_obs, paste(country, "_ch_obs.csv", sep=""), row.names = FALSE)

  #...................................    
  ## Write survey metadata with results of re-analysis
    # Format estimates
    x1 <- c("gamz_est", "gamz_lci", "gamz_uci", "samz_est", "samz_lci", "samz_uci",
      "gamm_est", "gamm_lci", "gamm_uci", "samm_est", "samm_lci", "samm_uci", "prop_flag" )
    surveys[, x1] <- round(surveys[, x1] * 100, 1)
    
    x1 <- c("wfhz_est", "wfhz_lci", "wfhz_uci", "mfaz_est", "mfaz_lci", "mfaz_uci")
    surveys[, x1] <- round(surveys[, x1], 3)

    x1 <- c("muac_est", "muac_lci", "muac_uci")
    surveys[, x1] <- round(surveys[, x1] * 10, 1)
    
    # Write
    write.csv(surveys, paste(country, "_survey_metadata_reanalysed.csv", sep = ""), row.names = FALSE)
    

#.........................................................................................
### Computing summary statistics from mortality surveys, overall and by year
#.........................................................................................
    
  #...................................  
  ## Create table that will hold output
  df <- subset(surveys, exclude == "N")  
  years <- sort(as.numeric(unique(df[, "year_survey"])))
  x1 <- c("overall", as.character(years) )
  x2 <- c("eligible surveys (N)", "percentage using cluster design", "sample size (n)", "GAM prevalence (WFH)", "SAM prevalence (WFH)", 
    "GAM prevalence (MUAC)", "SAM prevalence (MUAC)", "percentage of flags")
  out <- data.frame(matrix(NA, nrow = length(x2), ncol = length(x1) ) )
  colnames(out) <- x1
  rownames(out) <- x2
  
  #...................................          
  ## Filling table
  out["eligible surveys (N)", ] <- c(nrow(df), table(df[, "year_survey"]) )
  out["percentage using cluster design", ] <- round(c(
    table(df[, "survey_design"])["multi-stage cluster"] / nrow(df),
    table(df[, "year_survey"], df[, "survey_design"])[, "multi-stage cluster"] / table(df[, "year_survey"])
    ) * 100, 1)
  
  out["sample size (n)", ] <- f_calc_svy("sample_size", df, 0, years)
  out["GAM prevalence (WFH)", ] <- f_calc_svy("gamz_est", df, 1, years)
  out["SAM prevalence (WFH)", ] <- f_calc_svy("samz_est", df, 1, years)
  out["GAM prevalence (MUAC)", ] <- f_calc_svy("gamm_est", df, 1, years)
  out["SAM prevalence (MUAC)", ] <- f_calc_svy("samm_est", df, 1, years)
  out["percentage of flags", ] <- f_calc_svy("prop_flag", df, 1, years)

  #...................................  
  ## Writing output to a file
  out <- out[, c("overall", sort(as.character(years)) )]
  write.csv(out, paste(country, "_out_crude_survey_statistics.csv", sep="") )    

    
#.........................................................................................
### Producing graphs and figures of interesting statistics
#.........................................................................................
 
  #...................................   
  ## Survey quality score by year  
  plot_quality <- ggplot(df, aes(y = quality_score, x = as.character(year_survey)) ) +
    geom_boxplot(colour = palette_cb[3], size = 1, fill = palette_cb[3], alpha = 0.5) +
    theme_bw() +
    theme(axis.title = element_text(size = 10, colour = "grey20")) +
    scale_y_continuous("survey quality score") +
    scale_x_discrete ("")
  plot_quality
      
  plot_flags <- ggplot(df, aes(y = prop_flag, x = as.character(year_survey)) ) +
    geom_boxplot(colour = palette_cb[2], size = 1, fill = palette_cb[2], alpha = 0.5) +
    theme_bw() +
    theme(axis.title = element_text(size = 10, colour = "grey20")) +
    scale_y_continuous("percentage of flagged observations") +
    scale_x_discrete ("")
  plot_flags  
  
  #...................................  
  ## Trends over time in key demographic indicators
    # List of indicators of interest
    indicators <- c("gamz_est", "samz_est", "gamm_est", "samm_est")
    names(indicators) <- c(
      "GAM prevalence (WFH + oedema), %",
      "SAM prevalence (WFH + oedema), %",
      "GAM prevalence (MUAC + oedema), %",
      "SAM prevalence (MUAC + oedema), %"
    )
    
    # Prepare dataset for plotting
      # create dates
      df[, "month_survey"] <- df[, "month_end"]
      df[, "day_survey"] <- ifelse(df[, "month_start"] != df[, "month_end"], 1, 15) 
      df[, "date_survey"] <- ymd(paste(df[, "year_survey"], df[, "month_survey"], df[, "day_survey"], sep = "-"))
      
      # merge with admin0
      if ("admin1" %in% colnames(df)) {df <- merge(df, strata, by = c("admin1", "stratum"), all.x = TRUE) }
      if (! "admin1" %in% colnames(df)) {df <- merge(df, strata, by = "stratum", all.x = TRUE) }
      
    # Plot for each indicator  
    for (i in 1:length(indicators) )  {
      # prepare
      x1 <- df
      colnames(x1)[colnames(x1) == indicators[i]] <- "var"
      colours <- palette_cb[c(4,7,6,8)][1:length(unique(df$admin0))]
      
      # plot
      plot <- ggplot(x1, aes(x = date_survey, y = var, group = admin0, colour = admin0, fill = admin0) ) +
        geom_point(alpha = 0.7, size= 3) + 
        theme_bw() + 
        # theme(plot.margin = unit(c(0.5, 0.5, 1, 0.5), "cm") ) +
        scale_colour_manual(values = colours ) +
        scale_fill_manual(values = colours ) +
        guides(colour = guide_legend(title = "region"), fill = guide_legend(title = "region") ) +
        scale_x_date("", date_breaks = "6 months", date_labels = "%b-%Y", 
          limits = c(ymd(paste(y_analysis_start, m_analysis_start, 1, sep = "-")), 
            ymd(paste(y_analysis_end, m_analysis_end, 1, sep = "-"))) ) +
        scale_y_continuous(names(indicators)[i] ) +
        theme(legend.title = element_text(color = "grey20", size = 10),
          legend.position = "bottom",
          axis.title = element_text(color = "grey20", size = 10),
          axis.text.x = element_text(color = "grey20", size = 10, angle = 30, hjust = 1),
          axis.text.y = element_text(color = "grey20", size = 10, hjust = 1),
          legend.text = element_text(color = "grey20", size = 10)
        )
      print(plot)
      
      # assign plot name and save
      # assign name
      assign(paste("plot_", indicators[i], sep = ""), plot)
    }


    # Composite plot showing all indicators
    plot <- ggarrange(
      get(paste("plot_", indicators[1], sep = "")) + theme(axis.text.x=element_blank() ), 
      get(paste("plot_", indicators[2], sep = "")) +  theme(axis.text.x=element_blank() ), 
      get(paste("plot_", indicators[3], sep = "")) ,
      get(paste("plot_", indicators[4], sep = "")) ,
      plot_quality , 
      plot_flags , 
      ncol = 2, nrow = 3, common.legend = TRUE, legend = "bottom", align = "v") 
    plot
    
    ggsave(paste(country, "_svy_trends.png", sep = ""), height = 32, width = 25, units = "cm", dpi = "print")
        
  
#.........................................................................................
### ENDS
#.........................................................................................
