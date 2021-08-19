#..........................................................................................
### +++++++++++++++ SMALL-AREA PREDICTION OF ACUTE MALNUTRITION BURDEN ++++++++++++++++ ###
#..........................................................................................

#..........................................................................................
## ------------ 'CONTROL' R CODE TO CALL OTHER SCRIPTS AND PERFORM ANALYSIS ------------ ##
#..........................................................................................

                                          # Written by Francesco Checchi, LSHTM (Nov 2020)
                                          # francesco.checchi@lshtm.ac.uk 

#..........................................................................................
### Preparatory steps
#.........................................................................................

  #...................................      
  ## Install or load required R packages

    # List of required packages
    x1 <- c("anthro", "boot", "broom", "broom.mixed", "data.table", "dplyr", "flextable", "GGally", "ggplot2", 
      "ggpubr", "glmmTMB", "gtools", "lme4", "lubridate", "MASS", "mice", "multiwayvcov", "randomForest", 
      "randomForestExplainer", "ranger", "RColorBrewer", "rcompanion", "readxl", "readr", "reshape2", "rlang",  
      "ROCR", "scales", "survey", "tidyr", "tmap", "zoo")
    
    # Install any packages not yet installed
    x2 <- x1 %in% row.names(installed.packages())
    if (any(x2 == FALSE)) { install.packages(x1[! x2]) }

    # Load all packages    
    lapply(x1, library, character.only = TRUE)
    
  #...................................      
  ## Starting setup

    # Clean up from previous code / runs
    rm(list=ls(all=TRUE) )
  
    # Set font
    windowsFonts(Arial=windowsFont("Arial"))

    # Set working directory to where this file is stored
    current_path = rstudioapi::getActiveDocumentContext()$path 
    setwd(dirname(current_path ))
    print( getwd() )
    
      # set sub-directory where survey datasets are stored
      dir_surveys <- paste(getwd( ), "/survey_datasets", sep = "")
    
    # Initialise random numbers
    set.seed(123)
    
    # Colour-blind palette for graphing
    palette_cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
    show_col(palette_cb)

#.........................................................................................
### Specifying parameters
#.........................................................................................
    
  # Select country of analysis
  country <- "som"  # "ssd" for S Sudan, "som" for Somalia
  
  # Substrings to be used for finding livelihood variables and standardising livelihood terms
  livelihood_substrings <- list( 
    "livelihood" <- c("livel", "lz", "LHZ", "lhz", "LZ", "Livel", "lhood", "Lhood"),
    "agriculturalists" <- c("gric", "cultur", "farm", "cultiv", "Farm", "Cultur", "Cultiv"),
    "pastoralists" <- c("past", "Past", "herd", "Herd", "cattle", "Cattle"),
    "agropastoralists" <- c("grop", "grip", "groP", "griP", "gro P", "gri P", "gro-p", "gri-p"),
    "riverine" <- c("river", "River", "riv", "Riv"),
    "fishing" <- c("fish", "Fish"),
    "urban" <- c("urb", "Urb", "city", "City", "town", "Town"),
    "displaced" <- c("displ", "IDP", "camp", "intern", "Displ", "idp", "Camp", "PoC", "poc", "POC"),
    "refugee" <- c("Ref", "ref", "gee")
  )
  
  names(livelihood_substrings) <- c("livelihood", "agriculturalists", "pastoralists", "agropastoralists", 
                                    "riverine", "fishing", "urban", "displaced", "refugee")
    # note: pastoralists must come before agropastoralists as the former have the more ambiguous substrings
  
  
#.........................................................................................
### Calling bespoke functions
#.........................................................................................
  
source("malnut_pred_0_functions.R", echo = TRUE)

 
#.........................................................................................    
### Reading in required files and parameters
#.........................................................................................
  
  #...................................      
  ## General parameters
  gen_pars <- read_excel(paste(country, "_analysis_parameters_nut.xlsx", sep=""), sheet = "general_parameters")
    # remove tibble
    gen_pars <- as.data.frame(gen_pars)
  
    # Assign parameter objects and their values
    gen_pars <- gen_pars[, c("parameter", "value")]
    for (i in 1:nrow(gen_pars)) { 
      if (suppressWarnings(! is.na(as.numeric(gen_pars[i, "value"]) ) ) )
      { assign(gen_pars[i, "parameter"], as.numeric(gen_pars[i, "value"])) } 
      else 
      { assign(gen_pars[i, "parameter"], gen_pars[i, "value"]) }
      
      }

  #...................................      
  ## Predictor-specific parameters            
  var_pars <- read_excel(paste(country, "_analysis_parameters_nut.xlsx", sep=""), sheet = "predictor_parameters")
    # remove tibble
    var_pars <- as.data.frame(var_pars)

  #...................................      
  ## Analysis strata
  strata <- read_excel(paste(country, "_analysis_strata_nut.xlsx", sep=""), sheet = "strata")
    # remove tibble
    strata <- as.data.frame(strata)

    # Rename country-specific geographic units
    colnames(strata)[colnames(strata) == admin2_name] <- "stratum"
    colnames(strata)[colnames(strata) == admin1_name] <- "admin1"
        
    
  #...................................    
  ## Metadata of anthropometric surveys
  surveys <- read_excel(paste(country, "_survey_metadata_nut.xlsx", sep=""), sheet = "survey_metadata")
    # remove tibble
    surveys <- as.data.frame(surveys)

    # Rename country-specific geographic units
    colnames(surveys)[colnames(surveys) == admin2_name] <- "stratum"
    colnames(surveys)[colnames(surveys) == admin1_name] <- "admin1"


  #...................................
  ## Demographic data
  pop <- read.csv(paste(country, "_population_denoms_nut.csv", sep = ""), sep = ",")

    # Rename country-specific geographic units
    colnames(pop)[colnames(pop) == admin2_name] <- "stratum"
    colnames(pop)[colnames(pop) == admin1_name] <- "admin1"    

        
  #...................................   
  ## Predictor data
    # File name that contains these
    filename <- paste(country, "_predictor_data_nut.xlsx", sep="")
    
    # First read predictors table tab, which contains meta-data on all the other predictors...
    predictors <- read_excel(filename, sheet = "predictors_table")
      # exclude any predictors that are not going to be used for analysis
      predictors <- subset(predictors, used_in_analysis == "Y" )
      #get rid of tibble
      predictors <- data.frame(predictors)

    # ...also read variable dictionary which contains variable on whether any variable is to be read (used in analysis)
    dictionary <- read_excel(filename, sheet = "dictionary")
      #get rid of tibble
      dictionary <- data.frame(dictionary)
        
    # ...also read table of manual imputations that may need to be done to fill in specific missing predictor values
    manual_imputations <- read_excel(filename, sheet = "manual_imputations")
      #get rid of tibble
      manual_imputations <- data.frame(manual_imputations)
    
    # ...then read all of the predictors and name them as per the predictors table
    for (i in 1:nrow(predictors) )  {
      # which predictor is being read?
      x1 <- paste(predictors[i, "worksheet"])
      
      # read file
      x2 <- read_excel(filename, sheet = x1 )
      
      # check that import has been successful
      print("+++++++++++++++++++++++++++++++++++++")
      print(paste("now importing predictor...", x1) )
      str(x2)
      
      # only keep variables that will be used for analysis
      x1 <- subset(dictionary, worksheet == x1 & used_in_analysis == "Y")[, "variable"]
      x2 <- x2[, c(unlist(x1)) ]
      
      # rename country-specific geographic units
        colnames(x2) <- gsub(admin2_name, "stratum", colnames(x2) )
        colnames(x2) <- gsub(admin1_name, "admin1", colnames(x2) )
        
      # name the predictor as per the name of the worksheet
      assign(predictors[i, "worksheet"], get("x2") )
    }


#.........................................................................................                            
### Generating time series of stratum-months and other required metadata
#.........................................................................................
  
  #...................................    
  ## Create a time series of stratum-time (including burn-in/-out periods of x years before/after analysis period)

    # create a time unit variable tm (from month 1 to month T of analysis period T) 
    tm <- seq(1, (( y_analysis_end + burn_out_period - y_analysis_start + burn_in_period ) * 12 + m_analysis_end - m_analysis_start + 1 ), 1)
    
    # create a time series of stratum-year-months
    ts <- expand.grid(unlist(strata[, "stratum"] ), tm)
    colnames(ts) <- c("stratum", "tm")
    
    # work out corresponding year and month values
    ts[, "y"] <- floor( (ts[, "tm"] + m_analysis_start - 2) / 12) + y_analysis_start - burn_in_period
    ts[, "m"] <- (ts[, "tm"] + m_analysis_start - 1) - (ts[, "y"] - y_analysis_start + burn_in_period) * 12
    
    # merge admin1 back in
    ts <- merge(ts, strata[, c("stratum", "admin1")], by=c("stratum"), sort=TRUE)
    
    # sort time series
    ts <- ts[order(ts[, "stratum"], ts[, "tm"]), ]
  
      
  #...................................    
  ## Define stratum names and time units
    # Stratum names
    stratum_names <- as.character(unique(strata[, "stratum"]))
    stratum_names <- sort(stratum_names)
    
    # Time units  
    t_units <- unique(ts[, c("tm", "m", "y")])
  
    # Period start and end points
      # overall period of analysis
      tm_analysis_start <- t_units[t_units$y == y_analysis_start & t_units$m == m_analysis_start, "tm"]
      tm_analysis_end <- t_units[t_units$y == y_analysis_end & t_units$m == m_analysis_end, "tm"]
    
   
#.........................................................................................    
### Managing survey data
#.........................................................................................
      
source("malnut_pred_1_manage_surveys.R", echo = TRUE)      
      

#.........................................................................................    
### Managing predictor data
#.........................................................................................
      
source("malnut_pred_2_manage_predictors.R", echo = TRUE)      
      
 
#.........................................................................................    
### Fitting and evaluating predictive model for different anthropometric indicators
#.........................................................................................

source("malnut_pred_3_predictive_model.R", echo = TRUE)       


#.........................................................................................    
### Growing and evaluating random forest for different anthropometric indicators
#.........................................................................................

source("malnut_pred_4_random_forest.R", echo = TRUE)       
      
      
#.........................................................................................
### ENDS
#.........................................................................................
      
