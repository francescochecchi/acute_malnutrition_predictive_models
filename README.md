## Can we predict the burden of acute malnutrition in crisis-affected countries? Findings from Somalia and South Sudan

# Notes on data and R analysis code
19 August 2021  
Francesco Checchi  
Department of Infectious Disease Epidemiology  
Faculty of Epidemiology and Population Health, Keppel St  
London School of Hygiene and Tropical Medicine  
Francesco.checchi@lshtm.ac.uk  

## Background on the study
This repository contains R scripts (and, if permission is received to upload them, datasets) required to replicate a study to develop models to predict the burden of acute malnutrition in crisis settings. The study was conducted by the London School of Hygiene and Tropical Medicine (www.lshtm.ac.uk) and funded by Unicef. We combined previously collected datasets of ground anthropometric surveys with a variety of 'predictor' data on factors (food security, conflict intensity, climate, health services, etc.) that are theoretically associated causally with acute malnutrition. We explored both generalised linear models and random forest regressions, and for either of these approaches and a range of anthropometric indicators, quantified various metrics of predictive accuracy.

## Datasets
For each of Somalia and South Sudan, the analysis requires the following data files (where 'xxx' = 'som' for Somalia and 'ssd' for South Sudan):
* <xxx_analysis_strata_nut.xlsx>, which contains a list of geographic units;
* <xxx_predictor_data_nut.xlsx>, which contains all the predictor datasets (oner per worksheet), along with a table of predictors with options for their data management;
* <xxx_population_denoms_nut.csv>, which contains population and displacement estimates for each geographic unit and month in the analysis;
* <xxx_survey_metadata_nut.xlsx>, which contains metadata for each of the anthropometric surveys used in the analysis (one per row);
* a sub-dictionary named </survey_datasets> in which all raw datasets of each survey are stored (each dataset is named as per the unique ID in the survey metadata file, e.g. <survey001.csv>;
* <xxx_analysis_parameters_nut.xlsx>, which contains options for general and variable-specific parameters needed at various stages of the analysis.
All data files also include a ‘dictionary’ worksheet, whose column ‘use’ determines which variables are read into R analysis (below). Only variables marked as ‘yes’ for this column are necessary to replicate the analysis.

## R code scripts
The script `malnut_pred_0_control_code.R` loads necessary R packages, loads the above datasets, reads parameters found in these, sets up the analysis and calls dependent scripts (below). Only this script needs to be run to replicate the analysis for any one country and anthropometric indicator. Computation time should be about 20-30 min on a standard laptop. The only modification the user needs to make to this script is to specify the country as 'som' or 'ssd' (line 61).

The dependent script files, called from the above script, are, in order of implementation:
* `malnut_pred_0_functions.R`: this script contains all the user-defined functions (not all are in fact used);
* `malnut_pred_1_manage_surveys.R`: this script cleans, re-analyses and prepares anthropometric survey data for analysis;
* `malnut_pred_2_manage_predictors.R`: this script manages each predictor dataset as desired (combination with population dataset, aggregation, smoothing, interpolation, imputation, lagging, etc.), and merges data with anthropometric observations;
* `malnut_pred_3_predictive_model.R`: this script categorises data, splits them into training and holdout, performs univariate analysis, fits all possible models, does further evaluation on the best-fitting of these, selects the best-fitting, introduces interaction terms and also fits a mixed model if desired; the script computes various metrics of performance accuracy.
* `malnut_pred_4_random_forest.R`: this script grows a random forest model and evaluates its performance. The script does not depend on script 3 having been run.
Scripts 3 and 4 need to be run for each anthropometric indicator of interest.

## Requirements to replicate the analysis
Datasets and R scripts must be placed in the same folder, where output files will be saved. The folder is identified automatically when the control code is run.
It is recommended to run the code from the open-source RStudio interface (https://rstudio.com/ ).
