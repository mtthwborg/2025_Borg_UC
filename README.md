
# Overview
This repository contains R code and datasets that together provide a reproducible example to simulate the main methodology and results from Borg et al. 2022, *Nat. Commun.*, Anomalous temperatures increase occupational injuries, illnesses and their associated cost burden in Australia. These files can be downloaded after clicking on the green button *Code* at this GitHub repository.


# System Requirements
## Hardware requirements
The code should only require a standard computer with enough RAM to suport the user-defined operations. R memory usage when running this code should not exceed 1 GB of RAM.

## Software requirements
This code was tested using the following systems:
* macOS: Monterey Version 12.5.1
* Windows 10: 21H2

This software uses R and was tested using R versions 4.2.0 and 4.2.1. The required packages can be installed using the code below within approximately 30 seconds:
`install.packages(c('readxl','data.table','devtools','lubridate','stringr','zoo',
                   'mgcv','statmod','tweedie','dlnm','mixmeta','FluMoDL'))
devtools::install_github("anacv/HeatStress")`

The versions of the attached non-base packages used during testing were:
` [1] pryr_0.1.5        tweedie_2.3.5     statmod_1.4.37    FluMoDL_0.0.3     mvmeta_1.0.3      mixmeta_1.2.0    
 [7] dlnm_2.4.7        mgcv_1.8-40       nlme_3.1-159      zoo_1.8-10        stringr_1.4.1     lubridate_1.8.0  
[13] HeatStress_1.0.7  magrittr_2.0.3    data.table_1.14.2 readxl_1.4.1`


# Files
This repository comes with the .Rproj (project) file, datasets and analysis files.

## Data
The datasets for use with this example:
  * *brambilla.max.rda* is a derived meteorological data set sourced externally from [Brambilla et al. 2022, *Data Br.*, Hygrothermal climate analysis: An Australian dataset](https://doi.org/10.1016/j.dib.2022.108291).
  * *pop.rda* is the study data set of Australian worker population counts in Adelaide, Brisbane, Darwin, Hobart, Melbourne, Perth and Sydney, sourced from the Australian Bureau of Statistics (ABS). This includes estimates of the number of indoor and outdoor workers as described in the methodology by Borg et al. 2022.
  * *public.holidays.rda* is the study list of Australian public holidays from 2004 to 2023. This dataset is sourced from [https://doi.org/10.25909/6311e7a0dcb3f67](https://adelaide.figshare.com/articles/dataset/Public_holidays_in_Australian_capital_cities_from_2004_to_2023/20732449).
  * *school.holidays.rda* is the study list of Australian school holidays from 2004 to 2023 (as a binary variable and a factor variable with levels for each school holiday period). This dataset is sourced from [https://doi.org/10.25909/6311e7b3bc760](https://adelaide.figshare.com/articles/dataset/School_holidays_in_Australian_capital_cities_from_2004_to_2023/20732173).
  * A simulated claims data set is created in file "01_DataPrep.R"

*pop.rda*, *public.holidays.rda* and *school.holidays.rda* were also used in the main study analysis by Borg et al. 2022.
  
## Analysis
The numbered code files from *01_DataPrep.R* to *04_AnalysisStage2.R* reproduce the study analysis using the aforementioned datasets. They are designed to be run in numerical order:
  * *01_DataPrep.R* sets up the packages and data sets. It creates a simulated (fake) dataset to represent the number of claims and their associated total costs.
  * *02_AnalysisPrep.R* prepares the parameters for statistical analysis. The default parameters replicate those used for the main analysis with total costs as the outcome variable and wet bulb globe temperature as the exposure variable. The parameters can be changed to simulate many of the supplementary analyses.
  * *03_AnalysisStage1.R* runs the Stage 1 statistical analysis, the individual models prior to multivariate meta-analysis.
  * *04_AnalysisStage2.R* runs the Stage 2 statistical analysis. This includes the multivariate meta-analysis, refitting the Stage 1 models with the best linear unbiased predictors, and using the refitted models to generate the main study results.

*02_AnalysisPrep.R* includes a brief loop that can run the code for both *03_AnalysisStage1.R* and *04_AnalysisStage2.R* with multiple outcome variables. Please run *01_DataPrep.R* and the preceding code in *02_AnalysisPrep.R* before running this loop.


# Results
*02_AnalysisPrep.R* will create the folder in the directory to store the results from these latter two files. Within this folder, *03_AnalysisStage1.R* will create an additional folder with the names of the outcome variable and exposure variable, as well as an additional folder called "Stage 1" inside this to replace the results that it generates. *04_AnalysisStage2.R* will create an additional folder named "Stage 2" to store its results.

*.gitignore* is set to exclude the folder *Results* from the working directory.

This example only includes results from 2011 to 2015 and for Melbourne Indoors, Melbourne Outdoors, Sydney Indoors and Sydney Outdoors (four models in total, instead of 14). **Because the example claims dataset was simulated and not representative of real data, this example's results are neither similar to those of the main study nor should be used to draw any formal conclusions.**

## Stage 1
The files generated are:
  * *Reduced coef and model fit.csv*. This contains, for each model, the number of days included for analysis, the reduced coefficients from the distributed lag non-linear models, the AIC and the dispersion parameters. The AIC is summed across models.
  * *Coefficients.csv*. This includes a range of descriptive statistics for the model coefficients. This was used to determine whether the included model predictors were useful predictors or not and to guide modelling choices.
  * *Tweedie shape parameters.csv*. This includes the Tweedie shape parameters from the outcome variables. It is only generated if the outcome has a Tweedie distribution. It is included in case one wishes to repeat the models without needing to re-estimate the shape parameters, which can be time-consuming with large datasets.

## Stage 2
A folder will be created to store results representing each model . Three graphs will be created here to describe the exposure-lag-response relationships prior to multivariate meta-analysis. These 
The files generated are:
  * *MM tests.csv*. This contains the multivariate extension of the Cochran Q test and the I^2^ statistic for the multivariate meta-analysis.
  * *RR (95% CI).csv*. Generates the relative risk and 95% confidence intervals used to generate the overall exposure-response relationship graph.
  * *Overall e-r.png*. Generates the national overall cumulative exposure-response relationship graph.
  * *Overall lag.png*. Generates the national overall cumulative lag-response relationship graph including curves for nine different exposure percentiles ranging from 1st to 99th percentiles.
  * *Attributable fractions.csv*. This includes the proportion of the outcome attributable to heat and cold.
  * *Attributable number.csv*. This includes the number of outcome values attributable to heat and cold.
  * *Location level*. This is a folder including overall cumulative exposure-response graphs for each of the individual city and indoor/outdoor strata combinations, including city-level and indoor/outdoor-levels. There are additional overall cumulative exposure-response graphs for the two-level strata combinations, with a "z" in front of their names, which include histograms representing the spread of the exposure values at that level.
