
### Updated R code and sample data from Borg et al. (Journal) (Year)

--------------------------------------------------------------------------------

Anomalous temperatures increase occupational injuries, illnesses and their associated cost burden in Australia

Borg, Matthew A; Xiang, Jianjun; Anikeeva, Olga; Ostendorf, Bertram; Varghese, Blesson; Dear, Keith; Pisaniello, Dino; Hansen, Alana; Zander, Kerstin; Sim, Malcolm R; Bi, Peng

--------------------------------------------------------------------------------

The numbered code files from *01_DataPrep.R* to *04_AnalysisStage2.R* reproduce the study analysis using the provided data:
  * *01_DataPrep.R* sets up the packages and data sets
  * *02_AnalysisPrep.R* prepares the parameters for statistical analysis
  * *03_AnalysisStage1.R* runs the Stage 1 statistical analysis
  * *04_AnalysisStage2.R* runs the Stage 2 statistical analysis

The data for use with this example 
  * *brambilla.max.rda* is an external meteorological data set sourced from Brambilla et al. 2022, Data in Brief
  * *pop.rda* is the study data set of Australian worker population counts in Adelaide, Brisbane, Darwin, Hobart, Melbourne, Perth and Sydney, sourced from the Australian Bureau of Statistics (ABS). This includes estimates of the number of indoor and outdoor workers based on the ABS Census data.
  * *public.holidays.rda* is the study list of Australian public holidays from 2004 to 2019
  * *school.holidays.rda* is the study list of Australian school holidays from 2004 to 2021 (as a binary variable and a factor variable with levels for each school holiday period)
  * A simulated claims data set is created in file "01_DataPrep.R"
  
Due to the claims data set not being representative of real data, this examples' results are neither similar to those of the main study nor should be used to draw any formal conclusions.

These files can be downloaded after clicking on the green button *Code* at this GitHub repository.
