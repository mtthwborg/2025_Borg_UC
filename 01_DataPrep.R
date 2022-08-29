##################################################
### Load packages
##################################################

# Data manipulation
library(readxl) # Read Excel files (read_excel)
library(data.table) # data.tables. This is used extensively in the subsequent files
library(HeatStress) # WBGT calculations # install.packages("devtools"); devtools::install_github("anacv/HeatStress")
library(lubridate) # Date manipulation
library(stringr) # String commands
library(zoo) # rollapply()

# Statistical analysis
library(mgcv) # Generalized additive models
library(dlnm) # Distributed lag non-linear models
library(splines) # Splines
library(mixmeta) # Multivariate meta-analysis
library(FluMoDL) # Compute attributable risk
library(statmod) # Tweedie distibution
library(tweedie) # Estimate Tweedie shape parameter

# Graphing
library(pryr) # %<a-%



##################################################
### Prepare meteorological data and save as .rda file
###   This step is already done for you and is commented out
###   However, the code is included to demonstrate how the original data was extracted and edited
##################################################

# # Source: Brambilla et al. 2022, Data in Brief, https://doi.org/10.1016/j.dib.2022.108291
# 
# ## Extract data
# brambilla.loc <- paste0("your_directory_with_the_files") # Location of Brambilla et al. 2022 Appendix 1 Climate_files 
# brambilla.data <- list() # Save data in list
# for(a in c('Melbourne','Sydney')) {
#   for(b in 21:25) { # 2021-2015
#     if (a=='Sydney') {brambilla.data[[paste0(a,b)]] <- cbind('City'=a, 'Year'=b+1990, read_excel(paste0(brambilla.loc,'5_Sydney_Climatedata.xlsx'), sheet=b)[c(1:6,8:10)])}
#     else if (a=='Melbourne') {brambilla.data[[paste0(a,b)]] <- cbind('City'=a, 'Year'=b+1990, read_excel(paste0(brambilla.loc,'6_Melbourne_Climatedata.xlsx'), sheet=b)[c(1:6,8:10)])}
#   } 
# }
# 
# brambilla <- do.call(rbind.data.frame, brambilla.data) # Combine list into a data frame
# brambilla <- as.data.table(brambilla) # Save as data.table. This will remove row names
# colnames(brambilla) <- c('City','Year','Month','Day','Hour','temp','rh','ws','ap','dnr','dhr')
# 
# ## Obtain daily metrics at time of maximum temperature using BoM classification
# brambilla[,predate:=ISOdate(Year, Month, Day, Hour)] # Date
# brambilla[,Date:=as.Date(predate-9*3600+1)] # Bureau of Meteorology (BoM) measuring period for maximum temperature is 9am on same day to 9am next day, date_bom treats 9am as "midnight" starting the measuring period
# brambilla.max <- setDT(brambilla)[Date!='2010-12-31', .SD[which.max(temp)], by=c('Date','City')] # Metrics at time of maximum temperature
# 
# ## Calculate WBGT. (Saturation) vapor pressure and dew point temperature calculations use Buck's 1996 equation and air P
# brambilla.max[,sr := dnr+dhr] # Solar radiation.  Brambilla et al. 2022 uses Zhang Huang solar model which accommodates for zenith angle
# brambilla.max[,ef := fifelse(temp > 0,  1 + 10^(-4) * (7.2 + ap*(0.0320+5.9*10^(-6)*temp^2)),
#                              1 + 10^(-4) * (2.2 + ap*(0.0383+6.4*10^(-6)*temp^2)))] # Enhancement factor, used to account for ratio of (saturation) vapor pressure of moist air to that of pure water vapor
# brambilla.max[,svp := ef * fifelse(temp > 0, 6.1121*exp((18.678-temp/234.5)*temp/(257.14+temp)),
#                                    6.1115*exp((23.036-temp/333.7)*temp/(279.82+temp)))] # Saturation vapor pressure
# brambilla.max[,vp:=rh*svp/100] # Vapor pressure calculating using relative humidity and SVP
# brambilla.max[,s:=fifelse(temp > 0,  log(vp/ef/6.1121), log(vp/ef/6.1115))] # Intermediary step for dew point temperature
# brambilla.max[,dewpt:= fifelse(temp > 0,  234.5/2 * (18.678 - s - sqrt((18.678-s)^2 - 4 * 257.14 * s / 234.5)),
#                                333.7/2 * (23.036 - s - sqrt((23.036-s)^2 - 4 * 279.82 * s / 333.7)))] # Dew point temperature
# brambilla.max$wbgtb <- with(brambilla.max, wbgt.Bernard(tas=temp, dewp=dewpt))[["data"]] # Indoor WBGT
# brambilla.max$wbgtl <- with(brambilla.max, ifelse(City=='Melbourne', wbgt.Liljegren(tas=temp, dewp=dewpt, wind=ws, radiation=sr, lon=144.84, lat=-37.67, dates=Date, hour=T)[["data"]],
#                                                   wbgt.Liljegren(tas=temp, dewp=dewpt, wind=ws, radiation=sr, lon=151.18, lat=-33.95, dates=Date, hour=T)[["data"]])) # Outdoor WBGT
# brambilla.max[,':='(predate=NULL,dnr=NULL,dhr=NULL,ap=NULL,ef=NULL,svp=NULL,s=NULL)] # Removing unneeded variables calculation, retaining humidity metrics for sensitivity analyses and components of WBGT
# 
# ## Save daily maximum temperature metrics
# save(brambilla.max, file='brambilla.max.rda')



##################################################
### Create simulated claims data
##################################################

## Create randomly-generated claims data with Poisson and Tweedie assumed distributions for number of OIIs (claims) and total costs, respectively
## This simulated data has no relationship to any predictor variables
no.models <- 4 # 2 cities (Melbourne and Sydney) and their indoor/outdoor combinations
set.seed(7)
dates <- rep(seq(as_date("2011-01-01"), as_date("2015-12-31"), 1),each=no.models)
ldates <- length(dates)/no.models
claims <- data.table(
  Date = dates,
  City = rep(c('Melbourne','Melbourne','Sydney','Sydney'),ldates),
  outin = rep(c('Indoors','Outdoors'),ldates*2),
  'Number of OIIs' = c(rpois(ldates,60), rpois(ldates,22), rpois(ldates,145), rpois(ldates,50)), # Number of OIIs
  'Total costs' = c(rtweedie(ldates, xi=1.5, mu=240, phi=8), rtweedie(ldates, xi=1.7, mu=85, phi=9),
                    rtweedie(ldates, xi=1.4, mu=515, phi=11), rtweedie(ldates, xi=1.6, mu=200, phi=7)) # Total costs
)

# By variables for commands
by.vars <- c("Date","City","outin")
by.vars2 <-by.vars[!by.vars == 'outin']
by.vars3 <- by.vars[!by.vars == 'Date']



##################################################
### Merge data sets together
##################################################

load(file='brambilla.max.rda') # Climate data
load(file='public.holidays.rda') 
load(file='school.holidays.rda')
load(file='pop.rda') # Worker's population using ABS data. Due to estimating indoor/outdoor proportions, n can have decimal places

daily <- merge(merge(merge(claims, brambilla.max, by=c('Date','City')),
               public.holidays[!is.na(phol)], by=c('Date','City'), all.x=T),
                school.holidays, by=c('Date','City'), all.x=T)
daily[is.na(phol), phol := 0] # Replace NA in public holidays with 0

# Merge worker's month.y population
daily.ds <- merge(daily, pop, by=c("Year","Month","City","outin"), all.x=T)

# Create WBGT based on indoors or outdoors
daily.ds[, 'Maximum WBGT':=fifelse(outin=='Outdoors',wbgtl,wbgtb)]

# Stratum variables
daily.ds[,stratum := do.call(paste, c(mget(by.vars3), sep=' '))] # combine by.vars3 into a string as stratum
ds.stratum <- sort(unique(daily.ds[,stratum]))
length.ds.stratum <- length(ds.stratum)



##################################################
### Time variables
##################################################

## Financial year (July to June in Australia)
daily.ds[, FYear := fifelse(Month <= 6, Year-1, Year)]

## Month as a factor variable
daily.ds[,month := factor(Month, labels=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"))]

## Day of the week variables (one with all categories, rest are binary variables per day
daily.ds[, dow:=weekdays(Date)]
daily.ds[, Mon:=fifelse(str_detect(dow,'Mon'), 1, 0)]
daily.ds[, Tue:=fifelse(str_detect(dow,'Tue'), 1, 0)]
daily.ds[, Wed:=fifelse(str_detect(dow,'Wed'), 1, 0)]
daily.ds[, Thu:=fifelse(str_detect(dow,'Thu'), 1, 0)]
daily.ds[, Fri:=fifelse(str_detect(dow,'Fri'), 1, 0)]
daily.ds[, Sat:=fifelse(str_detect(dow,'Sat'), 1, 0)]
daily.ds[, Sun:=fifelse(str_detect(dow,'Sun'), 1, 0)]

## Format public holidays into a factor
daily.ds[, public.hol := factor(phol, levels=c(0,1), labels=c('No','Yes'))]

## Special holidays
daily.ds[, shol:=0] # Default value for special holidays
daily.ds[Month==12 & Day %in% c(23:30), shol:=1] # Christmas break
daily.ds[Month==12 & Day==31, shol:=2] # NYE
daily.ds[Month==1 & Day==1, shol:=3] # NYD
daily.ds[Month==1 & Day %in% c(2:4), shol:=4] # 2nd to 4th January
# daily.ds[City=='Brisbane' & str_detect(`Public holiday`, 'Royal Queensland Show') | str_detect(`Public holiday`, 'G20 Leaders'), shol:=5] # Royal Queensland Show and G20 Summit, only affect Brisbane.Not relevant for this example
daily.ds[City=='Sydney' & str_detect(`Public holiday`, 'Australia Day'), shol:=6] # Australia Day has a public celebration at the Sydney Opera House
# daily.ds[City=='Adelaide' & Year==2008 & Month==6 & Day %in% c(24:30), shol:=7] # Adelaide period of 24th Tue - 30th Mon Jun 08. Not relevant for this example
daily.ds[, shol := factor(shol, labels=c('No','Xmas period','NYE','NYD','2-4 Jan','Australia Day'))] # Make sphol an unordered factor

## 1st day of the month, excluding New Year's Day
daily.ds[, day1 := as.factor(fifelse(Day==1 & Month!=1, 'Yes','No'))]

## Numeric date for gam()
daily.ds[, date := as.numeric(Date)]



######################### END ############################
######################### END ############################
