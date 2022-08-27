

##################################################
### Prepare meteorological data from 
##################################################

# Source: Brambilla et al. 2022, Data in Brief, https://doi.org/10.1016/j.dib.2022.108291

# Extract data
brambilla.loc <- paste0(data.original.loc,'Brambilla 2022/Appendix 1/Climate_files/') # Location of Brambilla et al. 2022 Climate files 
brambilla.data <- list()
for(a in c('Melbourne','Sydney')) {
  for(b in 21:25) { # 2021-2015
    brambilla.data[[paste0(a,b)]] <- as.data.table(cbind(a,b+1990,read_excel(paste0(brambilla.loc,'5_Sydney_Climatedata.xlsx'), sheet=b)[c(1:6,8:10)])) 
  } # Extract required variables, sav in list
}
brambilla <- do.call(rbind.data.frame, brambilla.data) # Combine lists into one data frame (data.table)
colnames(brambilla) <- c('City','Year','Month','Day','Hour','temp','rh','ws','ap','dnr','dhr')

# Obtain metrics at time of maximum temperature using BoM classification
brambilla[,predate:=ISOdate(Year, Month, Day, Hour)] # Date
brambilla[,Date:=as.Date(predate-9*3600+1)] # Bureau of Meteorology (BoM) measuring period for maximum temperature is 9am on same day to 9am next day, date_bom treats 9am as "midnight" starting the measuring period
brambilla.max <- setDT(brambilla)[Date!='2010-12-31', .SD[which.max(temp)], by=c('Date','City')] # Metrics at time of maximum temperature

# Calculate WBGT. (Saturation) vapor pressure and dew point temperature calculations use Buck's 1996 equation and air P
brambilla.max[,sr := dnr+dhr] # Solar radiation.  Brambilla et al. 2022 uses Zhang Huang solar model which accommodates for zenith angle
brambilla.max[,ef := fifelse(temp > 0,  1 + 10^(-4) * (7.2 + ap*(0.0320+5.9*10^(-6)*temp^2)),
                             1 + 10^(-4) * (2.2 + ap*(0.0383+6.4*10^(-6)*temp^2)))] # Enhancement factor, used to account for ratio of (saturation) vapor pressure of moist air to that of pure water vapor
brambilla.max[,svp := ef * fifelse(temp > 0, 6.1121*exp((18.678-temp/234.5)*temp/(257.14+temp)),
                                   6.1115*exp((23.036-temp/333.7)*temp/(279.82+temp)))] # Saturation vapor pressure
brambilla.max[,vp:=rh*svp/100] # Vapor pressure calculating using relative humidity and SVP
brambilla.max[,s:=fifelse(temp > 0,  log(vp/ef/6.1121), log(vp/ef/6.1115))] # Intermediary step for dew point temperature
brambilla.max[,dewpt:= fifelse(temp > 0,  234.5/2 * (18.678 - s - sqrt((18.678-s)^2 - 4 * 257.14 * s / 234.5)),
                               333.7/2 * (23.036 - s - sqrt((23.036-s)^2 - 4 * 279.82 * s / 333.7)))] # Dew point temperature
brambilla.max$wbgtb <- with(brambilla.max, wbgt.Bernard(tas=temp, dewp=dewpt))[["data"]] # Indoor WBGT
brambilla.max$wbgtl <- with(brambilla.max, ifelse(City=='Melbourne', wbgt.Liljegren(tas=temp, dewp=dewpt, wind=ws, radiation=sr, lon=144.84, lat=-37.67, dates=Date, hour=T)[["data"]],
                                                  wbgt.Liljegren(tas=temp, dewp=dewpt, wind=ws, radiation=sr, lon=151.18, lat=-33.95, dates=Date, hour=T)[["data"]])) # Outdoor WBGT
brambilla.max[,':='(predate=NULL,dnr=NULL,dhr=NULL,ap=NULL,ef=NULL,svp=NULL,s=NULL)] # Removing unneeded variables calculation, retaining humidity metrics for sensitivity analyses and components of WBGT



##################################################
### Create simulated claims data
##################################################

# Create randomly-generated claims data with Poisson and Tweedie assumed distributions for number of OIIs (claims) and total costs, respectively

set.seed(3)
dates <- rep(seq(as_date("2011-01-01"), as_date("2015-12-31"), 1),each=4)
ldates <- length(dates)
claims <- data.table(
  Date = dates,
  City = rep(c('Melbourne','Melbourne','Sydney','Sydney'),ldates),
  outin = rep(c('Indoors','Outdoors'),ldates*2),
  no.claims = c(rpois(ldates,60), rpois(ldates,22), rpois(ldates,145), rpois(ldates,50)), # Number of OIIs
  total = c(rTweedie(ldates, p=1.5, phi=8), rTweedie(ldates, p=1.7, phi=9), rTweedie(ldates, p=1.4, phi=11), rTweedie(ldates, p=1.6, phi=7)) # Total costs
) # in rTweedie, p is shape parameter, phi is dispersion parameter

daily <- merge(claims, brambilla.max, by=c('Date','City'))



##################################################
### Public holidays
##################################################

# load(file=paste0(barra.loc,'barra_all.rda'))

# Upload public holidays
public.hols <- as.data.table(read_excel(paste0(data.noriginal.loc,'Public holidays.xlsx'), col_names=T))

# Remove unneeded variables
public.hols.vars.unneeded <- c(
  'Year', 'Month', 'Day', 'DOW', # already have split date variables
  'ACT', 'Regional_Tas', 'Regional_location', # these locations not investigated
  'Number', 'Check_ascending_date', 'Check_date' # not needed
)
public.hols[, (public.hols.vars.unneeded) := NULL] # delete vector of multiple columns

# Verticalise dataset to merge by state. 1 column for Date, City and public holiday
public.hols1 <- melt(public.hols, id.vars=c('Date','Public holiday'),
                     measure.vars=c('SA','QLD','NT','TAS','WA','VIC','NSW'),
                     variable.name='state', value.name='phol')

# Change state to City
public.hols1$City <- state.to.city.fn(public.hols1$state)
public.hols1[, ':='(state=NULL, Date=as.Date(Date))]

# Merge public holidays 
daily1 <- merge(daily, public.hols1[!is.na(phol)], by=c('Date','City'), all.x=T) # merge days without a public holiday
daily1[is.na(phol), phol := 0] # replace NA with 0



##################################################
### School holidays
##################################################

# schoolhols is a binary marker of whether there is school holidays or not (0 = no holiday, 1 = holiday)
# school.hols takes a value of 1 to 4 based on the school holiday periods (for years with 3 terms, the Easter break is longer and counts as a school holiday period), 0 during school term
load(paste0(data.noriginal.loc, 'school_holidays.rda')) # load condensed SWA dataset
daily1 <- merge(daily1, school.hols, by=c('Date','City'), all.x=T) # use date in daily1
rm(school.hols, public.hols1)



#################################################
### Public holiday variables
##################################################

## daily1[str_detect(get('Public holiday'), 'Bank Holiday') & City=='Sydney', phol:=0] # not an actual public holiday in NSW or ACT, but Picnic Day in NT is # REMOVED FROM FILE, NO LONGER NEEDED
daily1[, phol := factor(phol, levels=c(0,1), labels=c('No','Yes'))] # make phol a factor
daily1[, public.hol := phol] # all public holidays, including those removed later

daily1[str_detect(get('Public holiday'), 'hristmas'), ':='(xmas='Yes', phol='No')] # model xmas separately from phol
daily1[is.na(xmas) | xmas!='Yes', xmas:='No'] # model xmas separately from phol
daily1[, xmas := factor(xmas, labels=c('No','Yes'))] # make xmas a factor
daily1[public.hol %in% c(0,'No'), public.hol.x := 0] # public.hol.x is public holiday with 3 values: one for no public holiday, one for public holiday, and another for Xmas. Set if not phol
daily1[public.hol %in% c(1,'Yes'), public.hol.x := fifelse(xmas %in% c(1,'Yes'),2,1)] # Set if phol or xmas
daily1[, public.hol.x := factor(public.hol.x, labels=c('No','Yes','Christmas'))] # make phol an unordered factor

daily1[Month==12 & Day==31, ':='(nye='Yes', phol='No')] # model NYD / 1st Jan (regardless of phol or when observance occurs). Don't model Eve
daily1[is.na(nye) | nye!='Yes', nye:='No'] # fill in No for NYD
daily1[, nye := factor(nye, labels=c('No','Yes'))] # make nyd a factor

daily1[Month==1 & Day==1, ':='(nyd='Yes', phol='No')] # model NYD / 1st Jan (regardless of phol or when observance occurs). Don't model Eve
daily1[is.na(nyd) | nyd!='Yes', nyd:='No'] # fill in No for NYD
daily1[, nyd := factor(nyd, labels=c('No','Yes'))] # make nyd a factor

daily1[, xmasbreak:=as.factor(fifelse(Month==12 & Day %in% c(23:30), 'Yes', 'No'))] # xmaxbreak, associated with less OIs. NYD associated with significantly more than Xmax period, so make it separate
# daily1[, summerbreak:=as.factor(fifelse(Month==12 & Day %in% c(23:31) | Month==1 & Day %in% c(1:6), 'Yes', 'No'))] # summerbreak, associated with less OIs. specific pubhol vars can separate xmas, NYE and NYD
# & (!(dow %in% c('Saturday','Sunday')) | Day==25)

# Special holidays
daily1[, shol:=0] # default value for special holidays
daily1[xmasbreak=='Yes', shol:=1] # summer break
daily1[nye=='Yes', shol:=2] # nye
daily1[nyd=='Yes', shol:=3] # nyd
daily1[Month==1 & Day %in% c(2:4), shol:=4] # after nyd
# daily1[str_detect(`Public holiday`, 'Australia Day'), shol:=4] # Aus day
# daily1[FYear==2017 & Month==6 & Day %in% c(25:30), shol:=5] # last week, from Monday to Saturday
daily1[City=='Brisbane' & str_detect(`Public holiday`, 'Royal Queensland Show') | str_detect(`Public holiday`, 'G20 Leaders'), shol:=5] # Royal Queensland Show and G20 Summit, only affect Brisbane and significant on residual testing
daily1[City=='Sydney' & str_detect(`Public holiday`, 'Australia Day'), shol:=6] # Australia Day has a public celebration at the Sydney Opera House
daily1[City=='Adelaide' & Year==2008 & Month==6 & Day %in% c(24:30), shol:=7] # Adelaide period of 24th Tue - 30th Mon Jun 08, very negative residuals
# daily1[City=='Sydney' &  Year==2018 & Month==6 & Day %in% c(29:30), shol:=8] # Sydney period of 29th Fri - 30th Sat Jun 18, very negative residuals
# daily1[City=='Sydney' &  Year==2018 & Month==6 & Day %in% c(27:30), shol:=8] # Sydney period of 25th Mon - 30th Sat Jun 18, very negative residuals
# daily1[, shol := factor(shol, labels=c('No','Summer break','NYE','NYD','Australia Day','27-30/6/18','Other'))] # make sphol an unordered factor # Other = 'Brisbane only PH', 'Adelaide 24-30/6/08'
daily1[, shol := factor(shol, labels=c('No','Xmas period','NYE','NYD','2-4 Jan','RQS or G20 Leaders','Australia Day','Adelaide 24-30/6/08'))] # make sphol an unordered factor # Other = 'Brisbane only PH', 'Adelaide 24-30/6/08', 'Sydney 27-30/6/18'

# daily1[!(Day==1 & Month!=1), day1:=0] # Not 1st day of the month, excluding New Year's Day
# daily1[Day==1 & Month!=1, day1:=fifelse(FYear<2013, 1, 2)] # 1st day of the month, excluding New Year's Day, subdivide by <2013 and >=2013
# daily1[, day1 := factor(day1, labels=c('No','<2013','>=2013'))] # make factor
daily1[, day1 := as.factor(fifelse(Day==1 & Month!=1, 'Yes','No'))] # 1st day of the month, excluding New Year's Day

daily1[, ldfy := as.factor(fifelse(Day==30 & Month==6, 'Yes','No'))] # last day of financial year

