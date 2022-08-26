

#################################################
### Create sample dataset
##################################################

set.seed(3)
dates <- rep(seq(as_date("2005-07-01"), as_date("2018-06-30"), 1),each=6)
ldates <- length(dates)
claims <- data.table(
  Date = dates,
  City = rep(c('Brisbane','Brisbane','Melbourne','Melbourne','Sydney','Sydney'),ldates),
  outin = rep(c('Indoors','Outdoors'),ldates*3),
  no.claims = c(rpois(ldates,60), rpois(ldates,25), rpois(ldates,60), rpois(ldates,22), rpois(ldates,145), rpois(ldates,50)),
  no.claims = c(rtweedie(ldates,60), rtweedie(ldates,25), rtweedie(ldates,60), rpois(ldates,22), rpois(ldates,145), rpois(ldates,50)),
)
# Means and dispersion p


daily[City %in% c('Brisbane','Melbourne','Sydney'), summary(.SD), by=mget(by.vars3), .SDcols = 'no.claims']


# OI refers to occupational illness, IO refers to indoor/outdoor

# Percentile values
predper <- c(seq(0,1,0.1),2,2.5,3:33,100/3,34:66,200/3,67:97,97.5,98,seq(99,100,0.1)) # percentiles (%). 100/3 and 200/3 work even without rounding
predper.short <- c(1,2.5,10,25,50,75,90,97.5,99) # c(1,10,90,99) # percentiles to report RRs for
rTweedie(ldates,60, 1.5)
rtweedie(ldates,60, 1.5)
rtweedie(ldates,60, 1.7)
rTweedie(mu,p=1.5,phi=1.3)

f2 <- function(x) 0.2 * x^11 * (10 * (1 - x))^6 + 10 *
  (10 * x)^3 * (1 - x)^10
n <- 300
x <- runif(n)
mu <- exp(f2(x)/3+.1);x <- x*10 - 4
y <- rTweedie(mu,p=1.5,phi=1.3)
b <- gam(y~s(x,k=20),family=Tweedie(p=1.5))
b
plot(b) 

rtweedie(20,xi=1.5,mu=1,phi=1)


##################################################
### MODEL, INDOOR/OUTDOOR AND DISTRIBUTION
##################################################

# Indoors/outdoors marker and variable
io.vars <- c('outorin','Outorin','outdoors','Outdoors','outdoor','Outdoor','indoors','Indoors','Indoor','Indoor','outin','outorini','outini','Outin','Outorini','Outini')
if (any(by.vars %in% io.vars)) { # if a by-variable is indoors or outdoors
  io.marker <- T # flag marker to stratify by indoor/outdoor
  by.vars.oi <- by.vars[by.vars %in% io.vars] # save name of outorin var
} else {
  io.marker <- F # don't stratify by indoor/outdoor
  by.vars.oi <- NULL # remove name of outorin var
}

# Heatwave marker
if(str_detect(exposure.var, 'xcess heat '))  { # if excess heat factor or heatwave
  ehf.marker <- T
} else {
  ehf.marker <- F
}

# Type of outcome
if (str_detect(outcome.var,'claims') | str_detect(outcome.var,'injur') | str_detect(outcome.var,'disease') | str_detect(outcome.var,'llness')) {
  type.outcome <- 'oi' # occupational ilnesses
} else if (any(str_detect(outcome.var, c('non-zero','non-0','no zero','no 0'))) | isTRUE(plus1)) {
  type.outcome <- 'non-0' # costs without 0 values
} else {
  type.outcome <- 'cost' # costs with 0 values
}

# Outcome distribution
if (distribution.choice %in% c('auto','auto2','automatic','default',NULL)) {
  if (type.outcome=='oi') { # if number of claims/injuries/diseases
    if (distribution.choice=='auto2') {distribution <- 'poisson'}
    else {distribution <- 'quasipoisson'}
  } else if (type.outcome=='non-0') { # if only include clams with costs
    distribution <- Gamma(link = log)
  } else if (type.outcome=='cost') { # if cost
    distribution <- tweedie(var.power=1.7, link.power=0) # parameters restimated later; initial choice makes no difference
    tweedie.profile. <- var.power <- list() # if tweedie distribution, create tweedie.profile
  } else {
    distribution <- NULL
    print('NO DISTRIBUTION ASSIGNED')
  }
} else {
  distribution <- distribution.choice
}

# Model type (GLM, GAM, GNM)
if(model.type %in% c('a','gam','GAM')) { # if gam, assign gam
  modeltype <- 'gam'
} else if(model.type=='auto') { # if auto
  if(type.outcome=='oi') { # if number of claims/injuries/diseases
    modeltype <- 'glm' # use glm
  } else {
    modeltype <- 'gam' # use gam
  }
} else {
  modeltype <- model.type # use specified model.type
}

# Select gam or non-gam trend as appropriate
if(modeltype %in% c('a','gam','GAM') ) {  
  trnd <- trend.gam
} else { # glm (Date or date)
  trnd <- trend
}

# Add trnd to equations
mformula <- paste0(m.formula, trnd)
mformula.io <- paste0(m.formula.io, trnd)

# If seasonal trend df is automatic. Only relevant if adf is used as knots
if (select.sdf %in% c('auto','automatic','default',NULL)) {
  if(type.outcome=='oi') { # if number of claims/injuries/diseases
    s.sdf <- 'oi'
  } else {
    s.sdf <- 'cost.glm'
  }
} else {
  s.sdf <- select.sdf
}

# If humidity needs to be adjusted for a gam
if (modeltype=='gam') {
  if(!is.null(humid)) { # if humid present
    if (str_detect(humid,'ns\\(')) { # if humid has ns
      hum <- str_replace_all(humid, 'ns\\(', 's\\(') # replace ns with s
      hum <- str_replace_all(hum,',',",bs='cr',fx=T,k=") # add gam ns
    } else if (str_detect(humid,'bs\\(')) {  # if humid has bs
      hum <- str_replace_all(humid, 'bs\\(', 's\\(') # replace bs with s
      hum <- str_replace_all(hum,',',",bs='bs',fx=T,k=") # add gam ns
    }
  }
}
if(!is.null(humid)) {hum <- humid # if humid exists, assign hum
} else {hum <- NULL} # else remove hum

# Add trnd and (if present) hum to equations
mformula <- paste0(m.formula, trnd, hum)
mformula.io <- paste0(m.formula.io, trnd, hum)


##################################################
### Yaxis for graphs as percentages. Will be adjusted to RR in code as required
##################################################

if (rr.yaxis %in% c('auto','automatic','default',NULL)) { # affect individual model plots
  if(type.outcome=='oi') { # if number of claims/injuries/diseases
    rryaxis <- c(-20,60)
    rryaxis.s <- c(rryaxis[1]-40,rryaxis[2]) # rryaxis*1.5
  } else {
    rryaxis <- c(-20,60)
    rryaxis.s <- c(rryaxis[1]-40,rryaxis[2]) # rryaxis*1.5
  }
} else {
  rryaxis <- rr.yaxis
  rryaxis.s <- c(rryaxis[1]-40,rryaxis[2]) # rryaxis*1.5
}

if (bigvar.yaxis[1] %in% c('auto','automatic','default',NULL)) { # affect bigvar plots
  bigvaryaxis <- rryaxis
} else {
  bigvaryaxis <- bigvar.yaxis
}
# c(-20,30) appropriate for OIs, c(-10,60) for costs. Want to keep axes identical

# overall plots determined per plot



##################################################
### Folder and subfolder destinations
##################################################

### Text description of model

# Shorten outcome text
outcome.var.short <- gsub(" ", "", str_remove_all(outcome.var, ' \\(000s\\)'), fixed=T) # remove (000s) and spaces to shorten folder name
outcome.var.short <- str_replace_all(outcome.var.short, 'up to Jun 14',' fy5') # if up to fyear5, shorten name
outcome.var.short <- str_replace_all(outcome.var.short, 'in same financial year','fy0') # if fyear0 only, shorten name
outcome.var.short <- str_replace_all(outcome.var.short, 'oodsandservices','S') # Goodsandservices to GS
outcome.var.short <- str_replace_all(outcome.var.short, 'ealthservices','S') # Healthsandservices to GS
outcome.var.short <- str_replace_all(outcome.var.short, '-compensation','comp') # Non-compensation to GS
outcome.var.short <- str_replace_all(outcome.var.short, 'ompensation','omp') # Compensation to GS
outcome.var.short <- str_replace_all(outcome.var.short, 'umberof','umber') # Remove "number of" before injuries/illnesses/diseases
outcome.var.short <- str_remove_all(outcome.var.short, '\\$')

# Shorten exposure text
exposure.var.short <- gsub(" ", "", exposure.var, fixed=T)
exposure.var.short <- str_replace_all(exposure.var.short, 'Average','Ave')
exposure.var.short <- str_replace_all(exposure.var.short, 'Mean','Ave')
exposure.var.short <- str_replace_all(exposure.var.short, 'Maximum','Max')
exposure.var.short <- str_replace_all(exposure.var.short, 'Minimum','Min')
exposure.var.short <- str_replace_all(exposure.var.short, "temperature",'T') # shortening
exposure.var.short <- str_replace_all(exposure.var.short, "Steadman'sapparent",'SteadmanA') # shortening

# Filter check
if(is.null(filter.conds) | is.na(filter.conds) | filter.conds=='') { # if no filter
  by.filter.text <- NULL
} else {
  by.filter.text <- paste(filter.conds, collapse = '')
  by.filter.text <- str_remove_all(by.filter.text, 'inj.or.dis==') # takes up quite a lot of room otherwise
  by.filter.text <- str_replace_all(by.filter.text, 'Diseases and conditions','DaC') # takes up quite a lot of room otherwise
  by.filter.text <- str_remove_all(by.filter.text, "\\'") # remove symbols that may appear
  by.filter.text <- str_remove_all(by.filter.text, '\\"')
  by.filter.text <- str_remove_all(by.filter.text, '\\<')
  by.filter.text <- str_remove_all(by.filter.text, '\\>')
  by.filter.text <- str_remove_all(by.filter.text, '\\=')
  by.filter.text <- str_remove_all(by.filter.text, '\\$')
  by.filter.text <- str_remove_all(by.filter.text, ' ') # space
}

if(!is.null(strata.vars)) {
  by.strata.text <-  paste(strata.vars, collapse = '')
} else {
  by.strata.text <- NULL
}

if(str_detect(trnd, "time.strata" )) {
  time.strata.text <- time.strata.length
} else {
  time.strata.text <- NULL
}

if(oy==0) {
  oytext <- " Fyear 0"
} else {
  oytext <- NULL
}

if(!is.null(eknots)) {
  eknots.text <- round(eknots,2)
} else {
  eknots.text <- NULL
}

if(lspline=='integer') {
  lagknot.text <- NULL # knots irrelevant to unconstrained, do not report them
} else {
  lagknot.text <- paste0(' knts',paste0(round(lknots,2), collapse = ','))
}

# Cannot include ; or exceed 255 characters. : is converted to / and best avoided. Combined folder in OneDrive cannot exceed 520 characters
trnd.no.text <- c('\\*','\\/','.no.years','Date,','date,') # text to remove from trnd in folder. ,'date,' required for gam but forgot to use it
model.text <- c(paste0(distribution[1],' ',modeltype, oytext), # space inserted in oytext manually, accounting for whether its NULL or not
                # paste("By vars-", paste(by.vars, collapse = ',')), # too long
                by.filter.text,
                by.strata.text,
                # paste("Formula -",mformula), # too long
                paste0("Exp ", espline, ' knts', paste0(eknots.text, collapse = ','), hum),
                paste0("Lag ", lspline, ' max', lmax, lagknot.text),
                paste0("Trend ", str_remove_all(trnd, paste(trnd.no.text, collapse='|')), time.strata.text)) # shortened to just have spline and df (minus * years) mention, * also not compatible with some paths
print(model.text)

### Overall

## Outcome
outcome.exposure.l <- paste0(results.loc,outcome.var.short) # folder name, exclude 000s in name
suppressWarnings(dir.create(outcome.exposure.l)) # Create folder if it doesn't exist already. Warning if exists (doesn't replace), suppress this

## Exposure
outcome.exposure.l0 <- paste0(outcome.exposure.l,'/',by.vars.oi,exposure.var.short) # folder name, exclude 000s in name
suppressWarnings(dir.create(outcome.exposure.l0)) # Create folder if it doesn't exist already. Warning if exists (doesn't replace), suppress this

## Model details
outcome.exposure.loc0 <- paste(outcome.exposure.l0, paste(model.text, collapse = '; '), sep='/') # folder name
suppressWarnings(dir.create(outcome.exposure.loc0)) # Create folder if it doesn't exist already. Warning if exists (doesn't replace), supress this
outcome.exposure.loc <- paste0(outcome.exposure.loc0,'/') # folder destination (add /)

## .txt file with model parameters
fileConn <- file(paste0(outcome.exposure.loc,'Model specifications.txt'))
writeLines(c(model.text,mformula,add.formula.io), fileConn)
close(fileConn)


### Descriptive statistics
suppressWarnings(dir.create(paste0(outcome.exposure.loc,'Descriptive'))) # Create folder if it doesn't exist already. Warning if exists (doesn't replace)
outcome.exposure.loc.des <- paste0(outcome.exposure.loc,'Descriptive/') # folder destination (add /)

# Density plots 
suppressWarnings(dir.create(paste0(outcome.exposure.loc,'Density'))) # Create folder if it doesn't exist already. Warning if exists (doesn't replace)
outcome.exposure.loc.density <- paste0(outcome.exposure.loc,'Density/') # folder destination (add /)

suppressWarnings(dir.create(paste0(outcome.exposure.loc.density, 'Outcome'))) # Create folder if it doesn't exist already. Warning if exists (doesn't replace)
outcome.exposure.loc.density.outcome <- paste0(outcome.exposure.loc.density,'Outcome/') # folder destination (add /)

suppressWarnings(dir.create(paste0(outcome.exposure.loc.density, 'Exposure'))) # Create folder if it doesn't exist already. Warning if exists (doesn't replace)
outcome.exposure.loc.density.exposure <- paste0(outcome.exposure.loc.density,'Exposure/') # folder destination (add /)

# Histograms 
suppressWarnings(dir.create(paste0(outcome.exposure.loc,'Histograms'))) # Create folder if it doesn't exist already. Warning if exists (doesn't replace)
outcome.exposure.loc.histogram <- paste0(outcome.exposure.loc,'Histograms/') # folder destination (add /)

suppressWarnings(dir.create(paste0(outcome.exposure.loc.histogram, 'Outcome'))) # Create folder if it doesn't exist already. Warning if exists (doesn't replace)
outcome.exposure.loc.histogram.outcome <- paste0(outcome.exposure.loc.histogram,'Outcome/') # folder destination (add /)

suppressWarnings(dir.create(paste0(outcome.exposure.loc.histogram, 'Exposure'))) # Create folder if it doesn't exist already. Warning if exists (doesn't replace)
outcome.exposure.loc.histogram.exposure <- paste0(outcome.exposure.loc.histogram,'Exposure/') # folder destination (add /)

# Scatter 
suppressWarnings(dir.create(paste0(outcome.exposure.loc,'Scatter'))) # Create folder if it doesn't exist already. Warning if exists (doesn't replace)
outcome.exposure.loc.scatter<- paste0(outcome.exposure.loc,'Scatter/') # folder destination (add /)

suppressWarnings(dir.create(paste0(outcome.exposure.loc.scatter, 'Outcome'))) # Create folder if it doesn't exist already. Warning if exists (doesn't replace)
outcome.exposure.loc.scatter.outcome <- paste0(outcome.exposure.loc.scatter,'Outcome/') # folder destination (add /)

suppressWarnings(dir.create(paste0(outcome.exposure.loc.scatter, 'Exposure'))) # Create folder if it doesn't exist already. Warning if exists (doesn't replace)
outcome.exposure.loc.scatter.exposure <- paste0(outcome.exposure.loc.scatter,'Exposure/') # folder destination (add /)

# Trend, seasonality and noise 
suppressWarnings(dir.create(paste0(outcome.exposure.loc,'Trend seasonality noise'))) # Create folder if it doesn't exist already. Warning if exists (doesn't replace)
outcome.exposure.loc.ts <- paste0(outcome.exposure.loc,'Trend seasonality noise/') # folder destination (add /)


### Stage 1
suppressWarnings(dir.create(paste0(outcome.exposure.loc,'S1'))) # Create folder if it doesn't exist already. Warning if exists (doesn't replace)
outcome.exposure.loc.s1 <- paste0(outcome.exposure.loc,'S1/') # folder destination (add /)

suppressWarnings(dir.create(paste0(outcome.exposure.loc.s1,'r'))) # Create folder if it doesn't exist already. Warning if exists (doesn't replace)
outcome.exposure.loc.s1.r <- paste0(outcome.exposure.loc.s1,'r/') # folder destination (add /)

suppressWarnings(dir.create(paste0(outcome.exposure.loc.s1,'Model checks'))) # Create folder if it doesn't exist already. Warning if exists (doesn't replace)
outcome.exposure.loc.s1.mcheck <- paste0(outcome.exposure.loc.s1,'Model checks/') # folder destination (add /)

suppressWarnings(dir.create(paste0(outcome.exposure.loc.s1,'O e-r'))) # for overall e-r plots similar to S2. Not to be in publications, but for checking difference made by BLUPS
outcome.exposure.loc.s1.oer <- paste0(outcome.exposure.loc.s1,'O e-r/') # folder destination (add /)

if(modeltype %in% c('a','gam','GAM')) {
  suppressWarnings(dir.create(paste0(outcome.exposure.loc.s1,'Smoothers'))) # Create folder if it doesn't exist already. Warning if exists (doesn't replace)
  outcome.exposure.loc.s1.smoother <- paste0(outcome.exposure.loc.s1,'Smoothers/') # folder destination (add /)
} 



##################################################
### Data manipulation
##################################################

# Define outcomes
outcomes <- c("no.claims", "total","cmp","gas", "hs","ncmp","total0","cmp0","gas0", "hs0","ncmp0","total_5","cmp_5","gas_5","hs_5","ncmp_5") # outcomes for analysis

# Derived by variable vectors
by.vars1 <- replace(by.vars, by.vars=="state", 'city')
by.vars2 <- replace(by.vars, by.vars=="Date", 'quarter')
by.vars2 <- by.vars2[by.vars2 %in% c('Date','date','City','city','GCCSA','gccsa','SA4','sa4','Quarter','quarter','outorin','outorini','outin','outini','Industry','industry','Gender','gender','Sex','sex','n','hours')]
by.vars3 <- by.vars[!by.vars == 'Date']
by.vars4 <- c('Year', replace(by.vars, by.vars=="Date", 'Month'))

# SWA4: Filter by outorin and conditions (if applicable) and stratify by by-variables
if(isTRUE(io.marker)) {
  swa4.io <- swa4[str_detect(swa4[,get(by.vars.oi)],'door')] # only keep if indoors or outdoors, and not a number
} else {
  swa4.io <- swa4
}

if(is.null(filter.conds) | is.na(filter.conds) | filter.conds=='') { # mget enables character variables to work
  daily <- swa4.io[, lapply(.SD, sum, na.rm=T), keyby=mget(by.vars), .SDcols = outcomes]
  # daily <- swa4[inj.or.dis=='Injuries', lapply(.SD, sum, na.rm=T), keyby=mget(by.vars), .SDcols = outcomes] 
} else {
  daily <- swa4.io[eval(parse(text=filter.conds)), lapply(.SD, sum, na.rm=T), keyby=mget(by.vars), .SDcols = outcomes] # mget enables character variables to work
}
# rm(swa4.io) # save space

# Complete time series by inserting rows for each by.vars. Values filled with NA. No effect for daily as already complete. !!!syms(by.vars) enables character variables to work, all required to include all rows
daily <- merge(daily, distinct(expand_grid(swa4[,mget(by.vars)])), by=by.vars, all=T) # make sure all records are included, even if remove by filter
daily <- as.data.table(complete(daily, !!!syms(by.vars), fill=list(total=0,cmp=0,gas=0,hs=0,ncmp=0,no.claims=0, total0=0,cmp0=0,gas0=0,hs0=0,ncmp0=0, total_5=0,cmp_5=0,gas_5=0,hs_5=0,ncmp_5=0))) # complete time period and replace NA with 0
# daily <- merge(daily, complete(daily, !!!syms(by.vars), fill=list(total=0,cmp=0,gas=0,hs=0,ncmp=0,no.claims=0, total0=0,cmp0=0,gas0=0,hs0=0,ncmp0=0, total_5=0,cmp_5=0,gas_5=0,hs_5=0,ncmp_5=0)), by=c(by.vars,outcomes), all=T) # if merge with expand_grid, results in duplicates. no longer needed
daily <- daily[!is.na(no.claims)] # FIX. Remove NA values, which result above alongside 0. After complete, only useulf after merge where duplicate rwos, one with NA and one with 0, results
daily <- na.omit(daily, cols=by.vars) # remove rows where by.vars is missing (industry, outorin)
daily <- daily[City!='Hobart' | City=='Hobart' & Date>='2006-07-01',]  # Remove Hobart financial 05-06 dates

# Outcomes by by.vars 
daily.outcomes <- daily[, lapply(.SD, sum, na.rm=T), keyby=mget(by.vars3), .SDcols = outcomes] # mget enables character variables to work

# Adjust outcomes to be in 1000s. Greatly improves model fitting time and increases likelihood of converge
if(oy==0) {
  daily[, ':='('Total costs (000s)'=total0/1000, 'Compensation costs (000s)'=cmp0/1000, 'Goods and services costs (000s)'=gas0/1000, 'Health services costs (000s)'=hs0/1000, 'Non-compensation costs (000s)'=ncmp0/1000,
               'Total costs'=total0, 'Compensation costs'=cmp0, 'Goods and services costs'=gas0, 'Health services costs'=hs0, 'Non-compensation costs'=ncmp0,)]
} else {
  daily[, ':='('Total costs'=total, 'Compensation costs'=cmp, 'Goods and services costs'=gas, 'Health services costs'=hs, 'Non-compensation costs'=ncmp,
               'Total costs (000s)'=total/1000, 'Compensation costs (000s)'=cmp/1000, 'Goods and services costs (000s)'=gas/1000, 'Health services costs (000s)'=hs/1000, 'Non-compensation costs (000s)'=ncmp/1000,
               'Costs in same financial year (000s)'=total0/1000, 'Compensation costs in same financial year (000s)'=cmp0/1000, 'Goods and services costs in same financial year (000s)'=gas0/1000, 'Health services costs in same financial year (000s)'=hs0/1000, 'Non-compensation costs in same financial year (000s)'=ncmp0/1000,
               'Costs up to Jun 14 (000s)'=total_5/1000, 'Compensation costs up to Jun 14 (000s)'=cmp_5/1000, 'Goods and services costs up to Jun 14 (000s)'=gas_5/1000, 'Health services costs up to Jun 14 (000s)'=hs_5/1000, 'Non-compensation costs up to Jun 14 (000s)'=ncmp_5/1000)] # Costs if FYear <= 13 (all claims have at least 5 financial years to record financial payments, up to June 19)
}

# If non-zero costs, add 1
if(type.outcome=='non-0') {
  daily[, ':='('Total costs (000s)'=`Total costs (000s)`+0.001, 'Compensation costs (000s)'=`Compensation costs (000s)`+0.001,
               'Goods and services costs (000s)'=`Goods and services costs (000s)`+0.001,
               'Health services costs (000s)'=`Health services costs (000s)`+0.001,
               'Non-compensation costs (000s)'=`Non-compensation costs (000s)`+0.001,
               'Total costs'=`Total costs`+1, 'Compensation costs'=`Compensation costs`+1,
               'Goods and services costs'=`Goods and services costs`+1,
               'Health services costs'=`Health services costs`+1,
               'Non-compensation costs'=`Non-compensation costs`+1)]
  print('Added $1')
}

# Calculate cost per illness
daily[, ':='('Cost per OII (000s)'=`Total costs (000s)`/no.claims, 'Cost per OII'=`Total costs`/no.claims)]

# If analysing costs with (quasi-)Poisson, round costs,
if(str_detect(distribution[1], 'oisson') & type.outcome == 'cost') {
  daily[,c(outcome.var) := round(get(outcome.var))] # c and get to refer to columns, respectively
}

# Obtain by variable combinations
by.var.ds <- unique(daily[,mget(by.vars3)])



#################################################
### Merge data, and adjusting for indoor and outdoor + Creating supplementary climate data for DLNM 
##################################################

# ### Load climate station dataset
# load(paste0(daily.station.data.loc, 'daily.bom3.rda')) # load climate data
# climate <- daily.bom3

### Load and merge BARRA climate dataset
load(paste0(barra.loc, 'barra_all.rda')) # load climate data

# Barra data from start date - lmax (to adjust for days lost with dlnm) to end date. Accounts for Hobart starting in 06
climate <- barra_all1[City!='Hobart' & Date %in% seq(as.Date('2005-07-01')-lmax, as.Date('2018-06-30'), by = "day") |
                        City=='Hobart' & Date %in% seq(as.Date('2006-07-01')-lmax, as.Date('2018-06-30'), by = "day"),]
climate <- climate[City != 'Canberra'] # not for stats analysis

### Heat metrics adjusted for indoor/outdoor
if(isTRUE(io.marker)) { # if indoor/outdoor used
  if (!(by.vars.oi %in% colnames(climate))) { # if have not already duplicated supplementary barra data for outorin
    climate[, get('by.vars.oi'):='Indoors'] # Create outorin variable and set to indoors
    climate <- bind_rows(climate, climate)  # duplicate dataset
    climate[duplicated(climate), get('by.vars.oi'):='Outdoors'] # for duplicate, set outorin variable to outdoors
  }
}
climate <- unique(climate) # additional safety measure to ensure climate_supp isn't duplicated more than once
climate[City!='Hobart' & Date %in% seq(as.Date('2005-07-01')-lmax, as.Date('2005-06-30'), by = "day") |
          City=='Hobart' & Date %in% seq(as.Date('2006-07-01')-lmax, as.Date('2006-06-30'), by = "day"), dlnm.dummy:=1] # remove 05-06 fyear from Hobart


### Merge climate dataset and rename exposure
daily1 <- merge(daily, climate, by=by.vars, all=T) # outcome data will be NA as expected

exposure.old <- c('no.claims',
                  'max_airt','ave_airt','min_airt',
                  'max_wbgtb','ave_wbgtb','min_wbgtb',
                  'max_wbgtl','ave_wbgtl','min_wbgtl',
                  'max_steadman.id','ave_steadman.id','min_steadman.id',
                  'max_steadman.od','ave_steadman.od','min_steadman.od',
                  'max_heat.index','ave_heat.index','min_heat.index',
                  'max_humidex','ave_humidex','min_humidex',
                  'ehf','ehff','ehf.hi','ehff.hi'
                  # ,'max_sh','ave_sh','min_sh',
                  # 'max_rh','ave_rh','min_rh'
) # new names
exposure.new <- c("Number of illnesses",
                  "Maximum temperature","Average temperature","Minimum temperature",
                  "Maximum indoor WBGT","Average indoor WBGT","Minimum indoor WBGT",
                  "Maximum outdoor WBGT","Average outdoor WBGT","Minimum outdoor WBGT",
                  "Maximum indoor Steadman's apparent temperature","Average indoor Steadman's apparent temperature","Minimum indoor Steadman's apparent temperature",
                  "Maximum outdoor Steadman's apparent temperature","Average outdoor Steadman's apparent temperature","Minimum outdoor Steadman's apparent temperature",
                  'Maximum heat index','Average heat index','Minimum heat index',
                  'Maximum humidex','Average humidex','Minimum humidex',
                  'Excess heat factor','Excess heat factor (forward)','Excess heat index factor','Excess heat index factor (forward)'
                  # ,"Maximum specific humidity","Average specific humidity","Minimum specific humidity",
                  # "Maximum relative humidity","Average relative humidity","Minimum relative humidity"
) # new names
setnames(daily1, old=exposure.old, new=exposure.new)  # rename exposures. non-T metrics don't work well in code if >1 word
daily1[,':='("Maximum specific humidity"=max_sh, "Maximum relative humidity"=max_rh, 
             "Average specific humidity"=ave_sh, "Average relative humidity"=ave_rh)]

### Load and merge monthly workforce data
pop <- copy(labour.m)
pop[, ':='(City=state.to.city.fn(GCCSA), GCCSA=NULL, Date=NULL)] # convert GCCSA to city

if(isTRUE(io.marker)) {
  if(by.vars.oi=='outin') { # as outorin for labour.m is based on proportions of outorin, it is fine to use it with outin
    names(pop) = gsub('outorin', 'outin', names(pop)) # rename outorin in pop to outin without changing labour.m
  } 
  setnames(labour.m, 'outin', 'outorin', skip_absent=T) # must be done after changing pop, as labour.m otherwise affected
} else {
  pop <- pop[, lapply(.SD, sum, na.rm=T), by=mget(by.vars4), .SDcols = c("n","fulltime","parttime")] 
}

if (any(by.vars %in% c('outorini','outini'))) {
  # Load and merge quarterly workforce data
  daily1 <- daily1[,quarter := quarter.fn(Date)] # add quarter
  pop <- labgccsa.qi[, lapply(.SD, sum, na.rm=T), keyby=c('quarter','City','outorin'), .SDcols = c("n","hours")] # aggregate labgccsa if needed # replace(x, x==0, 1)
  pop[, quarter:=as.Date(quarter)]
  setnames(pop, 'outorin', by.vars.oi)
  daily.ds <- merge(daily1, pop, by=by.vars2, all.x=T)
} else {
  daily.ds <- merge(daily1, pop, by=by.vars4, all.x=T)
}

## Time strata
daily.ds[, ref.date :=fifelse(City=='Hobart','2006-07-01','2005-07-01')] # first Monday before study period adjusted for Tas (starts on Fri, Sat for Tas)
daily.ds[, time.strata := as.factor(floor(difftime(Date, ref.date, units='days')/time.strata.length)+0)] # time strata length from reference date, used by Blesson for time-stratified case-crossover
daily.ds[, my := as.factor(paste0(month,fyear))] # time strata length from reference date, used by Blesson for time-stratified case-crossover

## Fourier term
# daily.ds[, harominc(Date, nfreq=4, )]



##################################################
###  Patching certain data and years, only included if needed
##################################################

if(isTRUE(ehf.marker))  { # if excess heat factor or heatwave
  daily.ds[Month %in% 4:9, (outcome.var):=NA] # limit results to warm season by removing outcome data during cold season. Keeps temperature data intact for lag
  # daily.ds[Month==12 & Day %in% c(23:30), shol:='No'] # exclude as shol for EHF only. is better with it as shol
  day.before.melbourne.cup <- as.Date(c("2004-11-01","2005-10-31","2006-11-06","2007-11-05","2008-11-03","2009-11-02","2010-11-01","2011-10-31","2012-11-05","2013-11-04","2014-11-03","2015-11-02","2016-10-31","2017-11-06","2018-11-05")) # do manually, as shift seems to miss 31st October
  daily.ds[City=='Melbourne' & Date %in%  day.before.melbourne.cup, shol:='Day before Melbourne Cup'] # limit results to warm season by removing outcome data during cold season. Keeps temperature data intact for lag
  # daily.ds[Month==12 & Day==23, shol:='23rd Dec'] # 23rd Dec seems to results in high res for Melbourne and Sydney contrasting with generally negative res in Xmas break
  daily.ds[Month==12 & Day %in% c(24,25), shol:='Xmas Day or Eve'] # 23rd Dec seems to results in high res for Melbourne and Sydney contrasting with generally negative res in Xmas break
}

if(type.outcome %in% c('cost','non-0')) { # specifically Perth FYear 2011 for health service costs - it has a huge peak
  daily.ds[, fyear_p := fifelse(City=='Perth' & FYear==2011, 1, 0)]
  fyearp <- T
} else {
  daily.ds[, fyear_p := 0] # useless variable
  fyearp <- F
}

if(str_detect(outcome.var, '14'))  { # if detect 14 for up to June 2014 (outcome excludes last 5 financial years)
  daily.ds <- daily.ds[FYear <= 2013] # exclude last 5 financial years (14-18)
}



##################################################
### Adjust heat metrics for outdoors and indoors
##################################################

if (isTRUE(io.marker)) { # if the by-variables include a variable for indoor, outdoor or both
  daily.ds[, "Maximum Steadman's apparent temperature" := fifelse(get(by.vars.oi) %in% c('Outdoor','Outdoors'),get("Maximum outdoor Steadman's apparent temperature"),get("Maximum indoor Steadman's apparent temperature"))]
  daily.ds[, "Average Steadman's apparent temperature" := fifelse(get(by.vars.oi) %in% c('Outdoor','Outdoors'),get("Average outdoor Steadman's apparent temperature"),get("Average indoor Steadman's apparent temperature"))]
  daily.ds[, "Maximum WBGT" := fifelse(get(by.vars.oi) %in% c('Outdoor','Outdoors'),get("Maximum outdoor WBGT"),get("Maximum indoor WBGT"))]
  daily.ds[, "Average WBGT" := fifelse(get(by.vars.oi) %in% c('Outdoor','Outdoors'),get("Average outdoor WBGT"),get("Average indoor WBGT"))]
}



##################################################
### Set up . Create objects to store results
##################################################

daily.ds[,stratum := do.call(paste, c(mget(by.vars3), sep=' '))] # combine by.vars3 into a string as stratum
# daily.ds <- daily.ds[Date < as.Date('2018-06-25')] # last week has less OIs (particularly last 2 days), and has poor residual testing
daily.ds0 <- daily.ds[dlnm.dummy!=1 | is.na(dlnm.dummy)]
ds.stratum <- sort(unique(daily.ds$stratum))
length.ds.stratum <- length(ds.stratum)

# Matrices with 1 column: Number of outcomes, formula, exposure at 0, AIC and dispersion (if poisson used)
no.outcome <- no.years <- form <- m.aic <- m.dispersion <- m.r2 <- m.devexp <- matrix(NA, length.ds.stratum, 1, dimnames=list(ds.stratum)) # matrix of NAs, with a column to denote total number of outcomes
colnames(m.aic) <- 'AIC'
colnames(m.dispersion) <- 'Dispersion parameter'
colnames(m.r2) <- 'R^2'
colnames(m.devexp) <- 'Deviance explained'

# Exposure. Matrix of NAs, with a column to denote total number of outcomes
exposure <- matrix(NA, length.ds.stratum, 4 + length(predper.short), dimnames=list(ds.stratum, c('mean','range','min','max',paste0(predper.short,'%'))))

# EHF. Heatwave if >0, though as a continuous variable, they are essentially identical
exposure_ehf <- matrix(NA, length.ds.stratum, 4, dimnames=list(ds.stratum, c('Heatwave','Severe','Extreme','Very extreme'))) # matrix of NAs, with a column to denote total number of outcomes
ehf_threshold <- matrix(NA, length.ds.stratum, 4, dimnames=list(ds.stratum, c('Heatwave','Severe','Extreme','Very extreme'))) # matrix of NAs, with a column to denote total number of outcomes

# Percentile exposure per by variables
exposure.by <- matrix(NA, length.ds.stratum, length(predper), dimnames=list(ds.stratum, predper))

# Coefficients for overall cumulative summary
if(espline=='bs') {
  espline.length <- length(eknots) + 3 # not sure if correct, based on testing. If wrong, error: Error in coef[i, ] <- coef(red[[i]]) : number of items to replace is not a multiple of replacement length
  coef <- matrix(NA, length.ds.stratum, espline.length, dimnames=list(ds.stratum, paste0('b',rep(1:espline.length)))) # matrix of NAs, with rows per by.vars
} else {
  coef <- matrix(NA, length.ds.stratum, edf, dimnames=list(ds.stratum, paste0('b',rep(1:edf)))) # matrix of NAs, with rows per by.vars
}

# Covariance for overall cumulative summary
vcov <- vector("list", length.ds.stratum)
names(vcov) <- ds.stratum

# Tweedie max shape parameter values (if Tweedie distribution used)
m.tweedie.shape <- mean_rhs <- matrix(NA, length.ds.stratum, 1, dimnames=list(ds.stratum)) # matrix of NAs, with a column to denote optimised shape parameter
# m.tweedie.shape <- matrix(NA, length.ds.stratum,3, dimnames=list(ds.stratum)) # matrix of NAs, with columns to dneote optimised shape parameters and their 95% CI

# Matrix for k.index
no.smooth <- str_count(trnd, ',bs=') # number of smooths
gam.k <- matrix(NA, length.ds.stratum*no.smooth, 4, dimnames=list(rep(ds.stratum, each=no.smooth), c("k'","edf","k-index","p-value"))) # matrix of NAs, with rows per by.vars

# List objects
temps <- rhs <- outcomes <- model <- model.omit <- model2 <- m.coef <- res <- model.checks <- collin <- gam_k <- gam.concurvity <- check <- model.checks.month <- model.checks.week <- check.res.coef <- check.res.coef.sig <- exposure.rr.s1 <- red <- list() #  Model, residuals, residual length + plots and earity



################################################################################
# Descriptive temperature matrix
################################################################################

# wbgt.desc.v <- c('Maximum temperature','max_rh','max_ws','max_dswr','Maximum indoor WBGT','Maximum outdoor WBGT') # 'Average temperature','ave_rh','Average indoor WBGT','Average outdoor WBGT'
# ds.stratum.city <- ds.stratum[!str_detect(ds.stratum, 'utdoor')]
# wbgt.desc <- matrix(NA, length(ds.stratum.city), length(wbgt.desc.v), dimnames=list(ds.stratum.city)) # matrix of NAs, with a column to denote total number of outcomes
# colnames(wbgt.desc) <- wbgt.desc.v
# n.desc <- matrix(NA, length.ds.stratum, 1, dimnames=list(ds.stratum)) # matrix of NAs, with a column to denote total number of outcomes
# wbgt.desc.r <- 1
# n.mean <- n.years <- n.sd <- list()
# 
# # Fill in descriptive table for temperture metrics
# for(i in ds.stratum[!str_detect(ds.stratum, 'utdoor')]) {
#   .ds <- daily.ds[stratum==i] # dataset for each unique value #
#   .d <- .ds[-(1:lmax)] # dataset for each unique value # daily.ds[stratum=='Adelaide, Outdoors']
# 
#   # temp.desc[i,] <- as.matrix(round(.d[, lapply(.SD, mean, na.rm=T), .SDcols = temp.desc.v],2))
# 
#   # .wbgt.mean <- as.matrix(round.fn(.d[, lapply(.SD, mean, na.rm=T), .SDcols = wbgt.desc.v],wbgt.desc.r)) # overall mean
#   # .wbgt.years <- round.fn(.d[, lapply(.SD, mean, na.rm=T), by=FYear, .SDcols = wbgt.desc.v],wbgt.desc.r)[,-1] # values for each fyear
#   # wbgt.desc[i,] <- paste0(.wbgt.mean, ' (', apply(.wbgt.years, 2, min),' to ',apply(.wbgt.years, 2, max), ')') # range across fyears
# 
#   .wbgt.mean <- .d[, lapply(.SD, mean, na.rm=T), .SDcols = wbgt.desc.v]# overall mean
#   .wbgt.sd <- .d[, lapply(.SD, sd, na.rm=T), .SDcols = wbgt.desc.v]
#   wbgt.desc[i,] <- paste0(round.fn(.wbgt.mean,wbgt.desc.r), ' (', round.fn(.wbgt.mean-.wbgt.sd,wbgt.desc.r),'-', round.fn(.wbgt.mean+.wbgt.sd,wbgt.desc.r), ')') # range across fyears
# }
# 
# rownames(wbgt.desc) <- word(rownames(wbgt.desc), 1) # Only state city
# colnames(wbgt.desc) <- c('Maximum temperature (Â°C)','Relative humidity (%)','Wind speed (m/s)','Solar radiation (W/m2)','Indoor WBGT (Â°C)','Outdoor WBGT (Â°C)')
# # Units are degrees for temp, RH is %, wind speed is m/s at 10m, solar radiation is (W/m2)
# 
# # # Fill in descriptive table for n
# for(i in ds.stratum) {
#   .ds <- daily.ds[stratum==i] # dataset for each unique value #
#   .d <- .ds[-(1:lmax)] # dataset for each unique value # daily.ds[stratum=='Adelaide, Outdoors']
#   n.mean[[i]] <- mean(.d[day1=='Yes',n]) # mean per month. This works as pop changes every month
#   n.sd[[i]] <- sd(.d[day1=='Yes',n]) # mean per month. This works as pop changes every month
#   # n.years[[i]] <- .d[day1=='Yes', lapply(.SD, mean), by=FYear, .SDcols = 'n'][,-1] # values for each fyear
#   # n.desc[i,] <- paste0(round.fn(n.mean[[i]]), ' (', round.fn(min(.n.years[[i]])),' to ', round.fn(max(.n.years[[i]])), ')') # range across fyears
#   n.desc[i,] <- paste0(round.fn(n.mean[[i]]/1000), ' (', round.fn((n.mean[[i]]-n.sd[[i]])/1000),'-', round.fn((n.mean[[i]]+n.sd[[i]])/1000), ')') # range across fyears
# }
# 
# n.desc.in <- n.desc[str_detect(rownames(n.desc), 'ndoor')]
# n.desc.out <- n.desc[str_detect(rownames(n.desc), 'utdoor')]
# 
# t.desc <- t(cbind(wbgt.desc, 'Indoor workers (000s)'=n.desc.in, 'Outdoor workers (000s)'=n.desc.out))
# # smallest and largest n for all strata are the first and last years, except for Darwin. Darwin indoors lowest is 2nd fyear, Darwin outdoors highest is 2nd to last year
# 
# # Output file
# write_xlsx(as.data.table(cbind('City'=rownames(wbgt.desc), wbgt.desc)), path=paste0(swa.heat.desc.results.loc,'WBGT.xlsx'), col_names=T) # create csv file
# write_xlsx(as.data.table(cbind('Parameter'=rownames(t.desc), t.desc)), path=paste0(swa.heat.desc.results.loc,'Desc.xlsx'), col_names=T) # create csv file
# save(wbgt.desc, t.desc, n.mean, file=paste0(swa.heat.desc.results.loc, 'Desc.rda')) # save dataset after all changes for easier access



##################################################
### DESCRIPTIVE STATISTICS, do not include dummy dlnm rows
##################################################

# Density plots
png(file = paste0(outcome.exposure.loc.density.outcome, 'National', '.png')) # plot location. File name based on heat metric
plot(density(daily.ds0[!is.na(get(outcome.var)),get(outcome.var)]), main=outcome.var, xlab=outcome.var)
dev.off() # Save image + clear settings

png(file = paste0(outcome.exposure.loc.density.exposure, 'National', '.png')) # plot location. File name based on heat metric
plot(density(daily.ds0[!is.na(get(exposure.var)),get(exposure.var)]), main=exposure.var, xlab=exposure.var)
dev.off() # Save image + clear settings

# Histograms. Better represents zero mass of 0s
png(file = paste0(outcome.exposure.loc.histogram.outcome, 'National', '.png')) # plot location. File name based on heat metric
hist(daily.ds0[,get(outcome.var)], main=outcome.var, xlab=outcome.var, breaks=50)
dev.off() # Save image + clear settings

png(file = paste0(outcome.exposure.loc.histogram.exposure, 'National', '.png')) # plot location. File name based on heat metric
hist(daily.ds0[,get(exposure.var)], main=exposure.var, xlab=exposure.var, breaks=50)
dev.off() # Save image + clear settings

# Scatter plots
png(file = paste0(outcome.exposure.loc.scatter.outcome, 'National', '.png')) # plot location. File name based on heat metric
plot(x=daily.ds0[,get(exposure.var)], y=daily.ds0[,get(outcome.var)], pch=16, cex=0.5, col=daily.ds0[,dow], main='Scatter plot', ylab=outcome.var, xlab=exposure.var)
legend("topright", legend=levels(daily.ds0[,dow]), col=unique(daily.ds0[,dow]))
# ggplot(daily.ds0, aes(x=exposure.var, y=outcome.var, color=dow))
dev.off() # Save image + clear settings

png(file = paste0(outcome.exposure.loc.scatter.exposure, 'National', '.png')) # plot location. File name based on heat metric
plot(x=daily.ds0[,Date], y=daily.ds0[,get(exposure.var)], pch=16, cex=0.1, col=daily.ds0[,dow], main='Scatter plot', ylab=exposure.var, xlab='Date')
dev.off() # Save image + clear settings

# Create descriptive table for each stratum and overall
.descrip.overall <- with(daily.ds0, psych::describe(get(outcome.var)))
.descrip.stratum <- with(daily.ds0, as.data.table(matrix(unlist(psych::describeBy(get(outcome.var),stratum)), byrow=T, ncol=13)))
colnames(.descrip.stratum) <- c('vars','n','mean','sd','median','trimmed','mad','min','max','range','skew','kurtosis','se')   
.descrip <- rbind(.descrip.stratum, .descrip.overall)
.descrip[, c('Stratum','variance','mean/variance') := .(c(sort(unique(daily.ds$stratum)),'Overall'), sd^2, mean/(sd^2))] # add variance and dispersion factor
setcolorder(.descrip, 'Stratum') # reorganise stratum to be first
write.csv(.descrip, file=paste0(outcome.exposure.loc.des,'Outcome desc.csv'), na='') # create csv file



##################################################
### Time series components, do not include dummy dlnm rows
##################################################

# Shared graphical parameters. Add to end of ggplots
theme.bs <- theme_bw(base_size=6) # theme_bw()
theme.ts <- theme(plot.title=element_text(size=6), legend.title=element_text(size=5), legend.text=element_text(size=4))
# scale.date <- scale_x_continuous(breaks = c(2005,2007,2009,2011,2013,2015,2017))
dow.colours <- scale_color_manual(values=c('red','orange','yellow2','green3','green4','cyan','blue'))
# month.colours <- scale_color_manual(values=c('red','orange','yellow2','green3','green4','cyan','blue','dark blue','purple','pink','brown','grey'))
pubhol.shapes <- scale_shape_manual(values=c(1,2,4)) # circle, triangle and cross
# theme.bs + theme.ts + dow.colours + pubhol.shapes # combining as object fails

gc()
tic()
for(i in ds.stratum) {
  # .d <- daily.ds0[stratum==i] # dataset for each unique value # daily.d[stratum=='Adelaide, Outdoors']
  .ds <- daily.ds[stratum==i] # dataset for each unique value #
  .d <- .ds[-(1:lmax)] # dataset for each unique value # daily.ds[stratum=='Adelaide, Outdoors']
  .name <- unique(.d$stratum) # names with all of a, b and c together
  .outcome <- .d[,get(outcome.var)] # outcome
  .exposure <- .d[,get(exposure.var)] # exposure
  
  # Density plot: outcome
  png(file = paste0(outcome.exposure.loc.density.outcome, .name, ' d.png')) # plot location. File name based on heat metric
  plot(density(.outcome[!is.na(.outcome)]), main=outcome.var, xlab=outcome.var) # no missing data should be present
  dev.off()
  
  # Density plot: exposure
  png(file = paste0(outcome.exposure.loc.density.exposure, .name, ' d.png'))
  plot(density(.exposure[!is.na(.exposure)]), main=exposure.var, xlab=exposure.var) # density fails with missing data
  dev.off()
  
  # Histograms: outcome
  png(file = paste0(outcome.exposure.loc.histogram.outcome, .name, ' h.png'))
  hist(.outcome, main=outcome.var, xlab=outcome.var, breaks=50)
  dev.off() # Save image + clear settings
  
  # Histograms: exposure
  png(file = paste0(outcome.exposure.loc.histogram.exposure, .name, ' h.png'))
  hist(.exposure, main=outcome.var, xlab=outcome.var, breaks=50)
  dev.off() # Save image + clear settings
  
  # Scatter plot: outcome vs exposure
  png(file = paste0(outcome.exposure.loc.scatter.outcome, .name, ' soe.png'))
  plot(x=.exposure, y=.outcome, pch=16, cex=0.5, col=.d[,dow], main='Scatter plot', ylab=outcome.var, xlab=exposure.var)
  dev.off() # Save image + clear settings
  
  # Scatter plot: outcome vs time
  png(file = paste0(outcome.exposure.loc.scatter.outcome, .name, ' s.png'))
  plot(x=.d[,Date], y=.outcome, pch=16, cex=0.5, col=.d[,dow], main='Scatter plot', ylab=outcome.var, xlab='Date')
  dev.off() # Save image + clear settings
  
  # Scatter plot: exposure vs time
  png(file = paste0(outcome.exposure.loc.scatter.exposure, .name, ' s.png'))
  plot(x=.d[,Date], y=.exposure, pch=16, cex=0.1, col=.d[,dow], main='Scatter plot', ylab=exposure.var, xlab='Date')
  dev.off() # Save image + clear settings
  
  # Trend. Assume seasonality period is a week
  tns.file.name <- paste0(outcome.exposure.loc.ts, .name, ',tsn.png')
  if(isTRUE(tns) & nchar(tns.file.name) <= 255) { # if elect to save tns and file isn't too long
    .d[,.moving.ave := rollmean(get(outcome.var), 7, fill=NA)] # over a week. align=centre (default) # frollmean from data.table()
    .d[,.centered.moving.ave := rowMeans(cbind(.moving.ave, shift(.moving.ave)))]
    g.t <- ggplot(.d, aes(x=Date, y=.centered.moving.ave, col=dow)) + geom_point(aes(shape=public.hol.x)) + labs(title='Trend + seasonality: weekly rolling averages', x="Date", y=outcome.var) + theme.bs + theme.ts + dow.colours + pubhol.shapes
    # ggplot(.d, aes(x=Date, y=.centered.moving.ave, col=.d[,get(exposure.var)])) + geom_point(aes(shape=dow, size=public.hol)) + theme_bw() + labs(title=paste('Centered moving averages:',outcome.var), x="Date", y=outcome.var) # temperature varies with time and no clear pattern
    
    # Seasonality * noise
    .d[,.seasonality.noise := get(outcome.var) / .centered.moving.ave] # seasonality * noise = time series value / trend
    g.s <- ggplot(.d, aes(x=Date, y=.seasonality.noise, col=dow)) + geom_point(aes(shape=public.hol.x)) + labs(title='Noise', x="Date", y=paste('Proportion:',outcome.var)) + theme.bs + theme.ts + dow.colours + pubhol.shapes
    
    # Day of the week. Daily seasonality effect
    .d[, .dow.effect :=lapply(.SD, mean, na.rm=T), by=DOW, .SDcols=outcome.var] # better control by DOW than by adding month and Fyear
    g.d <- ggplot(.d, aes(x=Date, y=.dow.effect, col=dow)) + geom_point(aes(shape=public.hol.x)) + labs(title='Day of the week effect', x="Date", y=outcome.var) + theme.bs + theme.ts + dow.colours + pubhol.shapes
    
    # Public holiday
    .d[, .phol.effect:=lapply(.SD, mean, na.rm=T), by=public.hol.x, .SDcols=outcome.var] # better control by DOW than by adding month and Fyear
    g.p <- ggplot(.d, aes(x=Date, y=.phol.effect, col=dow)) + geom_point(aes(shape=public.hol.x)) + labs(title='Public holiday effect', x="Date", y=outcome.var) + theme.bs + theme.ts + dow.colours + pubhol.shapes
    
    # Day of the week * public holiday (seasonality) component, extracted from seasonality + noise
    # if use outcome instead of seasonality + noise, represents mean effect of seasonality on outcome
    .d[, .seasonality := lapply(.SD, mean, na.rm=T), by=.(DOW,public.hol.x), .SDcols='.seasonality.noise'] # less outcome on Xmas regardless of phol
    g.dp <- ggplot(.d, aes(x=Date, y=.seasonality, col=dow)) + geom_point(aes(shape=public.hol.x)) + labs(title='DoW and PH effect', x="Date", y=paste('Proportion:',outcome.var)) + theme.bs + theme.ts + dow.colours + pubhol.shapes
    
    # Noise
    .d[,.noise := .seasonality.noise / .dow.effect / .phol.effect * mean(.dow.effect) * mean(.phol.effect)] # Modelling dow without stratification by phol, more noise seen in Sat and especially Sun.
    g.n <- ggplot(.d, aes(x=Date, y=.noise, colour=dow)) + geom_point(aes(shape=public.hol.x)) + labs(title='Noise', x="Date", y=paste('Proportion:',outcome.var)) + theme.bs + theme.ts + dow.colours + pubhol.shapes
    
    # Noise, with interaction
    .d[,.noisei := .seasonality.noise / .seasonality * mean(.seasonality)] # Modelling dow without stratification by phol, more noise seen in Sat and especially Sun.
    g.ni <- ggplot(.d, aes(x=Date, y=.noise, colour=dow)) + geom_point(aes(shape=public.hol.x)) + labs(title='Noise', x="Date", y=paste('Proportion:',outcome.var)) + theme.bs + theme.ts + dow.colours + pubhol.shapes
    
    # # Noise adjusted for log(labour force) # negligible effect # Very likely absorbed into trend
    # .d[,.noise.lab := .noise / log(n) * mean(log(n))]
    # g5 <- ggplot(.d, aes(x=Date, y=.noise.lab, colour=dow)) + geom_point(aes(shape=public.hol.x)) + labs(title='Noise adjusted for log(labour force size)', x="Date", y=paste('Proportion:',outcome.var)) + theme.bs + theme.ts + dow.colours + pubhol.shapes
    # # unique(daily.d[City=='Melbourne' & outorin=='Indoors', n])
    
    # Adjust for exposure (residual noise not explained by trend + seasonality, dow, phol and temperature)
    .d[,.fit.temp := .noise / get(exposure.var) * mean(get(exposure.var))]
    g.e <- ggplot(.d, aes(x=Date, y=.fit.temp, colour=dow)) + geom_point(aes(shape=public.hol.x)) + labs(title='Fitting temperature to noise', x="Date", y=paste('Proportion:',outcome.var)) + theme.bs + theme.ts + dow.colours + pubhol.shapes
    
    # Adjust for exposure (residual noise not explained by trend + seasonality, dow, phol and temperature)
    .d[,.fit.tempi := .noisei / get(exposure.var) * mean(get(exposure.var))]
    g.ei <- ggplot(.d, aes(x=Date, y=.fit.tempi, colour=dow)) + geom_point(aes(shape=public.hol.x)) + labs(title='Fitting temperature to noise', x="Date", y=paste('Proportion:',outcome.var)) + theme.bs + theme.ts + dow.colours + pubhol.shapes
    # View(.d[stratum=="Hobart, Indoors" & .fit.temp > 100]) # 2015-08-04, a very cold day, with very high costs, on Tuesday, with 11 illnesses. Big outlier for costs but not injuries
    
    # Combine plots
    g0 <- arrangeGrob(g.t,g.s, g.d,g.p,g.n,g.e, g.dp,g.ni,g.ei, ncol=3, nrow=3, top=textGrob(i, gp=gpar(fontsize=6,font=3)))
    ggsave(file = paste0(outcome.exposure.loc.ts, .name, ',tsn.png'), g0)
    dev.off()
  }
}
toc()

# summary(.d[, .noise])
# .d[.noise>quantile(.noise, 0.95, na.rm=T)]
# table(.d[.noise>quantile(.noise, 0.95, na.rm=T), dow])
## prop.table(table(.d[.noise>quantile(.noise, 0.95, na.rm=T), dow])) # mostly Sunday (42%) and Sat (31%)
# table(.d[.noise>quantile(.noise, 0.95, na.rm=T), dow], .d[.noise>quantile(.noise, 0.95, na.rm=T), Month]) # no major pattern



################################################################################
# FIRST-STAGE ANALYSIS: MODEL FOR EACH BY-VARIABLE COMBINATION, REDUCE AND SAVE, NO POOLING, include dummy dlnm rows
################################################################################

set.seed(3) # likely not needed, but just in case
gc()
tic()
for(i in ds.stratum) {
  # i <- 'Hobart' # 'Darwin, Outdoors'
  print(paste('Stage 1 model:',i))
  .ds <- daily.ds[stratum==i] # dataset for each unique value #
  if(any(is.na(.ds[, get(outcome.var)]))) {
    .d <- .ds[!is.na(get(outcome.var))] # no missing outcome. SHOULD NOT, BUT BE CAUTIOUS, AS IT COULD UPSET WBGT ANALYSIS
    # .d <- .ds[-(1:lmax)] # dataset for each unique value # daily.ds[stratum=='Adelaide, Outdoors']
  }
  .name <- unique(.ds$stratum) # names with all of a, b and c together
  adf <- asdf[stratum==i, get(s.sdf)] # select appropriate sdf, if it's to be used
  
  # Dependent and independent variables using .ds
  .outcome <- .ds[,get(outcome.var)]
  outcomes[[i]] <- .d[,get(outcome.var)]
  no.outcome[i,] <- sum(outcomes[[i]] != 0) # number of non-zero (and non-missing) outcomes
  temps[[i]] <- .ds[,get(exposure.var)] # exposure. Include all temperature values used for dlnm including lagged values (early ones are used less, but so are day 0 values towards end of study)
  # if (no.outcome[i,] < 100) {next} # skip iteration if number of outcomes is less than 100, which can lead to non-convergence. Still results in model output
  
  # Define exposure
  if(isTRUE(ehf.marker)) {
    .cen <- 0 # centred on 0
  } else {
    .cen <- mean(temps[[i]], na.rm=T) # centre (reference value) on mean for crosspred
  }
  
  # Calculate crossbasis. Add group if restrict analysis to summer
  .cb <- crossbasis(temps[[i]], argvar=list(fun=espline, knots=quantile(temps[[i]], eknots, na.rm=T)), lag=lmax, arglag=li_arglag)
  
  # Mean exposure and range, for meta-predictors. Note, this is only dependent on location and indoor/outdoor
  exposure[i,1] <- mean(temps[[i]], na.rm=T) # mean
  exposure[i,2] <- diff(range(temps[[i]], na.rm=T)) # range
  exposure[i,3] <- min(temps[[i]], na.rm=T) # min
  exposure[i,4] <- max(temps[[i]], na.rm=T) # max
  exposure[i,5:(length(predper.short)+4)] <- quantile(temps[[i]], predper.short/100, na.rm=T)
  exposure.by[i,] <- quantile(temps[[i]], predper/100,na.rm=T) # all percentiles and related average exposure
  exposure_ehf[i,1] <- 0
  exposure_ehf[i,2] <- quantile(temps[[i]][temps[[i]]>=0],0.85)
  exposure_ehf[i,3] <- quantile(temps[[i]][temps[[i]]>=0],0.85)*2
  exposure_ehf[i,4] <- quantile(temps[[i]][temps[[i]]>=0],0.85)*3
  ehf_threshold[i,1] <- ecdf(temps[[i]])(0) # percentile corresponding with 0. Needed for EHF. If does not exist, it becomes 1
  ehf_threshold[i,2] <- ecdf(temps[[i]])(exposure_ehf[i,2]) # percentile corresponding with severe ehf heatwave
  ehf_threshold[i,3] <- ecdf(temps[[i]])(exposure_ehf[i,3]) # percentile corresponding with severe ehf heatwave * 2
  ehf_threshold[i,4] <- ecdf(temps[[i]])(exposure_ehf[i,4]) # percentile corresponding with severe ehf heatwave * 3
  rhs[[i]] <- .d[,ave_rh] # relative humidity
  mean_rhs[i,] <- mean(.d[,ave_rh], na.rm=T) # relative humidity
  
  # Calculate second cross-basis if requested
  if(str_detect(mformula,'.cb2')) {
    if(!is.null(exposure.var2)) {
      .exposure2 <- .ds[,get(exposure.var2)]
      .cen2 <- mean(.exposure2, na.rm=T) # centre (reference value) on mean for crosspred
      .cb2 <- crossbasis(.exposure2, argvar=list(fun=espline, knots=quantile(.exposure2, eknots, na.rm=T)),
                         lag=lmax, arglag=list(fun=lspline, knots=lknots)) # add group if restrict analysis to summer
    } else stop('If include a second crossbasis, need to include a second exposure variable (humidity)')
  }
  
  # Number of years in model: 13, 12 for TAS
  no.years[i,] <- max(.d[,FYear], na.rm=T) - min(.d[,FYear], na.rm=T) + 1
  .no.years <- no.years[i,]
  # .date.knots <- seq(range(.ds[,date])[1], range(.ds[,date])[2], length=.no.years*6)
  # .date.knots <-  unique(floor_date(.ds[,Date], "month"))[c(T,F)]
  # .date.knotsl <- list(.date.knots) # evenly spaced knots
  # .no.k <- length(.date.knots)
  
  # Adjust formula as needed
  if (isTRUE(io.marker) & str_detect(.name, "utdoor")) { # if analysing by indoor/outdoor and this is an outdoor category
    form[i,] <- mformula.io # use outdoor formula
  } else {
    form[i,] <- mformula # use indoor formula
  }
  if(!is.null(perth.patch) & isTRUE(fyearp) & str_detect(i, 'Perth')) { # if perth.patch is not NULL, insert formula patchin formula. In this case, only if relevant (costs and Perth)
    form[i,] <- paste0(form[i,], perth.patch)
  }
  
  # Generate model
  if(modeltype %in% c('a','gam','GAM')) {
    if (distribution[1]=='Tweedie') {
      model[[i]] <- gam(as.formula(form[i,]), family=tw(), .ds, na.action="na.exclude", method=gam.convergence, maxit=max.iter) #  family=twlss()
      m.tweedie.shape[i,] <- as.numeric(str_remove_all(model[[i]]$family[["family"]], "[Twediep=()]")) # save shape parameter
    } else {
      model[[i]] <- gam(as.formula(form[i,]), family=distribution, data=.ds, na.action="na.exclude", method=gam.convergence, maxit=max.iter)
    }
    if(no.smooth==1) {
      tryCatch({gam.k[i,] <- k.check.MAB(model[[i]], subsample=nrow(.ds)-lmax)}, warning=function(w) print(i)) # check k' and edf. MAB version removes NA residuals, which is required with DLNM. Default subsample is 5000, which exceeds data limit. Lower data limit to be no more than the number of rows in dataset excluding the dlnm lag NA rows
    } else if(no.smooth==2) {
      .no <- which(ds.stratum==i) 
      gam.k[c(2*.no-1,2*.no),] <- k.check.MAB(model[[i]], subsample=nrow(.ds)-lmax) # use if have >1 smooth
    } 
    # GAM SE, t-values, P-values and R^2
    .m.se <- summary(model[[i]])[["se"]]
    .m.t <- summary(model[[i]])[["p.t"]] # values not produced for gam smoothers
    .m.p <- summary(model[[i]])[["p.pv"]] # values not produced for gam smoothers, that is summary(model[[i]])$s.table[,4]
    m.r2[i,] <- summary(model[[i]])[["r.sq"]] # R2
    
    # GAM specific model parameters
    gam.concurvity[[i]] <- concurvity(model[[i]], full=T) # check concurvity
    
    print(paste0(i,'. Convergence: ',model[[i]]$converged,'. Smoothing parameter: ',model[[i]]$sp,'. Total df: ',sum(model[[i]]$edf)))
    png(file = paste0(outcome.exposure.loc.s1.smoother, .name,'- Date smoother', '.png')) # save plot for smoother
    plot(model[[i]], pages=0, main=paste0("Smoother for date - ",.name))
    dev.off()
  } else {
    if (distribution[1]=='Tweedie') { 
      var.power[[i]] <- tweedie.profile(formula=form[i,], data=.ds, p.vec=var.power.values, method='series', do.ci=F, do.smooth=tp.smooth) # estimate optimal shape parameter. Turn off CI for speed
      m.tweedie.shape[i,] <- var.power[[i]]$p.max # save shape parameter
      model[[i]] <- glm(as.formula(form[i,]), family=tweedie(var.power=m.tweedie.shape[i,], link.power=0), .ds, na.action="na.omit", control=list(maxit=max.iter)) # "na.exclude". na.omit works better for Tweedie e.g. AIC calculation
      # model.omit[[i]] <- glm(as.formula(form[i,]), family=tweedie(var.power=m.tweedie.shape[i,], link.power=0), .ds, na.action="na.omit", control=list(maxit=max.iter)) # NAs will cause tweedie AIC calculation to fail 
    } else {
      model[[i]] <- glm(as.formula(form[i,]), family=distribution, .ds, na.action="na.exclude", control=list(maxit=max.iter))
    }
    .collin <- suppressMessages(vif(model[[i]], terms='high-order')) # car::vif(model[[i]], terms='high-order') # suppress "there are higher-order terms (interactions) in this model consider setting type = 'predictor'; see ?vif"
    collin[[i]] <- data.table(Location=i, Coefficient=rownames(.collin), .collin)
    collin[[i]][!1, Location:=NA] # only keep location for 1st row
    # vif.gam.MAB(model[[i]]) # creating collinearity function for GAM
    # mgcv.helper::vif.gam(model[[i]])
    
    # GLM SE, t-values, P-values and R^2
    .glm.coef <- summary(model[[i]])$coefficients
    .m.se <- .glm.coef[,2]
    .m.t <- .glm.coef[,3]
    .m.p <- .glm.coef[,4]
    m.r2[i,] <- 1 - model[[i]]$deviance/model[[i]]$null.deviance # ?Nagelekre R&2
  }
  
  # Save coefficients
  .m.coef <- coef(model[[i]])
  .m.row <- names(.m.coef)
  .m.lci <- .m.coef - qnorm(0.975)*.m.se
  .m.uci <- .m.coef + qnorm(0.975)*.m.se
  suppressWarnings(m.coef[[i]] <- data.table(location=i, Parameter=.m.row, Coefficient=.m.coef, SE=.m.se, LCI=.m.lci, UCI=.m.uci, RR=exp(.m.coef), 'RR LCI'=exp(.m.lci), 'RR UCI'=exp(.m.uci), 't-value'=.m.t, 'p-value'=round(.m.p,8))) # combine coefficients into a data.table
  m.coef[[i]][!1, location:=NA] # only keep location for 1st row
  m.coef[[i]][`p-value` <0.1, sig:='.']
  m.coef[[i]][`p-value` <0.05, sig:='*']
  m.coef[[i]][`p-value` <0.01, sig:='**']
  m.coef[[i]][`p-value` <0.001, sig:='***']
  if(modeltype %in% c('a','gam','GAM')) {
    m.coef[[i]][str_detect(Parameter,'\\('), ':='(`t-value`=NA, `p-value`=NA)] # t-values and p-values not calculated for gam smoothers, so set NA instead of repeat for them
  }
  
  # Residuals. Default residuals are deviance, but changed to quantile if Tweedie
  if (distribution[1]=='Tweedie' | isTRUE(res.quantile)) {
    res[[i]] <- qresiduals(model[[i]]) # quantile residuals for glm, recommended when deviance and Pearson residuals are grossly non-normal. Doesn't work with nb
  } else {
    res[[i]] <- residuals(model[[i]]) # default is deviance for gam and glm
  }
  
  # Variables for model checking
  # if see Warning in as.data.table.list(x, keep.rownames = keep.rownames, check.names = check.names,  : Item 1 has <4748 rows but longest item has 4748, it ooccurs in .model checks
  suppressWarnings(.cooks <- cooks.distance(model[[i]])[!is.na(res[[i]])]) # as dlnm with missing residuals, will return warning In res/(1 - hat) : longer object length is not a multiple of shorter object length
  .model.checks <- model.checks[[i]] <- data.table(res=res[[i]][!is.na(res[[i]])], fitted=fitted.values(model[[i]])[!is.na(res[[i]])], cooksd=.cooks,
                                                   Date=.d[,Date], DoW=.d[,DoW], Week=as.factor(.d[,week]), month=.d[,month], Year=.d[,Year], FYear=.d[,FYear], public.hol=.d[,public.hol],
                                                   shol=.d[,shol], PH=.d[,get('Public holiday')], day1=.d[,day1], school.hols=.d[,school.hols],
                                                   outcome=outcomes[[i]], temp=.d[,get(exposure.var)])
  
  # Aggregated model checks
  model.checks.month[[i]] <- .model.checks[, lapply(.SD, mean, na.rm=T), by=.(month,Year), .SDcols='res']
  model.checks.month[[i]][,Date:=as.Date(paste(1, month, Year, sep="."), format = "%e.%b.%Y")] # create date for 1st month
  
  model.checks.week[[i]] <- merge(model.checks[[i]][, lapply(.SD, mean, na.rm=T), by=Week, .SDcols='res'], .ds[DoW=='Mon',.(Week,Date)], all.x=T) # add Date for easier checking later
  
  # Residual plots
  # TO UPDATE: From left to right, Sun=yellow, Sat=purple, Mon=black, Tues=red, Wed=green, Thurs=blue, Fri=teal
  plot.h %<a-% {hist(.model.checks[,res], main='Residual histogram', xlab='Residuals', breaks=50); abline(v=0, lty=2, lwd=2, col='blue')} # residual distribution
  plot.qq %<a-% {qqnorm(.model.checks[,res], main="QQ-plot", las=1, col=.model.checks[,DoW]); qqline(.model.checks[,res])} # QQplot
  # car::qqPlot(.model.checks[,res], main="QQ-plot",las=1) # QQplot with 95% confidence envelope, also prints most extreme positive and negative residual numbers
  plot.of %<a-% {plot(x=.model.checks[,outcome], y=.model.checks[,fitted], pch=16, cex=0.5, col=.model.checks[,DoW], main='Outcome vs fitted values', ylab="Outcome", xlab="Fitted values"); abline(a=0, b=1, lty=2, lwd=2, col='blue')} # Outcome vs fitted values, with straight ilne
  plot.rf %<a-% {plot(x=.model.checks[,fitted], y=.model.checks[,res], pch=16, cex=0.5, col=.model.checks[,DoW], main='Residuals vs fitted values', ylab="Residuals", xlab="Fitted values"); abline(h=0, lty=2, lwd=2, col='blue')} # residual variance, must remove NA residuals to plot against date
  # ggplot(.model.checks) + geom_point(aes(x=fitted, y=res, colour=dow)) + theme_bw() + labs(title='Residuals vs fitted values', x="Fitted values", y="Residuals") + geom_hline(yintercept=0, linetype='dashed', color='blue')
  plot.rd %<a-% {plot(x=.model.checks[,Date], .model.checks[,res], pch=16, cex=0.5, col=.model.checks[,DoW], main='Residuals vs date', ylab="Residuals", xlab="Date"); abline(h=0, lty=2, lwd=2, col='blue')} # residual variance, must remove NA residuals to plot against date
  plot.re %<a-% {plot(x=.model.checks[,temp], .model.checks[,res], pch=16, cex=0.5, col=.model.checks[,DoW], main='Residuals vs exposure', ylab="Residuals", xlab=exposure.var); abline(h=0, lty=2, lwd=2, col='blue')} # residual variance, must remove NA residuals to plot against date
  bplot.rd %<a-% boxplot(.model.checks[,res] ~ .model.checks[,DoW], main='Residuals vs DoW', ylab="Residuals", xlab="Day of the week")
  bplot.rm %<a-% boxplot(.model.checks[,res] ~ .model.checks[,month], main='Residuals vs month', ylab="Residuals", xlab="Day of the week")
  bplot.ry %<a-% boxplot(.model.checks[,res] ~ .model.checks[,FYear], main='Residuals vs Fyear', ylab="Residuals", xlab="Financial year")
  bplot.rp %<a-% boxplot(.model.checks[,res] ~ .model.checks[,public.hol], main='Residuals vs PH', ylab="Residuals", xlab="Public holiday?")
  bplot.ro %<a-% boxplot(.model.checks[,res] ~ .model.checks[,shol], main='Residuals vs specific PH', ylab="Residuals", xlab="Special PH")
  bplot.r1 %<a-% boxplot(.model.checks[,res] ~ .model.checks[,day1], main='Residuals vs Day 1 of month', ylab="Residuals", xlab="Day 1 of month")
  bplot.rs %<a-% boxplot(.model.checks[,res] ~ .model.checks[,school.hols], main='Residuals vs school holidays', ylab="Residuals", xlab="School holidays")
  plot.acf %<a-% {Acf(.model.checks[,res], na.action=na.omit, ylab='ACF', xlab="Lag (days)", main=''); title(main='Autocorrelation function', line=1)} # autocorrelation, without day 0. title sep to move closer to plot, as default is further away than normal
  plot.pacf %<a-% {Pacf(.model.checks[,res], na.action=na.omit, ylab='PACF', xlab="Lag (days)", main=''); title(main='Partial autocorrelation function', line=1)} # partial autocorrelation
  # plot(model[[i]], xlab='Date', ylab='Trend', se=T, main='Trend vs date') # se=T is default, plotting +/-2 SE 
  # plot(x=.model.checks[,Date], y=.model.checks[,cooksd], pch="*", cex=2, main="Cook's distance", xlab='Date', ylab="Cook's distance'"); abline(h = 4/nrow(.model.checks), col="red")  # Cook's distance
  # text(x=1:.model.checks[,.N], y=.model.checks[,cooksd], labels=.model.checks[cooksd > 4/.N, Date], col="blue")  # add labels
  plot.mrd %<a-% {plot(x=model.checks.month[[i]][,Date], model.checks.month[[i]][,res], pch=16, cex=0.5, main='Monthly residuals vs date', ylab="Residuals", xlab="Date"); abline(h=0, lty=2, lwd=2, col='blue')} # residual variance, must remove NA residuals to plot against date
  plot.wrd %<a-% {plot(x=model.checks.week[[i]][,Week], model.checks.week[[i]][,res], pch=16, cex=0.5, main='Weekly residuals vs Week', ylab="Residuals", xlab="Week"); abline(h=0, lty=2, lwd=2, col='blue')} # residual variance, must remove NA residuals to plot against date
  
  # Combine plots (large)
  png(file = paste0(outcome.exposure.loc.s1.mcheck,'zL ',.name,' mc.png')) # plot location
  par(mfrow=c(4,4), mar=rep(2.8,4), mgp=c(1.5,0.5,0)) # lower margin sizes
  plot.qq; plot.h; plot.of; plot.rf
  plot.rd; plot.re; bplot.rd; bplot.rm
  bplot.ry; bplot.rp ; bplot.ro; bplot.r1
  bplot.rs; plot.acf; plot.pacf; plot.wrd 
  dev.off()
  
  # Combine plots (small)
  png(file = paste0(outcome.exposure.loc.s1.mcheck, .name,' mc.png')) # plot location
  par(mfrow=c(3,3), mar=rep(3,4), mgp=c(1.5,0.5,0)) # lower margin sizes
  plot.qq; plot.h; plot.of
  plot.rf; plot.rd; plot.re
  bplot.rd; plot.acf; plot.pacf
  dev.off() 
  
  # Plotting residuals against covariates, to determine if modelling patterns missed 
  .check.fourier <- harmonic(.d[,Date], nfreq=6, period = nrow(.d)/.no.years)
  .check.res.coef <- coef(summary(lm(model.checks[[i]][,res] ~ log(.d[,n]) + .d[,Tue]+.d[,Wed]+.d[,Thu]+.d[,Fri]+.d[,Sat]+.d[,Sun] +
                                       .d[,public.hol] *.d[,DoW]*+ .d[,month]*.d[,DoW] + .d[,fyear]*.d[,DoW] + .d[,shol]*.d[,DoW] + .d[,day1]*.d[,DoW] + .d[,school.hols]*.d[,DoW] + 
                                       .d[,month]:.d[,fyear] + .d[,month]:.d[,day1] + .d[,month]:.d[,school.hols] +
                                       .d[,fyear]:.d[,day1] + .d[,fyear]:.d[,school.hols] + .d[,day1]:.d[,school.hols] + .check.fourier)))
  .check.res.coef.names <- rownames(.check.res.coef)
  .check.res.coef.sig <-  case_when(
    .check.res.coef[,4] < 0.001 ~ '***',
    .check.res.coef[,4] < 0.01 ~ '**',
    .check.res.coef[,4] < 0.05 ~ '*',
    .check.res.coef[,4] < 0.1 ~ '.',
    TRUE ~ as.character(NA))
  check.res.coef[[i]] <- data.table(Location=i, Coefficient=.check.res.coef.names, .check.res.coef, Significance=.check.res.coef.sig)
  check.res.coef[[i]][!1, Location:=NA] # only keep location for 1st row
  check.res.coef.sig[[i]] <- check.res.coef[[i]][,.(Coefficient, Significance)] # condensed for side-to-side comparison
  setnames(check.res.coef.sig[[i]], old=c('Coefficient','Significance'), new=c('Coefficient',paste(i,'Sig'))) # rename with i for merging
  
  if(distribution[1]=='quasipoisson') { #[1] is to avoid warnings from using multiple elements
    if(modeltype %in% c('a','gam','GAM')) {
      m.aic[i,] <- fqaic.fn(model[[i]], gam=T) # if gam, use gam option. CURRENTLY REULSTS IN INF. Not sure (I suspect not) if it applies corrected AIC, but also not sure if GAM AIC works on it
    } else {
      m.aic[i,] <- fqaic.fn(model[[i]])
    } 
  } else if(distribution[1]=='Tweedie') {
    if(modeltype %in% c('a','gam','GAM')) {
      m.aic[i,] <- AIC(model[[i]]) # conditional AIC, which is better for GAM # model[[i]]$aic is uncorrected AIC
    } else {
      m.aic[i,] <- AICtweedie(model[[i]])
    }
  } else {
    m.aic[i,] <- AIC(model[[i]]) # model[[i]]$aic is identical
  }
  m.devexp[i,] <- (model[[i]]$null.deviance - model[[i]]$deviance) / model[[i]]$null.deviance # deviance explained
  m.dispersion[i,] <- sum(residuals(model[[i]], type="pearson")^2, na.rm=T)/df.residual(model[[i]]) # summary(model[[i]])[["dispersion"]]
  # for(i in ds.stratum) { # chi-square test for over-dispersion
  #   print(c(i, round(sum(residuals(model[[i]], type="pearson")^2, na.rm=T)/df.residual(model[[i]]),3),
  #           round(pchisq(sum(residuals(model[[i]], type="pearson")^2, na.rm=T), df.residual(model[[i]]), lower.tail=F),3)))
  # }
  
  # Sum (accumulate) effects of all lags in order to eliminate one dimension of the association
  # Predicted effects: extract parameters from model corresponding to .cb variables through functions coef and vcov
  .pred1 <- crosspred(basis=.cb, model=model[[i]], cen=.cen, model.link='log') # must set log for tweedie package 
  if(str_detect(form[i,],'.cb2')) { # if above When used with .cb2, Error in crosspred(basis = .cb, model = model[[i]], cen = .cen) : coef/vcov not consistent with basis matrix. See help(crosspred). Does work with pred2 though
    .pred2 <- crosspred(basis=.cb2, model=model[[i]], cen=.cen2, model.link='log')
  }
  # .pred1 <- crosspred(basis=.cb, cen=.cen, coef=model[[i]]$coefficients, vcov=model[[i]]$coefficients, model.link='log')
  
  # RR is y-axis for (quasi-)poisson, but not for Tweedie model. So what is y-axis for Tweedie?
  # Plot exposure-lag-response relationship
  png(file = paste0(outcome.exposure.loc.s1.r, .name,', e-l-r relationship.png')) # plot location
  plot(.pred1, "3d", xlab=exposure.var, ylab="Lag", zlab="Relative risk",
       main=paste0("e-l-r relationship: ",.name))
  dev.off()
  
  # Plot exposure-lag-response relationship as contour plot
  png(file = paste0(outcome.exposure.loc.s1.r, .name,', e-l-r contour plot.png'))
  plot(.pred1, "contour", xlab=exposure.var, ylab="Lag", main=paste0("e-l-r relationship: ",.name))
  dev.off()
  
  # plot(.pred1, "3d", xlab=exposure.var, ylab="Lag", zlab="Relative risk", main=paste0("e-l-r relationship: ",.name))
  
  # Plot exposure-response relationship at day 0 (no lag)
  # : allow in title names but not file names
  # add ylim for claims, but not for costs
  png(file = paste0(outcome.exposure.loc.s1.r,.name, ', e-r relationship day 0', '.png'))
  plot(.pred1, "slices", lag=0, lwd=1, col=2,  xlab=paste(exposure.var,'?(C)'), # ylim=c(0.85,1.10),
       ylab="Relative risk", main= paste0("Exposure-response relationship at day 0: ", .name)) 
  dev.off() # Save image + clear settings
  
  # Plot exposure-response relationship at all lag days
  png(file = paste0(outcome.exposure.loc.s1.r,.name, ', e-r relationships', '.png'))
  plot(.pred1, "slices", lag=0, lwd=1, col=2, xlab=paste(exposure.var,'?(C)'), # ylim=c(0.85,1.10), 
       ylab="Relative risk", main= paste0("Exposure-response relationships: ", .name))
  for(l in seq(1:lmax)) {
    lines(.pred1, "slices", lag=l, col=l+2) # for each lag day excluding lag 0 (main plot)
  }
  legend("bottomright", paste("Lag",seq(0:lmax)-1), col=(seq(0:lmax)+1), lwd=1) # ?topleft
  dev.off()
  
  # Overall cumulative summary for main model
  red[[i]] <- crossreduce(.cb, model[[i]], cen=.cen, model.link='log') # reduce exposure-lag-response association to overall cumulative exposure-response association, using mean
  coef[i,] <- coef(red[[i]])
  vcov[[i]] <- vcov(red[[i]])
  
  exposure.rr.s1[[i]] <- data.frame(matrix(nrow=length(red[[i]]$predvar), ncol=4))
  colnames(exposure.rr.s1[[i]]) <- c("temp", "RRfit", "RRlow", "RRhigh")
  exposure.rr.s1[[i]]$temp <- red[[i]]$predvar # = .pred1$predvar
  
  if(any(str_detect(names(red[[i]]), 'RRfit'))) { # default dlnm behaviour for log-link
    exposure.rr.s1[[i]]$RRfit  <- red[[i]]$RRfit # already %, do not need to convert like in supp. Replaced RRfit/low/high with fit/low/high. Crosspred does output se, ?can manully calculate robust SE
    exposure.rr.s1[[i]]$RRlow <- red[[i]]$RRlow
    exposure.rr.s1[[i]]$RRhigh <- red[[i]]$RRhigh
  }  else if(any(str_detect(names(red[[i]]), 'fit'))) { # if Tweedie glm, which outputs a log-link value not recognised by dlnm
    exposure.rr.s1[[i]]$RRfit  <- exp(red[[i]]$fit) # already %, do not need to convert like in supp. Replaced RRfit/low/high with fit/low/high. Crosspred does output se, ?can manully calculate robust SE
    exposure.rr.s1[[i]]$RRlow <- exp(red[[i]]$low)
    exposure.rr.s1[[i]]$RRhigh <- exp(red[[i]]$high)
  } 
  
  # Overall cumulative exposure-response relationship plot. We're not using these plots for results, but for checking
  # Convert RR into %
  .exposure.rr.s1.rrfit <- (exposure.rr.s1[[i]]$RRfit-1)*100 
  .exposure.rr.s1.rrlow <- (exposure.rr.s1[[i]]$RRlow-1)*100
  .exposure.rr.s1.rrhigh <- (exposure.rr.s1[[i]]$RRhigh-1)*100
  
  .ind1 <- .pred1$predvar<=.cen # values below centering value (blue)
  .ind2 <- .pred1$predvar>=.cen # values above centering value (red)
  
  png(file = paste0(outcome.exposure.loc.s1.oer, .name,', Overall e-r.png')) # plot location. File name based on heat metric
  plot(exposure.rr.s1[[i]]$temp, .exposure.rr.s1.rrfit, type="n", ylim=c(rryaxis[1], rryaxis[2]), lwd=2, col="white",
       main=i, ylab="Percent change (%)", xlab=paste(exposure.var,'(Â°C)'), 
       cex.main=cexmain, cex.lab=cexlab, cex.axis=cexaxis, lab=c(6,5,7))
  .erplot %<a-% { # save main plot additions.
    polygon(c(exposure.rr.s1[[i]]$temp,rev(exposure.rr.s1[[i]]$temp)),c(.exposure.rr.s1.rrlow,rev(.exposure.rr.s1.rrhigh)), col="grey89", border=F) # 95% CI envelope
    lines(exposure.rr.s1[[i]]$temp[.ind1],.exposure.rr.s1.rrfit[.ind1],col=4,lwd=2); # cold, left of cen
    lines(exposure.rr.s1[[i]]$temp[.ind2],.exposure.rr.s1.rrfit[.ind2],col=2,lwd=2); # hot, right of cen
    abline(h=0) # horizontal line
  }
  .lines %<a-% { #  vertical lines
    abline(v=.cen,lty=2); # centering value
    abline(v=c(exposure[i,c("1%","99%")]),lty=3) # 1st and 99th percentiles
  }
  .erplot # insert main plot additions
  .lines # insert lines
  dev.off() # Save image + clear settings
} 
toc()
# traceback()



################################################################################
# Comvine and review results of Stage 1 Analysis
################################################################################

# Combine number of outcomes, DLNM coefficients, AIC and dispersion parameter, including totals for outcomes and AIC

m1.results <- rbind(cbind(no.outcome, coef, m.aic, m.r2, m.devexp, m.dispersion), c(sum(no.outcome),rep(NA, ncol(coef)),sum(m.aic), rep(NA,3)))
colnames(m1.results) <- c('Number of outcomes',colnames(coef),'AIC','R^2','Deviance explained','Dispersion parameter')
write.csv(m1.results, file=paste0(outcome.exposure.loc.s1,'Outcome, coefficients and model fit.csv'), na='', row.names=T) # create csv file
save(m1.results, file=paste0(outcome.exposure.loc.s1,'Outcome, coefficients and model fit.rda')) # save dataset for easier access
write.csv(sum(m.aic), file=paste0(outcome.exposure.loc.s1,'Total AIC.csv'), na='', row.names=T)

# Output model coefficients and RRs
write.csv(do.call(rbind, m.coef), file=paste0(outcome.exposure.loc.s1,'Coefficients.csv'), na='') # Model coefficients, SE, 95% CI and t-tests
# do.call(rbind, m.coef)[str_detect(parameter,'day1')] 

# Modelling against residuals. Not output varies
check.res.coef.s <- check.res.coef.sig %>% purrr::reduce(full_join, by = "Coefficient") # p-values only, next to each other
check.res.coef.f <- do.call(rbind, check.res.coef) # all measurements, arranged vertically
write.csv(check.res.coef.s, file=paste0(outcome.exposure.loc.s1.mcheck,'zLM residuals with covariates.csv'), na='', row.names=F)
write.csv(check.res.coef.f, file=paste0(outcome.exposure.loc.s1.mcheck,'zLM residuals with covariatesf.csv'), na='', row.names=F) # f short for full
# View(check.res.coef.s[str_detect(Coefficient,'Sat') | str_detect(Coefficient,'Sun')])

# Collinearity or GAM specific checks
if(modeltype %in% c('a','gam','GAM')) {
  write.csv(gam.k, file=paste0(outcome.exposure.loc.s1,'k-index.csv'), na='', row.names=T) # create csv file
  write.csv(do.call(rbind, gam.concurvity), file=paste0(outcome.exposure.loc.s1,'Concurvity.csv'), na='', row.names=F)
} else {
  write.csv(do.call(rbind, collin), file=paste0(outcome.exposure.loc.s1,'Collinearity.csv'), na='', row.names=F)
}

# Tweedie shape parameters
if (distribution[1]=='Tweedie') { 
  write.csv(m.tweedie.shape, file=paste0(outcome.exposure.loc.s1,'Tweedie shape parameters.csv'), na='', row.names=T)
} 

# # Check results or residuals to determine patterns
# for(i in ds.stratum) {
# # for(i in c('Sydney, Indoors','Sydney, Outdoors')) {
#   print(i)
#   # print(summary(model[[i]]))
#   # print(model.checks.week[[i]][abs(res)>2])
#   # print(model.checks[[i]][Week %in% model.checks.week[[i]][abs(res)>2, Week] ]) # print weeks with high mean absolute residuals])
#   # print(model.checks[[i]][str_detect(PH,'Australia Day')])
#   # print(cor(model.checks[[i]][,fitted], model.checks[[i]][,outcome])) # Hobart & especially Darwin models are poor predictors of outcome. All others have correlation > 0.5
#   check[[i]] <- model.checks[[i]][abs(res) > 3.1] # model.checks[[i]][abs(res) > quantile(abs(res),0.995)]
#   # check[[i]] <- model.checks[[i]][res > 2.5] # model.checks[[i]][abs(res) > quantile(abs(res),0.995)]
#   print(head(check[[i]][,-'day1'][order(res)],100))
#   # print(model.checks[[i]][month=='Jun' & Year=='2008' & day(Date) %in% c(24:30)]) # 25-30 June 18
#   # print(model.checks[[i]][month=='Jun' & fyear=='20' & Week ==678]) # 25-30 June 18
#   # print(model.checks[[i]][month=='Jun' & day(Date)==29]) # 25-30 June 18
#   # print(model.checks[[i]][month=='Jun', lapply(.SD, mean, na.rm=T), by=day(Date), .SDcols = 'res']) # days in Jun, mean reisdual across all fyears. Sydney both, Melbourne both and Adelaide indoors only day in June for abs(res) > 1 (Melbourne indoors > 2). Rarely other days exceeded 1
#   # print(model.checks[[i]][month=='Jan', lapply(.SD, mean, na.rm=T), by=day(Date), .SDcols = 'res']) # no noticeable pattern for Jul. Jan inconsistenly from 2nd to 7th Jan (2nd Jan more likely), probably 2nd to 4th Jan
#   # print(model.checks[[i]][month=='Dec', lapply(.SD, mean, na.rm=T), by=day(Date), .SDcols = 'res']) # a few days before Xmas, Xmas or Boxing Day. 24th and 25th is safe (negative everywhere except Darwin indoors where 24 is negligbly positive). 23rd often positive (probably not in xmas period)
#   # print(model.checks[[i]][month=='Apr', lapply(.SD, mean, na.rm=T), by=day(Date), .SDcols = 'res']) # a few days before Xmas, Xmas or Boxing Day. 24th and 25th is safe (negative everywhere except Darwin indoors where 24 is negligbly positive). 23rd often positive (probably not in xmas period)
#   # print(model.checks[[i]][, lapply(.SD, mean, na.rm=T), by=month, .SDcols = 'res']) # Darwin and Hobart have all months negative. Negative common when excessive seasonality, but degree of negative is higher
#   # print(table(Month.day=check[[i]][,.(month,day(Date))]))
#   # print(table(Month.DoW=check[[i]][,.(month,DoW)]))
#   # print(table(Month.fyear=check[[i]][,.(month,fyear)]))
#   # print(table(public.hol=check[[i]][,public.hol]))
#   # print(table(shol=check[[i]][,shol]))
#   # print(table(school.hols=check[[i]][,school.hols]))
#   # # print(table(check[[i]][month=='Dec',day(Date)]))
#   # # print(table(check[[i]][month=='Jan',day(Date)]))
#   # # print(table(check[[i]][month=='Jun' & fyear=='2017',day(Date)]))
#   #
#   # # print(check[[i]][month=='Jan' & day(Date) %in% c(1:7) & abs(res) > 3,][order(res)])
#   # print(check[[i]][month=='Jun' & fyear=='2017',])
#   #
#   # # .ds <- daily.ds[stratum==i]
#   # # print(doBy::summaryBy(total ~ DoW, data=.ds, FUN=c(min, max, mean, sd), na.rm=T))
# }
# plot(y=daily.ds[City=='Perth' & Year %in% c(2011,2012), get(outcome.var)], x=daily.ds[City=='Perth' & Year %in% c(2011,2012), Date])
# plot(y=daily.ds[City=='Perth' & FYear==2011, get(outcome.var)], x=daily.ds[City=='Perth' & FYear==2011, Date])
# form
# predict(model[[i]], list=(month='Jan'))
# predict(model[[i]], type='response')
# 
# jtools::effect_plot(model[[i]], pred = bs(Date,6*.no.years))
# plot(y=daily.ds[City=='Perth', get(outcome.var)], x=daily.ds[City=='Perth', Date])


# for(i in ds.stratum) {
#   print(i)
#   .ds <- daily.ds[stratum==i] # dataset for each unique value #
#   # .ds <- daily.ds[stratum=='Sydney, Indoors']
#   .d <- .ds[-(1:lmax)] # dataset for each unique value # daily.ds[stratum=='Adelaide, Outdoors']
#   .name <- unique(.ds$stratum) # names with all of a, b and c together
#   adf <- asdf[stratum==i, get(s.sdf)] # select appropriate sdf, if it's to be used
# 
#   # Dependent and independent variables using .ds
#   .outcome <- .ds[,get(outcome.var)]
#   no.outcome[i,] <- sum(.outcome != 0) # number of non-zero outcomes
#   temps[[i]] <- .ds[,get(exposure.var)] # exposure
#   # if (no.outcome[i,] < 100) {next} # skip iteration if number of outcomes is less than 100, which can lead to non-convergence. Still results in model output
# 
#   # Define exposure
#   .cen <- mean(temps[[i]], na.rm=T) # centre (reference value) on mean for crosspred
# 
#   print(cor(model.checks[[i]][,res], .d[,get(outcome.var)]))
# }



######################### END ############################
######################### END ############################


