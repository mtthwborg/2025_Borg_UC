
########################################################################################################
# SECOND-STAGE ANALYSIS: MULTIVARIATE META-ANALYSIS FOR OVERALL RESULTS
########################################################################################################
##################################################
### Multivariate meta-analysis functions
##################################################

## Perform Wald test with mixmeta or mvmeta object
## fwald2 from Gasparrini et al. 2015 Lancet with additional "short" option (if T, also outputs Wald statistic and df)
fwald2.fn <- function(model, var, short=F) {
  ind <- grep(var,names(coef(model)))
  coef <- coef(model)[ind]
  vcov <- vcov(model)[ind,ind]
  waldstat <- coef%*%solve(vcov)%*%coef
  df <- length(coef)
  p <- 1-pchisq(waldstat,df)
  if(isTRUE(short)) { # return p-value only
    return(p)
  } else { # return statistic, df and p-value
    results <- c(waldstat,df,p)
    names(results) <- c('Wald statistic','df','p-value')
    return(results)
  }
}


## Report mixmeta Cochran Q test, I^2 and AIC statistics from a mixmeta or mvmeta object (mv.)
mv.results.fn <- function(mv.) {
  sum.mv.q <- summary(mv.)[["qstat"]] # Q-stat values. Can also get with qtest()
  .vars <- attr(mv.[['terms']],"term.labels") # fixed variables
  if(identical(.vars,character(0))) { # if no fixed predictors
    .vars.wald <- NULL
  } else {
    .vars.wald <- c(sapply(.vars, function(w) fwald2.fn(mv., w))) # results from multivariate Wald test
    names(.vars.wald) <- paste(rep(.vars, each=3), c('W','df','p'), sep = ".")
  }
  .mv.results <- c('Q-statistic'=sum.mv.q$Q[1], 'df'=sum.mv.q$df[1], 'P-value'=sum.mv.q$pvalue[1], 'I^2'=summary(mv.)[["i2stat"]][1], 'AIC'=AIC(mv.), .vars.wald)
  names(.mv.results) <- str_remove_all(names(.mv.results), '..all') # remove ..all  from mixmeta code
  return(.mv.results) # return combined results
}


## Wrap format() around round(). Enables more control over presentation of rounded results. Converts results to character
# round.fn <- function(x, digits=0) {format(round(x, digits), nsmall=digits, trim=T)}
round.fn <- function(x, round=0,trim=T,big.mark=",",significant=NULL, justify=c("left","right","centre","none"),width=NULL,na.encode=T,scientific=NA,big.interval=3L,small.mark="",small.interval=5L,decimal.mark=getOption("OutDec"),zero.print=NULL,drop0trailing=F,...) {
  .x <- round(x, digits=round)
  .x <- format(.x, nsmall=round, trim=trim, big.mark=big.mark, digits=significant,justify=justify,width=width,na.encode=na.encode,scientific=scientific,big.interval=big.interval,small.mark=small.mark,small.interval=small.interval,decimal.mark=decimal.mark,zero.print=zero.print,drop0trailing=drop0trailing,...)
  return(.x)
}



##################################################
### Meta-analysis: results of reduced coef and compute BLUPS (best linear unbiased prediction)
##################################################

exposure.country.mean <- colMeans(exposure.by) # could arguably weight by area number of years (study frame), but Gasaparrini didn't do this either and impact likely negligible

## Calculate potential meta-analysis predictors 
exposure.mean <- exposure[,1] # mean exposure
mm.pred.m <- matrix(rep(1,length.ds.stratum), nrow=length.ds.stratum, dimnames=list(ds.stratum,'1'))
by.vars.ds <- cbind(unique(daily.ds[,mget(by.vars3)]), mm.pred.m) # data frame with categories and mixmeta.predictors

## Meta-analysis. Default methods is REML
mv <- mixmeta(coef~1, vcov, data=by.vars.ds)
mv.results <- mv.results.fn(mv) 
write.csv(mv.results, file=paste0(s2results, 'MM tests.csv'), na='', row.names=T) # create csv file


## Pool estimated overall cumulative exposure-response associations to obtain overall estimation at each location level for main model
blups <- blup(mv, vcov=T) # BLUPs. Element for each stratum



##################################################
### Define minimum occupational injuries values (not default, but still outputted)
##################################################

# Generate matrix for storing re-centred results
minpercity <- mintempcity <- rep(NA,length.ds.stratum)
names(mintempcity) <- names(minpercity) <- ds.stratum

# Define minimum injuries values, excluding very low and hot Ts (<1st and >99th percentiles)
for(i in ds.stratum) {
  .ds <- daily.ds[stratum==i] # dataset for each unique value # daily.ds[stratum=='Hobart: Indoor']
  
  # Use same exposure-response relationship, but for exposure percentiles 1-99
  .predvar <- quantile(temps[[i]], 1:99/100,na.rm=T)
  .argvar.bound <- list(x=.predvar, fun=espline, knots=quantile(temps[[i]], eknots, na.rm=T), Bound=range(temps[[i]], na.rm=T)) # list of arguments
  
  # Redefine function using all arguments, boundary knots included
  .blups <- blups[[which(ds.stratum==i)]] 
  .bvar.bound <- do.call(onebasis, .argvar.bound) # basis for exposure-response
  .bvar1 <- .bvar.bound%*%.blups$blup # 
  minpercity[i] <- (1:99)[which.min(.bvar1)] # if stop using 1:99, may need to change code so that position links to actual percentile
  mintempcity[i] <- quantile(temps[[i]], minpercity[i]/100, na.rm=T)
}



################################################################################
# Predict pooled overall cumulative associations
################################################################################

## Mean values of mixmeta predictors as a data frame. In this case none, so no change
datanew <- as.data.frame(t(colMeans(mm.pred.m))) # names of datanew and mvpred must match
mvpred <- predict.mixmeta(mv, datanew, vcov=T, format="list") # Results are identical to that of mixmeta because there are no meta-predictors

## Define exposure percentile with lowest incidence of injuries
bvar <- crossbasis(exposure.country.mean, argvar=list(fun=espline, knots=exposure.country.mean[paste0(eknots*100)]), arglag=li_arglag) # crossbasis
bvar1 <- bvar%*%mvpred$fit # multiply coefficients for combined effect
bvar1 <- bvar1[order(bvar1[,1]),, drop=F] # sort in ascending order (not required, but easier for visualisation
bvar1 <- bvar1[between(rownames(bvar1), cenpen.minmax[1], cenpen.minmax[2]),] # exclude percentiles outside range for selecting optimal percentile, converts to vector
cenindcountry <- as.numeric(names(bvar1[which.min(bvar1)])) # exposure percentile with lowest incidence of outcome (optimal percentile)

# Define centering percentile for country, using either median or above exposure percentile 
if (is.null(cenpen)) {
  cenpercountry <- pmin(pmax(predper[cenindcountry],cenpen.minmax[1]),cenpen.minmax[2]) # min/max is 10th/90th. chooses 10, as incidence lowest at 10
} else if(is.numeric(cenpen)) {
  cenpercountry <- cenpen # use defined centering value
} else if(cenpen %in% c('mean','Mean','average','Average')) {
  cenpercountry <- 'mean' # to be used as code for mean
} 

if(is.numeric(cenpercountry)) {
  centre <- exposure.country.mean[paste0(cenpercountry)] # centered prediction (from corresponding percentile)
} else {
  centre <- sum(sapply(temps, sum)) / sum(sapply(temps, length)) # mean exposure across all strata
}

# Predict pooled overall cumulative associations (RR) for main model
cp <- crosspred(bvar, coef=mvpred$fit, vcov=mvpred$vcov, model.link="log", at=exposure.country.mean, cen=centre)



################################################################################
# LAG-RESPONSE ASSOCIATIONS, USING SELECTED CENTERING VALUE
################################################################################

# Objects to store overall cumulative exposure results (lag-response association at set percentiles vs centering value)
cvlag.length <- length(crossreduce(.cb, model=model[[i]], "var", value=exposure.by[1,'1'], cen=.cen, model.link='log')$coefficients) # only interested in length, uses last cb from S1 but that is irrelevant

coeflag1 <- coeflag2<- coeflag10 <- coeflag25 <- coeflag50 <- coeflag75 <- coeflag90<- coeflag975 <- coeflag99 <- matrix(NA, length.ds.stratum, cvlag.length, dimnames=list(ds.stratum)) # not really sure why columns is 3. I though 8 was for columns for lag 0 + each lag day. CHANGED WHEN I INCREASED PARAMETES FROM edf TO edf + 1
cpmodel <- vcovlag1 <- vcovlag2 <- vcovlag10 <- vcovlag25 <- vcovlag50 <- vcovlag75 <- vcovlag90 <- vcovlag975 <- vcovlag99 <- vector("list", length.ds.stratum)
names(cpmodel)<-names(vcovlag1)<-names(vcovlag2)<-names(vcovlag10)<-names(vcovlag25)<-names(vcovlag50)<-names(vcovlag75)<-names(vcovlag90)<-names(vcovlag975)<-names(vcovlag99) <- ds.stratum

# Run model for each by.vars with re-centered values to obtain lag-reponse relationship at the 1st and 99th percentile

for(i in ds.stratum) {
  # print(paste('Stage 2 individual models:',i))
  .ds <- daily.ds[stratum==i] # dataset for each unique value
  .name <- unique(.ds$stratum) # names with all of a, b and c together
  
  .outcome <- .ds[,get(outcome.var)] # outcome
  
  # Calculate crossbasis (comments in Gasparrini and Martinez-Solanas code state it is centred, although this isn't the case). Model appears same, with changes made for crosspred and crossreduce with cen 
  .cb <- crossbasis(temps[[i]], argvar=list(fun=espline, knots=quantile(temps[[i]], eknots, na.rm=T)), lag=lmax, arglag=li_arglag) # lag=lag2
  
  # Model using crossbasis
  .no.years <- no.years[i,]
  
  # Predictions and reduction to lag-response at 1st (extreme cold) and 99th (extreme hot) percentiles with new centering (changes coef and vcov)
  if (is.numeric(cenpercountry)) {
    .cen <- quantile(temps[[i]], cenpercountry/100, na.rm=T) # quantile
  }  else {
    .cen <- exposure.mean[i] # mean
  }
  
  # Predictions and reduction to lag-response at set percentiles. Requires centering required as it changes coef-vcov
  .perc <-  exposure.by[i,] # all relevant percentiles per i
  # 1st percentile
  redlag1 <- crossreduce(.cb, model=model[[i]], "var", value=.perc['1'], cen=.cen, model.link='log')
  coeflag1[i,] <- coef(redlag1)
  vcovlag1[[i]] <- vcov(redlag1)
  # 2.5th percentile
  redlag2 <- crossreduce(.cb, model=model[[i]], "var", value=.perc['2.5'], cen=.cen, model.link='log')
  coeflag2[i,] <- coef(redlag2)
  vcovlag2[[i]] <- vcov(redlag2)
  # 10th percentile
  redlag10 <- crossreduce(.cb, model=model[[i]], "var", value=.perc['10'], cen=.cen, model.link='log')
  coeflag10[i,] <- coef(redlag10)
  vcovlag10[[i]] <- vcov(redlag10)
  # 25th percentile
  redlag25 <- crossreduce(.cb, model=model[[i]], "var", value=.perc['25'], cen=.cen, model.link='log')
  coeflag25[i,] <- coef(redlag25)
  vcovlag25[[i]] <- vcov(redlag25)
  # 50th percentile
  redlag50 <- crossreduce(.cb, model=model[[i]], "var", value=.perc['50'], cen=.cen, model.link='log')
  coeflag50[i,] <- coef(redlag50)
  vcovlag50[[i]] <- vcov(redlag50)
  # 75th percentile
  redlag75 <- crossreduce(.cb, model=model[[i]], "var", value=.perc['75'], cen=.cen, model.link='log')
  coeflag75[i,] <- coef(redlag75)
  vcovlag75[[i]] <- vcov(redlag75)
  # 90th percentile
  redlag90 <- crossreduce(.cb, model=model[[i]], "var", value=.perc['90'], cen=.cen, model.link='log')
  coeflag90[i,] <- coef(redlag90)
  vcovlag90[[i]] <- vcov(redlag90)
  # 97.5th percentile
  redlag975 <- crossreduce(.cb, model=model[[i]], "var", value=.perc['97.5'], cen=.cen, model.link='log') 
  coeflag975[i,] <- coef(redlag975)
  vcovlag975[[i]] <- vcov(redlag975)
  # 99th percentile
  redlag99 <- crossreduce(.cb, model=model[[i]], "var", value=.perc['99'], cen=.cen, model.link='log') 
  coeflag99[i,] <- coef(redlag99)
  vcovlag99[[i]] <- vcov(redlag99)
}



# Run meta-analysis with lag models 
mvlag1 <- mixmeta(coeflag1~1, vcovlag1, data=list(ds.stratum)) 
mvlag2 <- mixmeta(coeflag2~1, vcovlag2, data=list(ds.stratum))
mvlag10 <- mixmeta(coeflag10~1, vcovlag10, data=list(ds.stratum)) 
mvlag25 <- mixmeta(coeflag25~1, vcovlag25, data=list(ds.stratum))
mvlag50 <- mixmeta(coeflag50~1, vcovlag50, data=list(ds.stratum))
mvlag75 <- mixmeta(coeflag75~1, vcovlag75, data=list(ds.stratum)) 
mvlag90 <- mixmeta(coeflag90~1, vcovlag90, data=list(ds.stratum))
mvlag975 <- mixmeta(coeflag975~1, vcovlag975, data=list(ds.stratum)) 
mvlag99 <- mixmeta(coeflag99~1, vcovlag99, data=list(ds.stratum))
# for(p in c(1,2.5,10,25,50,75,90,97.5,99))

# Predict pooled coefficients
mvpredlag1 <- predict(mvlag1,datanew,vcov=T,format="list")
mvpredlag2 <- predict(mvlag2,datanew,vcov=T,format="list")
mvpredlag10 <- predict(mvlag10,datanew,vcov=T,format="list")
mvpredlag25 <- predict(mvlag25,datanew,vcov=T,format="list")
mvpredlag50 <- predict(mvlag50,datanew,vcov=T,format="list")
mvpredlag75 <- predict(mvlag75,datanew,vcov=T,format="list")
mvpredlag90 <- predict(mvlag90,datanew,vcov=T,format="list")
mvpredlag975 <- predict(mvlag975,datanew,vcov=T,format="list")
mvpredlag99 <- predict(mvlag99,datanew,vcov=T,format="list")

# Obtain predictions for lag 0 to "lmax"
blag <- do.call(onebasis,c(list(x=seq(0,lmax)),attr(.cb,"arglag"))) # uses CB attributes (same for all by.vars)

# Predict pooled lag-response associations for heat and cold
cplag1 <- crosspred(blag,coef=mvpredlag1$fit,vcov=mvpredlag1$vcov, model.link="log", at=0:lmax)
cplag2 <- crosspred(blag,coef=mvpredlag2$fit,vcov=mvpredlag2$vcov, model.link="log", at=0:lmax)
cplag10 <- crosspred(blag,coef=mvpredlag10$fit,vcov=mvpredlag10$vcov, model.link="log", at=0:lmax)
cplag25 <- crosspred(blag,coef=mvpredlag25$fit,vcov=mvpredlag25$vcov, model.link="log", at=0:lmax)
cplag50 <- crosspred(blag,coef=mvpredlag50$fit,vcov=mvpredlag50$vcov, model.link="log", at=0:lmax)
cplag75 <- crosspred(blag,coef=mvpredlag75$fit,vcov=mvpredlag75$vcov, model.link="log", at=0:lmax)
cplag90 <- crosspred(blag,coef=mvpredlag90$fit,vcov=mvpredlag90$vcov, model.link="log", at=0:lmax)
cplag975 <- crosspred(blag,coef=mvpredlag975$fit,vcov=mvpredlag975$vcov, model.link="log", at=0:lmax)
cplag99 <- crosspred(blag,coef=mvpredlag99$fit,vcov=mvpredlag99$vcov, model.link="log", at=0:lmax)



################################################################################
# Main results: RRs (95% CI) of set exposure percentiles compared to minimum occupational-injuries percentile (MOIP)
################################################################################

# Short percentile
predprer.short.rep <- 1:length(predper.short) # number of percentiles
results.short <- matrix(NA, nrow=length(predper.short), ncol=4) # create empty matrix of results
results.short[,1] <- predper.short # 1st column is percentile values
results.short[predprer.short.rep,2:4] <- c(cp$allRRfit[predper %in% predper.short], # column 2 is RR corresponding to each percentile (row), 
                                           cp$allRRlow[predper %in% predper.short], # column 3 is RRlow
                                           cp$allRRhigh[predper %in% predper.short]) # column 4 is RRhigh
results.short <- round(results.short, digits=3)

# All percentiles
predprer.rep <- 1:length(predper) # number of percentiles
results <- matrix(NA, nrow=length(predper), ncol=4) # create empty matrix of results
results[,1] <- predper # 1st column is percentile values
results[predprer.rep,2:4] <- c(cp$allRRfit,cp$allRRlow,cp$allRRhigh) # column 2/3/4 is RR/RRlow/RRhigh corresponding to each percentile (row), 
results <- round(results, digits=3)


# Export results
colnames(results) <- c('Percentile','RR','RRlow','RRhigh') # column names
write.csv(results, file=paste0(s2results,'RR (95% CI).csv'), na='', row.names=F) # create csv file



###############################################################################
# COMPUTE ATTRIBUTABLE INJURIES FOR EACH CITY, WITH EMPIRICAL CI ESTIMATED USING RE-CENTERED BASES
################################################################################

# Create vectors to store total injuries, accounting for missing
totclaims <- rep(NA,length.ds.stratum)
names(totclaims) <- ds.stratum

# Objects to store attributable injuries (simulations)
sim.names <- c("Total","Cold","Heat","Extreme cold","Moderate cold","Moderate heat","Extreme heat")
length.sim <- length(sim.names)
matsim <- matrix(NA, length.ds.stratum, length.sim, dimnames=list(ds.stratum, sim.names)) # matrix: attributable injuries
arraysim <- array(NA, dim=c(length.ds.stratum, length.sim, nsim), dimnames=list(ds.stratum, sim.names)) # array: attributable injuries CI

# Run loop
gc()
time <- proc.time()[3]
for(i in ds.stratum){
  print(paste('Attributable risk:',i))
  
  .ds <- daily.ds[stratum==i] # dataset for each unique value
  .outcome <- .ds[,get(outcome.var)] # .define outcome
  
  # Predictions and reduction to lag-response at 1st (extreme cold) and 99th (extreme hot) percentiles with new centering (changes coef and vcov)
  if (is.numeric(cenpercountry)) {
    .cen <- quantile(temps[[i]], cenpercountry/100, na.rm=T)
  } else {
    .cen <- exposure.mean[i] # mean
  }
  
  # Derive cross-basis. NB: centering point different than original choice of 75th (i see no difference)
  .cb <- crossbasis(temps[[i]], argvar=list(fun=espline, knots=quantile(temps[[i]], eknots, na.rm=T)), lag=lmax, arglag=li_arglag)
  
  # Compute attributable numbers with reduced coefficients, based on centering value
  .perc <-  exposure.by[i,] # all percentiles and related average exposure
  .blups <- blups[[which(ds.stratum==i)]]
  
  matsim[i,"Total"] <- attrdl(temps[[i]], .cb, .outcome, coef=.blups$blup, vcov=.blups$vcov, type="an", dir='forw', cen=.cen)
  matsim[i,"Cold"] <- attrdl(temps[[i]], .cb, .outcome, coef=.blups$blup, vcov=.blups$vcov, type="an", dir='forw', cen=.cen,
                             range=c(-100, .cen))
  matsim[i,"Heat"] <- attrdl(temps[[i]], .cb, .outcome, coef=.blups$blup, vcov=.blups$vcov, type="an", dir='forw', cen=.cen,
                             range=c(.cen, 100))
  matsim[i,"Extreme cold"] <- attrdl(temps[[i]], .cb, .outcome, coef=.blups$blup, vcov=.blups$vcov, type="an", dir='forw', cen=.cen,
                                     range=c(-100,.perc["2.5"]))
  matsim[i,"Moderate cold"] <- attrdl(temps[[i]], .cb, .outcome, coef=.blups$blup, vcov=.blups$vcov, type="an", dir='forw', cen=.cen,
                                         range=c(.perc["2.5"],.cen))
  
  matsim[i,"Moderate heat"] <- attrdl(temps[[i]], .cb, .outcome, coef=.blups$blup, vcov=.blups$vcov, type="an", dir='forw', cen=.cen,
                                         range=c(.cen,.perc["97.5"]))
  matsim[i,"Extreme heat"] <- attrdl(temps[[i]], .cb, .outcome, coef=.blups$blup, vcov=.blups$vcov, type="an", dir='forw', cen=.cen,
                                     range=c(.perc["97.5"],100))
  
  # Compute empirical occurrences of attributable numbers used to derive CIs
  # Based on overall curve, patterns seems consistent from mean till generally about 10th/90th percentile, and then a sharper change at 2.5th/97.5th percentile
  arraysim[i,"Total",] <- attrdl(temps[[i]], .cb, .outcome, coef=.blups$blup, vcov=.blups$vcov, type="an", dir='forw', cen=.cen,
                                 sim=T, nsim=nsim)
  arraysim[i,"Cold",] <- attrdl(temps[[i]], .cb, .outcome, coef=.blups$blup, vcov=.blups$vcov, type="an", dir='forw', cen=.cen,
                                range=c(-100, .cen), sim=T, nsim=nsim)
  arraysim[i,"Heat",] <-  attrdl(temps[[i]], .cb, .outcome, coef=.blups$blup, vcov=.blups$vcov, type="an", dir='forw', cen=.cen,
                                 range=c(.cen,100), sim=T, nsim=nsim)
  arraysim[i,"Extreme cold",] <- attrdl(temps[[i]], .cb, .outcome, coef=.blups$blup, vcov=.blups$vcov, type="an", dir='forw', cen=.cen,
                                        range=c(-100,.perc["2.5"]), sim=T, nsim=nsim)
  arraysim[i,"Moderate cold",] <- attrdl(temps[[i]], .cb, .outcome, coef=.blups$blup, vcov=.blups$vcov, type="an", dir='forw', cen=.cen,
                                            range=c(.perc["2.5"],.cen), sim=T, nsim=nsim)
  arraysim[i,"Moderate heat",] <- attrdl(temps[[i]], .cb, .outcome, coef=.blups$blup, vcov=.blups$vcov, type="an", dir='forw', cen=.cen,
                                            range=c(.cen,.perc["97.5"]), sim=T, nsim=nsim)
  arraysim[i,"Extreme heat",] <- attrdl(temps[[i]], .cb, .outcome, coef=.blups$blup, vcov=.blups$vcov, type="an", dir='forw', cen=.cen,
                                        range=c(.perc["97.5"], 100), sim=T, nsim=nsim)

  totclaims[i] <- sum(.outcome, na.rm=T) # store total injuries (account for missing)
}
proc.time()[3]-time


### Attributable numbers ###
# City-specific attributable numbers. Already have row names
ancitylow <- apply(arraysim,c(1,2),quantile,0.025)
ancityhigh <- apply(arraysim,c(1,2),quantile,0.975)
ancity <- matrix(paste0(round.fn(matsim),' (',round.fn(ancitylow),' to ',round.fn(ancityhigh),')'), ncol=length.sim) # insert comma every 3 digits, no rounding
rownames(ancity) <- ds.stratum

# Total attributable numbers
antot <- colSums(matsim) # sum through strata
antotlow <- apply(apply(arraysim,c(2,3),sum),1,quantile,0.025)
antothigh <- apply(apply(arraysim,c(2,3),sum),1,quantile,0.975)
antotal <- matrix(paste0(round.fn(antot),' (',round.fn(antotlow),' to ',round.fn(antothigh),')'), ncol=length.sim)
rownames(antotal) <- 'Total'
colnames(ancity) <- sim.names # colnames

# City and IO attributable numbers
an_city <- .an_citylow <- .an_cityhigh <- matrix(NA, length.ds.stratum/2, length.sim, dimnames=list(sort(unique(daily.ds[,City])))) # matrix of NAs, AN for cities
colnames(an_city) <- sim.names # colnames
.an_city <- rowsum(matsim, as.integer(gl(nrow(matsim), 2, nrow(matsim)))) # sum every 2 rows for city AN
for(i in 1:(length.ds.stratum/2)) { # for each city
  .an_cityv <- colSums(arraysim[(i*2-1):(i*2),,]) # sum results within same city, still keeping all repetitions
  .an_citylow[i,] <- apply(.an_cityv,1,quantile,0.025) # lower quantile
  .an_cityhigh[i,] <- apply(.an_cityv,1,quantile,0.975) # upper quantile
}
an_city <- matrix(paste0(round.fn(.an_city),' (',round.fn(.an_citylow),' to ',round.fn(.an_cityhigh),')'), ncol=length.sim) # combine for presentation
rownames(an_city) <- rownames(.an_citylow) <- rownames(.an_cityhigh) <- sort(unique(daily.ds[,City])) # rownames

an_io <- .an_iolow <- .an_iohigh <- matrix(NA, length.ds.stratum/length(unique(daily.ds[,City])), length.sim, dimnames=list(sort(unique(daily.ds[,outin])))) # matrix of NAs, AN for indoor/outdoor
.an_io <- rowsum(matsim, rep(1:2,length.ds.stratum/2)) # sum every 2 rows for city AN
for(i in 1:2) { # for each indoor/outdoor combination
  .an_iov <- colSums(arraysim[seq(0,length.ds.stratum-2,by=2)+i,,]) # sum results for each indoor or outdoor, still keeping all repetitions
  .an_iolow[i,] <- apply(.an_iov,1,quantile,0.025) # lower quantile
  .an_iohigh[i,] <- apply(.an_iov,1,quantile,0.975) # upper quantile
}

an_io <- matrix(paste0(round.fn(.an_io),' (',round.fn(.an_iolow),' to ',round.fn(.an_iohigh),')'), ncol=length.sim) # combine for presentation
rownames(an_io) <- rownames(.an_iolow) <- rownames(.an_iohigh) <- sort(unique(daily.ds[,outin])) # rownames
colnames(an_io) <- colnames(antotal) <- sim.names # colnames
ancountry <- rbind(antotal, an_io, an_city, ancity) # all ANs
  

# Export attributable numbers
write.csv(ancountry, file=paste0(s2results, 'Attributable numbers.csv'), na='', row.names=T) # create csv file


### Attributable fractions ###

totclaimtot <- sum(totclaims) # total injuries

# City-specific
afcit <- matsim/totclaims*100
afcitylow <- ancitylow/totclaims*100
afcityhigh <- ancityhigh/totclaims*100
afcity <- matrix(paste0(round.fn(afcit, af.round),' (',round.fn(afcitylow, af.round),' to ',round.fn(afcityhigh, af.round),')'), ncol=length.sim)
rownames(afcity) <- ds.stratum

# Total
aftot <- antot/totclaimtot*100
aftotlow <- antotlow/totclaimtot*100
aftothigh <- antothigh/totclaimtot*100
aftotal <- matrix(paste0(round.fn(aftot, af.round),' (',round.fn(aftotlow, af.round),' to ',round.fn(aftothigh, af.round),')'), ncol=length.sim)
rownames(aftotal) <- 'Total'
colnames(afcity) <- colnames(aftotal) <- sim.names # colnames

# City and IO attributable numbers
totclaimscity <- rollapply(totclaims, 2, by=2, sum)
afcitcity <- .an_city/totclaimscity*100
afcitycitylow <- .an_citylow/totclaimscity*100
afcitycityhigh <- .an_cityhigh/totclaimscity*100
afcitycity <- matrix(paste0(round.fn(afcitcity, af.round),' (',round.fn(afcitycitylow, af.round),' to ',round.fn(afcitycityhigh, af.round),')'), ncol=length.sim)
rownames(afcitycity) <- sort(unique(daily.ds$City)) # rownames

totclaimsio <- c(sum(totclaims[c(TRUE, FALSE)]), sum(totclaims[c(FALSE, TRUE)]))
afcitio <- .an_io/totclaimsio*100
afcityiolow <- .an_iolow/totclaimsio*100
afcityiohigh <- .an_iohigh/totclaimsio*100
afcityio <- matrix(paste0(round.fn(afcitio, af.round),' (',round.fn(afcityiolow, af.round),' to ',round.fn(afcityiohigh, af.round),')'), ncol=length.sim)
rownames(afcityio) <- sort(unique(daily.ds[,outin])) # rownames
colnames(afcitycity) <- colnames(afcityio) <- sim.names # colnames
  
afcountry <- rbind(aftotal, afcityio, afcitycity, afcity) # all AFs


# Export attributable fractions
write.csv(afcountry, file=paste0(s2results, 'Attributable fractions.csv'), na='', row.names=T) # create csv file


################################################################################
# Plots: overall relationship
################################################################################

indlab <- predper %in% c(0,2.5,10,25,50,75,90,75,97.5,100) # requires predper, so must be in S2

# Plot: Overall cumulative exposure-response association (all exposure values, lag reduced)
oer.yaxis <- seq(0.9,1.5,by=0.1)

png(file = paste0(s2results, 'Overall e-r.png'), res=gdpi, width=glength, height=glength) # plot location. File name based on heat metric
par(mar=c(4.1,3,1.6,0)) # inner graph margins, as much whitespace removed as possible

plot(cp,lwd=2,col="white",axes=F, ylim=c(min(oer.yaxis), max(oer.yaxis)), xlab='', ylab='') # str_remove_all(paste0('Percent change in ',tolower(outcome.var)," (%)"),'\\(000s\\)') # ylim=c(floor(min(cp[["allRRfit"]])*10)/10-0.05, ceiling(max(cp[["allRRfit"]])*10)/10+0.05)
ind1 <- cp$predvar<=cp$cen
ind2 <- cp$predvar>=cp$cen
lines(cp$predvar[ind1],cp$allRRfit[ind1],col=4,lwd=2)
lines(cp$predvar[ind2],cp$allRRfit[ind2],col=2,lwd=2)
axis(1, at=exposure.country.mean[indlab], labels=predper[indlab], cex.axis=0.8, gap.axis=0.1, mgp=c(2.5,0.3,0))
mtext('Percentile', 1, at=7.8, outer=F, cex=0.8, adj=1)
axis(1, at=seq(ceiling(min(cp$predvar)),max(cp$predvar),by=1), cex.axis=0.8, line=xline2, mgp=c(2.5,0.3,0)) # round min inwards with ceiling, otherwise starting point is decimals. From there on, proceed by 2
mtext('Value', 1, at=7.8, outer=F, cex=0.8, line=xline2, adj=1)
axis(2, at=oer.yaxis, labels=(oer.yaxis-1)*100, cex.axis=0.8, las=1, mgp=c(2.5,0.8,0))
title(xlab=paste0(exposure.var," (°C)"), line=3.1) # xtitle, moved away from graph
abline(v=c(exposure.country.mean[c("2.5","97.5")]), lty=c(3,3)) # choosing these in part because previous studies highlighted 1&99, in part because these corresponding to extreme, and in part because there's extreme heat results
title(main=str_replace(str_remove_all(outcome.var,'\\(000s\\)'),'Number of OIIs','Number of injuries and illnesses'), line=nline) # title, move closer to graph
title(ylab='Percent change (%)', line=2.1) # title, move closer to graph
abline(v=cp$cen,lty=2) # centre line

dev.off() # Save image + clear settings


# Plot: overall cumulative lag-response associations
.cexaxis <- 1

lag.yaxis <- seq(0.98,1.02,by=0.01)
.axis %<a-% {axis(2,at=lag.yaxis, labels=(lag.yaxis-1)*100, cex.axis=.cexaxis, las=1)} # shared axis as percentage change, tick labels are horizontal
lag.ylim <- c(min(lag.yaxis),max(lag.yaxis))

png(file = paste0(s2results, 'Overall lag.png'), res=gdpi, width=glength.3by3, height=glength.3by3) # plot location. File name based on heat metric
par(mfrow=c(3,3), mar=c(3,2.5,1,1), oma=c(0,0,0,0), mgp=c(1.5,0.5,0)) # 3*3, space between borders and text

plot(cplag1,ylab="Percent change (%)",xlab="Lag (days)",lwd=2,cex.axis=.cexaxis,
     ylim=lag.ylim,yaxt='n',col="steelblue3", ci.arg=list(density=20,col="steelblue3"))
title(main='1st percentile', line=tline) # move title closer to graph
.axis
plot(cplag2,ylab="Percent change (%)",xlab="Lag (days)",lwd=2,cex.axis=.cexaxis,
     ylim=lag.ylim,yaxt='n',col="steelblue2", ci.arg=list(density=20,col="steelblue2"))
title(main='2.5th percentile', line=tline)
.axis
plot(cplag10,ylab="Percent change (%)",xlab="Lag (days)",lwd=2,cex.axis=.cexaxis,
     ylim=lag.ylim,yaxt='n',col="steelblue1", ci.arg=list(density=20,col="steelblue1"))
title(main='10th percentile', line=tline)
.axis
plot(cplag25,ylab="Percent change (%)",xlab="Lag (days)",lwd=2,cex.axis=.cexaxis,
     ylim=lag.ylim,yaxt='n',col="skyblue1", ci.arg=list(density=20,col="skyblue1"))
title(main='25th percentile', line=tline)
.axis
plot(cplag50,ylab="Percent change (%)",xlab="Lag (days)",lwd=2,cex.axis=.cexaxis,
     ylim=lag.ylim,yaxt='n',col="green", ci.arg=list(density=20,col="green")) # practically identical to mean, thus no relationship apparent, but included for comparison and c(3,3) instead of c(2,4)
title(main='Median', line=tline)
.axis
plot(cplag75,ylab="Percent change (%)",xlab="Lag (days)",lwd=2,cex.axis=.cexaxis,
     ylim=lag.ylim,yaxt='n',col="orangered", ci.arg=list(density=20,col="orangered"))
title(main='75th percentile', line=tline)
.axis
plot(cplag90,ylab="Percent change (%)",xlab="Lag (days)",lwd=2,cex.axis=.cexaxis,
     ylim=lag.ylim,yaxt='n',col="red1", ci.arg=list(density=20,col="red1"))
title(main='90th percentile', line=tline)
.axis
plot(cplag975,ylab="Percent change (%)",xlab="Lag (days)",lwd=2,cex.axis=.cexaxis,
     ylim=lag.ylim,yaxt='n',col="red3", ci.arg=list(density=20,col="red3"))
title(main='97.5th percentile', line=tline)
.axis
plot(cplag99,ylab="Percent change (%)",xlab="Lag (days)",lwd=2,cex.axis=.cexaxis,
     ylim=lag.ylim,yaxt='n',col="red4", ci.arg=list(density=20,col="red4"))
title(main='99th percentile', line=tline)
.axis

dev.off() # Save image + clear settings




################################################################################
#### Plots: overall cumulative exposure-response associations by by-variables
################################################################################

.exposure.rr <- list()
s2locations <- paste0(s2results,'Location level/')
dir.create(paste0(s2results,'Location level')) # Warning if exists (doesn't replace)

# Loop over each strata. Plots are individualised
for(i in ds.stratum) {
  .ds <- daily.ds[stratum==i] # dataset for each unique value
  .name <- unique(.ds$stratum) # names with all of a, b and c together
  
  if (is.numeric(cenpercountry)) {
    .cen <- quantile(temps[[i]], cenpercountry/100, na.rm=T)
  } else {
    .cen <- exposure.mean[i] # mean
  }
  
  .argvar <- list(x=temps[[i]], fun=espline, knots=quantile(temps[[i]], eknots, na.rm=T))
  
  .bvar <- do.call(onebasis, .argvar)
  .blups <- blups[[which(ds.stratum==i)]] # get correct blups
  .pred2 <- crosspred(.bvar, coef=.blups$blup, vcov=.blups$vcov, model.link="log", by=0.1, cen=.cen) #  overall cumulative exposure-response relationship, for each percentile. Specific humidity fails at this step because of Error in seq.default(from = min(pretty), to = to, by = by) : 'from' must be a finite number
  
  .exposure.rr[[i]] <- data.frame(matrix(nrow=length(.pred2$predvar), ncol=4))
  colnames(.exposure.rr[[i]]) <- c("temp", "RRfit", "RRlow", "RRhigh")
  .exposure.rr[[i]]$temp <- .pred2$predvar
  .exposure.rr[[i]]$RRfit <- (.pred2$allRRfit-1)*100
  .exposure.rr[[i]]$RRlow <- (.pred2$allRRlow-1)*100
  .exposure.rr[[i]]$RRhigh <- (.pred2$allRRhigh-1)*100
  
  # Values above or below centering value (for red/blue colour)
  .ind1 <- .pred2$predvar<=.cen
  .ind2 <- .pred2$predvar>=.cen
  
  # Overall cumulative exposure-response relationship plot
  # Could technically extract specific RRs values from RRfit/RRlow/RRhigh, but AF more useful as it covers a range of RRs
  png(file = paste0(s2locations, .name,', Overall e-r.png'), res=gdpi, width=glength.3by3, height=glength.3by3) # plot location. File name based on heat metric
  plot(.exposure.rr[[i]]$temp, .exposure.rr[[i]]$RRfit, type="n", ylim=rryaxis, lwd=2, col="white",
       main=i, ylab="Percent change (%)", xlab=paste(exposure.var,'(°C)'),
       cex.main=cexmain, cex.lab=cexlab, cex.axis=cexaxis, lab=c(6,5,7)) # main plot
  .erplot %<a-% { # save main plot additions
    polygon(c(.exposure.rr[[i]]$temp,rev(.exposure.rr[[i]]$temp)),c(.exposure.rr[[i]]$RRlow,rev(.exposure.rr[[i]]$RRhigh)), col="grey89", border=F) # 95% CI envelope
    lines(.exposure.rr[[i]]$temp[.ind1],.exposure.rr[[i]]$RRfit[.ind1],col=4,lwd=2); # cold, left of cen
    lines(.exposure.rr[[i]]$temp[.ind2],.exposure.rr[[i]]$RRfit[.ind2],col=2,lwd=2); # hot, right of cen
    abline(h=0) # horizontal line
  }
  .lines %<a-% { #  vertical lines
    abline(v=.cen,lty=2); # centering value
    abline(v=c(exposure[i,c("2.5%","97.5%")]),lty=3) # extreme percentiles
  }
  .erplot # insert main plot additions
  .lines # insert lines
  dev.off() # Save image + clear settings
  
  # Plot + histogram (slice plot)
  png(file = paste0(s2locations, 'z', .name,', Overall e-r.png'), res=gdpi, width=glength.3by3, height=glength.3by3) # plot location. File name based on heat metric
  plot(.exposure.rr[[i]]$temp, .exposure.rr[[i]]$RRfit, type="n", ylim=c(rryaxis[1]-50,rryaxis[2]), lwd=2, col="white", ylab="Percent change (%)", xlab='',
       mgp=c(1.3,0.4,0), cex.lab=cexlab*0.7, cex.axis=cexaxis*0.7, yaxt='n', axes=F) # main plot, but no axes # cexlab*.75, cexaxis*.75
  .erplot
  title(main=i, line=tline*2, cex.main=cexmain*0.7) # closer to plot
  title(xlab=paste(exposure.var,'(°C)'), line=1.1, cex.lab=cexlab*0.7) # xtitle, moved away from graph
  axis(1, col.axis="black", cex.axis=cexaxis*0.7, tck=tck.length, mgp=c(2.5,0.3,0)) # x-axis
  axis(2, at=seq(rryaxis[1],rryaxis[2],by=20), col.axis="black", cex.axis=cexaxis*0.7, las=1, tck=tck.length, mgp=c(2.5,0.4,0)) # y-axis on L
  
  .breaks <- c(min(temps[[i]],na.rm=T)-1, seq(.pred2$predvar[1], .pred2$predvar[length(.pred2$predvar)],length=30), max(temps[[i]],na.rm=T)+1)
  hist <- hist(temps[[i]],breaks=.breaks,plot=F) # histogram
  hist$density <- hist$density/max(hist$density)*0.7 # density value
  prop <- max(hist$density)/max(hist$counts) # convert proportion so that it fits on same scale as other plot
  counts <- pretty(hist$count,3) # counts per bin
  par(new=TRUE) # add to graph
  plot(hist,ylim=c(0,max(hist$density)*2.0),axes=F,ann=F,col=grey(0.95),freq=F) # plot histogram
  axis(4, at=counts*prop, labels=counts, cex.axis=0.7, las=1, tck=tck.length, mgp=c(2.5,0.4,0)) # axis on R indicating bin counts
  .lines
  dev.off() # Save image + clear settings
}

# Loop over each strata. Plots are combined
png(file = paste0(s2locations, 'zOverall e-rs.png'), res=350, width=1700, height=2400) # plot location # paste0('/Users/MatthewBorg/', 'zOverall e-rs.png')
layout(matrix(c(1:14,0),ncol=3,byrow=T))
par(mfrow=c(5,3), mar=c(2.1,2.2,1,1.5), oma=c(0,0,0,0), mgp=c(1.5,0.4,0), las=1) # 4*4, space between borders and text

for(i in ds.stratum) {
  .ds <- daily.ds[stratum==i] # dataset for each unique value
  .name <- unique(.ds$stratum) # names with all of a, b and c together
  # .name <- str_replace(str_replace(str_replace(str_replace(str_replace(str_replace(str_replace(str_replace(str_replace(.name, ',',''), 'Adelaide','Ade'), 'Brisbane','Bri'), 'Canberra','Can'), 'Darwin','Dar'), 'Hobart','Hob'), 'Melbourne','Mel'), 'Perth','Per'), 'Sydney','Syd')
  
  if (is.numeric(cenpercountry)) {
    .cen <- quantile(temps[[i]], cenpercountry/100, na.rm=T)
  } else {
    .cen <- exposure.mean[i] # mean
  }
  
  .argvar <- list(x=temps[[i]], fun=espline, knots=quantile(temps[[i]], eknots, na.rm=T))
  
  .bvar <- do.call(onebasis, .argvar)
  .blups <- blups[[which(ds.stratum==i)]] # get correct blups
  .pred2 <- crosspred(.bvar, coef=.blups$blup, vcov=.blups$vcov, model.link="log", by=0.1, cen=.cen) #  overall cumulative exposure-response relationship, for each percentile. Specific humidity fails at this step because of Error in seq.default(from = min(pretty), to = to, by = by) : 'from' must be a finite number
  
  # Values above or below centering value (for red/blue colour)
  .ind1 <- .pred2$predvar<=.cen
  .ind2 <- .pred2$predvar>=.cen
  
  # Overall cumulative exposure-response relationship plot
  .erplot %<a-% { # save main plot additions
    polygon(c(.exposure.rr[[i]]$temp,rev(.exposure.rr[[i]]$temp)),c(.exposure.rr[[i]]$RRlow,rev(.exposure.rr[[i]]$RRhigh)), col="grey89", border=F) # 95% CI envelope
    lines(.exposure.rr[[i]]$temp[.ind1],.exposure.rr[[i]]$RRfit[.ind1],col=4,lwd=2); # cold, left of cen
    lines(.exposure.rr[[i]]$temp[.ind2],.exposure.rr[[i]]$RRfit[.ind2],col=2,lwd=2); # hot, right of cen
    abline(h=0) # horizontal line
  }
  .lines %<a-% { #  vertical lines
    abline(v=.cen,lty=2); # centering value
    abline(v=c(exposure[i,c("2.5%","97.5%")]),lty=3) # extreme percentiles
  }
  
  # Plot + histogram (slice plot)
  plot(.exposure.rr[[i]]$temp, .exposure.rr[[i]]$RRfit, type="n", ylim=c(rryaxis[1]-50,rryaxis[2]), lwd=2, col="white", ylab="Percent change (%)", xlab='',
       mgp=c(1.3,0.4,0), cex.lab=cexlab*0.7, cex.axis=cexaxis*0.7, yaxt='n', axes=F) # main plot, but no axes # cexlab*.75, cexaxis*.75
  .erplot
  title(main=i, line=tline*2, cex.main=cexmain*0.7) # closer to plot
  title(xlab=paste(exposure.var,'(°C)'), line=1.1, cex.lab=cexlab*0.7) # xtitle, moved away from graph
  axis(1, col.axis="black", cex.axis=cexaxis*0.7, tck=tck.length, mgp=c(2.5,0.3,0)) # x-axis
  axis(2, at=seq(rryaxis[1],rryaxis[2],by=20), col.axis="black", cex.axis=cexaxis*0.7, las=1, tck=tck.length, mgp=c(2.5,0.4,0)) # y-axis on L
  
  .breaks <- c(min(temps[[i]],na.rm=T)-1, seq(.pred2$predvar[1], .pred2$predvar[length(.pred2$predvar)],length=30), max(temps[[i]],na.rm=T)+1)
  hist <- hist(temps[[i]],breaks=.breaks,plot=F) # histogram
  hist$density <- hist$density/max(hist$density)*0.7 # density value
  prop <- max(hist$density)/max(hist$counts) # convert proportion so that it fits on same scale as other plot
  counts <- pretty(hist$count,3) # counts per bin
  par(new=TRUE) # add to graph
  plot(hist,ylim=c(0,max(hist$density)*2.0),axes=F,ann=F,col=grey(0.95),freq=F) # plot histogram
  axis(4, at=counts*prop, labels=counts, cex.axis=0.7, las=1, tck=tck.length, mgp=c(2.5,0.4,0)) # axis on R indicating bin counts
  .lines
  # dev.off() # Save image + clear settings
}
dev.off()

# Export exposure RRs
exposure.rr <- do.call(rbind,.exposure.rr) # convert to df
exposure.rr <- cbind(gsub("\\..*","",rownames(exposure.rr)), exposure.rr) # combine with rownames formatted to remove .
rownames(exposure.rr) <- NULL
colnames(exposure.rr) <- c('City',exposure.var,'RRfit','RRlow','RRhigh')
write.csv(exposure.rr, file=paste0(s2locations, 'Exposure values RR.csv'), na='', row.names=F) # create csv file



################################################################################
# CITY AND OUTORIN OVERALL E-R CURVES
################################################################################

### Full meta-analysis with random-predictor. BLUPs will be conditional expectations given the random effects (Sera et al 2019)

  for(j in by.vars3) { # City and outorin
    # j <- 'City'
    mixv_j <- mixmeta(coef~1, vcov, random=~1|by.vars.ds[,get(j)], data=by.vars.ds) # BLUP will be conditionally dependent on random predictor. Usage of predictor like this identical to City (city results identical for both indoor and outdoor) or outorin (outorin results identical for all cities)
    if(j=='City') {
      jstrat_no <- seq(0,length.ds.stratum-2,by=2)+1 # every odd number
      ds.stratumj <- word(unique(ds.stratum))[jstrat_no] # new ds.stratum (technically
    } else {
      jstrat_no <- c(1:2)
      ds.stratumj <- word(unique(ds.stratum), start=-1L)[jstrat_no]
    }
    blups_j <-  blup(mixv_j, vcov=T)[jstrat_no] # BLUPs for each j
    
    for(l in jstrat_no) { # repeat for each strata (city or out/in)
      m <- which(jstrat_no==l) # ordering of l, instead of l itself, which is needed for extracting blups
      .name <- ds.stratumj[m] # get name for plot
      mvpred_jl <- list('fit'=blups_j[[m]][["blup"]], 'vcov'=blups_j[[m]][["vcov"]]) # predictions are blups, instead of entire mixmeta. give same names as predict.mixmeta
      exposure.jl.mean <- exposure.by[l,] # mean for relevant strata already prepared
      
      # Define exposure percentile with lowest incidence of injuries
      bvarjl <- crossbasis(exposure.jl.mean, argvar=list(fun=espline, knots=exposure.jl.mean[paste0(eknots*100)]), arglag=li_arglag) # crossbasis
      bvar1jl <- bvarjl%*%mvpred_jl$fit # multiply coefficients for combined effect
      bvar1jl <- bvar1jl[order(bvar1jl[,1]),, drop=F] # sort in ascending order (not required, but easier for visualisation
      bvar1jl <- bvar1jl[between(rownames(bvar1jl), cenpen.minmax[1], cenpen.minmax[2]),] # exclude percentiles outside range for selecting optimal percentile, converts to vector
      cenindjl <- as.numeric(names(bvar1jl[which.min(bvar1jl)])) # exposure percentile with lowest incidence of outcome (optimal percentile)
      
      # Define centering percentile for country
      if (is.null(cenpen)) {
        cenperjl <- pmin(pmax(predper[cenindjl],cenpen.minmax[1]),cenpen.minmax[2]) # min/max is 10th/90th. chooses 10, as incidence lowest at 10
      } else if(is.numeric(cenpen)) {
        cenperjl <- cenpen # use defined centering value
      } else if(cenpen %in% c('mean','Mean','average','Average')) {
        cenperjl <- 'mean' # to be used as code for mean
      }
      
      if (is.numeric(cenperjl)) {
        centrejl <- exposure.j.mean[paste0(cenperjl)] # centered prediction (from corresponding percentile)
      } else {
        centrejl <- sum(sapply(temps[l], sum)) / sum(sapply(temps[l], length)) # mean exposure across all strata
      }
      cpjl <- crosspred(bvarjl, coef=mvpred_jl$fit, vcov=mvpred_jl$vcov, model.link="log", at=exposure.jl.mean, cen=centrejl)
      
      # Plot: Overall cumulative exposure-response association (all exposure values, lag reduced)
      .yaxis.jl <- bigvaryaxis/100+1 # adjusting axis from % to RR
      
      png(file = paste0(s2locations, .name, ' Overall e-r.png'), res=gdpi, width=glength, height=glength) # plot location. File name based on heat metric
      plot(cpjl,lwd=2,col="white",axes=F, ylab=str_remove_all(paste0('Percent change in ',tolower(outcome.var)," (%)"),'\\(000s\\)'),
           xlab='', ylim=.yaxis.jl, mgp=c(2.5,0.5,0))
      title(main=.name, line=nline) # closer to plot
      title(xlab=paste0(exposure.var," (°C)"), line=3.2) # xtitle, moved away from graph
      ind1 <- cpjl$predvar<=cpjl$cen
      ind2 <- cpjl$predvar>=cpjl$cen
      lines(cpjl$predvar[ind1],cpjl$allRRfit[ind1],col=4,lwd=2)
      lines(cpjl$predvar[ind2],cpjl$allRRfit[ind2],col=2,lwd=2)
      axis(1, at=exposure.jl.mean[indlab], labels=predper[indlab], cex.axis=0.8, gap.axis=0.1, mgp=c(2.5,0.3,0))
      mtext('Percentile', 1, at=min(cpjl$predvar)-0.5, outer=F, cex=0.8, adj=1)
      axis(1, at=seq(ceiling(min(cpjl$predvar)),max(cpjl$predvar),by=1), cex.axis=0.8, line=xline2, mgp=c(2.5,0.3,0))
      mtext('Value', 1, at=min(cpjl$predvar)-0.5, outer=F, cex=0.8, line=xline2, adj=1)
      axis(2,cex.axis=0.8, at=seq(.yaxis.jl[1],.yaxis.jl[2], by=bigvarby/100), labels=seq(from=bigvaryaxis[1], to=bigvaryaxis[2], by=bigvarby)) 
      abline(v=cpjl$cen,lty=2) # centre line
      abline(v=c(exposure.jl.mean[c("2.5","97.5")]), lty=c(3,3)) # correspond to mean and extreme lines
      dev.off() # Save image + clear settings
    }
  }
  
  
  ### Full meta-analysis with random-predictor. BLUPs will be conditional expectations given the random effects (Sera et al 2019)
  # outcome.exposure.loc.bv
  png(file = paste0(s2locations, 'Overall e-rs.png'), res=gdpi, width=glength.3by3, height=glength.3by3) # plot location. File name based on heat metric
  par(mfrow=c(3,3), mar=c(3,2.5,1,1), oma=c(0,0,0,0), mgp=c(1.5,0.5,0)) # 3*3, space between borders and text
    bigvarc <- 0 # default marker for recordPlot()
    for(j in rev(by.vars3)) { # City and outorin
      # j <- 'City'
      mixv_j <- mixmeta(coef~1, vcov, random=~1|by.vars.ds[,get(j)], data=by.vars.ds) # BLUP will be conditionally dependent on random predictor. Usage of predictor like this identical to City (city results identical for both indoor and outdoor) or outorin (outorin results identical for all cities)
      if(j=='City') {
        jstrat_no <- seq(0,length.ds.stratum-2,by=2)+1 # every odd number
        ds.stratumj <- word(unique(ds.stratum))[jstrat_no] # new ds.stratum (technically
      } else {
        jstrat_no <- c(1:2)
        ds.stratumj <- word(unique(ds.stratum), start=-1L)[jstrat_no]
      }
      blups_j <-  blup(mixv_j, vcov=T)[jstrat_no] # BLUPs for each j
      
      for(l in jstrat_no) { # repeat for each strata (city or out/in)
        m <- which(jstrat_no==l) # ordering of l, instead of l itself, which is needed for extracting blups
        .name <- ds.stratumj[m] # get name for plot
        mvpred_jl <- list('fit'=blups_j[[m]][["blup"]], 'vcov'=blups_j[[m]][["vcov"]]) # predictions are blups, instead of entire mixmeta. give same names as predict.mixmeta
        exposure.jl.mean <- exposure.by[l,] # mean for relevant strata already prepared
        
        # Define exposure percentile with lowest incidence of injuries
        bvarjl <- crossbasis(exposure.jl.mean, argvar=list(fun=espline, knots=exposure.jl.mean[paste0(eknots*100)]), arglag=li_arglag) # crossbasis
        bvar1jl <- bvarjl%*%mvpred_jl$fit # multiply coefficients for combined effect
        bvar1jl <- bvar1jl[order(bvar1jl[,1]),, drop=F] # sort in ascending order (not required, but easier for visualisation
        bvar1jl <- bvar1jl[between(rownames(bvar1jl), cenpen.minmax[1], cenpen.minmax[2]),] # exclude percentiles outside range for selecting optimal percentile, converts to vector
        cenindjl <- as.numeric(names(bvar1jl[which.min(bvar1jl)])) # exposure percentile with lowest incidence of outcome (optimal percentile)
        
        # Define centering percentile for country
        if (is.null(cenpen)) {
          cenperjl <- pmin(pmax(predper[cenindjl],cenpen.minmax[1]),cenpen.minmax[2]) # min/max is 10th/90th. chooses 10, as incidence lowest at 10
        } else if(is.numeric(cenpen)) {
          cenperjl <- cenpen # use defined centering value
        } else if(cenpen %in% c('mean','Mean','average','Average')) {
          cenperjl <- 'mean' # to be used as code for mean
        } 
        
        if (is.numeric(cenperjl)) {
          centrejl <- exposure.j.mean[paste0(cenperjl)] # centered prediction (from corresponding percentile)
        } else {
          centrejl <- sum(sapply(temps[l], sum)) / sum(sapply(temps[l], length)) # mean exposure across all strata
        }
        
        cpjl <- crosspred(bvarjl, coef=mvpred_jl$fit, vcov=mvpred_jl$vcov, model.link="log", at=exposure.jl.mean, cen=centrejl)
        .yaxis.jl <- bigvaryaxis/100+1 # adjusting axis from % to RR
        
        bigvarc <- bigvarc + 1 # change recordPlot marker
        plot(cpjl,lwd=2,col="white",ylab="Percent change (%)", # str_remove_all(paste0('RR for ',tolower(outcome.var)),'\\(000s\\)') # 'Relative risk'
             xlab='', ylim=.yaxis.jl, yaxt='n')
        title(main=.name, line=0.2) # closer to plot
        title(xlab=paste0(exposure.var," (°C)"), line=1.5) # xtitle, moved away from graph
        ind1 <- cpjl$predvar<=cpjl$cen
        ind2 <- cpjl$predvar>=cpjl$cen
        lines(cpjl$predvar[ind1],cpjl$allRRfit[ind1],col=4,lwd=2)
        lines(cpjl$predvar[ind2],cpjl$allRRfit[ind2],col=2,lwd=2)
        axis(2,cex.axis=0.8, at=seq(.yaxis.jl[1],.yaxis.jl[2], by=bigvarby/100), labels=seq(from=bigvaryaxis[1], to=bigvaryaxis[2], by=bigvarby), las=1)
        abline(v=cpjl$cen,lty=2) # centre line
        abline(v=c(exposure.jl.mean[c("2.5","97.5")]), lty=c(3,3)) # correspond to mean and extreme lines
        assign(paste0('bigvar',bigvarc), recordPlot()) # assign recordPlot based on marker
      }
    }
  dev.off() 



######################### END ############################
######################### END ############################
  