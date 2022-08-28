
########################################################################################################
# SECOND-STAGE ANALYSIS: MULTIVARIATE META-ANALYSIS FOR OVERALL RESULTS
########################################################################################################
##################################################
### Folder and subfolder destinations
##################################################

### Stage 2
suppressWarnings(dir.create(paste0(outcome.exposure.loc,'S2 Cen',cenpen,' ',mmpred))) # Create folder if it doesn't exist already. Warning if exists (doesn't replace)
outcome.exposure.loc.s2 <- paste0(outcome.exposure.loc,'S2 Cen',cenpen,' ',mmpred,'/') # folder destination (add /)

# Meta-analysis
suppressWarnings(dir.create(paste0(outcome.exposure.loc.s2,'Meta-analysis'))) # Create folder if it doesn't exist already. Warning if exists (doesn't replace)
outcome.exposure.loc.mm <- paste0(outcome.exposure.loc.s2,'Meta-analysis/') # folder destination (add /)

# Attributable
suppressWarnings(dir.create(paste0(outcome.exposure.loc.s2,'Attributable'))) # Create folder if it doesn't exist already. Warning if exists (doesn't replace)
outcome.exposure.loc.attr <- paste0(outcome.exposure.loc.s2,'Attributable/') # folder destination (add /)

# National
suppressWarnings(dir.create(paste0(outcome.exposure.loc.s2,'National'))) # Create folder if it doesn't exist already. Warning if exists (doesn't replace)
outcome.exposure.loc.s2.national <- paste0(outcome.exposure.loc.s2,'National/') # folder destination (add /)

# By-variables (formerly supplementary)
suppressWarnings(dir.create(paste0(outcome.exposure.loc.s2,'By-vars'))) # Create folder if it doesn't exist already. Warning if exists (doesn't replace)
outcome.exposure.loc.s2.byvar <- paste0(outcome.exposure.loc.s2,'By-vars/') # folder destination (add /)



##################################################
### Meta-analysis: results of reduced coef and compute BLUPS (best linear unbiased prediction)
##################################################

exposure.country.mean <- colMeans(exposure.by) # could arguably weight by area number of years (study frame), but Gasaparrini didn't do this either and impact likely negligible

# Calculate potential meta-analysis predictors 
exposure.mean <- exposure[,1] # mean and range of exposure as per Gasparrani et al 2015 Lancet
exposure.range <- exposure[,2]
exposure.ci <- exposure[,"97.5%"] - exposure[,"2.5%"]
exposure.iqr <- exposure[,"75%"] - exposure[,"25%"]
exposure.min <- exposure[,3]
exposure.max <- exposure[,4]
exposure.025 <- exposure[,"2.5%"]
exposure.05 <- exposure.by[,"5"]
exposure.1 <- exposure[,"10%"]
exposure.25 <- exposure[,"25%"]
exposure.5 <- exposure[,"50%"]
exposure.75 <- exposure[,"75%"]
exposure.9 <- exposure[,"90%"]
exposure.95 <- exposure.by[,"95"]
exposure.975 <- exposure[,"97.5%"]
# cor(exposure.mean, exposure.95) # check if predictors are highly correlated  # single T values collinear with one another e.g. mean, high low/Ts (max and min are about 74%, rest are high 90s), and range values (range/ci/iqr) collinear with one another but not with single T values
# exposure.test <- exposure.mean + exposure.95 # testing combining variables, as they are collinear # performs halfway between the two vars, doesn't help

# Meta-analysis. Default methods is reml
mv <- mixmeta(coef~1, vcov, data=by.vars.ds)
mv.results <- mv.results.fn(mv) 
mv.results1 <- c(mv.results[6:8], mv.results[1:4])
names(mv.results1) <- c('Wald statistic', 'Wald df', 'Wald p-value', names(mv.results[1:4]))
save(mv.results1, file=paste0(outcome.exposure.loc.mm, 'mv.results1.rda'))
write.csv(mv.results1, file=paste0(outcome.exposure.loc.mm, 'MM tests.csv'), na='', row.names=T) # create csv file

# Meta-analysis, with exclusion of Darwin and/or Hobart. No random prediction
mv.all <- mixmeta(as.formula(paste0('coef~',mm.pred.eq)), vcov, data=by.vars.ds)
mv.nod <- mixmeta(as.formula(paste0("coef[!str_detect(rownames(coef),'Darwin'),]~",mm.pred.eq)), vcov[!str_detect(names(vcov),'Darwin')], data=by.vars.ds[!str_detect(City,'Darwin')])
mv.noh <- mixmeta(as.formula(paste0("coef[!str_detect(rownames(coef),'Hobart'),]~",mm.pred.eq)), vcov[!str_detect(names(vcov),'Hobart')], data=by.vars.ds[!str_detect(City,'Hobart')])
coef.nodh <- coef[!(str_detect(rownames(coef),'Darwin') | str_detect(rownames(coef),'Hobart')),] # needed to be recognise by mixmeta
mv.nodh <-  mixmeta(as.formula(paste("coef.nodh~",mm.pred.eq,collapse='')), vcov[!str_detect(names(vcov),paste(c('Darwin','Hobart'),collapse='|'))], data=by.vars.ds[!str_detect(City,paste(c('Darwin','Hobart'),collapse='|'))])
coef.nodha <- coef[!(str_detect(rownames(coef),'Darwin') | str_detect(rownames(coef),'Hobart') | str_detect(rownames(coef),'Adelaide')),] # needed to be recognise by mixmeta. Adelaide shouldn't be a sore thumb: this determines difference
mv.nodha <-  mixmeta(as.formula(paste("coef.nodha~",mm.pred.eq,collapse='')), vcov[!str_detect(names(vcov),paste(c('Darwin','Hobart','Adelaide'),collapse='|'))], data=by.vars.ds[!str_detect(City,paste(c('Darwin','Hobart','Adelaide'),collapse='|'))])
mv.results.ex <- rbind(mv.results.fn(mv.all), mv.results.fn(mv.nod), mv.results.fn(mv.noh), mv.results.fn(mv.nodh), mv.results.fn(mv.nodha))
rownames(mv.results.ex) <- c('All cities','Exclude Darwin','Exclude Hobart','Exclude Darwin and Hobart','Exclude Adelaide, Darwin and Hobart')
# colnames(mv.results.ex) <- c('Wald statistic', 'Wald df', 'Wald p-value', names(mv.results[1:4]))
save(mv.results.ex, file=paste0(outcome.exposure.loc.mm, 'mv.results.ex.rda'))
write.csv(mv.results.ex, file=paste0(outcome.exposure.loc.mm, 'MM tests ex.csv'), na='', row.names=T) # create csv file


## Pool estimated overall cumulative exposure-response associations to obtain overall estimation at each location level for main model
blups0 <- blup(mv, vcov=T) # BLUPS. Element for each stratum
# mv$dim$q

# Add city-level and indoor-outdoor results to blups. Created by meta-analysing the BLUP results for the indoor and outdoor populations of each city
city_results <- .city_results <- io_results <- .io_results <- list()
# if (isTRUE(io.marker)) { # if stratify by indoors/outdoors
#   for(i in 1:(length.ds.stratum/2)) { # for each city
#     .blups.coef <- rbind(blups0[[i*2-1]][["blup"]], blups0[[i*2]][["blup"]]) # blup coefficient for both indoors and outdoors
#     .blups.vcov <- list(blups0[[i*2-1]][["vcov"]], blups0[[i*2]][["vcov"]]) # blup variance-covariance for both indoors and outdoors
#     names(.blups.coef) <- names(.blups.vcov) <- c(names(vcov)[i*2-1], names(vcov)[i*2]) # name both the above
#     .city_results[[i]] <- mixmeta(.blups.coef, .blups.vcov) # combine indoor and outdoor results, no meta predictors (fails if included)
#     # city_results[[i]] <- blup(.city_results[[i]], vcov=T)[["1"]] # blup here is identical to mixmeta results; blup used only to keep the format consistent with the overall blup. "1" and "2" are identical (indoors and outdoors), so only include one of them
#     city_results[[i]] <-  list('blup'=.city_results[[i]][["coefficients"]], 'vcov'=.city_results[[i]][["vcov"]]) # extract mixmeta results and format identical to blups
#   }
#   for(i in 1:2) { # for each indoor/outdoor combination
#     .blups.io.coef <- rbind(blups0[[i]][["blup"]], blups0[[i+2]][["blup"]], blups0[[i+4]][["blup"]], blups0[[i+6]][["blup"]], blups0[[i+8]][["blup"]], blups0[[i+10]][["blup"]], blups0[[i+12]][["blup"]]) # blup coefficient for each city
#     .blups.io.vcov <- list(blups0[[i]][["vcov"]], blups0[[i+2]][["vcov"]], blups0[[i+4]][["vcov"]], blups0[[i+6]][["vcov"]], blups0[[i+8]][["vcov"]], blups0[[i+10]][["vcov"]], blups0[[i+12]][["vcov"]]) # blup variance-covariance for each city
#     names(.blups.io.coef) <- names(.blups.io.vcov) <- c(names(vcov)[i], names(vcov)[i+2], names(vcov)[i+4], names(vcov)[i+6], names(vcov)[i+8], names(vcov)[i+10], names(vcov)[i+12]) # name both the above
#     .io_results[[i]] <- mixmeta(.blups.io.coef, .blups.io.vcov) # combine indoor and outdoor results, no meta predictors (fails if included)
#     # io_results[[i]] <- blup(.io_results[[i]], vcov=T) # blup here is NOT identical to mixmeta results: results vary for each city
#     io_results[[i]] <- list('blup'=.io_results[[i]][["coefficients"]], 'vcov'=.io_results[[i]][["vcov"]])
#   }
#   blups <- c(blups0, city_results, io_results) # order: individual results, city results, indoor/outdoor results
# } else {
#   blups <- blups0 # nothing to add
# }
blups <- blups0

### Heat metrics adjusted for indoor/outdoor, if not already done (helps deals with potential bugs)
if(isTRUE(io.marker)) { # if indoor/outdoor used
  if (!(by.vars.oi %in% colnames(climate))) { # if have not already duplicated supplementary barra data for outorin
    climate[, get('by.vars.oi'):='Indoors'] # Create outorin variable and set to indoors
    climate <- bind_rows(climate, climate)  # duplicate dataset
    climate[duplicated(climate), get('by.vars.oi'):='Outdoors'] # for duplicate, set outorin variable to outdoors
  }
}



##################################################
### Define minimum occupational injuries values (no longer used in analysis, but still outputted)
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

# Mean values of mixmeta predictors as a data frame
# datanew <- data.frame(exposure.mean=mean(tapply(exposure[,1],ds.stratum,mean)),
#                       exposure.max=mean(tapply(exposure[,4],ds.stratum,mean))) # formerly exposure[,2] for range
# datanew <- as.data.frame(t(colMeans(mm.pred))) # names of datanew and mvpred must match
# as.data.frame(t(colMeans(cbind(exposure.mean, exposure.range))))
datanew <- as.data.frame(t(colMeans(mm.pred.m))) # names of datanew and mvpred must match

mvpred <- predict.mixmeta(mv, datanew, vcov=T, format="list") # With mixed model, results are the same
# predict(mv, by.vars.ds[get(by.vars.oi=='Indoors'], vcov=T, format="list")
# colMeans(mv[["fitted.values"]]) # pooled coefficients equal to colmeans


# Define exposure percentile with lowest incidence of injuries
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
save(cp, file=paste0(outcome.exposure.loc.s2, 'cp.rda')) # save for fwald1

# Cumulative exposure-response RRs (associations) at set percentiles
cp.rrfit <- cp$allRRfit[predper %in% c(1,2.5,10,25,50,75,90,97.5,99)] # cp$allRRfit[names(cp$allRRfit) %in% exposure.country.mean[c('1','2.5','10','25','50','75','90','97.5','99')]]
cp.rrlow <- cp$allRRlow[predper %in% c(1,2.5,10,25,50,75,90,97.5,99)] # cp$allRRlow[names(cp$allRRlow) %in% exposure.country.mean[c('1','2.5','10','25','50','75','90','97.5','99')]]
cp.rrhigh <- cp$allRRhigh[predper %in% c(1,2.5,10,25,50,75,90,97.5,99)] # cp$allRRhigh[names(cp$allRRhigh) %in% exposure.country.mean[c('1','2.5','10','25','50','75','90','97.5','99')]]
cp.rr2 <- cbind(outcome.var, matrix(paste0(round.fn(cp.rrfit, 2),' (',round.fn(cp.rrlow, 2),' to ',round.fn(cp.rrhigh, 2),')'), ncol=length(cp.rrfit)))
cp.rr3 <- cbind(outcome.var, matrix(paste0(round.fn(cp.rrfit, 3),' (',round.fn(cp.rrlow, 3),' to ',round.fn(cp.rrhigh, 3),')'), ncol=length(cp.rrfit)))
colnames(cp.rr2) <- colnames(cp.rr3) <- c('Outcome','1%','2.5%','10%','25%','50%','75%','90%','97.5%','99%')
write.csv(cp.rr2, file=paste0(outcome.exposure.loc.s2, 'CP RR round2.csv'), na='', row.names=T) # create csv file
write.csv(cp.rr3, file=paste0(outcome.exposure.loc.s2, 'CP RR round3.csv'), na='', row.names=T) # create csv file


## Cumulative exposure-response RRs (associations)
# Originally to be at set values (25,26,28,30) matching Ma et al. 2019 and the ISO WBGT thresholds for acclimatised workers (33 is resting), but this information is in overall expoure response-curve
cp.vpct <- predper # doesn't have clinical relevance as it's pooled from many sites, but is handy to know for reading graphs
cp.vrrfit <- cp$allRRfit# [cp$predvar > 25] # c(25,26,28,30,33)
cp.vrrlow <- cp$allRRlow
cp.vrrhigh <- cp$allRRhigh
cp.vrr <- cbind('Exposure percentile'=cp.vpct, 'Exposure'=as.numeric(names(cp.vrrfit)), 'RR'=cp.vrrfit, 'RRlow'=cp.vrrlow, 'RRhigh'=cp.vrrhigh)
# cp.vrr <- cbind(outcome.var, matrix(paste0(round.fn(cp.vrrfit, 4),' (',round.fn(cp.vrrlow, 4),' to ',round.fn(cp.vrrhigh, 4),')'), ncol=length(cp.vrrfit)))
write.csv(cp.vrr, file=paste0(outcome.exposure.loc.s2, 'CP Exposure.csv'), na='', row.names=T) # create csv file

# Table format
cp.table <- data.table(cp.vpct, round.fn(as.numeric(names(cp.vrrfit)),2), 'RR'=paste0(round.fn(cp.vrr[,'RR'],2),' (',round.fn(cp.vrr[,'RRlow'],2),' to ',round.fn(cp.vrr[,'RRhigh'],2),')'))
colnames(cp.table) <- c(paste(exposure.var, 'percentile'), exposure.var, paste(str_remove_all(outcome.var,'\\ \\(000s\\)'), 'RR (95% CI)'))
write.csv(cp.table, file=paste0(outcome.exposure.loc.s2, 'CP Exposure table.csv'), na='', row.names=T) # create csv file


## City-specific and national exposure percentiles of minimum occupational OII (PMOI) and corresponding exposure
mincity <- rbind(cbind(minpercity,mintempcity),c(as.numeric(names(centre)), centre)) # combine PMOI, including selected overall (within range)
colnames(mincity) <- c('Percentile',exposure.var)
rownames(mincity) <- c(ds.stratum,'Overall') # replaces moip as the overall name
write.csv(mincity, file=paste0(outcome.exposure.loc.s2, 'Exposure PMOIs.csv'), na='', row.names=T) # create csv file



################################################################################
# LAG-RESPONSE ASSOCIATIONS, USING SELECTED CENTERING VALUE
################################################################################

# Objects to store overall cumulative exposure results (lag-response association at set percentiles vs centering value)
cvlag.length <- length(crossreduce(.cb, model=model[[i]], "var", value=exposure.by[1,'1'], cen=.cen, model.link='log')$coefficients) # only interested in length, uses last cb from S1 but that is irrelevant

coeflag1 <- coeflag2<- coeflag10 <- coeflag25 <- coeflag50 <- coeflag75 <- coeflag90<- coeflag975 <- coeflag99 <- matrix(NA, length.ds.stratum, cvlag.length, dimnames=list(ds.stratum)) # not really sure why columns is 3. I though 8 was for columns for lag 0 + each lag day. CHANGED WHEN I INCREASED PARAMETES FROM edf TO edf + 1
cpmodel <- vcovlag1 <- vcovlag2 <- vcovlag10 <- vcovlag25 <- vcovlag50 <- vcovlag75 <- vcovlag90 <- vcovlag975 <- vcovlag99 <- vector("list", length.ds.stratum)
names(cpmodel)<-names(vcovlag1)<-names(vcovlag2)<-names(vcovlag10)<-names(vcovlag25)<-names(vcovlag50)<-names(vcovlag75)<-names(vcovlag90)<-names(vcovlag975)<-names(vcovlag99) <- ds.stratum

# Run model for each by.vars with re-centered values to obtain lag-reponse relationship at the 1st and 99th percentile
tic()
for(i in ds.stratum) {
  # print(paste('Stage 2 individual models:',i))
  .ds <- daily.ds[stratum==i] # dataset for each unique value
  .d <- .ds[-(1:lmax)] # dataset for each unique value # daily.ds[stratum=='Adelaide, Outdoors']
  .name <- unique(.ds$stratum) # names with all of a, b and c together
  
  .outcome <- .ds[,get(outcome.var)] # outcome
  
  # Calculate crossbasis (comments in Gasparrini and Martinez-Solanas code state it is centred, although this isn't the case). Model appears same, with changes made for crosspred and crossreduce with cen 
  .cb <- crossbasis(temps[[i]], argvar=list(fun=espline, knots=quantile(temps[[i]], eknots, na.rm=T)), lag=lmax, arglag=li_arglag) # lag=lag2
  
  # Model using crossbasis
  .no.years <- no.years[i,]
  
  # Contrary to comments, model is unchanged, except partially from GAM shape parameter being rounded to 3 dp. Quicker, easier and more precise just to use old model
  # if(model.type %in% c('a','gam','GAM')) {
  #   if (distribution[1]=='Tweedie') { # use same Tweedie power function as before, to save time with re-estimating. As gam only returns a rounded (to 3 dp) shape parameter, model will be negligbly less precise
  #     model2[[i]] <- gam(as.formula(form[i,]), family=Tweedie(m.tweedie.shape[i,], link='log'), .ds, na.action="na.exclude", method=gam.convergence, maxit=max.iter)
  #   } else {
  #     model2[[i]] <- gam(as.formula(form[i,]), family=distribution, .ds, na.action="na.exclude", method=gam.convergence, maxit=max.iter)
  #   }
  # }
  # else {
  #   if (distribution[1]=='Tweedie') { # use same Tweedie power function as before, to save time with re-estimating. Exact value produced from glm (automatically produced by to 3 dp)
  #     model2[[i]] <- glm(as.formula(form[i,]), family=tweedie(var.power=m.tweedie.shape[i,], link.power=0), .ds, na.action="na.exclude", control=list(maxit=max.iter))
  #   } else {
  #     model2[[i]] <- glm(as.formula(form[i,]), family=distribution, .ds, na.action="na.exclude", control=list(maxit=max.iter))
  #   }
  # }
  model2[[i]] <- model[[i]]
  
  # Predictions and reduction to lag-response at 1st (extreme cold) and 99th (extreme hot) percentiles with new centering (changes coef and vcov)
  if (is.numeric(cenpercountry)) {
    .cen <- quantile(temps[[i]], cenpercountry/100, na.rm=T) # quantile
  }  else {
    .cen <- exposure.mean[i] # mean
  }
  
  # Predictions and reduction to lag-response at set percentiles. Requires centering required as it changes coef-vcov
  .perc <-  exposure.by[i,] # all relevant percentiles per i
  # 1st percentile
  redlag1 <- crossreduce(.cb, model=model2[[i]], "var", value=.perc['1'], cen=.cen, model.link='log')
  coeflag1[i,] <- coef(redlag1)
  vcovlag1[[i]] <- vcov(redlag1)
  # 2.5th percentile
  redlag2 <- crossreduce(.cb, model=model2[[i]], "var", value=.perc['2.5'], cen=.cen, model.link='log')
  coeflag2[i,] <- coef(redlag2)
  vcovlag2[[i]] <- vcov(redlag2)
  # 10th percentile
  redlag10 <- crossreduce(.cb, model=model2[[i]], "var", value=.perc['10'], cen=.cen, model.link='log')
  coeflag10[i,] <- coef(redlag10)
  vcovlag10[[i]] <- vcov(redlag10)
  # 25th percentile
  redlag25 <- crossreduce(.cb, model=model2[[i]], "var", value=.perc['25'], cen=.cen, model.link='log')
  coeflag25[i,] <- coef(redlag25)
  vcovlag25[[i]] <- vcov(redlag25)
  # 50th percentile
  redlag50 <- crossreduce(.cb, model=model2[[i]], "var", value=.perc['50'], cen=.cen, model.link='log')
  coeflag50[i,] <- coef(redlag50)
  vcovlag50[[i]] <- vcov(redlag50)
  # 75th percentile
  redlag75 <- crossreduce(.cb, model=model2[[i]], "var", value=.perc['75'], cen=.cen, model.link='log')
  coeflag75[i,] <- coef(redlag75)
  vcovlag75[[i]] <- vcov(redlag75)
  # 90th percentile
  redlag90 <- crossreduce(.cb, model=model2[[i]], "var", value=.perc['90'], cen=.cen, model.link='log')
  coeflag90[i,] <- coef(redlag90)
  vcovlag90[[i]] <- vcov(redlag90)
  # 97.5th percentile
  redlag975 <- crossreduce(.cb, model=model2[[i]], "var", value=.perc['97.5'], cen=.cen, model.link='log') 
  coeflag975[i,] <- coef(redlag975)
  vcovlag975[[i]] <- vcov(redlag975)
  # 99th percentile
  redlag99 <- crossreduce(.cb, model=model2[[i]], "var", value=.perc['99'], cen=.cen, model.link='log') 
  coeflag99[i,] <- coef(redlag99)
  vcovlag99[[i]] <- vcov(redlag99)
}
toc()


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

# # Checking values
# cenpercountry # minimum occupational-injuries percentile (MOIP)
# predper[cp$allRRfit==1] # confirming MOIP. This is at a WGBT of 11.0407549172094

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
colnames(results.short) <- colnames(results) <- c('Percentile','RR','RRlow','RRhigh') # column names
write.csv(results.short, file=paste0(outcome.exposure.loc.s2.national, 'RR (95% CI) short.csv'), na='', row.names=F) # create csv file
write.csv(results, file=paste0(outcome.exposure.loc.s2.national,'RR (95% CI).csv'), na='', row.names=F) # create csv file



###############################################################################
# COMPUTE ATTRIBUTABLE INJURIES FOR EACH CITY, WITH EMPIRICAL CI ESTIMATED USING RE-CENTERED BASES
################################################################################

# Create vectors to store total injuries, accounting for missing
totclaims <- rep(NA,length.ds.stratum)
names(totclaims) <- ds.stratum

# Objects to store attributable injuries (simulations)
if(isTRUE(a.full)) {
  sim.names <- c("Total","Cold","Heat","Extreme cold","Moderate cold",'Mild cold','Non-extreme cold','Mild heat',"Moderate heat","Non-extreme heat","Extreme heat")
} else {
  im.names <- c("Total","Cold","Heat","Extreme cold","Non-extreme cold","Non-extreme heat","Extreme heat")
}


length.sim <- length(sim.names)
matsim <- matrix(NA, length.ds.stratum, length.sim, dimnames=list(ds.stratum, sim.names)) # matrix: attributable injuries
arraysim <- array(NA, dim=c(length.ds.stratum, length.sim, nsim), dimnames=list(ds.stratum, sim.names)) # array: attributable injuries CI

# Run loop
gc()
print('Compute attributable numbers')
tic()
for(i in ds.stratum){
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
  matsim[i,"Non-extreme cold"] <- attrdl(temps[[i]], .cb, .outcome, coef=.blups$blup, vcov=.blups$vcov, type="an", dir='forw', cen=.cen,
                                         range=c(.perc["2.5"],.cen))
  
  matsim[i,"Non-extreme heat"] <- attrdl(temps[[i]], .cb, .outcome, coef=.blups$blup, vcov=.blups$vcov, type="an", dir='forw', cen=.cen,
                                         range=c(.cen,.perc["97.5"]))
  matsim[i,"Extreme heat"] <- attrdl(temps[[i]], .cb, .outcome, coef=.blups$blup, vcov=.blups$vcov, type="an", dir='forw', cen=.cen,
                                     range=c(.perc["97.5"],100))
  if(isTRUE(a.full)) {
    matsim[i,"Moderate cold"] <- attrdl(temps[[i]], .cb, .outcome, coef=.blups$blup, vcov=.blups$vcov, type="an", dir='forw', cen=.cen,
                                        range=c(.perc["2.5"],.perc["10"])) # range=c(.perc[1], .cen)
    matsim[i,"Mild cold"] <- attrdl(temps[[i]], .cb, .outcome, coef=.blups$blup, vcov=.blups$vcov, type="an", dir='forw', cen=.cen,
                                    range=c(.perc["10"],.cen))
    matsim[i,"Mild heat"] <- attrdl(temps[[i]], .cb, .outcome, coef=.blups$blup, vcov=.blups$vcov, type="an", dir='forw', cen=.cen,
                                    range=c(.cen,.perc["90"]))
    matsim[i,"Moderate heat"] <- attrdl(temps[[i]], .cb, .outcome, coef=.blups$blup, vcov=.blups$vcov, type="an", dir='forw', cen=.cen,
                                        range=c(.perc["90"],.perc["97.5"])) # range=c(.cen, .perc[2])
  }
  
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
  arraysim[i,"Non-extreme cold",] <- attrdl(temps[[i]], .cb, .outcome, coef=.blups$blup, vcov=.blups$vcov, type="an", dir='forw', cen=.cen,
                                            range=c(.perc["2.5"],.cen), sim=T, nsim=nsim)
  arraysim[i,"Non-extreme heat",] <- attrdl(temps[[i]], .cb, .outcome, coef=.blups$blup, vcov=.blups$vcov, type="an", dir='forw', cen=.cen,
                                            range=c(.cen,.perc["97.5"]), sim=T, nsim=nsim)
  arraysim[i,"Extreme heat",] <- attrdl(temps[[i]], .cb, .outcome, coef=.blups$blup, vcov=.blups$vcov, type="an", dir='forw', cen=.cen,
                                        range=c(.perc["97.5"], 100), sim=T, nsim=nsim)
  if(isTRUE(a.full)) {
    arraysim[i,"Moderate cold",] <- attrdl(temps[[i]], .cb, .outcome, coef=.blups$blup, vcov=.blups$vcov, type="an", dir='forw', cen=.cen,
                                           range=c(.perc["2.5"],.perc["10"]), sim=T, nsim=nsim)
    arraysim[i,"Mild cold",] <- attrdl(temps[[i]], .cb, .outcome, coef=.blups$blup, vcov=.blups$vcov, type="an", dir='forw', cen=.cen,
                                       range=c(.perc["10"],.cen), sim=T, nsim=nsim)
    arraysim[i,"Mild heat",] <- attrdl(temps[[i]], .cb, .outcome, coef=.blups$blup, vcov=.blups$vcov, type="an", dir='forw', cen=.cen,
                                       range=c(.cen,.perc["90"]), sim=T, nsim=nsim)
    arraysim[i,"Moderate heat",] <- attrdl(temps[[i]], .cb, .outcome, coef=.blups$blup, vcov=.blups$vcov, type="an", dir='forw', cen=.cen,
                                           range=c(.perc["90"], .perc["97.5"]), sim=T, nsim=nsim)
  }
  totclaims[i] <- sum(.outcome, na.rm=T) # store total injuries (account for missing)
}
toc()


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
if (isTRUE(io.marker)) { # if stratify by indoors/outdoors
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
  
  an_io <- .an_iolow <- .an_iohigh <- matrix(NA, length.ds.stratum/length(unique(daily.ds[,City])), length.sim, dimnames=list(sort(unique(daily.ds[,get(by.vars.oi)])))) # matrix of NAs, AN for indoor/outdoor
  .an_io <- rowsum(matsim, rep(1:2,length.ds.stratum/2)) # sum every 2 rows for city AN
  for(i in 1:2) { # for each indoor/outdoor combination
    .an_iov <- colSums(arraysim[seq(0,12,by=2)+i,,]) # sum results for each indoor or outdoor, still keeping all repetitions
    .an_iolow[i,] <- apply(.an_iov,1,quantile,0.025) # lower quantile
    .an_iohigh[i,] <- apply(.an_iov,1,quantile,0.975) # upper quantile
  }
  an_io <- matrix(paste0(round.fn(.an_io),' (',round.fn(.an_iolow),' to ',round.fn(.an_iohigh),')'), ncol=length.sim) # combine for presentation
  rownames(an_io) <- rownames(.an_iolow) <- rownames(.an_iohigh) <- sort(unique(daily.ds[,get(by.vars.oi)])) # rownames
  colnames(an_io) <- colnames(antotal) <- sim.names # colnames
  ancountry <- rbind(antotal, an_io, an_city, ancity) # all ANs
} else {
  ancountry <- rbind(antotal, ancity) # ANs without extra stratification
}

# Export attributable numbers
write.csv(ancountry, file=paste0(outcome.exposure.loc.attr, 'Numbers.csv'), na='', row.names=T) # create csv file


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
if (isTRUE(io.marker)) { # if stratify by indoors/outdoors
  totclaimscity <- c(totclaims[1]+totclaims[2],totclaims[3]+totclaims[4],totclaims[5]+totclaims[6],totclaims[7]+totclaims[8],totclaims[9]+totclaims[10],totclaims[11]+totclaims[12],totclaims[13]+totclaims[14])
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
  rownames(afcityio) <- sort(unique(daily.ds[,get(by.vars.oi)])) # rownames
  colnames(afcitycity) <- colnames(afcityio) <- sim.names # colnames
  
  afcountry <- rbind(aftotal, afcityio, afcitycity, afcity) # all AFs
} else {
  afcountry <- rbind(aftotal, afcity) # AFs without extra stratification
}

# Export attributable fractions
write.csv(afcountry, file=paste0(outcome.exposure.loc.attr, 'Fractions.csv'), na='', row.names=T) # create csv file

# Export AF and AN together
acountry <- `dim<-`(str_replace(sprintf('%s,%s', ancountry, afcountry),',',', '), dim(ancountry)) # Combine AF and AN. sprintf pastes elements together with , as unchangeable sep. Replace it with ', '. Then reformat into matrix 
dimnames(acountry) <- dimnames(ancountry)
write.csv(acountry, file=paste0(outcome.exposure.loc.attr, 'Both.csv'), na='', row.names=T) # create csv file

save(ancountry, afcountry, acountry, file=paste0(outcome.exposure.loc.attr, 'acountry.rda')) # save dataset after all changes for easier access



################################################################################
# Plots: overall relationship
################################################################################

indlab <- predper %in% c(0,2.5,10,25,50,75,90,75,97.5,100) # requires predper, so must be in S2

# Plot: Overall cumulative exposure-response association (all exposure values, lag reduced)
oer.yaxis <- seq(0.9,1.5,by=0.1)

png(file = paste0(outcome.exposure.loc.s2.national, 'Overall e-r.png'), res=gdpi, width=glength, height=glength) # plot location. File name based on heat metric
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

png(file = paste0(outcome.exposure.loc.s2.national, 'Overall lag.png'), res=gdpi, width=glength.3by3, height=glength.3by3) # plot location. File name based on heat metric
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
  png(file = paste0(outcome.exposure.loc.s2.byvar, .name,', Overall e-r.png'), res=gdpi, width=glength.3by3, height=glength.3by3) # plot location. File name based on heat metric
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
  png(file = paste0(outcome.exposure.loc.s2.byvar, 'z', .name,', Overall e-r.png'), res=gdpi, width=glength.3by3, height=glength.3by3) # plot location. File name based on heat metric
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
png(file = paste0(outcome.exposure.loc.s2.byvar, 'zOverall e-rs.png'), res=350, width=1700, height=2400) # plot location # paste0('/Users/MatthewBorg/', 'zOverall e-rs.png')
if (isTRUE(io.marker)) { 
  layout(matrix(c(1:14,0),ncol=3,byrow=T))
  par(mfrow=c(5,3), mar=c(2.1,2.2,1,1.5), oma=c(0,0,0,0), mgp=c(1.5,0.4,0), las=1) # 4*4, space between borders and text
} else {
  layout(matrix(c(1:7,0),ncol=2,byrow=T))
  par(mfrow=c(4,2), mar=c(2.1,2.2,1,1.5), oma=c(0,0,0,0), mgp=c(1.5,0.4,0), las=1) # 4*4, space between borders and text
}

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
write.csv(exposure.rr, file=paste0(outcome.exposure.loc.s2.national, 'Exposure values RR.csv'), na='', row.names=F) # create csv file



################################################################################
# CITY AND OUTORIN OVERALL E-R CURVES
################################################################################

### Full meta-analysis with random-predictor. BLUPs will be conditional expectations given the random effects (Sera et al 2019)
if (isTRUE(io.marker)) { 
  
  # Create folder
  suppressWarnings(dir.create(paste0(outcome.exposure.loc.s2,'Big-vars'))) # Create folder if it doesn't exist already. Warning if exists (doesn't replace)
  outcome.exposure.loc.bv <- paste0(outcome.exposure.loc.s2,'Big-vars/') # folder destination (add /)
  
  for(j in by.vars3) { # City and outorin
    # j <- 'City'
    mixv_j <- mixmeta(as.formula(paste0('coef~',mm.pred.eq)), vcov, random=~1|by.vars.ds[,get(j)], data=by.vars.ds) # BLUP will be conditionally dependent on random predictor. Usage of predictor like this identical to City (city results identical for both indoor and outdoor) or outorin (outorin results identical for all cities)
    if(j=='City') {
      jstrat_no <- c(1,3,5,7,9,11,13)
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
      
      png(file = paste0(outcome.exposure.loc.bv, .name, ' Overall e-r.png'), res=gdpi, width=glength, height=glength) # plot location. File name based on heat metric
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
  png(file = paste0(outcome.exposure.loc.bv, 'Overall e-rs.png'), res=gdpi, width=glength.3by3, height=glength.3by3) # plot location. File name based on heat metric
  par(mfrow=c(3,3), mar=c(3,2.5,1,1), oma=c(0,0,0,0), mgp=c(1.5,0.5,0)) # 3*3, space between borders and text
  if (isTRUE(io.marker)) {
    bigvarc <- 0 # default marker for recordPlot()
    for(j in rev(by.vars3)) { # City and outorin
      # j <- 'City'
      mixv_j <- mixmeta(as.formula(paste0('coef~',mm.pred.eq)), vcov, random=~1|by.vars.ds[,get(j)], data=by.vars.ds) # BLUP will be conditionally dependent on random predictor. Usage of predictor like this identical to City (city results identical for both indoor and outdoor) or outorin (outorin results identical for all cities)
      if(j=='City') {
        jstrat_no <- c(1,3,5,7,9,11,13)
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
  }
  dev.off() 
}



######################### END ############################
######################### END ############################