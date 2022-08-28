


predper <- c(seq(0,1,0.1),2,2.5,3:33,100/3,34:66,200/3,67:97,97.5,98,seq(99,100,0.1)) # percentiles (%). 100/3 and 200/3 work even without rounding
predper.short <- c(1,2.5,10,25,50,75,90,97.5,99) # c(1,10,90,99) # percentiles to report RRs for



##################################################
### MODEL, INDOOR/OUTDOOR AND DISTRIBUTION
##################################################

# Type of outcome
if (str_detect(outcome.var,'claims') | str_detect(outcome.var,'injur') | str_detect(outcome.var,'disease') | str_detect(outcome.var,'llness') | str_detect(outcome.var,'OII')) {
  type.outcome <- 'oi' # Injury and/or illness (disease or condition)
} else {
  type.outcome <- 'cost' # Costs
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



##################################################
### Yaxis for graphs as percentages. Will be adjusted to RR in code as required
##################################################

rryaxis <- c(-20,60) # rr.yaxis
rryaxis.s <- c(rryaxis[1]-40,rryaxis[2]) # rryaxis*1.5
bigvaryaxis <-c(-15,60)  # bigvar.yaxis




##################################################
### Set up . Create objects to store results
##################################################


# Matrices with 1 column: Number of outcomes, formula, exposure at 0, AIC and dispersion (if poisson used)
no.outcome <- no.years <- form <- m.aic <- m.dispersion <- m.r2 <- m.devexp <- matrix(NA, length.ds.stratum, 1, dimnames=list(ds.stratum)) # matrix of NAs, with a column to denote total number of outcomes
colnames(m.aic) <- 'AIC'
colnames(m.dispersion) <- 'Dispersion parameter'
colnames(m.r2) <- 'R^2'
colnames(m.devexp) <- 'Deviance explained'

# Exposure. Matrix of NAs, with a column to denote total number of outcomes
exposure <- matrix(NA, length.ds.stratum, 4 + length(predper.short), dimnames=list(ds.stratum, c('mean','range','min','max',paste0(predper.short,'%'))))

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
m.tweedie.shape <- matrix(NA, length.ds.stratum, 1, dimnames=list(ds.stratum)) # matrix of NAs, with a column to denote optimised shape parameter
# m.tweedie.shape <- matrix(NA, length.ds.stratum,3, dimnames=list(ds.stratum)) # matrix of NAs, with columns to dneote optimised shape parameters and their 95% CI

# Matrix for k.index
no.smooth <- str_count(trnd, ',bs=') # number of smooths
gam.k <- matrix(NA, length.ds.stratum*no.smooth, 4, dimnames=list(rep(ds.stratum, each=no.smooth), c("k'","edf","k-index","p-value"))) # matrix of NAs, with rows per by.vars

# List objects
temps <- outcomes <- model <- model.omit <- model2 <- m.coef <- res <- model.checks <- collin <- gam_k <- gam.concurvity <- check <- model.checks.month <- model.checks.week <- check.res.coef <- check.res.coef.sig <- exposure.rr.s1 <- red <- list() #  Model, residuals, residual length + plots and earity



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
      var.power[[i]] <- tweedie.profile(formula=form[i,], data=.ds, p.vec=seq(1.02,1.98,(1.98-1.02)/5), method='series', do.ci=F, do.smooth=T) # estimate optimal shape parameter for 10 values ranging from 1.02 to 1.98
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


