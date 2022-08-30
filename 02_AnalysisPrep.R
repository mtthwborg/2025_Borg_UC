##################################################
### ADJUST THESE PARAMETERS TO ADJUST STAGE 1 (AND 2 INDIRECTLY)
##################################################

## Outcome variable
outcome.var <- 'Total costs' # "Number of illnesses" or "Total costs"

## Exposure variable
exposure.var <- 'Maximum WBGT' 


### Model parameters

## Exposure parameters
espline <- 'ns' # For use with splines package as part of dlnm
eknots <- 0.5 # Knot locations for exposure variable as quantiles. Can also be a vector such as c(0.1,0.9)
edf <- length(eknots) + 1 # Number of coefficients

## Lag parameters
lspline <- 'ns' # For use with splines package as part of dlnm
lmax <- 20 # Duration of lag period (days)
no.lknots <- 1 # Ignored if lspline = 'integer'(unconstrained)

if(lspline=='integer') { # unconstrained dlnm, 1 parameter per lag
  li_arglag <- list(fun=lspline) # knots not useable
} else { # constrained dlnm
  lknots <- seq(from=lmax/(no.lknots+1),to=lmax*no.lknots/(no.lknots+1),length.out=no.lknots)  # equally spaced knots, automated based on lagmax and no.lagknots
  # lknots <- exp(seq(from=log(lmax)/(no.lknots+1),to=log(lmax)*no.lknots/(no.lknots+1),length.out=no.lknots))  # equally spaced knots on lag scale, automated based on lagmax and no.lagknots
  li_arglag <- list(fun=lspline, knots=lknots)
}

## Other
max.iter <- 200 # max number of model iterations # 25 default for glm, 200 for gam (gam requires more iterations)
gam.convergence <- 'REML' # model convergence method for GAM (only used if GAM selected)


### Model specifications

## Type of (generalised) model
model.type <- 'auto' # l/glm = linear, a/gam = additive, auto=auto

## Seasonality / long-term trends for GLM and GAM
mdf <- '4' # Seasonality df
trend <- paste0("ns(Date,",mdf,"*.no.years)")
trend.gam <- paste0("s(date,bs='cr',k=",mdf,"*.no.years,fx=T)") # If use "s(", MUST include "bs=" even if default tp, as later code requires detecting "bs="

## Base formula. Must include ".outcome ~ .cb" which are coded to be the outcome variable and crossbasis (with exposure variable and lag) respectively. if crossbasis for humidity, included as ,cb2
m.formula <- paste0('.outcome ~ .cb + Tue+Wed+Thu+Fri+Sat+Sun +Sat*public.hol + Sun:public.hol + month + offset(log(n)) + Sat*school.hols + Sun:school.hols + Sat*day1 + Sun:day1 + shol +') # Feb+Mar+Apr+May+Jun+Jul+Aug+Sep+Oct+Nov+Dec # month + offset(log(n)) +


## Distribution choice. auto2 is Poisson for OIIs, Tweedie for costs. auto uses quasipoisson instead of Poisson. Can be set to a distribution of your choice to force said distirbution (if it exists)
distribution.choice <- 'auto2' # 'auto' 'quasipoisson' 'Tweedie' # if use costs with quasipoisson, must use non-000 results



##################################################
### ADJUST THESE PARAMETERS TO ADJUST STAGE 2
##################################################

## Centering percentile
# cenpen <- NULL # Default, use minimum OII percentile
# cenpen <- 50 # Median
cenpen <- 'mean' # 'mean'
cenpen.minmax <- c(10,90) # max and min value for cenpen. Ignored if cenpen set manually # c(10,90) # c(25,75)


## Attributable numbers and fractions
nsim <- 5000 # Number of simulation runs for computing empirical CI
attrdl.dir <- 'forw' # Direction, either forw or back. back not allowed for reduced estimates
a.full <- F # if F, don't do mild and moderate heat and cold
af.round <- 2 # Round attributable fraction to number of decimal figures


## Graphical parameters
nline <- 0.8 # distance from title to plot for plots
tline <- 0.1 # distance from title to plot for condensed plots
xline2 <- 1.6 # location of second xaxis
tck.length <- -0.02 # tick length for condensed plots

gdpi <- 300 # usual journal quality is 300
glength <- 1600 # to maintain default scale, use 480*gdpi/72. 1600 if 300, 2100 if 400
glength.3by3 <- 1900 # slightly larger to fit 3*3 lag text. 1900 if 300, 2500 of 400

cexmain <- 1.25
cexlab <- 1.25
cexaxis <- 1

# y-axes (affect S2 graphs)
rr.yaxis <- 'auto'
bigvar.yaxis <- c(-15,60)
bigvarby <- 15
oer.yaxis <- seq(0.9,1.5,by=0.1)



##################################################
### Loop to run 03_AnalysisStage1.r and 04_AnalysisStage2.r for both outcome variables
### Alternatively, can change the outcome variable (outcome.var) manually above and re-run analysis
##################################################

for (outcome.var in c('Number of OIIs',"Total costs")) { 
  source('03_AnalysisStage1.r')
  source('04_AnalysisStage2.r')
  gc()
}



######################### END ############################
######################### END ############################
