library(deSolve)
library(CollocInfer)
library(fda)

#Set seed
set.seed(1234)

#Read in data
data("ChemoRMData")

#Store time and data
time = ChemoRMTime
data = log(ChemoRMData)

#Rosenzweig-MacArthur model (1969)
RosMac = function(t, x, p, more) {
  p = exp(p); x = exp(x)
  dx = x
  dx[,"C"] = p["rho"] * (1 - x[, "C"] / p["kappaC"]) - p["gamma"] * x[, "B"] / (p["kappaB"] + x[, "C"])
  dx[, "B"] <- p["chi"] * p["gamma"] * x[, "C"] /
    (p["kappaB"] + x[, "C"]) - p["delta"]
  return(dx)
}

#Define initial parameters
varnames = RMvarnames
varnames
parnames = RMparnames
parnames

#Set up basis functions
rr = range(time)
breaks = seq(rr[1],rr[2],by=1)

ChemoRMbasis = create.bspline.basis(rr,norder=4,breaks=breaks)

#Set up initial coefficients and parameters
#Coefficients
fd0 = smooth.basis(time, data, fdPar(ChemoRMbasis, int2Lfd(2), 10))
coef0 = fd0$fd$coef
colnames(coef0) = varnames
#Parameters
pars = log(ChemoRMPars)
names(pars) = parnames

#There are obtained estimates of rho and kappaC from experiments without rotifer predation, so we will fix these by defining active pars are the indices of the parameters in p that we wish to estimate. We can call activepars in the active argument in the ParsMatchOpt and outeropt, as well as Profile.LS.
activepars = 3:6

#Set up profiling objects
out = LS.setup(pars = pars, coefs=coef0, basisvals=ChemoRMbasis, fn=RosMac, lambda=c(5e4,5e2), times=time)
lik = out$lik
proc = out$proc

#Obtain initial parameter estimates. The "active" elements tell ParsMatchOpt which elements of pars to estimate.
res1 = ParsMatchOpt(pars,coef0,proc,active=activepars)

#Run profiling procedure
res2 = outeropt(data, time, res1$pars, coef0, lik, proc, active=activepars)

#Construct diagnostic plots
out2 = CollocInferPlots(res2$coefs, res2$pars, lik, proc, times=time, data=data)

#Obtain confidence intervals
covar = Profile.covariance(res2$pars, times = time, data = data, coefs = res2$coefs, lik = lik, proc = proc, active = activepars)
CIs = cbind(res2$pars[activepars] - 2 * sqrt(diag(covar)),res2$pars[activepars] + 2 * sqrt(diag(covar)))
rownames(CIs) = parnames[activepars]
exp(CIs)

#FPE procedure with 5-day look-ahead window
whichtimes = cbind(1:102,6:107)
lambdafac = c(0.1,0.5,1,5,10)
FPEs = 0*lambdafac
temp.pars = res2$pars
temp.coefs = res2$coefs

for (ilam in 1:length(lambdafac)){
  t.res = Profile.LS(RosMac,data,time,temp.pars,temp.coefs,ChemoRMbasis,lambdafac[ilam]*c(5e4,5e2),active=activepars, out.meth="nlminb")
  FPEs[ilam] = forward.prediction.error(time,data,t.res$coefs,lik,proc,t.res$pars,whichtimes)
  temp.pars = t.res$pars
  temp.coefs = t.res$coefs
}
FPEs

####################
#### Extensions ####
####################

#Code the model for use with CollocInfer:
RosMac2 = function(t, x, p, more) {
  p = exp(p)
  dx = x
  dx[, "C1"] = p["rho1"] * x[, "C1"] * (1 - x[, "C1"] / p["kappaC1"] - x[, "C2"] / p["kappaC2"]) - p["pi"] * p["gamma"] * x[, "C1"] * x[, "B"] /
    (p["kappaB"] + p["pi"] * x[, "C1"] + x[, "C2"])
  dx[, "C2"] = p["rho2"] * x[, "C2"] *
    (1 - x[, "C1"] / p["kappaC1"] - x[, "C2"] / p["kappaC2"]) - p["gamma"] * x[, "C2"] * x[, "B"] / (p["kappaB"] + p["pi"] * x[, "C1"] + x[, "C2"])
  dx[, "B"] = p["chi"] * p["gamma"] * (p["pi"] * x[, "C1"] + x[, "C2"]) * x[, "B"] / (p["kappaB"] + p["pi"] * x[, "C1"] + x[, "C2"]) - p["delta"] * x[, "B"]
  return(dx)
}

#We've hardcoded that we will estimate log parameters, since all parameters should be positive. To generate data:
RMpars = c(0.2,0.025,0.125,2.2e4,1e5,5e6,1,1e9,0.3)
RMParnames = c("pi","rho1","rho2","kappaC1","kappaC2","gamma","chi","kappaB","delta")
logpars = log(RMpars)
names(logpars) = RMParnames

#Initial conditions
RMVarnames = c("C1","C2","B")
x0 = c(50,50,2)
names(x0) = RMVarnames


##############################
#### Log transformed data ####
##############################
#All state variables are strictly positive, so it may be useful to model z = log x. The model is implemented by incorporating the log transform in a form suitable for use with "lsoda":
RosMac2ODE = function(t,z,p){
  p = exp(p)
  x = exp(z)
  dx = x
  dx["C1"] = p["rho1"] * x["C1"] * (1 - x["C1"] / p["kappaC1"] - x["C2"] / p["kappaC2"]) - p["pi"] * p["gamma"] * x["C1"] * x["B"]/(p["kappaB"] + p["pi"] * x["C1"] + x["C2"])
  dx["C2"] = p["rho2"] * x["C2"] * (1 - x["C2"] / p["kappaC2"] - x["C1"] / p["kappaC1"]) - p["gamma"] * x["C2"] * x["B"] / (p["kappaB"] + p["pi"] * x["C1"] + x["C2"])
  dx["B"] = p["chi"] * p["gamma"] * (p["pi"] * x["C1"] + x["C2"]) * x["B"] / (p["kappaB"] + p["pi"] * x["C1"] + x["C2"]) - p["delta"] * x["B"]
  return(list(dx / x))
}

#Now we generate data by solving the ODE and adding noise. We'll pretend we can measure each of the three species and that the experiment lasts 200 days and is sampled daily.
time = 0:200
res0 = lsoda(y=log(x0),times=time,func=RosMac2ODE,parms=logpars)
head(res0)
matplot(exp(res0[,2:4]), type="l")
#Add noise
data = res0[,2:4] + 0.2*matrix(rnorm(length(time)*3),201,3)
matplot(data)

#Set up basis expansions, smooth data
rr = range(time)
breaks = seq(rr[1],rr[2],by=1)
RMbasis = create.bspline.basis(rr,breaks=breaks)

#Obtain initial set of coefficients
coef0 = smooth.basis(time, data, fdPar(RMbasis,int2Lfd(2),10))$fd$coef
colnames(coef0) = RMVarnames

#To set up our equations on the log scale, we could transform RosMac2ODE into a form appropriate for CollocInfer, but this requires writing out a new function. So we'll use posproc = TRUE option in LS.setup, which makes the transformation automatically. 
#Our data is already measured on the log-scale, so we don't need the corresponding poslik = TRUE option
out = LS.setup(pars=logpars,coefs=coef0,basisvals=RMbasis,fn=RosMac2,lambda=1e5,times=time,posproc=TRUE)
lik = out$lik
proc = out$proc

#Get initial parameter estimates
res1 = ParsMatchOpt(logpars,coef0,proc)

#Profile
res3 = outeropt(data, time, res1$pars, coef0, lik, proc)
exp(res3$pars)

#####################################
#### Indirectly observed systems ####
#####################################
#We'll mimic not being able to distinguish C1 and C2
data2 = cbind(log(exp(data[, "C1"]) + exp(data[, "C2"])), data[, "B"])

#We need a function to go from states to observations. The structure of this function should look exactly like the functions we employ for the right hand side of an ODE. This function is added to the likfn argument to LS.setup, smooth.LS or profile.LS
RMobsfn = function(t, x, p, more) {
  x = exp(x)
  y = cbind(x[,"C1"] + x[,"C2"], x[,"B"])
  return(log(y))
}

#Get proc and lik objects
out = LS.setup(pars=logpars,coefs=coef0,basisvals=RMbasis,fn=RosMac2,times=time,lambda=1e5,posproc=TRUE,likfn=RMobsfn)
lik2 = out$lik
proc2 = out$proc

#We can't start from coef0 as a pre-smooth since we do not have separate data from C1 and C2, so we set the corresponding coefficients to zero
coef02 = coef0
coef02[,1:2] = 0

#Get initial parameter estimates and then profile
Fres3 = FitMatchOpt(coef02,1:2,res1$pars,proc2)
res32 = outeropt(data2,time,res1$pars,Fres3$coefs,lik2,proc2)
exp(res32$pars)

