library(deSolve)
library(CollocInfer)
library(fda)
#################################
#### Example data generation ####
#################################
#Define vectors of names for the indices of arrays containing parameters and values of the trajectories:
FhNvarnames = c("V","R")
FhNparnames = c("a","b","c")

#To use lsoda, we need to specify initial conditions, respectively:
x0 = c(-1,1)
names(x0) = FhNvarnames

#Specify the initial values for the parameter vector theta = c(a,b,c)
FhNpars = c(0.2,0.2,3)
names(FhNpars) = FhNparnames

#We need a function to evalute the FitzHugh-Nagumo equations. The function should return a names list, with a member names fn and containing the function evaluation
fhn.ode = function(times,x,p){
  dx = x
  dimnames(dx) = dimnames(x)
  dx["V"] = p["c"] * (x["V"] - x["V"]^3 / 3 + x["R"])
  dx["R"] = - (x["V"] - p["a"] + p["b"] * x["R"]) / p["c"]
  return(list(dx))
}

#Now define a time vector from 0 to 20
FhNplottimes = seq(0,20,0.05)

#Call function lsoda to compute the approximate solution at the time values specified above. The output should be an nx3 matrix where the first column is the evaluation times.
out = lsoda(x0,times=FhNplottimes,fhn.ode,FhNpars)

#Plotting out should obtain a wave:
matplot(out[,1],out[,2:3], type="l")
legend("bottomleft",c("V","R"))

#Now we will generate some data with noise drawn from a normal distribution with mean zero and standard deviation 0.1.
FhNtimes = seq(0,20,0.5)
FhNn = length(FhNtimes)
out = lsoda(x0,time=FhNtimes,fhn.ode,FhNpars)
FhNdata = out[,2:3] + 0.1*matrix(rnorm(2*FhNn), FhNn, 2)

#######################################################
#### Basic functions for representing trajectories ####
#######################################################
#First define the range of the basis and the breaks
FhNrange = c(0,20)
breaks = seq(0,20,0.5)

#Then specify the order of the polynomial segments, which is one more than the highest power in the polynomial. We will use cubic polynomial splines. 
FhNbasis = create.bspline.basis(FhNrange,norder=4,breaks=breaks)

#Obtain initial estimates of matrix of coefficients C by employing the smooth.fd function in fda. We assume a lambda of 1.
FhNfdPar = fdPar(FhNbasis, int2Lfd(2), 1)

#Smooth the data
DEfd0 = smooth.basis(FhNtimes,FhNdata,FhNfdPar)$fd

#Examine the fit of the smoothed data to the data
plotfit.fd(FhNdata, FhNtimes, DEfd0)

#To obtain the coefficients of this smooth, which will be used later:
coefs0 = DEfd0$coefs

###########################
#### Using CollocInfer ####
###########################
#To obtain parameter estimates from CollocInfer, we need 1) a basis system, 2) a function to evaluate and 3) initial values to start off the iterative estimation of theta and matrix C

#First set up the function to evaluate f(x;t,theta). The package requires that f is a function with arguments: t (a vector of evaluation times), x (a matrix of evaluations of x(t) at the times in t with rows corresponding to time points), p (a vector of parameters) and more (a list of additional inputs f might require).

#For the FitzHugh-Nagumo equations, the function is:
fhn.fun = function(times,x,p,more){
  dx = x
  dx[,"V"] = p["c"] * (x[,"V"] - x[,"V"]^3 / 3 + x[,"R"])
  dx[,"R"] = - (x[,"V"] - p["a"] + p["b"] * x[,"R"]) / p["c"]
  return(dx)
}

#We require a basis system, which we'll recreate here
FhNrange = c(0,20)
breaks = seq(0,20,0.5)
FhNbasis = create.bspline.basis(FhNrange,norder=4,breaks=breaks)

#We require the following input: data (an nxd matrix of data values, with missing values listed as NA), times (a vector of n observation times)

#Set up some required structures. There is a function for the common choice of sum of squared errors criteron. LS.setup returns a single list object with two names members, each of the list classes lik and proc that define the map from x to y and the map from x to Dx respectively.
lambda = 1000 #just chosen for now, but can be chosen using methods
profile.obj = LS.setup(pars=FhNpars,fn=fhn.fun,lambda=lambda,times=FhNtimes,coefs=coefs0,basisvals=FhNbasis)

proc = profile.obj$proc
lik = profile.obj$lik

#We obtained initial coefficients using the smoothing commands from the fda package. Next we will use gradient matching to get an initial estimate of parameters. The function ParsMatchOpt minimises ISSE(theta,x-hat) and returns a named list containing the optimized parameter values
Pres0 = ParsMatchOpt(FhNpars,coefs0,proc)
pars1 = Pres0$pars
pars1
#compare to true values for a,b,c=c(0.2,0.2,3)
#####################################################################
#### Running generalized profiling to obtain parameter estimates ####
#####################################################################
#A first step in generalized profiling routines is to find x that minimises the inner criterion PENSSE(x,theta). This is achieved by inneropt:
Ires1 = inneropt(FhNdata,times=FhNtimes,pars1,coefs0,lik,proc)
coefs1 = Ires1$coefs

#Now carry out profiling itself:
Ores2 = outeropt(FhNdata,FhNtimes,pars1,coefs1,lik,proc)

#Output is a list with members pars, coefs and res giving the optimum parameters and coefficients. The estimated parameter vector is:
Ores2$pars

#############################################################
#### Using FPE to choose the smoothing parameters lambda ####
#############################################################
#For use forwards prediction error, we use the function forward.prediction.error. Besides the results of profiling, we need to define starting and ending times. These are given in a matrix where it is the index of the given observation times. Here we will start from each time point and predict 10 ahead:
whichtimes = cbind(1:31,11:41)

#Now we call the FPE function:
FPE = forward.prediction.error(FhNtimes, FhNdata, Ores2$coefs, lik, proc, whichtimes)

#Now we'll cycle through values of lambda increasing as powers of 10, report the resulting parameter estimates and overlay the estimated trajectory on a plot of the data. To make things more computationally efficient, we'll update the starting value for the parameters and coefficients to be the estimate of the previous lambda.
matplot(FhNtimes,FhNdata,pch=FhNvarnames)
lambdas = 10^(0:7)
FPEs = 0*lambdas
temp.pars = FhNpars
temp.coefs = coefs0
for (ilam in 1:length(lambdas)){
  print(paste("lambda = ", lambdas[ilam]))
  t.Ores = Profile.LS(fhn.fun, FhNdata, FhNtimes, temp.pars, temp.coefs, FhNbasis, lambdas[ilam])
  print(t.Ores$pars)
  temp.pars = t.Ores$pars
  temp.coefs = t.Ores$coefs
  temp.lik = t.Ores$lik
  temp.proc = t.Ores$proc
  t.fd = fd(t.Ores$coefs,FhNbasis)
  lines(t.fd,lwd=2,col=ceiling(ilam/2),lty=1)
  FPEs[ilam] = forward.prediction.error(FhNtimes,FhNdata,t.Ores$coefs,temp.lik,temp.proc,t.Ores$pars,whichtimes)
}
FPEs

################################################
#### Using CollocInfer for diagnostic plots ####
################################################
#There are two things we can look at: how well the estimated trajectory x(t,theta) fits the data and how well Dx(t,theta) matches f(x,theta)
out1 = CollocInferPlots(Ores2$coefs, Ores2$pars,lik,proc,times=FhNtimes,data=FhNdata)

####################################
#### Estimating sample variance ####
####################################
#To evaluate statistical uncertainty, we use a Newey-West estimate of the covariance of the score function for the outer optimisation.
covar = Profile.covariance(Ores2$pars, times=FhNtimes, data=FhNdata, coefs=Ores2$coefs, lik=lik, proc=proc)
covar
#To construct confidence intervals, we add and subtract two standard deviations
CIs = cbind(Ores2$pars - 2 * sqrt(diag(covar)),Ores2$pars + 2*sqrt(diag(covar)))
rownames(CIs) = FhNparnames
CIs
