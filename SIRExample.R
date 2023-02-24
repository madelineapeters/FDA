library(deSolve)
library(CollocInfer)
library(fda)

set.seed(1234)
#################################
#### Example data generation ####
#################################
#Define vectors of names for the indices of arrays containing parameters and values of the trajectories:
SIRvarnames = c("S","I","R")
SIRparnames = c("alpha","beta","kappa")

#To use lsoda, we need to specify initial conditions, respectively:
x0 = c(10000,1,0)
names(x0) = SIRvarnames

#Specify the initial values for the parameter vector theta = c(a,b,c)
SIRpars = c(0.4,0.001,0.5)
names(SIRpars) = SIRparnames

#We need a function to evalute the SIR equations. The function should return a names list, with a member names fn and containing the function evaluation
sir.ode = function(times,x,p){
  dx = x
  dimnames(dx) = dimnames(x)
  dx["S"] = p["alpha"] * x["R"] - p["beta"] * x["S"] * x["I"]
  dx["I"] = p["beta"] * x["S"] * x["I"] - p["kappa"] * x["I"]
  dx["R"] = p["kappa"] * x["I"] - p["alpha"] * x["R"]
  return(list(dx))
}

#Now define a time vector from 0 to 20
sirplottimes = seq(0,10,0.25)

#Call function lsoda to compute the approximate solution at the time values specified above. The output should be an nx4 matrix where the first column is the evaluation times.
out = lsoda(x0,times=sirplottimes,sir.ode,SIRpars)

#Plotting out should show SIR dynamics:
matplot(out[,1],out[,2:4], type="l")
legend("topright",c("S","I","R"))

#Now we will generate some data with noise drawn from a normal distribution with mean zero and standard deviation 200.
SIRtimes = seq(0,10,0.25)
SIRn = length(SIRtimes)
out = lsoda(x0,time=SIRtimes,sir.ode,SIRpars)
SIRdata = out[,2:4] + 200*matrix(rnorm(3*SIRn), SIRn, 3)
SIRdata[SIRdata < 0] = 0
SIRdata = round(SIRdata)
matplot(SIRdata[,1:3], type="l",xlab="Time",ylab="S,I,R")
#######################################################
#### Basis functions for representing trajectories ####
#######################################################
#First define the range of the basis and the breaks
SIRrange = c(0,10)
breaks = seq(0,10,0.25)

#Then specify the order of the polynomial segments, which is one more than the highest power in the polynomial. We will use cubic polynomial splines. 
SIRbasis = create.bspline.basis(SIRrange,norder=4,breaks=breaks)

#Obtain initial estimates of matrix of coefficients C by employing the smooth.fd function in fda. We assume a lambda of 1.
SIRfdPar = fdPar(SIRbasis, int2Lfd(2), 1)

#Can smooth the data without enforcing positive values
DEfd0 = smooth.basis(SIRtimes,SIRdata,SIRfdPar)$fd

#To obtain the coefficients of this smooth:
coefs0 = DEfd0$coefs

#Examine the fit of the smoothed, positive data to the data
plotfit.fd(SIRdata, SIRtimes, DEfd0)

#Can smooth data but maintain positive values and obtain coefficients (to be used) (doesn't seem to work)
SIRfd.Pos = smooth.pos(SIRtimes,SIRdata,SIRfdPar)
SIRWfd = SIRfd.Pos$Wfdobj
#coefs0 = SIRWfd

###########################
#### Using CollocInfer ####
###########################
#To obtain parameter estimates from CollocInfer, we need 1) a basis system, 2) a function to evaluate and 3) initial values to start off the iterative estimation of theta and matrix C

#First set up the function to evaluate f(x;t,theta). The package requires that f is a function with arguments: t (a vector of evaluation times), x (a matrix of evaluations of x(t) at the times in t with rows corresponding to time points), p (a vector of parameters) and more (a list of additional inputs f might require).

#For the SIR equations, the function is:
sir.fun = function(times,x,p,more){
  dx = x
  dx[,"S"] = p["alpha"] * x[,"R"] - p["beta"] * x[,"S"] * x[,"I"]
  dx[,"I"] = p["beta"] * x[,"S"] * x[,"I"] - p["kappa"] * x[,"I"]
  dx[,"R"] = p["kappa"] * x[,"I"] - p["alpha"] * x[,"R"]
  return(dx)
}

#We require a basis system, which we created above as SIRbasis

#Set up some required structures. There is a function for the common choice of sum of squared errors criteron. LS.setup returns a single list object with two names members, each of the list classes lik and proc that define the map from x to y and the map from x to Dx respectively.
lambda = 100 #just chosen for now, but can be chosen using methods
profile.obj = LS.setup(pars=SIRpars,fn=sir.fun,lambda=lambda,times=SIRtimes,coefs=coefs0,basisvals=SIRbasis)

proc = profile.obj$proc
lik = profile.obj$lik

#We obtained initial coefficients using the smoothing commands from the fda package. Next we will use gradient matching to get an initial estimate of parameters. The function ParsMatchOpt minimises ISSE(theta,x-hat) and returns a named list containing the optimized parameter values
Pres0 = ParsMatchOpt(SIRpars,coefs0,proc)
pars1 = Pres0$pars
pars1

#####################################################################
#### Running generalized profiling to obtain parameter estimates ####
#####################################################################
#A first step in generalized profiling routines is to find x that minimises the inner criterion PENSSE(x,theta). This is achieved by inneropt:
Ires1 = inneropt(SIRdata,times=SIRtimes,pars1,coefs0,lik,proc)
coefs1 = Ires1$coefs

#Now carry out profiling itself:
Ores2 = outeropt(SIRdata,SIRtimes,pars1,coefs1,lik,proc)

#Output is a list with members pars, coefs and res giving the optimum parameters and coefficients. The estimated parameter vector is:
Ores2$pars

#############################################################
#### Using FPE to choose the smoothing parameters lambda ####
#############################################################
#For use forwards prediction error, we use the function forward.prediction.error. Besides the results of profiling, we need to define starting and ending times. These are given in a matrix where it is the index of the given observation times. Here we will start from each time point and predict 5 ahead:
whichtimes = cbind(1:5,6:10)

#Now we'll cycle through values of lambda increasing as powers of 10, report the resulting parameter estimates and overlay the estimated trajectory on a plot of the data. To make things more computationally efficient, we'll update the starting value for the parameters and coefficients to be the estimate of the previous lambda.
matplot(SIRtimes,SIRdata,pch=SIRvarnames)
lambdas = 10^(0:7)
FPEs = 0*lambdas
temp.pars = SIRpars
temp.coefs = coefs0
#NOTE: Get error with forward.prediction.error for lambda = 100 or higher, but plotting looks like lambda = 100 is smallest lambda that produces smooth fit to data and accurate parameter estimates
for (ilam in 1:length(lambdas)){
  print(paste("lambda = ", lambdas[ilam]))
  t.Ores = Profile.LS(sir.fun, SIRdata, SIRtimes, temp.pars, temp.coefs, SIRbasis, lambdas[ilam])
  print(t.Ores$pars)
  temp.pars = t.Ores$pars
  temp.coefs = t.Ores$coefs
  temp.lik = t.Ores$lik
  temp.proc = t.Ores$proc
  t.fd = fd(t.Ores$coefs,SIRbasis)
  lines(t.fd,lwd=2,col=ceiling(ilam/2),lty=1)
  FPEs[ilam] = forward.prediction.error(SIRtimes,SIRdata,t.Ores$coefs,temp.lik,temp.proc,t.Ores$pars,whichtimes)
}
FPEs

##################################################################
#### Getting final coefficients and final parameter estimates ####
##################################################################
#Select lambda based on FPEs
lambda = 10^2
final.profile = Profile.LS(sir.fun, SIRdata, SIRtimes, temp.pars, temp.coefs, SIRbasis,lambda)
final.pars = final.profile$pars
final.coefs = final.profile$coefs
final.lik = final.profile$lik
final.proc = final.profile$proc

########################################
#### Calculate confidence intervals ####
########################################
covar = Profile.covariance(final.pars, times = SIRtimes, data = SIRdata, coefs = final.coefs, lik = final.lik, proc = final.proc)
covar
CIs = cbind(final.pars - 2 * sqrt(diag(covar)),final.pars + 2 * sqrt(diag(covar)) )
rownames(CIs) = SIRparnames
CIs

#####################################################
#### Making data and initial coefficients for missing R ####
#####################################################
#Set all values for R to 0
data2 = SIRdata
data2[,3] = 0

#Set the coefficients for R to 0
coefs0.2 = coefs0
coefs0.2[,3] = 0

#Choose initial parameter values
par.guess = pars1 #not a perfect guess

#Get proc and lik (lik not needed)?
lambda = 100
profile.obj = LS.setup(pars=par.guess, fn=sir.fun, lambda=lambda, times=SIRtimes, coefs=coefs0.2, basisvals=SIRbasis)

proc = profile.obj$proc
lik = profile.obj$lik

#Try to estimate the third column of coefs0.2 so that the whole trajectory matches our function as well as possible with the first two columns being held fixed. In other words, optimized ISSE(theta,x) but over some components of x rather than over theta.
Fres3 = FitMatchOpt(coefs0.2,3,par.guess,proc) #how do we have a proc?

#Proceed with profiling
Ores4 = outeropt(SIRdata,SIRtimes,par.guess,Fres3$coefs,lik,proc)
Ores4$pars #parameter estimates are not good

#Try Profile.LS instead
lambda = 10^2
miss.profile = Profile.LS(sir.fun, SIRdata, SIRtimes, par.guess, Fres3$coefs, SIRbasis,lambda)
miss.pars = miss.profile$pars
miss.coefs = miss.profile$coefs
miss.proc = miss.profile$proc
miss.lik = miss.profile$lik

#Get CIs
covar = Profile.covariance(miss.pars, times = SIRtimes, data = SIRdata, coefs = miss.coefs, lik = miss.lik, proc = miss.proc)
covar
CIs = cbind(miss.pars - 2 * sqrt(diag(covar)),miss.pars + 2 * sqrt(diag(covar)) )
rownames(CIs) = SIRparnames
CIs


