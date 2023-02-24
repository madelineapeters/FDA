#Load required packages
library(deSolve)
library(CollocInfer)
library(fda)

#Set seed
set.seed(1234)

#Define vectors of names for the state variable and for the parameters:
LOGvarnames = c("N")
LOGparnames = c("r","K")

#To use lsoda, we need to specify initial conditions:
x0 = c(2500)
names(x0) = LOGvarnames

#Specify the values for the parameter vector theta = c(r,K):
LOGpars = c(1.0,10000)
names(LOGpars) = LOGparnames

#We need a function to evaluate the logistic equation:
log.ode = function(times,x,p){
  dx = x
  dimnames(dx) = dimnames(x)
  dx["N"] = (p["r"] * x["N"] * (p["K"] - x["N"])) / p["K"]
  return(list(dx))
}

#Set time vector
LOGtimes = seq(0,10,0.25)
LOGn = length(LOGtimes)

#Call function "lsoda" to compute the approximate solution at the time values specified above. The output should be an nx2 matrix where the first column is the evaluation times.
out = lsoda(x0,times=LOGtimes,log.ode,LOGpars)

#Add noise
data1 = out[,2] + 200*matrix(rnorm(LOGn), LOGn, 1)
colnames(data1) = "N"
data2 = out[,2] + 200*matrix(rnorm(LOGn), LOGn, 1)
colnames(data2) = "N"

#Create array and put both data sets (replicates) into the array
alldat = array(0,c(LOGn,2,1))
alldat[,1,] = data1
alldat[,2,] = data2

#Obtain three dimensional coefficient array, which we can get by smoothing the both data sets and putting the resulting coefficients together. First define the range of the basis and the breaks:
LOGrange = c(0,10)
breaks = seq(0,10,0.25)

#Then specify the order of the polynomial segments, which is one more than the highest power in the polynomial. We will use cubic polynomial splines.
LOGbasis = create.bspline.basis(LOGrange,norder=4,breaks=breaks)

#Obtain initial estimates of matrix of coefficients C by employing the "smooth.fd" function in fda
LOGfdPar = fdPar(LOGbasis, int2Lfd(2), 1)
LOGfd1 = smooth.basis(LOGtimes,data1,LOGfdPar)$fd
coefs1 = LOGfd1$coefs
colnames(coefs1) = "N"
LOGfd2 = smooth.basis(LOGtimes,data2,LOGfdPar)$fd
coefs2 = LOGfd2$coefs
colnames(coefs2) = "N"
coefs = array(0,c(dim(coefs1)[1],2,1))
coefs[,1,] = coefs1
coefs[,2,] = coefs2

#Set up function for use with CollocInfer
log.fun = function(times,x,p,more){
  dx = x
  dx[,"N"] = (p["r"] * x[,"N"] * (p["K"] - x[,"N"])) / p["K"]
  return(dx)
}

#Use LS.setup with the coefficients to get lik and proc objects
out = LS.setup(pars=LOGpars,coefs=coefs,basisvals=LOGbasis,fn=log.fun,lambda=100,times=LOGtimes,data=alldat,names=LOGvarnames)
lik = out$lik
proc = out$proc

#Get initial parameter estimates (doesn't work)
LOGParsMatchOpt = ParsMatchOpt(LOGpars,coefs,proc)
pars0 = LOGParsMatchOpt$pars
pars0

#Use our outputs above to call inneropt and outeropt
par.guess = c(1,10000)
names(par.guess) = LOGparnames
inner.res = inneropt(data=out$data,times=out$times,pars=par.guess,coefs=out$coefs,lik=lik,proc=proc)
outer.res = outeropt(data=out$data,times=out$times,pars=par.guess,coefs=inner.res$coefs,lik=lik,proc=proc)
outer.res$pars

outer.res$coefs

