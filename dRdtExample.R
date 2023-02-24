################################
#### Load required packages ####
################################
library(tidyverse)
library(deSolve)
library(fda)
########################################
#### Simulate data with step-wise B ####
########################################
#Define function that descibes ODE of changes in R as a function of time
dRdt = function(t, x, parms) {
  with(as.list(c(parms, x)), {
    if (t<5){B=10000}else{B=500000}
    dR = B - R*mu
    res = c(dR)
    list(res,B)
  })
}
# Specify initial conditions and parameter values
y0 = c(R = 5000000)
pars = c(mu=0.025)
#Specify timescale of interest
tout=seq(from=0,to=10,by=0.5)
#ODE integration and view head of data
out=ode(y=y0, times=tout, func=dRdt, parms=pars)
head(out)
#Convert output to a dataframe
out = as.data.frame(out)
#Create measured R values from true R values and normal distribution
out$Rmeas = out$R
#out$Rmeas = sapply(out$R,FUN=function(x){rnorm(1,x,200000)})
#Plot measured R values (points) and true R values (line)
plot(out$time,out$Rmeas,)
lines(out$time,out$R)

###########################################
#### Create FDA object for step-wise B ####
###########################################
Routput = as.matrix(out)
rMeasData = as.matrix(out[,4])  #transpose here if necessary (should be long)

#create vector of times (with one for each day)
daytime  = seq(from=0,to=10,by=0.5)

#Plot the number of cases over time
plot(daytime, rMeasData, "b", lwd=2,
     xlab="Time", ylab="R measured", cex=1.2)

#Create a B-spline basis, where the time interval is 1 to 15, there are 8 basis functions and the order is 4 (implicit)
rMeasbasis    = create.bspline.basis(c(0,10),10)

#Evaluate the basis at the vector of times (days) specified above
rMeasbasismat = eval.basis(rMeasbasis, daytime)

#Compute coefficients for the functional data object by the usual equations for regression coefficients, b = (X′X)^(−1) X′y.
rMeascoef = solve(crossprod(rMeasbasismat),
                 crossprod(rMeasbasismat,rMeasData))
#Can also fit data by least squares
rMeascoef = lsfit(rMeasbasismat, rMeasData, intercept=FALSE)$coef

#Create functional data object using the coefficients, the basis and a list of names
rMeasfd   = fd(rMeascoef, rMeasbasis,
              list("Time", "Individual", "R measured"))
plot(rMeasfd, lty=1, lwd=2, col=1)

#Plot the fit of the FD object to the actual data
plotfit.fd(rMeasData, daytime, rMeasfd,
           lty=1, lwd=2, main='')

##############################################
#### Apply linear operator for step-wise B####
##############################################
#Define linear operator that's just the first derivative
velLfd = int2Lfd(1)
#Evaluate the FD object at vector daytime and apply linear operator defined above
rMeasVelfd = eval.fd(daytime,rMeasfd,velLfd)
#Calculate first derivative of FD object
D2tempfd = deriv.fd(rMeasfd, 1)
#Compare output of linear operator evaluation versus first derivative function
plot(daytime,rMeasVelfd)
lines(D2tempfd)

#Create FD object multiplied by constant mu
mufd = 0.025*rMeasfd
#Plot what should be B, i.e., dRdt + mu*R; lay over step function that is the true B
plot(mufd+D2tempfd,ylab="Estimated B",xlab="Day")
lines(c(0,4.99,5,10),c(10000,10000,500000,500000))

#######################################
#### Simulate data with constant B ####
#######################################
#Define function that descibes ODE of changes in R as a function of time
dRdt = function(t, x, parms) {
  with(as.list(c(parms, x)), {
    dR = B - R*mu
    res = c(dR)
    list(res,B)
  })
}
# Specify initial conditions and parameter values
y0 = c(R = 5000000)
pars = c(mu=0.025,B=300000)
#Specify timescale of interest
tout=seq(from=0,to=10,by=0.5)
#ODE integration and view head of data
out=ode(y=y0, times=tout, func=dRdt, parms=pars)
head(out)
#Convert output to a dataframe
out = as.data.frame(out)
#Create measured R values from true R values and normal distribution
out$Rmeas = out$R
#out$Rmeas = sapply(out$R,FUN=function(x){rnorm(1,x,200000)})
#Plot measured R values (points) and true R values (line)
plot(out$time,out$Rmeas,)
lines(out$time,out$R)

###########################################
#### Create FDA object for constant B ####
###########################################
Routput = as.matrix(out)
rMeasData = as.matrix(out[,4])  #transpose here if necessary (should be long)

#create vector of times (with one for each day)
daytime  = seq(from=0,to=10,by=0.5)

#Plot the number of cases over time
plot(daytime, rMeasData, "b", lwd=2,
     xlab="Time", ylab="R measured", cex=1.2)

#Create a B-spline basis, where the time interval is 1 to 15, there are 8 basis functions and the order is 4 (implicit)
rMeasbasis    = create.bspline.basis(c(0,10),10)

#Evaluate the basis at the vector of times (days) specified above
rMeasbasismat = eval.basis(rMeasbasis, daytime)

#Compute coefficients for the functional data object by the usual equations for regression coefficients, b = (X′X)^(−1) X′y.
rMeascoef = solve(crossprod(rMeasbasismat),
                  crossprod(rMeasbasismat,rMeasData))
#Can also fit data by least squares
rMeascoef = lsfit(rMeasbasismat, rMeasData, intercept=FALSE)$coef

#Create functional data object using the coefficients, the basis and a list of names
rMeasfd   = fd(rMeascoef, rMeasbasis,
               list("Time", "Individual", "R measured"))
plot(rMeasfd, lty=1, lwd=2, col=1)

#Plot the fit of the FD object to the actual data
plotfit.fd(rMeasData, daytime, rMeasfd,
           lty=1, lwd=2, main='')

##############################################
#### Apply linear operator for constant B####
##############################################
#Define linear operator that's just the first derivative
velLfd = int2Lfd(1)
#Evaluate the FD object at vector daytime and apply linear operator defined above
rMeasVelfd = eval.fd(daytime,rMeasfd,velLfd)
#Calculate first derivative of FD object
D2tempfd = deriv.fd(rMeasfd, 1)
#Compare output of linear operator evaluation versus first derivative function
plot(daytime,rMeasVelfd)
lines(D2tempfd)

#Create FD object multiplied by constant mu
mufd = 0.025*rMeasfd
#Plot what should be B, i.e., dRdt + mu*R; lay over step function that is the true B
plot(mufd+D2tempfd,ylab="Estimated B",xlab="Day",ylim=c(200000,400000))
abline(a=300000,b=0,col="red")



########################################
#### Simulate data with parabolic B ####
########################################
#Define function that descibes ODE of changes in R as a function of time
dRdt = function(t, x, parms) {
  with(as.list(c(parms, x)), {
    B=-(5000)*(t-5)^2+300000
    dR = B - R*mu
    res = c(dR)
    list(res,B)
  })
}
# Specify initial conditions and parameter values
y0 = c(R = 5000000)
pars = c(mu=0.025)
#Specify timescale of interest
tout=seq(from=0,to=10,by=0.5)
#ODE integration and view head of data
out=ode(y=y0, times=tout, func=dRdt, parms=pars)
head(out)
#Convert output to a dataframe
out = as.data.frame(out)
#Create measured R values from true R values and normal distribution
out$Rmeas = out$R
#out$Rmeas = sapply(out$R,FUN=function(x){rnorm(1,x,200000)})
#Plot measured R values (points) and true R values (line)
plot(out$time,out$Rmeas,)
lines(out$time,out$R)

###########################################
#### Create FDA object for parabolic B ####
###########################################
Routput = as.matrix(out)
rMeasData = as.matrix(out[,4])  #transpose here if necessary (should be long)

#create vector of times (with one for each day)
daytime  = seq(from=0,to=10,by=0.5)

#Plot the number of cases over time
plot(daytime, rMeasData, "b", lwd=2,
     xlab="Time", ylab="R measured", cex=1.2)

#Create a B-spline basis, where the time interval is 1 to 15, there are 8 basis functions and the order is 4 (implicit)
rMeasbasis    = create.bspline.basis(c(0,10),10)

#Evaluate the basis at the vector of times (days) specified above
rMeasbasismat = eval.basis(rMeasbasis, daytime)

#Compute coefficients for the functional data object by the usual equations for regression coefficients, b = (X′X)^(−1) X′y.
rMeascoef = solve(crossprod(rMeasbasismat),
                  crossprod(rMeasbasismat,rMeasData))
#Can also fit data by least squares
rMeascoef = lsfit(rMeasbasismat, rMeasData, intercept=FALSE)$coef

#Create functional data object using the coefficients, the basis and a list of names
rMeasfd   = fd(rMeascoef, rMeasbasis,
               list("Time", "Individual", "R measured"))
plot(rMeasfd, lty=1, lwd=2, col=1)

#Plot the fit of the FD object to the actual data
plotfit.fd(rMeasData, daytime, rMeasfd,
           lty=1, lwd=2, main='')
##############################################
#### Apply linear operator for parabolic B####
##############################################
#Define linear operator that's just the first derivative
velLfd = int2Lfd(1)
#Evaluate the FD object at vector daytime and apply linear operator defined above
rMeasVelfd = eval.fd(daytime,rMeasfd,velLfd)
#Calculate first derivative of FD object
Rmeasfd.D1 = deriv.fd(rMeasfd, 1)
#Compare output of linear operator evaluation versus first derivative function
plot(daytime,rMeasVelfd)
lines(Rmeasfd.D1)

#Create FD object multiplied by constant mu
mufd = 0.025*rMeasfd
#Plot what should be B, i.e., dRdt + mu*R; lay over step function that is the true B
plot(mufd+Rmeasfd.D1,ylab="Estimated B",xlab="Day")
lines(seq(0,10,0.5),sapply(seq(0,10,0.5),FUN=function(x){-(5000)*(x-5)^2+300000}),col="red")

###################################################
#### Simulate data with parabolic B with error ####
###################################################
#Define function that descibes ODE of changes in R as a function of time
dRdt = function(t, x, parms) {
  with(as.list(c(parms, x)), {
    B=-(5000)*(t-5)^2+300000
    dR = B - R*mu
    res = c(dR)
    list(res,B)
  })
}
# Specify initial conditions and parameter values
y0 = c(R = 5000000)
pars = c(mu=0.025)
#Specify timescale of interest
tout=seq(from=0,to=10,by=0.5)
#ODE integration and view head of data
out=ode(y=y0, times=tout, func=dRdt, parms=pars)
head(out)
#Convert output to a dataframe
out = as.data.frame(out)
#Create measured R values from true R values and normal distribution
#out$Rmeas = out$R
out$Rmeas = sapply(out$R,FUN=function(x){rnorm(1,x,100000)})
#Plot measured R values (points) and true R values (line)
plot(out$time,out$Rmeas)
lines(out$time,out$R)

######################################################
#### Create FDA object for parabolic B with error ####
######################################################
Routput = as.matrix(out)
rMeasData = as.matrix(out[,4])  #transpose here if necessary (should be long)

#create vector of times (with one for each day)
daytime  = seq(from=0,to=10,by=0.5)

#Plot the number of cases over time
plot(daytime, rMeasData, "b", lwd=2,
     xlab="Time", ylab="R measured", cex=1.2)

#Create a B-spline basis, where the time interval is 1 to 15, there are 8 basis functions and the order is 4 (implicit)
rMeasbasis    = create.bspline.basis(c(0,10),4)

#Evaluate the basis at the vector of times (days) specified above
rMeasbasismat = eval.basis(rMeasbasis, daytime)

#Compute coefficients for the functional data object by the usual equations for regression coefficients, b = (X′X)^(−1) X′y.
rMeascoef = solve(crossprod(rMeasbasismat),
                  crossprod(rMeasbasismat,rMeasData))
#Can also fit data by least squares
rMeascoef = lsfit(rMeasbasismat, rMeasData, intercept=FALSE)$coef

#Create functional data object using the coefficients, the basis and a list of names
rMeasfd   = fd(rMeascoef, rMeasbasis,
               list("Time", "Individual", "R measured"))
plot(rMeasfd, lty=1, lwd=2, col=1)

#Plot the fit of the FD object to the actual data
plotfit.fd(rMeasData, daytime, rMeasfd,
           lty=1, lwd=2, main='')
##########################################################
#### Apply linear operator for parabolic B with error ####
##########################################################
#Define linear operator that's just the first derivative
velLfd = int2Lfd(1)
#Evaluate the FD object at vector daytime and apply linear operator defined above
rMeasVelfd = eval.fd(daytime,rMeasfd,velLfd)
#Calculate first derivative of FD object
Rmeasfd.D1 = deriv.fd(rMeasfd, 1)
#Compare output of linear operator evaluation versus first derivative function
plot(daytime,rMeasVelfd)
lines(Rmeasfd.D1)

#Create FD object multiplied by constant mu
mufd = 0.025*rMeasfd
#Plot what should be B, i.e., dRdt + mu*R; lay over step function that is the true B
plot(mufd+Rmeasfd.D1,ylab="Estimated B",xlab="Day")
lines(seq(0,10,0.5),sapply(seq(0,10,0.5),FUN=function(x){-(5000)*(x-5)^2+300000}),col="red")

##################################################
#### Simulate data with constant B with error ####
##################################################
#Define function that descibes ODE of changes in R as a function of time
dRdt = function(t, x, parms) {
  with(as.list(c(parms, x)), {
    dR = B - R*mu
    res = c(dR)
    list(res,B)
  })
}
# Specify initial conditions and parameter values
y0 = c(R = 5000000)
pars = c(mu=0.025,B=300000)
#Specify timescale of interest
tout=seq(from=0,to=10,by=0.5)
#ODE integration and view head of data
out=ode(y=y0, times=tout, func=dRdt, parms=pars)
head(out)
#Convert output to a dataframe
out = as.data.frame(out)
#Create measured R values from true R values and normal distribution
#out$Rmeas = out$R
out$Rmeas = sapply(out$R,FUN=function(x){rnorm(1,x,100000)})
#Plot measured R values (points) and true R values (line)
plot(out$time,out$Rmeas,)
lines(out$time,out$R)

#####################################################
#### Create FDA object for constant B with error ####
#####################################################
Routput = as.matrix(out)
rMeasData = as.matrix(out[,4])  #transpose here if necessary (should be long)

#create vector of times (with one for each day)
daytime  = seq(from=0,to=10,by=0.5)

#Plot the number of cases over time
plot(daytime, rMeasData, "b", lwd=2,
     xlab="Time", ylab="R measured", cex=1.2)

#Create a B-spline basis, where the time interval is 1 to 15, there are 8 basis functions and the order is 4 (implicit)
rMeasbasis    = create.bspline.basis(c(0,10),4)

#Evaluate the basis at the vector of times (days) specified above
rMeasbasismat = eval.basis(rMeasbasis, daytime)

#Compute coefficients for the functional data object by the usual equations for regression coefficients, b = (X′X)^(−1) X′y.
rMeascoef = solve(crossprod(rMeasbasismat),
                  crossprod(rMeasbasismat,rMeasData))
#Can also fit data by least squares
rMeascoef = lsfit(rMeasbasismat, rMeasData, intercept=FALSE)$coef

#Create functional data object using the coefficients, the basis and a list of names
rMeasfd   = fd(rMeascoef, rMeasbasis,
               list("Time", "Individual", "R measured"))
plot(rMeasfd, lty=1, lwd=2, col=1)

#Plot the fit of the FD object to the actual data
plotfit.fd(rMeasData, daytime, rMeasfd,
           lty=1, lwd=2, main='')

#########################################################
#### Apply linear operator for constant B with error ####
#########################################################
#Define linear operator that's just the first derivative
velLfd = int2Lfd(1)
#Evaluate the FD object at vector daytime and apply linear operator defined above
rMeasVelfd = eval.fd(daytime,rMeasfd,velLfd)
#Calculate first derivative of FD object
D2tempfd = deriv.fd(rMeasfd, 1)
#Compare output of linear operator evaluation versus first derivative function
plot(daytime,rMeasVelfd)
lines(D2tempfd)

#Create FD object multiplied by constant mu
mufd = 0.025*rMeasfd
#Plot what should be B, i.e., dRdt + mu*R; lay over step function that is the true B
plot(mufd+D2tempfd,ylab="Estimated B",xlab="Day",ylim=c(0,400000))
abline(a=300000,b=0,col="red")


