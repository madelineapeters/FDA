library(CollocInfer)

function (times, x, p, more) 
{
  m2 = 0
  betabasis = more$betabasis
  betanbasis = betabasis$nbasis
  m1 = m2 + 1
  m2 = m2 + betanbasis
  betacoef = p[m1:m2]
  betamat = eval.basis(times, betabasis)
  betavec = betamat %*% betacoef
  alphabasis = more$alphabasis
  alphanbasis = alphabasis$nbasis
  m1 = m2 + 1
  m2 = m2 + alphanbasis
  alphacoef = p[m1:m2]
  alphamat = eval.basis(times, alphabasis)
  alphavec = alphamat %*% alphacoef
  rain = eval.fd(times, more$rainfd)
  f = -betavec * x + alphavec * rain
  return(f)
}


data(NSdata)

yobs = NSgroundwater
zobs = NSrainfall
tobs = NStimes

N = length(tobs)
rangeval = c(0,N-1)

#Plot groundwater data
par(mfrow = c(2,1))
plot(tobs, yobs, "p",  xlab="Time (hours)", 
     ylab="Groundwater level (metres)")
plot(tobs, zobs, "p",  xlab="Time (hours)", 
     ylab="Hourly rainfall (millimetres)")

#Set up functional data object for rainfall
norder = 1
nbasis = N + norder - 2
rainbasis = create.bspline.basis(rangeval, nbasis, norder)
rainfd = smooth.basis(tobs,zobs,rainbasis)$fd

#Plot fit of functional data object to data
plotfit.fd(zobs, tobs, rainfd, nfine=N)

#Set up the basis for groundwater, with a knot placed at each mid-hour point

knots = c(rangeval[1],seq(rangeval[1]+0.5,rangeval[2]-0.5,len=N-1),rangeval[2])
norder = 3
nbasis = length(knots) + norder - 2
basisobj = create.bspline.basis(rangeval,nbasis,norder,knots)

#Create constant coefficient basis for constant beta and alpha
conbasis = create.constant.basis(rangeval)
nbetabasis = 1
betabasis = conbasis
nalphabsis = 1
alphabasis = conbasis

#Store the bases and rainfall object in a list "more"
more = vector("list",0)
more$betabasis = betabasis
more$alphabasis = alphabasis
more$rainfd = rainfd

#Set up functional parameter object that defines the roughness penalty
lambdaDE = 1
penorder = 1
GfdPar = fdPar(basisobj,penorder,lambdaDE)

#Smooth the data with this roughness penalty to get initial values for coefficients defining G(t)
DEfd0 = smooth.basis(tobs,yobs,GfdPar)$fd

#Plot the data and the fit
plotfit.fd(yobs,tobs,DEfd0,xlab="Time (hours)",ylab="Groundwater level (metres)",title="")

#Extract the coefficients
coefs0 = DEfd0$coefs

#Set up initial parameter values for constant bases
pars0 = matrix(0, nbetabasis+nalphabsis, 1)

#Profile estimation step
control.out = list()
control.out$trace = 6
control.out$tol = 1e-6

lambda = 1e0

#Estimate the parameters using differencing
res0 = Profile.LS(make.NS()$fn, yobs, tobs, pars0, coefs0, basisobj, lambda, more = more, out.meth="ProfileGN",control.out=control.out)

#Estimate the parameters using the derivative evaluation functions
lambda = 1e0
NSfun = make.NS()
res1 = Profile.LS(NSfun, yobs, tobs, pars0, coefs0, basisobj, lambda, more=more, out.meth='ProfileGN',control.out=control.out)
#Set up the functional data object for the fit
DEfd1 = fd(res1$coefs, basisobj)

#Try increasing lambda
lambda = 1e2
res2 = Profile.LS(NSfun, yobs, tobs, res1$pars, res1$coefs, basisobj, lambda,more = more, out.meth='ProfileGN',control.out=control.out)
res2$pars
DEfd2 = fd(res2$coefs, basisobj)

lambda = 1e4
res3 = Profile.LS(NSfun, yobs, tobs, res2$pars, res2$coefs, basisobj, lambda,
                  more = more, out.meth='ProfileGN', 
                  control.out=control.out)
res3$pars  #  display the parameter estimates

DEfd3 = fd(res3$coefs, basisobj)

#Compare fits to the data
par(mfrow=c(1,1))
plotfit.fd(yobs, tobs, DEfd1, 
           xlab="Time (hours)", ylab="Groundwater level (metres)", 
           title="")
lines(DEfd2)
lines(DEfd3)

