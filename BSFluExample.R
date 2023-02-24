library(fda)
library(pomp)
data(package='pomp')
bsflu

###################
#### Bedridden ####
###################
#Read in data as matrix
fluBDaily = as.matrix(bsflu[,2])
fluBdata = fluBDaily #transpose here if necessary - should be 1 column with 14 rows, one row for each day

#create vector of times (with one for each day)
daytime  = ((bsflu$day)+0.5)  

#Plot the number of cases over time
plot(daytime, fluBdata, "b", lwd=2,
     xlab="Day", ylab="Number bedridden", cex=1.2)

#Create a B-spline basis, where the time interval is 1 to 15, there are 8 basis functions and the order is 4 (implicit)
fluBbasis    = create.bspline.basis(c(1,15),8)

#Evaluate the basis at the vector of times (days) specified above
fluBbasismat = eval.basis(fluBbasis, daytime)

#Compute coefficients for the functional data object by the usual equations for regression coefficients, b = (X′X)^(−1) X′y. X is the fluBbasismat (14 x 8 matrix of our basis evaluated at vector daytime) and y is the data fluBdata (14 x 1 matrix). crossprod(x,y) is equivalent to t(x) %*% y. The first crossprod term is thus t(fluBbasismat) %*% fluBbasismat, so 8x14 %*% 14x8 = 8x8. The second crossprod term is t(fluBbasismat) %*% fluBdata, so 8x14 %*% 14x1 = 8x1. Output should then be 8x8 * 8x1 = 8x1.
fluBcoef = solve(crossprod(fluBbasismat),
                 crossprod(fluBbasismat,fluBdata))
dim(fluBcoef)

#Create functional data object using the coefficients, the basis and a list of names
fluBfd   = fd(fluBcoef, fluBbasis,
              list("Day","Year", "Number bedridden"))
plot(fluBfd, lty=1, lwd=2, col=1)

#Plot the fit of the FD object to the actual data
plotfit.fd(fluBdata, daytime, fluBfd,
           lty=1, lwd=2, main='')

######################
#### Convalescent ####
######################
fluCDaily = as.matrix(bsflu[,3])
fluCdata = fluCDaily #transpose here if necessary

#create vector of times (with one for each day)
daytime  = ((1:14)+0.5)  

#Plot the number of cases over time
plot(daytime, fluCdata, "b", lwd=2,
     xlab="Day", ylab="Number convalescent", cex=1.2)

#Create a B-spline basis, where the time interval is 1 to 15, there are 8 basis functions and the order is 4 (implicit)
fluCbasis    = create.bspline.basis(c(1,15),8)

#Evaluate the basis at the vector of times (days) specified above
fluCbasismat = eval.basis(fluCbasis, daytime)

#Compute coefficients for the functional data object by the usual equations for regression coefficients, b = (X′X)^(−1) X′y. X is the fluBbasismat (14 x 8 matrix of our basis evaluated at vector daytime) and y is the data fluBdata (14 x 1 matrix). crossprod(x,y) is equivalent to t(x) %*% y. The first crossprod term is thus t(fluBbasismat) %*% fluBbasismat, so 8x14 %*% 14x8 = 8x8. The second crossprod term is t(fluBbasismat) %*% fluBdata, so 8x14 %*% 14x1 = 8x1. Output should then be 8x8 * 8x1 = 8x1.
fluCcoef = solve(crossprod(fluCbasismat),
                 crossprod(fluCbasismat,fluCdata))
#Can also fit data by least squares
fluCcoef = lsfit(fluCbasismat, fluCdata, intercept=FALSE)$coef

#Create functional data object using the coefficients, the basis and a list of names
fluCfd   = fd(fluCcoef, fluCbasis,
              list("Day", "Year", "Number convalescent"))
plot(fluCfd, lty=1, lwd=2, col=1)

#Plot the fit of the FD object to the actual data
plotfit.fd(fluCdata, daytime, fluCfd,
           lty=1, lwd=2, main='')

#Create vector of times for which to evaluate FD object or derivative of FD object
dayvec = seq(1,14,0.25)
#Evaluate FD object and 1st and 2nd derivatives of FD object at vector of times specified above
fluCVec = eval.fd(dayvec,fluCfd)
DfluCVec = eval.fd(dayvec, fluCfd, 1)
D2fluCVec = eval.fd(dayvec, fluCfd, 2)
#Plot evaluated FD, D FD and D^2 FD
plot(dayvec,fluCVec,type="l",ylim=c(-150,175))
lines(dayvec, DfluCVec, lty=2)
lines(dayvec, D2fluCVec, lty=3)
#Plot first derivative versus second derivative, label points with time
plot(DfluCVec,D2fluCVec,type="l")
text(DfluCVec,D2fluCVec,dayvec)

#Define linear differential operator, in this case just the second derivative,  and evaluate FD object with LD operator applied at times in dayvec. Plot the evaluation and compare to second derivative of FD object evaluated above.
accelLfd = int2Lfd(2)
LfluCmat = eval.fd(dayvec, fluCfd, accelLfd)
plot(LfluCmat)
lines(D2fluCVec)

#The function smooth.basis is provided to produce the same results as the fd function above (as well as much more) without the need to explicitly evaluate the basis functions. Need to provide the vector of times, the data and the basis.
fluCfd.smooth = smooth.basis(daytime,fluCdata,fluCbasis)
plot(fluCfd.smooth$fd)
plotfit.fd(fluCdata, daytime, fluCfd.smooth$fd,
           lty=1, lwd=2, main='')
#To obtain the degrees of freedom used to obtain the fitted curve
fluCfd.smooth$df
#To obtain the the value of the generalized cross-validation criterion: a measure of lack of fit discounted for degrees of freedom. If there are multiple curves, a vector is returned containing gcv values for each curve.
fluCfd.smooth$gcv

#The coefficient estimates are obtained from the data in the vector fluCdata (y) by multiplying this vector by a matrix (from evaluating flucCbasis at specified time points). That is, the coefficients map to y via y2cmap
y2cMap = solve(crossprod(fluCbasismat),
               t(fluCbasismat))

#Compute the penalty matrix (should be 8x8)
Rmat = eval.penalty(fluCbasis, 2)
dim(Rmat)

#Define a functional parameter object that penalizes the roughness of acceleration by using the second derivative in the roughness penalty
fluCfdPar = fdPar(fluCbasis, 2, 0.1)

fluCfd.smoothPar = smooth.basis(daytime,fluCdata,fluCfdPar)
plot(fluCfd.smoothPar$fd)
plotfit.fd(fluCdata, daytime, fluCfd.smoothPar$fd,
           lty=1, lwd=2, main='')
lines(fluCfd.smooth$fd,lwd=2,lty=2)

#Log lambda values
loglam         = seq(-6, 0, 0.25)
#Vector to output GCV to
Gcvsave        = rep(NA, length(loglam))
names(Gcvsave) = loglam
Dfsave         = Gcvsave
for(i in 1:length(loglam)){
  fluCfdPari  = fdPar(fluCbasis, Lfdobj=2, 10^loglam[i])
  fluCSm.i    = smooth.basis(daytime, fluCdata, fluCfdPari)
  Gcvsave[i] = sum(fluCSm.i$gcv)
  Dfsave[i]  = fluCSm.i$df
}
Gcvsave

fluCfdPar = fdPar(fluCbasis, 2, 10^(-0.75))

fluCfd.smoothPar = smooth.basis(daytime,fluCdata,fluCfdPar)
plotfit.fd(fluCdata, daytime, fluCfd.smoothPar$fd,
           lty=1, lwd=2, main='')
abline(a=0,b=0)
fluCfd.Pos = smooth.pos(daytime,fluCdata,fluCfdPar)
fluCWfd = fluCfd.Pos$Wfdobj
fluCWfd$fdnames = list("Day",
                   "Set",
                   "Count")
plot(fluCWfd)

precfit = exp(eval.fd(dayvec, fluCWfd))
plot(daytime, fluCdata, type="p", cex=1.2,
     xlab="Day",
     ylab="Count")
lines(dayvec, precfit,lwd=2)
