library(tidyverse)
library(fda)
library(pomp)

head(blowflies)
ggplot(data=blowflies %>% filter(.,day%in%250:300))+geom_point(aes(x=day,y=count,col=factor(set)))+geom_line(aes(x=day,y=count,group=factor(set)))+xlab("Day")+ylab("Count")+theme_classic()


bfDaily = blowflies %>% filter(.,day%in%250:300) %>% spread(.,set,count) %>% na.omit()
bfdata = bfDaily %>% select(.,-day) %>% as.matrix #transpose here if necessary

#create vector of times (with one for each day)
daytime  = bfDaily$day  

#Plot the average number of blowflies as a function of days
plot(daytime, apply(bfdata,1,mean), "b", lwd=2,
     xlab="Day", ylab="Count", cex=1.2)

#Create a B-spline basis, where the time interval is 1 to 15, there are 8 basis functions and the order is 4 (implicit)
bfbasis    = create.bspline.basis(c(min(daytime),max(daytime)),8)

#Evaluate the basis at the vector of times (days) specified above
bfbasismat = eval.basis(bfbasis, daytime)

#Compute coefficients for the functional data object by the usual equations for regression coefficients, b = (X′X)^(−1) X′y.
bfcoef = solve(crossprod(bfbasismat),
                 crossprod(bfbasismat,bfdata))
#Can also fit data by least squares (should be identical)
bfcoef = lsfit(bfbasismat, bfdata, intercept=FALSE)$coef

#Create functional data object using the coefficients, the basis and a list of names
bffd   = fd(bfcoef, bfbasis,
              list("Day", "Set", "Count"))
plot(bffd, lty=1, lwd=2, col=1)

#Plot the fit of the FD object to the actual data (choose set)
plotfit.fd(bfdata[,1], daytime, bffd[1],
          lty=1, lwd=2, main='')

#Create vector of times for which to evaluate FD object or derivative of FD object
dayvec = seq(250,300,1)
#Evaluate FD object and 1st and 2nd derivatives of FD object at vector of times specified above
bfVec = eval.fd(dayvec,bffd)
DbfVec = eval.fd(dayvec, bffd, 1)
D2bfVec = eval.fd(dayvec, bffd, 2)
#Plot evaluated FD
plot(dayvec,bfVec[,1],type="l")
lines(dayvec,bfVec[,2],lty=2)
lines(dayvec,bfVec[,3],lty=3)
lines(dayvec,bfVec[,4],lty=4)
#Plot first derivative versus second derivative, label points with time
plot(DbfVec[,1],D2bfVec[,1],type="l")
text(DbfVec[,1],D2bfVec[,1],dayvec)

#Define linear differential operator, in this case just the second derivative,  and evaluate FD object with LD operator applied at times in dayvec. Plot the evaluation and compare to second derivative of FD object evaluated above.
accelLfd = int2Lfd(2)
Lbfmat = eval.fd(dayvec, bffd, accelLfd)
plot(Lbfmat[,1])
lines(D2bfVec[,1])

#The function smooth.basis is provided to produce the same results as the fd function above (as well as much more) without the need to explicitly evaluate the basis functions. Need to provide the vector of times, the data and the basis.
bffd.smooth = smooth.basis(daytime,bfdata,bfbasis)
plot(bffd.smooth$fd)
plotfit.fd(bfdata, daytime, bffd.smooth$fd,
           lty=1, lwd=2, main='')
#To obtain the degrees of freedom used to obtain the fitted curve
bffd.smooth$df
#To obtain the the value of the generalized cross-validation criterion: a measure of lack of fit discounted for degrees of freedom. If there are multiple curves, a vector is returned containing gcv values for each curve.
bffd.smooth$gcv

#The coefficient estimates are obtained from the data in the vector fluCdata (y) by multiplying this vector by a matrix (from evaluating flucCbasis at specified time points). That is, the coefficients map to y via y2cmap
y2cMap = solve(crossprod(bfbasismat),
               t(bfbasismat))

#Compute the penalty matrix (should be 8x8)
Rmat = eval.penalty(bfbasis, 2)
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
  bffdPari  = fdPar(bfbasis, Lfdobj=2, 10^loglam[i])
  bfSm.i    = smooth.basis(daytime, bfdata, bffdPari)
  Gcvsave[i] = sum(bfSm.i$gcv)
  Dfsave[i]  = bfSm.i$df
}
Gcvsave
plot(Gcvsave)
bffdPar = fdPar(bfbasis, 2, 10^(-0.75))

bffd.smoothPar = smooth.basis(daytime,bfdata,bffdPar)
plotfit.fd(bfdata, daytime, bffd.smoothPar$fd,
           lty=1, lwd=2, main='')

