#Read in full data and assign to MontrealTemp
cat(MontrealTemp, file='MtlDaily.txt')

#Read in MtlDaily.txt, covert to vector and put into matrix of 34 rows (years) and 365 columns (days)
MtlDaily = matrix(scan("MtlDaily.txt",0),34,365)
#For thaw data, get days 16 to 47 and transpose, so it's days 16 to 47 in rows and 34 columns (years 1961 to 1994)
thawdata = t(MtlDaily[,16:47])

#Create vector of days
daytime  = ((16:47)+0.5)
#Plot the average temperature as a function of days
plot(daytime, apply(thawdata,1,mean), "b", lwd=2,
     xlab="Day", ylab="Temperature (deg C)", cex=1.2)

#Fit these data by regression analysis by using a matrix of values of a basis system taken at the times in vector daytime. Here we construct a basis system over the interval [16,48] using seven cubic B-splines, and evaluate this basis at these points to produce a 32 by 7 matrix. By default the knots are equally spaced over this interval.

thawbasis    = create.bspline.basis(c(16,48),7)
thawbasismat = eval.basis(thawbasis, daytime)

#Compute coefficients for our functional data object by the usual equa-tions for regression coefficients, b = (X′X)^(−1) X′y. X is the thatbasismat (matrix of our basis at vector daytime) and y is the data thawdata (32 x 34 matrix)
thawcoef = solve(crossprod(thawbasismat),
                 crossprod(thawbasismat,thawdata))

#Construct a functional data object by combining coefficients with our basis object
thawfd   = fd(thawcoef, thawbasis,
              list("Day", "Year", "Temperature (deg C)"))
plot(thawfd, lty=1, lwd=2, col=1)

# Figure 4.5

plotfit.fd(thawdata[,1], daytime, thawfd[1],
           lty=1, lwd=2, main='')
plotfit.fd(thawdata[,2], daytime, thawfd[2],
           lty=1, lwd=2, main='')
