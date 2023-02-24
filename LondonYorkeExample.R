library(tidyverse)
library(pomp)
data(package='pomp')

#Prepare data to include only measles cases in Baltimore and recode month to integer; then plot cases as a function of month
LYdata = filter(LondonYorke,disease=="measles",town=="Baltimore",year!="1972") #remove 1972 because only 6 months of data
LYdata$month.int = as.integer(LYdata$month)
head(LYdata)

ggplot()+
  geom_point(data=LYdata,aes(x=month.int,y=cases))+
  geom_line(data=LYdata,aes(x=month.int,y=cases,group=year))+
  xlab("Month")+ylab("Cases")+theme_classic()+
  scale_x_continuous(breaks=c(1:12))

#Convert data into matrix of 33 rows (years) and 12 columns (months)
BltMonthly = LYdata %>% select(.,year,month.int,cases) %>%  spread(.,month.int,cases) %>% select(.,-year) %>% as.matrix
dim(BltMonthly)

#Keep all days and transpose so it's months in rows (1 to 12) and 33 columns (years)
casedata = t(BltMonthly)

#Create vector of months
monthtime  = ((1:12)+0.5)

#Plot the average number of cases as a function of month
plot(monthtime, apply(casedata,1,mean), "b", lwd=2,
     xlab="Month", ylab="Number of cases", cex=1.2)

#Fit these data by regression analysis by using a matrix of values of a basis system taken at the times in vector monthtime. Here we construct a basis system over the interval [1,13] using seven cubic B-splines, and evaluate this basis at these points to produce a 12 by 7 matrix. By default the knots are equally spaced over this interval.

casebasis    = create.bspline.basis(c(1,13),10)
casebasismat = eval.basis(casebasis, monthtime)

#Compute coefficients for our functional data object by the usual equations for regression coefficients, b = (X′X)^(−1) X′y. X is the casebasismat (matrix of our basis at vector monthtime) and y is the data casedata (12 x 33 matrix)
casecoef = solve(crossprod(casebasismat),
                 crossprod(casebasismat,casedata))

#Construct a functional data object by combining coefficients with our basis object
casefd   = fd(casecoef, casebasis,
              list("Month", "Year", "Number of cases"))
plot(casefd, lty=1, lwd=2, col=1)

#Plot an individual year's data with its function
plotfit.fd(casedata[,1], monthtime, casefd[1],
           lty=1, lwd=2, main='')
plotfit.fd(casedata[,3], monthtime, casefd[3],
           lty=1, lwd=2, main='')

meancasefd = mean(casefd)
sumcasefd = sum(casefd)
plot(casefd)
lines(sumcasefd*(1/33),lwd=3)
