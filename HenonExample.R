library(CollocInfer)

#x(i+1)=1-ax_i^2+y_i
#y(i+1)=bx_i

#Set parameters for a and b that yield chaotic behaviour
hpars = c(1.4,0.3)

#Generate some data, leaving off the first 20 time points as transients
ntimes = 200
x = c(-1,1) #initial values
X = matrix(0,ntimes+20,2)
X[1,] = x
for (i in 2:(ntimes+20)){
  X[i,] = make.Henon()$ode(i, X[i-1,], hpars, NULL)
}
X = X[20+1:ntimes,]
dim(X)
Y = X + 0.5*matrix(rnorm(ntimes*2),ntimes,2)
t = 1:ntimes

#Start generalized profiling
hpars2 = c(1.3,0.4)
lambda = 1000
coefs = as.matrix(Y)

profile.obj = LS.setup(pars=hpars2,coefs=coefs,fn=make.Henon(),basisvals=NULL,lambda=lambda,times=t,discrete=TRUE)
lik = profile.obj$lik
proc = profile.obj$proc
Ores1 = outeropt(data=Y,times=t,pars=hpars2,coefs=coefs,lik=lik,proc=proc,in.meth="nlminb",out.meth="nlminb")
parest = Ores1$pars
parest
