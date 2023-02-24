library(CollocInfer)
make.SIR = function(){
  fn = function (times, y, p, more) {
    r = y
    r[, 1] = y[,1] - (p[1]*y[,1]*y[,2])
    r[, 2] = y[,2] + (p[1]*y[,1]*y[,2]) - p[2]*y[,2]
    r[, 3] = y[,3] + p[2]*y[,2]
    
    return(r)
  }
  ode = function (times, y, p, more) {
    r = y
    r[1] = y[1] - (p[1]*y[1]*y[2])
    r[2] = y[2] + (p[1]*y[1]*y[2]) - p[2]*y[2]
    r[3] = y[3] + p[2]*y[2]
    return(r)
  }
  dfdx = function (times, y, p, more) {
    r = array(0, c(length(times), ncol(y), ncol(y)))
    r[, 1, 1] = 1 - (p[1]*y[,2])
    r[, 1, 2] = -(p[1]*y[,1])
    r[, 2, 1] = (p[1]*y[,2])
    r[, 2, 2] = 1 + p[1]*y[,1] - p[2]
    r[, 3, 2] = p[2]
    r[, 3, 3] = 1
    return(r)
  }
  dfdp = function (times, y, p, more) {
    r = array(0, c(length(times), ncol(y), length(p)))
    r[, 1, 1] = -y[, 1]*y[, 2]
    r[, 2, 1] = y[, 1]*y[, 2]
    r[, 2, 2] = -y[, 2]
    r[, 3, 2] = y[, 2]
    return(r)
  }
  d2fdx2 = function (times, y, p, more) {
    r = array(0, c(length(times), rep(ncol(y), 3)))
    r[, 1, 1, 2] = -p[1]
    r[, 1, 2, 1] = -p[1]
    r[, 2, 1, 2] = p[1]
    r[, 2, 2, 1] = p[1]
    return(r)
  }
  d2fdxdp = function (times, y, p, more) {
    r = array(0, c(length(times), ncol(y), ncol(y), length(p)))
    r[, 1, 1, 1] = -y[, 2]
    r[, 1, 2, 1] = -y[, 1]
    r[, 2, 1, 1] = y[, 2]
    r[, 2, 2, 1] = y[, 1]
    r[, 2, 2, 2] = -1
    r[, 3, 2, 2] = 1
    return(r)
  }
  list.obj = list(fn,ode,dfdx,dfdp,d2fdx2,d2fdxdp)
  names(list.obj) = c("fn","ode","dfdx","dfdp","d2fdx2","d2fdxdp")
  return(list.obj)
}
make.SIR()

#Set parameters for a, b and c
testpars = c(4/160000,3.5)

#Generate some data
ntimes = 48
x = c(160000-3,3,0) #initial values
X = matrix(0,ntimes,3)
X[1,] = x
for (i in 2:ntimes){
  X[i,] = make.SIR()$ode(i, X[i-1,], testpars, NULL)
}
head(X)
tail(X)
matplot(X)
dim(X)
#Y = X + 100*matrix(rnorm(ntimes*3),ntimes,3)
#Y[Y<0] = 0
#Y = round(Y)
#matplot(Y,type="l")
t = 1:ntimes
#plot(Y)

#Start generalized profiling
testpars2 = c(3/160000,3)
lambda = 1000
coefs = as.matrix(X)

profile.obj = LS.setup(pars=testpars2,coefs=coefs,fn=make.SIR(),basisvals=NULL,lambda=lambda,times=t,discrete=TRUE)
lik = profile.obj$lik
proc = profile.obj$proc

Ores1 = outeropt(data=X,times=t,pars=testpars2,coefs=coefs,lik=lik,proc=proc,in.meth="nlminb",out.meth="nlminb")
parest = Ores1$pars
parest
