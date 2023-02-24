make.test = function(){
  fn = function (times, y, p, more) {
    r = y
    r[, 1] = 1 - y[, 1]^2 + p[1] * y[, 2]
    r[, 2] = p[2] * y[,1] + p[3]
    return(r)
  }
  ode = function (times, y, p, more) {
    r = y
    r[1] = 1 - y[1]^2 + p[1] * y[2]
    r[2] = p[2] * y[1] + p[3]
    return(r)
  }
  dfdx = function (times, y, p, more) {
    r = array(0, c(length(times), ncol(y), ncol(y)))
    r[, 1, 1] = - 2 * y[, 1]
    r[, 1, 2] = p[1]
    r[, 2, 1] = p[2]
    return(r)
  }
  dfdp = function (times, y, p, more) {
    r = array(0, c(length(times), ncol(y), length(p)))
    r[, 1, 1] = y[, 2]
    r[, 2, 2] = y[, 1]
    r[, 2, 3] = 1
    return(r)
  }
  d2fdx2 = function (times, y, p, more) {
    r = array(0, c(length(times), rep(ncol(y), 3)))
    r[, 1, 1, 1] = -2
    return(r)
  }
  d2fdxdp = function (times, y, p, more) {
    r = array(0, c(length(times), ncol(y), ncol(y), length(p)))
    r[, 1, 2, 1] = 1
    r[, 2, 1, 2] = 1
    
    return(r)
  }
  list.obj = list(fn,ode,dfdx,dfdp,d2fdx2,d2fdxdp)
  names(list.obj) = c("fn","ode","dfdx","dfdp","d2fdx2","d2fdxdp")
  return(list.obj)
}
make.test()

#Set parameters for a, b and c
testpars = c(1.4,0.3,0.01)

#Generate some data, leaving off the first 20 time points as transients
ntimes = 200
x = c(-1,1) #initial values
X = matrix(0,ntimes+20,2)
X[1,] = x
for (i in 2:(ntimes+20)){
  X[i,] = make.test()$ode(i, X[i-1,], testpars, NULL)
}
head(X)
tail(X)
plot(X[,1],X[,2])
X = X[20+1:ntimes,]
dim(X)
Y = X + 0.025*matrix(rnorm(ntimes*2),ntimes,2)
t = 1:ntimes
plot(Y)

#Start generalized profiling
testpars2 = c(1.5,0.2,0.02)
lambda = 1000
coefs = as.matrix(Y)

profile.obj = LS.setup(pars=testpars2,coefs=coefs,fn=make.test(),basisvals=NULL,lambda=lambda,times=t,discrete=TRUE)
lik = profile.obj$lik
proc = profile.obj$proc

Ores1 = outeropt(data=Y,times=t,pars=testpars2,coefs=coefs,lik=lik,proc=proc,in.meth="nlminb",out.meth="nlminb")
parest = Ores1$pars
parest
