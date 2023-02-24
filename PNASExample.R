library(tidyverse)
library(CollocInfer)

data.full = read_csv("data.csv",
         col_types="iiinnnn"
) %>%
  mutate(
    mouseid=sprintf("%02d-%02d",box,mouse),
    box=sprintf("%02d",box),
    mouse=sprintf("%02d",mouse),
    paba=as.factor(box),
    rbc_density=rbc_density/1000
  ) %>%
  mutate(
    paba=recode(
      paba,
      "01"="0.05","02"="0.005","03"="0.0005","04"="0","05"="control"
    )
  ) %>%
  select(
    day,
    Pd=ama_density,
    RBC=rbc_density,
    Ter119=ter119_density,
    CD71=cd71_density,
    mouseid,
    paba,box,mouse
  ) %>%
  arrange(mouseid,day) %>%
  mutate(
    paba=as.character(paba),
    Ter119=ifelse(Ter119==0,NA,Ter119),
    CD71=ifelse(CD71==0,NA,CD71)
  ) %>%
  mutate(
    Eryth=(1-CD71/Ter119)*RBC,
    Retic=CD71/Ter119*RBC
  )
#Select mouse and box
data = data.full %>% filter(.,mouse=="01",box=="01") %>% select(.,day,Pd,Eryth,Retic)

#Get matrix of data and time vector
data.mat = data %>% select(.,-day) %>% as.matrix
data.mat[is.na(data.mat)] = 0
times = data$day

#Create bases for R(t), W(t) and N(t)
N = length(times)
rangeval = c(0,N-1)
breaks = times
Rbasis = create.bspline.basis(rangeval, breaks=breaks)
Wbasis = Rbasis
Nbasis = Rbasis

#Store the bases in a list "more"
more = vector("list",0)
more$Rbasis = Rbasis
more$Wbasis = Wbasis
more$Nbasis = Nbasis

#Set up functional parameter object that defines the roughness penalty
lambdaDE = 1
penorder = 2
GfdPar = fdPar(basisobj,penorder,lambdaDE)

#Smooth the data with this roughness penalty to get initial values for coefficients defining Eryths(t), Retics(t) and Pd(t)
DEfd0 = smooth.basis(times,data.mat,GfdPar)$fd

#Plot the data and the fit
plotfit.fd(data.mat,times,DEfd0,xlab="Time (hours)",ylab="Eryths,Retics,Parasites",title="")

#Extract the coefficients
coefs0 = DEfd0$coefs

#Set up initial parameter values for  bases
pars0 = matrix(0, Rbasis$nbasis+Nbasis$nbasis+Wbasis$nbasis, 1)

#variables:
#parameters: 
make.ME = function(){
  fn = function (times, y, p, more) {
    
    m2 = 0
    Rbasis = more$Rbasis
    Rnbasis = Rbasis$nbasis
    m1 = m2 + 1
    m2 = m2 + Rnbasis
    Rcoef = p[m1:m2]
    Rmat = eval.basis(times, Rbasis)
    Rvec = Rmat %*% Rcoef
    
    Nbasis = more$Nbasis
    Nnbasis = Nbasis$nbasis
    m1 = m2 + 1
    m2 = m2 + Nnbasis
    Ncoef = p[m1:m2]
    Nmat = eval.basis(times, Nbasis)
    Nvec = Nmat %*% Ncoef
    
    Wbasis = more$Wbasis
    Wnbasis = Wbasis$nbasis
    m1 = m2 + 1
    m2 = m2 + Wnbasis
    Wcoef = p[m1:m2]
    Wmat = eval.basis(times, Wbasis)
    Wvec = Wmat %*% Wcoef
    
    r = y
    Ktemp = 
    r[, "M"] = 
    r[, "E"] = 
    
    return(r)
  }
  list.obj = list(fn)
  names(list.obj) = c("fn")
  return(list.obj)
}

obs.fun = function(times,x,p,more){
  y = cbind(x[,"S"]+x[,"R"], x[,"I"])
  return(y)
}


