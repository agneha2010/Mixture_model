library(dplyr)
library(invgamma)

param.val = list(n = 50,
                 a = 0.1, ### allele frequency
                 mus = c(0,1),
                 sds = rep(1,2), ### keep the sds same
                 p = 0.2,
                 true.beta = 1
)

param.val$a1 = exp(param.val$mus[1]*param.val$true.beta+
                     (param.val$true.beta^2)*(param.val$sds[1]^2)/2)
param.val$a2 = exp(param.val$mus[2]*param.val$true.beta+
                     (param.val$true.beta^2)*(param.val$sds[2]^2)/2)

param.val$p.star = param.val$p*param.val$a1/(param.val$p*param.val$a1+(1-param.val$p)*param.val$a2)
param.val$true.alpha = -log(param.val$p*param.val$a1+(1-param.val$p)*param.val$a2)
param.val$mus.star = c(param.val$mus[1]+param.val$true.beta*(param.val$sds[1]^2),
                       param.val$mus[2]+param.val$true.beta*(param.val$sds[2]^2))

param.val$alpha.ini.lb = param.val$true.alpha-2
param.val$alpha.ini.ub = param.val$true.alpha+2
param.val$beta.ini.lb = param.val$true.beta-2
param.val$beta.ini.ub = param.val$true.beta+2
param.val$true.mu0 = sum(c(param.val$p,1-param.val$p)*param.val$mus)
param.val$mu0.ini.lb = param.val$true.mu0-2
param.val$mu0.ini.ub = param.val$true.mu0+2
param.val$true.mu1 = sum(c(param.val$p.star,1-param.val$p.star)*param.val$mus.star)
param.val$mu1.ini.lb = param.val$true.mu1-2
param.val$mu1.ini.ub = param.val$true.mu1+2

data.generation = function(seed,param.val){
  #seed=2
  a = param.val$a
  n = param.val$n
  
  ## theta is probability proband genotype is 1
  theta = 2*a-a^2
  # fix lambda for Gp=1
  lambda1_11 = (-a^3-4*a^2+7*a+2)/(4*(2-a))
  lambda1_10 = (2+a)*(1-a)^2/(4*(2-a))
  lambda1_01 = (1+a)*(1-a)^2/(2*(2-a))
  lambda1_00 = (1-a)^3/(2*(2-a))
  cut10 = c(lambda1_00,lambda1_01,lambda1_10,lambda1_11)
  
  # fix lambda for Gp=0
  lambda0_11 = a*(-a^2-a+2)/(4*(1-a))
  lambda0_10 = a*(a^2-3*a+2)/(4*(1-a))
  lambda0_01 = a*(1-a)/2
  lambda0_00 = (1-a)*(2-a)/2
  cut00 = c(lambda0_00,lambda0_01,lambda0_10,lambda0_11)
  
  cut0=cumsum(cut00)
  cut1=cumsum(cut10)
  
  k = 3*(seed-1)+1
  
  #set.seed(k)
  Gp = rbinom(n,1,theta)
  n1 = length(Gp[Gp==1])
  n0 = length(Gp[Gp==0])
  
  ### data for proband 
  components <- sample(1:2,prob=c(param.val$p,1-param.val$p),size=n0,replace=TRUE)
  z0 = rnorm(n=n0,mean=param.val$mus[components],sd=param.val$sds[components])
  components <- sample(1:2,prob=c(param.val$p.star,1-param.val$p.star),size=n1,replace=TRUE)
  z1 = rnorm(n=n1,mean=param.val$mus.star[components],sd=param.val$sds[components])
  
  u1=runif(n1)
  u0=runif(n0)
  
  #simulate t for Gp=1
  x1=rep(0,n1)
  y1=rep(0,n1)
  
  n11=length(u1[u1<=cut1[1]])
  n12=length(u1[u1>cut1[1] & u1<=cut1[2]])
  n13=length(u1[u1>cut1[2] & u1<=cut1[3]])
  n14=length(u1[u1>cut1[3]])
  
  components <- sample(1:2,prob=c(param.val$p,1-param.val$p),size=n11,replace=TRUE)
  x1[u1<cut1[1]]=rnorm(n=n11,mean=param.val$mus[components],sd=param.val$sds[components])
  components <- sample(1:2,prob=c(param.val$p,1-param.val$p),size=n11,replace=TRUE)
  y1[u1<cut1[1]]=rnorm(n=n11,mean=param.val$mus[components],sd=param.val$sds[components])
  
  components <- sample(1:2,prob=c(param.val$p,1-param.val$p),size=n12,replace=TRUE)
  x1[u1>cut1[1] & u1<=cut1[2]]=rnorm(n=n12,mean=param.val$mus[components],sd=param.val$sds[components])
  components <- sample(1:2,prob=c(param.val$p.star,1-param.val$p.star),size=n12,replace=TRUE)
  y1[u1>cut1[1] & u1<=cut1[2]]=rnorm(n=n12,mean=param.val$mus.star[components],sd=param.val$sds[components])
  
  components <- sample(1:2,prob=c(param.val$p.star,1-param.val$p.star),size=n13,replace=TRUE)
  x1[u1>cut1[2] & u1<=cut1[3]]=rnorm(n=n13,mean=param.val$mus.star[components],sd=param.val$sds[components])
  components <- sample(1:2,prob=c(param.val$p,1-param.val$p),size=n13,replace=TRUE)
  y1[u1>cut1[2] & u1<=cut1[3]]=rnorm(n=n13,mean=param.val$mus[components],sd=param.val$sds[components])
  
  components <- sample(1:2,prob=c(param.val$p.star,1-param.val$p.star),size=n14,replace=TRUE)
  x1[u1>cut1[3]]=rnorm(n=n14,mean=param.val$mus.star[components],sd=param.val$sds[components])
  components <- sample(1:2,prob=c(param.val$p.star,1-param.val$p.star),size=n14,replace=TRUE)
  y1[u1>cut1[3]]=rnorm(n=n14,mean=param.val$mus.star[components],sd=param.val$sds[components])
  
  #simulate t for Gp=0
  x0=rep(0,n0)
  y0=rep(0,n0)
  
  n01=length(u0[u0<=cut0[1]])
  n02=length(u0[u0>cut0[1] & u0<=cut0[2]])
  n03=length(u0[u0>cut0[2] & u0<=cut0[3]])
  n04=length(u0[u0>cut0[3]])
  
  components <- sample(1:2,prob=c(param.val$p,1-param.val$p),size=n01,replace=TRUE)
  x0[u0<cut0[1]]=rnorm(n=n01,mean=param.val$mus[components],sd=param.val$sds[components])
  components <- sample(1:2,prob=c(param.val$p,1-param.val$p),size=n01,replace=TRUE)
  y0[u0<cut0[1]]=rnorm(n=n01,mean=param.val$mus[components],sd=param.val$sds[components])
  
  components <- sample(1:2,prob=c(param.val$p,1-param.val$p),size=n02,replace=TRUE)
  x0[u0>cut0[1] & u0<=cut0[2]]=rnorm(n=n02,mean=param.val$mus[components],sd=param.val$sds[components])
  components <- sample(1:2,prob=c(param.val$p.star,1-param.val$p.star),size=n02,replace=TRUE)
  y0[u0>cut0[1] & u0<=cut0[2]]=rnorm(n=n02,mean=param.val$mus.star[components],sd=param.val$sds[components])
  
  components <- sample(1:2,prob=c(param.val$p.star,1-param.val$p.star),size=n03,replace=TRUE)
  x0[u0>cut0[2] & u0<=cut0[3]]=rnorm(n=n03,mean=param.val$mus.star[components],sd=param.val$sds[components])
  components <- sample(1:2,prob=c(param.val$p,1-param.val$p),size=n03,replace=TRUE)
  y0[u0>cut0[2] & u0<=cut0[3]]=rnorm(n=n03,mean=param.val$mus[components],sd=param.val$sds[components])
  
  components <- sample(1:2,prob=c(param.val$p.star,1-param.val$p.star),size=n04,replace=TRUE)
  x0[u0>cut0[3]]=rnorm(n=n04,mean=param.val$mus.star[components],sd=param.val$sds[components])
  components <- sample(1:2,prob=c(param.val$p.star,1-param.val$p.star),size=n04,replace=TRUE)
  y0[u0>cut0[3]]=rnorm(n=n04,mean=param.val$mus.star[components],sd=param.val$sds[components])
  
  t=c(x0,x1,y0,y1,z0,z1)
  t0 = c(x0,y0)
  t1 = c(x1,y1)
  
  return(list(t,x0,y0,x1,y1,c(n0,n1),z0,z1))
}

# library(dplyr)
# library(ggplot2)
# df = data.frame(data=c(x0,x1),cat=as.factor(c(rep(0,length(x0)),rep(1,length(x1)))))
# ggplot(df, aes(data, fill = cat)) + geom_histogram(alpha = 0.2)

#### semi-parametric approach
obj = function(par,t.list,n){
  #par = c(1.99,-1.99,0.3)
  #t.list = t.values[[1]]
  t = t.list[[1]]
  x0 = t.list[[2]]
  y0 = t.list[[3]]
  x1 = t.list[[4]]
  y1 = t.list[[5]]
  n0 = t.list[[6]][1]
  n1 = t.list[[6]][2]
  z0 = t.list[[7]]
  z1 = t.list[[8]]
  
  alp=par[1]
  bet=-B+2*B/(1+par[2]^2)
  a.par = par[3]^2/(1+par[3]^2)
  te=exp(alp+bet*t)-1
  
  tau=function(a)
  {
    val= mean(te/(1+a*te))
    return(val)
  }
  
  lef=-1/max(te)+0.00000001
  rig=-1/min(te)-0.00000001
  tt=uniroot(tau,c(lef,rig))$root
  
  theta.par = 2*a.par-a.par^2
  #lambda for Gp=1
  lambda1_11.par = (-a.par^3-4*a.par^2+7*a.par+2)/(4*(2-a.par))
  lambda1_10.par = (2+a.par)*(1-a.par)^2/(4*(2-a.par))
  lambda1_01.par = (1+a.par)*(1-a.par)^2/(2*(2-a.par))
  lambda1_00.par = (1-a.par)^3/(2*(2-a.par))
  cut1.par = c(lambda1_00.par,lambda1_01.par,lambda1_10.par,lambda1_11.par)
  
  #lambda for Gp=0
  lambda0_11.par = a.par*(-a.par^2-a.par+2)/(4*(1-a.par))
  lambda0_10.par = a.par*(a.par^2-3*a.par+2)/(4*(1-a.par))
  lambda0_01.par = a.par*(1-a.par)/2
  lambda0_00.par = (1-a.par)*(2-a.par)/2
  cut0.par = c(lambda0_00.par,lambda0_01.par,lambda0_10.par,lambda0_11.par)
  
  #log likelihood part for Gp=0
  tm0=cut0.par[1]+cut0.par[2]*exp(alp+y0*bet)+cut0.par[3]*exp(alp+bet*x0)+cut0.par[4]*exp(2*alp+bet*(x0+y0))
  #log likelihood part for Gp=1
  tm1=cut1.par[1]+cut1.par[2]*exp(alp+y1*bet)+cut1.par[3]*exp(alp+bet*x1)+cut1.par[4]*exp(2*alp+bet*(x1+y1))
  tp=alp+z1*bet
  
  value= n1*log(theta.par)+n0*log(1-theta.par)+sum(log(tm1))+sum(log(tm0))+
    sum(tp)-sum(log(1+tt*te))
  
  return(-value)
}

### wrapper function for optim to ignore runs by producing NA that doesn't converge
f.tryCatch = function(t.par,init){
  tryCatch(expr = {
    test = optim(c(init[1],init[2],init[3]),obj,control=list(maxit=500),t=t.par,n=param.val$n)
    tt=test$par
    te=-B+2*B/(1+tt[2]^2)
    ta = tt[3]^2/(1+tt[3]^2)
    
    c(tt[1],te,ta,test$counts,test$convergence,test$value)
    #c(round(test$estimate,3),test$code,test$iterations)
  },
  error = function(e){
    c(rep(NA,7))
  })
}

################# parametric part approach
obj.parametric.trio = function(par,t.list,n){
  t = t.list[[1]]
  x0 = t.list[[2]]
  y0 = t.list[[3]]
  x1 = t.list[[4]]
  y1 = t.list[[5]]
  n0 = t.list[[6]][1]
  n1 = t.list[[6]][2]
  z0 = t.list[[7]]
  z1 = t.list[[8]]
  
  mu0.par = par[1]
  mu1.par = par[2]
  sigma.par = par[3]
  a.par = par[4]^2/(1+par[4]^2)
  
  theta.par = 2*a.par-a.par^2
  #lambda for Gp=1
  lambda1_11.par = (-a.par^3-4*a.par^2+7*a.par+2)/(4*(2-a.par))
  lambda1_10.par = (2+a.par)*(1-a.par)^2/(4*(2-a.par))
  lambda1_01.par = (1+a.par)*(1-a.par)^2/(2*(2-a.par))
  lambda1_00.par = (1-a.par)^3/(2*(2-a.par))
  cut1.par = c(lambda1_00.par,lambda1_01.par,lambda1_10.par,lambda1_11.par)
  
  #lambda for Gp=0
  lambda0_11.par = a.par*(-a.par^2-a.par+2)/(4*(1-a.par))
  lambda0_10.par = a.par*(a.par^2-3*a.par+2)/(4*(1-a.par))
  lambda0_01.par = a.par*(1-a.par)/2
  lambda0_00.par = (1-a.par)*(2-a.par)/2
  cut0.par = c(lambda0_00.par,lambda0_01.par,lambda0_10.par,lambda0_11.par)
  
  #likelihood part for Gp=0
  tm0=sum(log(cut0.par[4]*dnorm(x0,mu1.par,sigma.par)*dnorm(y0,mu1.par,sigma.par)+
                cut0.par[3]*dnorm(x0,mu1.par,sigma.par)*dnorm(y0,mu0.par,sigma.par)+
                cut0.par[2]*dnorm(x0,mu0.par,sigma.par)*dnorm(y0,mu1.par,sigma.par)+
                cut0.par[1]*dnorm(x0,mu0.par,sigma.par)*dnorm(y0,mu0.par,sigma.par)))
  
  #likelihood part for Gp=1
  tm1=sum(log(cut1.par[4]*dnorm(x1,mu1.par,sigma.par)*dnorm(y1,mu1.par,sigma.par)+
                cut1.par[3]*dnorm(x1,mu1.par,sigma.par)*dnorm(y1,mu0.par,sigma.par)+
                cut1.par[2]*dnorm(x1,mu0.par,sigma.par)*dnorm(y1,mu1.par,sigma.par)+
                cut1.par[1]*dnorm(x1,mu0.par,sigma.par)*dnorm(y1,mu0.par,sigma.par)))
  
  tp1=sum(dnorm(z1,mu1.par,sigma.par,log=T))
  tp0=sum(dnorm(z0,mu0.par,sigma.par,log=T))
  
  value= n1*log(theta.par)+n0*log(1-theta.par)+tm0+tm1+tp0+tp1
  
  return(-value)
}

### wrapper function for optim to ignore runs by producing NA that doesn't converge
f.tryCatch.parametric.trio = function(t.par,init){
  tryCatch(expr = {
    test = optim(c(init[1],init[2],init[3],init[4]),obj.parametric.trio,control=list(maxit=500),t=t.par,n=param.val$n)
    tt=test$par
    ta = tt[4]^2/(1+tt[4]^2)
    c(tt[1:3],ta,test$counts,test$convergence,test$value)
  },
  error = function(e){
    c(rep(NA,8))
  })
}

################# parametric pairs
obj.parametric.pairs = function(par,t.list,n){
  t = t.list[[1]]
  x0 = t.list[[2]]
  y0 = t.list[[3]]
  x1 = t.list[[4]]
  y1 = t.list[[5]]
  n0 = t.list[[6]][1]
  n1 = t.list[[6]][2]
  z0 = t.list[[7]]
  z1 = t.list[[8]]
  
  mu0.par = par[1]
  mu1.par = par[2]
  sigma.par = par[3]
  a.par = par[4]^2/(1+par[4]^2)
  
  theta.par = 2*a.par-a.par^2
  #lambda for Gp=1
  lambda1_1.par = (-a.par^3+a.par^2+a.par)/(a.par*(2-a.par))
  lambda1_0.par = 1-lambda1_1.par
  cut1.par = c(lambda1_0.par,lambda1_1.par)
  
  #lambda for Gp=0
  lambda0_1.par = a.par
  lambda0_0.par = 1-lambda0_1.par
  cut0.par = c(lambda0_0.par,lambda0_1.par)
  
  #likelihood part for Gp=0
  tm0=sum(log(cut0.par[2]*dnorm(x0,mu1.par,sigma.par)+
                cut0.par[1]*dnorm(x0,mu0.par,sigma.par)))
  
  #likelihood part for Gp=1
  tm1=sum(log(cut1.par[2]*dnorm(x1,mu1.par,sigma.par)+
                cut1.par[1]*dnorm(x1,mu0.par,sigma.par)))
  
  lambda1_1s.par = (a.par^4-6*a.par^3+5*a.par^2+4*a.par)/(4*a.par*(2-a.par))
  lambda1_0s.par = 1-lambda1_1s.par
  cut1s.par = c(lambda1_0s.par,lambda1_1s.par)
  
  #lambda for Gp=0
  lambda0_1s.par = (-a.par^4+6*a.par^3-9*a.par^2+4*a.par)/(4*(1-a.par)^2)
  lambda0_0s.par = 1-lambda0_1s.par
  cut0s.par = c(lambda0_0s.par,lambda0_1s.par)
  
  ts0=sum(log(cut0s.par[2]*dnorm(y0,mu1.par,sigma.par)+
                cut0s.par[1]*dnorm(y0,mu0.par,sigma.par)))
  
  #likelihood part for Gp=1
  ts1=sum(log(cut1s.par[2]*dnorm(y1,mu1.par,sigma.par)+
                cut1s.par[1]*dnorm(y1,mu0.par,sigma.par)))
  
  tp1=sum(dnorm(z1,mu1.par,sigma.par,log=T))
  tp0=sum(dnorm(z0,mu0.par,sigma.par,log=T))
  
  value = n1*log(theta.par)+n0*log(1-theta.par)+tm0+tm1+ts0+ts1+tp0+tp1
  
  return(-value)
}

### wrapper function for optim to ignore runs by producing NA that doesn't converge
f.tryCatch.parametric.pairs = function(t.par,init){
  tryCatch(expr = {
    test = optim(c(init[1],init[2],init[3],init[4]),obj.parametric.pairs,control=list(maxit=500),t=t.par,n=param.val$n)
    tt=test$par
    ta = tt[4]^2/(1+tt[4]^2)
    c(tt[1:3],ta,test$counts,test$convergence,test$value)
  },
  error = function(e){
    c(rep(NA,8))
  })
}

################# parametric mother only
obj.parametric.just.moth = function(par,t.list,n){
  t = t.list[[1]]
  x0 = t.list[[2]]
  y0 = t.list[[3]]
  x1 = t.list[[4]]
  y1 = t.list[[5]]
  n0 = t.list[[6]][1]
  n1 = t.list[[6]][2]
  z0 = t.list[[7]]
  z1 = t.list[[8]]
  
  mu0.par = par[1]
  mu1.par = par[2]
  sigma.par = par[3]
  a.par = par[4]^2/(1+par[4]^2)
  
  theta.par = 2*a.par-a.par^2
  # lambda1_1.par = (-a.par^3+a.par^2+a.par)/(a.par*(2-a.par))
  # lambda1_0.par = 1-lambda1_1.par
  # cut1.par = c(lambda1_0.par,lambda1_1.par)
  # 
  # #lambda for Gp=0
  # lambda0_1.par = a.par
  # lambda0_0.par = 1-lambda0_1.par
  # cut0.par = c(lambda0_0.par,lambda0_1.par)
  # 
  # #likelihood part for Gp=0
  # tm0=sum(log(cut0.par[2]*dnorm(x0,mu1.par,sigma.par)+
  #               cut0.par[1]*dnorm(x0,mu0.par,sigma.par)))
  # 
  # #likelihood part for Gp=1
  # tm1=sum(log(cut1.par[2]*dnorm(x1,mu1.par,sigma.par)+
  #               cut1.par[1]*dnorm(x1,mu0.par,sigma.par)))
  
  lambda1_1s.par = (a.par^4-6*a.par^3+5*a.par^2+4*a.par)/(4*a.par*(2-a.par))
  lambda1_0s.par = 1-lambda1_1s.par
  cut1s.par = c(lambda1_0s.par,lambda1_1s.par)
  
  #lambda for Gp=0
  lambda0_1s.par = (-a.par^4+6*a.par^3-9*a.par^2+4*a.par)/(4*(1-a.par)^2)
  lambda0_0s.par = 1-lambda0_1s.par  
  cut0s.par = c(lambda0_0s.par,lambda0_1s.par)
  
  ts0=sum(log(cut0s.par[2]*dnorm(y0,mu1.par,sigma.par)+
                cut0s.par[1]*dnorm(y0,mu0.par,sigma.par)))
  
  #likelihood part for Gp=1
  ts1=sum(log(cut1s.par[2]*dnorm(y1,mu1.par,sigma.par)+
                cut1s.par[1]*dnorm(y1,mu0.par,sigma.par)))
  
  tp1=sum(dnorm(z1,mu1.par,sigma.par,log=T))
  tp0=sum(dnorm(z0,mu0.par,sigma.par,log=T))
  
  value= n1*log(theta.par)+n0*log(1-theta.par)+ts0+ts1+tp0+tp1
  
  return(-value)
}

### wrapper function for optim to ignore runs by producing NA that doesn't converge
f.tryCatch.parametric.just.moth = function(t.par,init){
  tryCatch(expr = {
    test = optim(c(init[1],init[2],init[3],init[4]),obj.parametric.just.moth,control=list(maxit=500),t=t.par,n=param.val$n)
    tt=test$par
    ta = tt[4]^2/(1+tt[4]^2)
    c(tt[1:3],ta,test$counts,test$convergence,test$value)
  },
  error = function(e){
    c(rep(NA,8))
  })
}

rep0 = 200 ## number of repetitions
B = 4 ## bound for beta
k = 20 ## Number of random initital values, 20+1

## generate data by calling the function
t.values = list()
for (j in 1:rep0){
  t.values[[j]] = data.generation(j,param.val)
}

start.time <- Sys.time()

est = matrix(0,rep0,10)
est.parametric.trio = matrix(0,rep0,12)
est.parametric.pairs = matrix(0,rep0,12)
est.parametric.just.moth = matrix(0,rep0,12)
for (j in 1:rep0){
  #j=1
  print(j)
  count = 0
  while(count==0){
    ini.val.ab = c(runif(1,param.val$alpha.ini.lb,param.val$alpha.ini.ub),
                   runif(1,param.val$beta.ini.lb,param.val$beta.ini.ub),rnorm(1))
    temp=f.tryCatch(t.values[[j]],ini.val.ab)
    val=temp[7]
    count = ifelse(is.na(val),0,1)
  }
  ini.val.ms = c(runif(1,param.val$mu0.ini.lb,param.val$mu0.ini.ub),
                 runif(1,param.val$mu1.ini.lb,param.val$mu1.ini.ub),
                 sqrt(rinvgamma(1,shape=1,rate=1)),rnorm(1))
  temp.parametric.trio = f.tryCatch.parametric.trio(t.values[[j]],ini.val.ms)
  temp.parametric.pairs = f.tryCatch.parametric.pairs(t.values[[j]],ini.val.ms)
  temp.parametric.just.moth = f.tryCatch.parametric.just.moth(t.values[[j]],ini.val.ms)
  val.parametric.trio = temp.parametric.trio[8]
  val.parametric.pairs = temp.parametric.pairs[8]
  val.parametric.just.moth = temp.parametric.just.moth[8]
  for(jj in 1:k)
  {
    #jj = 1
    count = 0
    while(count==0){
      ini.val.ab = c(runif(1,param.val$alpha.ini.lb,param.val$alpha.ini.ub),
                     runif(1,param.val$beta.ini.lb,param.val$beta.ini.ub),rnorm(1))
      alp=ini.val.ab[1]
      bet=-B+2*B/(1+ini.val.ab[2]^2)
      te=exp(alp+bet*t.values[[j]][[1]])-1
      if(min(te)*max(te)<0)
      {
        temp1 = f.tryCatch(t.values[[j]],ini.val.ab)
        val1 = temp1[7]
        count = ifelse(is.na(val1),0,1)
        if(count == 1 & val1 < val)
        {
          temp = temp1
          val = val1
          ini.val.ab1 = ini.val.ab
        }
      }
    }
    ini.val.ms = c(runif(1,param.val$mu0.ini.lb,param.val$mu0.ini.ub),
                   runif(1,param.val$mu1.ini.lb,param.val$mu1.ini.ub),
                   sqrt(rinvgamma(1,shape=1,rate=1)),rnorm(1))
    temp1.parametric.trio = f.tryCatch.parametric.trio(t.values[[j]],ini.val.ms)
    temp1.parametric.pairs = f.tryCatch.parametric.pairs(t.values[[j]],ini.val.ms)
    temp1.parametric.just.moth = f.tryCatch.parametric.just.moth(t.values[[j]],ini.val.ms)
    val1.parametric.trio = temp1.parametric.trio[8]
    val1.parametric.pairs = temp1.parametric.pairs[8]
    val1.parametric.just.moth = temp1.parametric.just.moth[8]
    if(val1.parametric.trio < val.parametric.trio)
    {
      temp.parametric.trio = temp1.parametric.trio
      val.parametric.trio = val1.parametric.trio
      ini.val.ms1 = ini.val.ms
    } 
    if(val1.parametric.pairs < val.parametric.pairs)
    {
      temp.parametric.pairs = temp1.parametric.pairs
      val.parametric.pairs = val1.parametric.pairs
      ini.val.ms1 = ini.val.ms
    } 
    if(val1.parametric.just.moth < val.parametric.just.moth)
    {
      temp.parametric.just.moth = temp1.parametric.just.moth
      val.parametric.just.moth = val1.parametric.just.moth
      ini.val.ms1 = ini.val.ms
    } 
  }
  est[j,] = c(ini.val.ab1,temp)
  est.parametric.trio[j,] = c(ini.val.ms1,temp.parametric.trio)
  est.parametric.pairs[j,] = c(ini.val.ms1,temp.parametric.pairs)
  est.parametric.just.moth[j,] = c(ini.val.ms1,temp.parametric.just.moth)
}

end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

colnames(est) = c("alpha.ini","beta.ini","a.ini","alpha.est","beta.est","a.est","iteration","nothing","convergence.code","min.value")
colnames(est.parametric.trio) = c("mu0.ini","mu1.ini","sigma.ini","a.ini","mu0.est","mu1.est","sigma.est","a.est","iteration","nothing","convergence.code","min.value")
colnames(est.parametric.pairs) = c("mu0.ini","mu1.ini","sigma.ini","a.ini","mu0.est","mu1.est","sigma.est","a.est","iteration","nothing","convergence.code","min.value")
colnames(est.parametric.just.moth) = c("mu0.ini","mu1.ini","sigma.ini","a.ini","mu0.est","mu1.est","sigma.est","a.est","iteration","nothing","convergence.code","min.value")

write.csv(est,"aUnknown_est.semiparametric_proband_genotype_n50_beta1_a0.1.csv")
write.csv(est.parametric.trio,"aUnknown_est.parametric.trio_proband_genotype_n50_beta1_a0.1.csv")
write.csv(est.parametric.pairs,"aUnknown_est.parametric.pairs_proband_genotype_n50_beta1_a0.1.csv")
write.csv(est.parametric.just.moth,"aUnknown_est.parametric.just.moth_proband_genotype_n50_beta1_a0.1.csv")
