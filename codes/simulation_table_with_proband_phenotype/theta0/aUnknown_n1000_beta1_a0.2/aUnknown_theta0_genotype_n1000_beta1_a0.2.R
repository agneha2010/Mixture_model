library(dplyr)
library(invgamma)

param.val = list(n = 1000,
                 a = 0.2,
                 true.mu0 = 0,
                 true.mu1 = 1,
                 true.sigma = 1,
                 alpha.ini.lb = -2.5,
                 alpha.ini.ub = 1.5,
                 beta.ini.lb = -1,
                 beta.ini.ub = 3,
                 mu0.ini.lb = -2,
                 mu0.ini.ub = 2,
                 mu1.ini.lb = -1,
                 mu1.ini.ub = 3
)


data.generation = function(seed,param.val){
  #seed=2
  a = param.val$a
  n = param.val$n
  u0 = param.val$true.mu0
  u1 = param.val$true.mu1
  s2 = param.val$true.sigma
  lambda_11 = a*(-a^2-a+2)/(4*(1-a))
  lambda_10 = a*(a^2-3*a+2)/(4*(1-a))
  lambda_01 = a*(1-a)/2
  lambda_00 = (1-a)*(2-a)/2
  cut0 = c(lambda_00,lambda_01,lambda_10,lambda_11)
  cut=cumsum(cut0)
  
  ### data for proband 
  z = rnorm(n,u0,s2)
  
  u=runif(n)
  
  x=rep(0,n)
  y=rep(0,n)
  
  n1=length(u[u<=cut[1]])
  n2=length(u[u>cut[1] & u<=cut[2]])
  n3=length(u[u>cut[2] & u<=cut[3]])
  n4=length(u[u>cut[3]])
  
  #set.seed(k+1)
  x[u<cut[1]]=rnorm(n1,u0,s2)
  #set.seed(k+2)
  y[u<cut[1]]=rnorm(n1,u0,s2)
  
  #set.seed(k+1)
  x[u>cut[1] & u<=cut[2]]=rnorm(n2,u0,s2)
  #set.seed(k+2)
  y[u>cut[1] & u<=cut[2]]=rnorm(n2,u1,s2)
  
  #set.seed(k+1)
  x[u>cut[2] & u<=cut[3]]=rnorm(n3,u1,s2)
  #set.seed(k+2)
  y[u>cut[2] & u<=cut[3]]=rnorm(n3,u0,s2)
  
  #set.seed(k+1)
  x[u>cut[3]]=rnorm(n4,u1,s2)
  #set.seed(k+2)
  y[u>cut[3]]=rnorm(n4,u1,s2)
  t=c(x,y,z)
  
  return(t)
}

#### semi-parametric
obj = function(par,t,n){
  alp = par[1]
  bet = -B+2*B/(1+par[2]^2)
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
  
  lambda_11.par = a.par*(-a.par^2-a.par+2)/(4*(1-a.par))
  lambda_10.par = a.par*(a.par^2-3*a.par+2)/(4*(1-a.par))
  lambda_01.par = a.par*(1-a.par)/2
  lambda_00.par = (1-a.par)*(2-a.par)/2
  cut0.par = c(lambda_00.par,lambda_01.par,lambda_10.par,lambda_11.par)
  x = t[1:n]
  y = t[(n+1):(2*n)]
  z = t[(2*n+1):(3*n)]
  
  tm=cut0.par[1]+cut0.par[2]*exp(alp+y*bet)+cut0.par[3]*exp(alp+bet*x)
  tm=tm+cut0.par[4]*exp(2*alp+bet*(x+y))
  
  value=sum(log(tm))-sum(log(1+tt*te))
  return(-value)
}


f.tryCatch = function(t.par,init){
  tryCatch(expr = {
    test = optim(c(init[1],init[2],init[3]),obj,control=list(maxit=500),t=t.par,n=param.val$n)
    tt = test$par
    te = -B+2*B/(1+tt[2]^2)
    ta = tt[3]^2/(1+tt[3]^2)
    
    c(tt[1],te,ta,test$counts,test$convergence,test$value)
    #c(round(test$estimate,3),test$code,test$iterations)
  },
  error = function(e){
    c(rep(NA,7))
  })
}

################# parametric part
obj.parametric.trio = function(par,t,n){
  mu0.par = par[1]
  mu1.par = par[2]
  sigma.par = par[3]
  a.par = par[4]^2/(1+par[4]^2)
  
  lambda_11.par = a.par*(-a.par^2-a.par+2)/(4*(1-a.par))
  lambda_10.par = a.par*(a.par^2-3*a.par+2)/(4*(1-a.par))
  lambda_01.par = a.par*(1-a.par)/2
  lambda_00.par = (1-a.par)*(2-a.par)/2
  cut0.par = c(lambda_00.par,lambda_01.par,lambda_10.par,lambda_11.par)
  x = t[1:n]
  y = t[(n+1):(2*n)]
  z = t[(2*n+1):(3*n)]
  
  value = sum(log(cut0.par[4]*dnorm(x,mu1.par,sigma.par)*dnorm(y,mu1.par,sigma.par)+
                    cut0.par[3]*dnorm(x,mu1.par,sigma.par)*dnorm(y,mu0.par,sigma.par)+
                    cut0.par[2]*dnorm(x,mu0.par,sigma.par)*dnorm(y,mu1.par,sigma.par)+
                    cut0.par[1]*dnorm(x,mu0.par,sigma.par)*dnorm(y,mu0.par,sigma.par)))+
    sum(dnorm(z,mu0.par,sigma.par,log=T))
  return(-value)
}

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
obj.parametric.pairs = function(par,t,n){
  #t = t.values[[j]]
  #par = ini.val.ms
  #n = param.val$n
  mu0.par = par[1]
  mu1.par = par[2]
  sigma.par = par[3]
  a.par = par[4]^2/(1+par[4]^2)
  
  x = t[1:n]
  y = t[(n+1):(2*n)]
  z = t[(2*n+1):(3*n)]
  
  theta.par = 2*a.par-a.par^2
  #lambda for Gp=0, mother
  lambda0_1.par = a.par
  lambda0_0.par = 1-lambda0_1.par
  cut1.par = c(lambda0_0.par,lambda0_1.par)
  
  #likelihood part for Gp=0
  tm1=sum(log(cut1.par[2]*dnorm(x,mu1.par,sigma.par)+
                cut1.par[1]*dnorm(x,mu0.par,sigma.par)))
  
  #lambda for Gp=0, sister
  lambda0_1s.par = (-a.par^4+6*a.par^3-9*a.par^2+4*a.par)/(4*(1-a.par)^2)
  lambda0_0s.par = 1-lambda0_1s.par
  cut1s.par = c(lambda0_0s.par,lambda0_1s.par)
  
  #likelihood part for Gp=0
  ts1=sum(log(cut1s.par[2]*dnorm(y,mu1.par,sigma.par)+
                cut1s.par[1]*dnorm(y,mu0.par,sigma.par)))
  
  value = tm1+ts1+sum(dnorm(z,mu0.par,sigma.par,log=T))
  
  return(-value)
}

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
obj.parametric.just.moth = function(par,t,n){
  mu0.par = par[1]
  mu1.par = par[2]
  sigma.par = par[3]
  a.par = par[4]^2/(1+par[4]^2)
  
  x = t[1:n]
  y = t[(n+1):(2*n)]
  z = t[(2*n+1):(3*n)]
  
  # lambda for Gp=0, mother
  # lambda0_1.par = a.par
  # lambda0_0.par = 1-lambda0_1.par
  # cut1.par = c(lambda0_0.par,lambda0_1.par)
  # # 
  # likelihood part for Gp=0
  # tm1=sum(log(cut1.par[2]*dnorm(x,mu1.par,sigma.par)+
  #               cut1.par[1]*dnorm(x,mu0.par,sigma.par)))
  
  #lambda for Gp=0, sister
  lambda0_1s.par = (-a.par^4+6*a.par^3-9*a.par^2+4*a.par)/(4*(1-a.par)^2)
  lambda0_0s.par = 1-lambda0_1s.par
  cut1s.par = c(lambda0_0s.par,lambda0_1s.par)
  
  #likelihood part for Gp=1
  ts1=sum(log(cut1s.par[2]*dnorm(y,mu1.par,sigma.par)+
                cut1s.par[1]*dnorm(y,mu0.par,sigma.par)))
  
  value = ts1+sum(dnorm(z,mu0.par,sigma.par,log=T))
  return(-ts1)
}

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

rep0 = 200
B = 4
k = 20
t.values = list()
for (j in 1:rep0){
  t.values[[j]] = data.generation(j,param.val)
}

mu2beta.f = function(par){
  mu0 = par[1]
  mu1 = par[2]
  sigma = par[3]
  a = par[4]
  beta = (mu1-mu0)/sigma^2
  alpha = -(beta*mu0+(beta^2)*(sigma^2)/2)
  return(c(alpha,beta,a))
}
#mu2beta.f(c(rnorm(2),1,runif(1)))

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
      te=exp(alp+bet*t.values[[j]])-1
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
    ini.val.ms=c(runif(1,param.val$mu0.ini.lb,param.val$mu0.ini.ub),
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

write.csv(est,"aUnknown_est.semiparametric_theta0_genotype_n1000_beta1_a0.2.csv")
write.csv(est.parametric.trio,"aUnknown_est.parametric.trio_theta0_genotype_n1000_beta1_a0.2.csv")
write.csv(est.parametric.pairs,"aUnknown_est.parametric.pairs_theta0_genotype_n1000_beta1_a0.2.csv")
write.csv(est.parametric.just.moth,"aUnknown_est.parametric.just.moth_theta0_genotype_n1000_beta1_a0.2.csv")
