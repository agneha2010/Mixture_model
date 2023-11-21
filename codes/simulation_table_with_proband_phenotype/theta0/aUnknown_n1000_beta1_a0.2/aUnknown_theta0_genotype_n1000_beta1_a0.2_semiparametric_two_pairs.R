library(dplyr)
library(invgamma)

param.val = list(n = 1000,
                 a = 0.2,
                 true.mu0 = 0,
                 true.mu1 = 1,
                 true.sigma = 1
)

mu2beta.f = function(par){
  mu0 = par[1]
  mu1 = par[2]
  sigma = par[3]
  a = par[4]
  beta = (mu1-mu0)/sigma^2
  alpha = -(beta*mu0+beta^2*sigma^2/2)
  return(c(alpha,beta,a))
}
tmp = mu2beta.f(c(param.val$true.mu0,param.val$true.mu1,param.val$true.sigma,param.val$a))

param.val$true.alpha = tmp[1]
param.val$true.beta = tmp[2]
param.val$alpha.ini.lb = param.val$true.alpha-2
param.val$alpha.ini.ub = param.val$true.alpha+2
param.val$beta.ini.lb = param.val$true.beta-2
param.val$beta.ini.ub = param.val$true.beta+2
param.val$mu0.ini.lb = param.val$true.mu0-2
param.val$mu0.ini.ub = param.val$true.mu0+2
param.val$mu1.ini.lb = param.val$true.mu1-2
param.val$mu1.ini.ub = param.val$true.mu1+2


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

#### semi-parametric approach
obj.semiparametric.pairs = function(par,t,n){
  #t.list = t.values[[1]]; par = ini.val.ab
  
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
  
  #### mother part
  
  #lambda for Gp=0
  lambda0_1m.par = a.par
  lambda0_0m.par = 1-lambda0_1m.par
  cut0m.par = c(lambda0_0m.par,lambda0_1m.par)
  
  #### sister part
  
  #lambda for Gp=0
  lambda0_1s.par = (-a.par^4+6*a.par^3-9*a.par^2+4*a.par)/(4*(1-a.par)^2)
  lambda0_0s.par = 1-lambda0_1s.par
  cut0s.par = c(lambda0_0s.par,lambda0_1s.par)
  
  x = t[1:n]
  y = t[(n+1):(2*n)]
  z = t[(2*n+1):(3*n)]
  
  #log likelihood part for Gp=0
  t0=(cut0m.par[1]+cut0m.par[2]*exp(alp+x*bet))*(cut0s.par[1]+cut0s.par[2]*exp(alp+y*bet))
  
  value = sum(log(t0))-sum(log(1+tt*te))
  
  return(-value)
}

### wrapper function for optim to ignore runs by producing NA that doesn't converge
f.tryCatch.semiparametric.pairs = function(t.par,init){
  #t.par = t.values[[1]]; init = ini.val.ab
  tryCatch(expr = {
    test = optim(c(init[1],init[2],init[3]),obj.semiparametric.pairs,control=list(maxit=500),t=t.par,n=param.val$n)
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

rep0 = 200 ## number of repetitions
B = 4 ## bound for beta
k = 20 ## Number of random initital values, 20+1

## generate data by calling the function
t.values = list()
for (j in 1:rep0){
  t.values[[j]] = data.generation(j,param.val)
}

start.time <- Sys.time()

est.semiparametric.pairs = matrix(0,rep0,10)
for (j in 1:rep0){
  #j=1
  print(j)
  count = 0
  while(count==0){
    ini.val.ab = c(runif(1,param.val$alpha.ini.lb,param.val$alpha.ini.ub),
                   runif(1,param.val$beta.ini.lb,param.val$beta.ini.ub),rnorm(1))
    temp=f.tryCatch.semiparametric.pairs(t.values[[j]],ini.val.ab)
    val=temp[7]
    count = ifelse(is.na(val),0,1)
  }
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
        temp1 = f.tryCatch.semiparametric.pairs(t.values[[j]],ini.val.ab)
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
  }
  est.semiparametric.pairs[j,] = c(ini.val.ab1,temp)
}

end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

colnames(est.semiparametric.pairs) = c("alpha.ini","beta.ini","a.ini","alpha.est","beta.est","a.est","iteration","nothing","convergence.code","min.value")
#c(mu2beta.f(c(param.val$true.mu0,param.val$true.mu1,param.val$true.sigma,param.val$a)),round(colMeans(est.semiparametric.pairs,na.rm = T),2)[4:6],round(apply(est.semiparametric.pairs,2, function(x) sd(x,na.rm=T)),3)[4:6])
write.csv(est.semiparametric.pairs,"aUnknown_est.semiparametric.pairs_proband_genotype_n1000_beta1_a0.2.csv")
