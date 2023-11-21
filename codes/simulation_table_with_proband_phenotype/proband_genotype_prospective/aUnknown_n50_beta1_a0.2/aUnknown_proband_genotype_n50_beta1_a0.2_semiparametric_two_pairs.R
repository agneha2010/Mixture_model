library(dplyr)
library(invgamma)

param.val = list(n = 50,
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
mu2beta.f = function(par){
  mu0 = par[1]
  mu1 = par[2]
  sigma = par[3]
  a = par[4]
  beta = (mu1-mu0)/sigma^2
  alpha = -(beta*mu0+beta^2*sigma^2/2)
  return(c(alpha,beta,a))
}
#mu2beta.f(c(2,1.5,1,0.3))


data.generation = function(seed,param.val){
  #seed=2
  a = param.val$a
  n = param.val$n
  mu0 = param.val$true.mu0
  mu1 = param.val$true.mu1
  s2 = param.val$true.sigma
  
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
  
  # n1 = n/2
  # n0 = n/2
  
  ### data for proband 
  z0 = rnorm(n0,mu0,s2)
  z1 = rnorm(n1,mu1,s2)
  
  u1=runif(n1)
  u0=runif(n0)
  
  #simulate t for Gp=1
  x1=rep(0,n1)
  y1=rep(0,n1)
  
  n11=length(u1[u1<=cut1[1]])
  n12=length(u1[u1>cut1[1] & u1<=cut1[2]])
  n13=length(u1[u1>cut1[2] & u1<=cut1[3]])
  n14=length(u1[u1>cut1[3]])
  
  #set.seed(k+1)
  x1[u1<cut1[1]]=rnorm(n11,mu0,s2)
  #set.seed(k+2)
  y1[u1<cut1[1]]=rnorm(n11,mu0,s2)
  
  #set.seed(k+1)
  x1[u1>cut1[1] & u1<=cut1[2]]=rnorm(n12,mu0,s2)
  #set.seed(k+2)
  y1[u1>cut1[1] & u1<=cut1[2]]=rnorm(n12,mu1,s2)
  
  #set.seed(k+1)
  x1[u1>cut1[2] & u1<=cut1[3]]=rnorm(n13,mu1,s2)
  #set.seed(k+2)
  y1[u1>cut1[2] & u1<=cut1[3]]=rnorm(n13,mu0,s2)
  
  #set.seed(k+1)
  x1[u1>cut1[3]]=rnorm(n14,mu1,s2)
  #set.seed(k+2)
  y1[u1>cut1[3]]=rnorm(n14,mu1,s2)
  
  #simulate t for Gp=0
  x0=rep(0,n0)
  y0=rep(0,n0)
  
  n01=length(u0[u0<=cut0[1]])
  n02=length(u0[u0>cut0[1] & u0<=cut0[2]])
  n03=length(u0[u0>cut0[2] & u0<=cut0[3]])
  n04=length(u0[u0>cut0[3]])
  
  #set.seed(k+1)
  x0[u0<cut0[1]]=rnorm(n01,mu0,s2)
  #set.seed(k+2)
  y0[u0<cut0[1]]=rnorm(n01,mu0,s2)
  
  #set.seed(k+1)
  x0[u0>cut0[1] & u0<=cut0[2]]=rnorm(n02,mu0,s2)
  #set.seed(k+2)
  y0[u0>cut0[1] & u0<=cut0[2]]=rnorm(n02,mu1,s2)
  
  #set.seed(k+1)
  x0[u0>cut0[2] & u0<=cut0[3]]=rnorm(n03,mu1,s2)
  #set.seed(k+2)
  y0[u0>cut0[2] & u0<=cut0[3]]=rnorm(n03,mu0,s2)
  
  #set.seed(k+1)
  x0[u0>cut0[3]]=rnorm(n04,mu1,s2)
  #set.seed(k+2)
  y0[u0>cut0[3]]=rnorm(n04,mu1,s2)
  
  t=c(x0,x1,y0,y1,z0,z1)
  t0 = c(x0,y0)
  t1 = c(x1,y1)
  
  return(list(t,x0,y0,x1,y1,c(n0,n1),z0,z1))
}

#### semi-parametric approach
obj.semiparametric.pairs = function(par,t.list,n){
  #t.list = t.values[[1]]; par = ini.val.ab
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
  #### mother part
  #lambda for Gp=1
  lambda1_1m.par = (-a.par^3+a.par^2+a.par)/(a.par*(2-a.par))
  lambda1_0m.par = 1-lambda1_1m.par
  cut1m.par = c(lambda1_0m.par,lambda1_1m.par)
  
  #lambda for Gp=0
  lambda0_1m.par = a.par
  lambda0_0m.par = 1-lambda0_1m.par
  cut0m.par = c(lambda0_0m.par,lambda0_1m.par)
  
  #### sister part
  #lambda for Gp=1
  lambda1_1s.par = (a.par^4-6*a.par^3+5*a.par^2+4*a.par)/(4*a.par*(2-a.par))
  lambda1_0s.par = 1-lambda1_1s.par
  cut1s.par = c(lambda1_0s.par,lambda1_1s.par)
  
  #lambda for Gp=0
  lambda0_1s.par = (-a.par^4+6*a.par^3-9*a.par^2+4*a.par)/(4*(1-a.par)^2)
  lambda0_0s.par = 1-lambda0_1s.par
  cut0s.par = c(lambda0_0s.par,lambda0_1s.par)
  
  #log likelihood part for Gp=0
  t0=(cut0m.par[1]+cut0m.par[2]*exp(alp+x0*bet))*(cut0s.par[1]+cut0s.par[2]*exp(alp+y0*bet))
  #log likelihood part for Gp=1
  t1=(cut1m.par[1]+cut1m.par[2]*exp(alp+x1*bet))*(cut1s.par[1]+cut1s.par[2]*exp(alp+y1*bet))
  tp=alp+z1*bet
  
  value= n1*log(theta.par)+n0*log(1-theta.par)+sum(log(t1))+sum(log(t0))+
    sum(tp)-sum(log(1+tt*te))
  
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
      te=exp(alp+bet*t.values[[j]][[1]])-1
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
c(mu2beta.f(c(param.val$true.mu0,param.val$true.mu1,param.val$true.sigma,param.val$a)),round(colMeans(est.semiparametric.pairs,na.rm = T),2)[4:6],round(apply(est.semiparametric.pairs,2, function(x) sd(x,na.rm=T)),3)[4:6])
write.csv(est.semiparametric.pairs,"aUnknown_est.semiparametric.pairs_proband_genotype_n50_beta1_a0.2.csv")
