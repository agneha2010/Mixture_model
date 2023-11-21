library(dplyr)
library(invgamma)

param.val = list(n = 1000,
                 a = 0.1,
                 true.mu0 = 1,
                 true.mu1 = 0,
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

data.generation = function(param.val){
  a = param.val$a
  n = param.val$n
  mu0 = param.val$true.mu0
  mu1 = param.val$true.mu1
  s2 = param.val$true.sigma
  
  ## theta is probability proband genotype is 1
  theta = 2*a-a^2
  
  #### mother part
  #lambda for Gp=1
  lambda1_1m = (-a^3+a^2+a)/(a*(2-a))
  lambda1_0m = 1-lambda1_1m
  cut1m = c(lambda1_0m,lambda1_1m)
  
  #lambda for Gp=0
  lambda0_1m = a
  lambda0_0m = 1-lambda0_1m
  cut0m = c(lambda0_0m,lambda0_1m)
  
  #### sister part
  #lambda for Gp=1
  lambda1_1s = (a^4-6*a^3+5*a^2+4*a)/(4*a*(2-a))
  lambda1_0s = 1-lambda1_1s
  cut1s = c(lambda1_0s,lambda1_1s)
  
  #lambda for Gp=0
  lambda0_1s = (-a^4+6*a^3-9*a^2+4*a)/(4*(1-a)^2)
  lambda0_0s = 1-lambda0_1s
  cut0s = c(lambda0_0s,lambda0_1s)
  
  #### daughter part
  #lambda for Gp=1
  lambda1_1d = 1-(a*(1-a)^2/(a*(2-a)))
  lambda1_0d = 1-lambda1_1d
  cut1d = c(lambda1_0d,lambda1_1d)
  
  #lambda for Gp=0
  lambda0_1d = (1-a)
  lambda0_0d = 1-lambda0_1d
  cut0d = c(lambda0_0d,lambda0_1d)
    
  ## genotype data for proband
  Gp = rbinom(n,1,theta)
  n1 = length(Gp[Gp==1])
  n0 = length(Gp[Gp==0])
  
  ### phenotype data for proband 
  p0 = rnorm(n0,mu0,s2)
  p1 = rnorm(n1,mu1,s2)
  
  ### data for proband-mother-sister
  u1=runif(n1)
  u0=runif(n0)
  
  #### relation proband-mother-sister
  ### mothers and sisters of probands carriers
  nm10 = length(u1[u1<=cut1m[1]])
  nm11 = n1-nm10
  ns10 = length(u1[u1<=cut1s[1]])
  ns11 = n1-ns10
  
  m1=rep(0,n1)
  s1=rep(0,n1)
  
  m1[u1<=cut1m[1]] = rnorm(nm10,mu0,s2)
  s1[u1<=cut1s[1]] = rnorm(ns10,mu0,s2)
  
  m1[u1>cut1m[1]] = rnorm(nm11,mu1,s2)
  s1[u1>cut1s[1]] = rnorm(ns11,mu1,s2)
  
  ### mothers and sisters of probands non-carriers
  nm00 = length(u0[u0<=cut0m[1]])
  nm01 = n0-nm00
  ns00 = length(u0[u0<=cut0s[1]])
  ns01 = n0-ns00
  
  m0=rep(0,n0)
  s0=rep(0,n0)
  
  m0[u0<=cut0m[1]] = rnorm(nm00,mu0,s2)
  s0[u0<=cut0s[1]] = rnorm(ns00,mu0,s2)
  
  m0[u0>cut0m[1]] = rnorm(nm01,mu1,s2)
  s0[u0>cut0s[1]] = rnorm(ns01,mu1,s2)
  
  r_pms = data.frame(r="pms",g=c(Gp[Gp==0],Gp[Gp==1]),p=c(p0,p1),r1=c(m0,m1),r2=c(s0,s1))
  
  ### data for proband-mother-sister
  u1=runif(n1)
  u0=runif(n0)
  
  #### relation proband-daughter-sister
  ### daughters and sisters of probands carriers
  nd10 = length(u1[u1<=cut1d[1]])
  nd11 = n1-nd10
  ns10 = length(u1[u1<=cut1s[1]])
  ns11 = n1-ns10
  
  d1=rep(0,n1)
  s1=rep(0,n1)
  
  d1[u1<=cut1d[1]] = rnorm(nd10,mu0,s2)
  s1[u1<=cut1s[1]] = rnorm(ns10,mu0,s2)
  
  d1[u1>cut1d[1]] = rnorm(nd11,mu1,s2)
  s1[u1>cut1s[1]] = rnorm(ns11,mu1,s2)
  
  ### daughters and sisters of probands non-carriers
  nd00 = length(u0[u0<=cut0d[1]])
  nd01 = n0-nd00
  ns00 = length(u0[u0<=cut0s[1]])
  ns01 = n0-ns00
  
  d0=rep(0,n0)
  s0=rep(0,n0)
  
  d0[u0<=cut0d[1]] = rnorm(nd00,mu0,s2)
  s0[u0<=cut0s[1]] = rnorm(ns00,mu0,s2)
  
  d0[u0>cut0d[1]] = rnorm(nd01,mu1,s2)
  s0[u0>cut0s[1]] = rnorm(ns01,mu1,s2)
  
  r_pds = data.frame(r="pds",g=c(Gp[Gp==0],Gp[Gp==1]),p=c(p0,p1),r1=c(d0,d1),r2=c(s0,s1))
  
  ### data for proband-sister1-sister2
  us1_1=runif(n1)
  us1_0=runif(n0)
  
  us2_1=runif(n1)
  us2_0=runif(n0)
  
  #### relation proband-sister-sister
  ### sister1 and sister2 of probands carriers
  ns1_10 = length(us1_1[us1_1<=cut1s[1]])
  ns1_11 = n1-ns1_10
  ns2_10 = length(us2_1[us2_1<=cut1s[1]])
  ns2_11 = n1-ns2_10
  
  s1_1=rep(0,n1)
  s2_1=rep(0,n1)
  
  s1_1[us1_1<=cut1s[1]] = rnorm(ns1_10,mu0,s2)
  s2_1[us2_1<=cut1s[1]] = rnorm(ns2_10,mu0,s2)
  
  s1_1[us1_1>cut1s[1]] = rnorm(ns1_11,mu1,s2)
  s2_1[us2_1>cut1s[1]] = rnorm(ns2_11,mu1,s2)
  
  ### mothers and sisters of probands non-carriers
  ns1_00 = length(us1_0[us1_0<=cut0s[1]])
  ns1_01 = n0-ns1_00
  ns2_00 = length(us2_0[us2_0<=cut0s[1]])
  ns2_01 = n0-ns2_00
  
  s1_0=rep(0,n0)
  s2_0=rep(0,n0)
  
  s1_0[us1_0<=cut0s[1]] = rnorm(ns1_00,mu0,s2)
  s2_0[us2_0<=cut0s[1]] = rnorm(ns2_00,mu0,s2)
  
  s1_0[us1_0>cut0s[1]] = rnorm(ns1_01,mu1,s2)
  s2_0[us2_0>cut0s[1]] = rnorm(ns2_01,mu1,s2)
  
  r_pss = data.frame(r="pss",g=c(Gp[Gp==0],Gp[Gp==1]),p=c(p0,p1),r1=c(s1_0,s1_1),r2=c(s2_0,s2_1))
  
  ## number of trios for proband-mother-sister, proband-daughter-sister, proband-sister1-sister2
  z0 = rmultinom(1,n0,c(0.4,0.3,0.4))
  z1 = rmultinom(1,n1,c(0.4,0.3,0.4))
  
  mydata0 = rbind(r_pms[1:z0[1],], r_pds[1:z0[2],], r_pss[1:z0[3],])
  mydata1 = rbind(r_pms[(n0+1):(n0+1+z1[1]-1),], r_pds[(n0+1):(n0+1+z1[2]-1),], r_pss[(n0+1):(n0+1+z1[3]-1),])
  mydata = rbind(mydata0,mydata1)
  
  m0_pms = mydata$r1[mydata$r=="pms" & mydata$g==0]
  m1_pms = mydata$r1[mydata$r=="pms" & mydata$g==1]
  s0_pms = mydata$r2[mydata$r=="pms" & mydata$g==0]
  s1_pms = mydata$r2[mydata$r=="pms" & mydata$g==1]
  
  d0_pds = mydata$r1[mydata$r=="pds" & mydata$g==0]
  d1_pds = mydata$r1[mydata$r=="pds" & mydata$g==1]
  s0_pds = mydata$r2[mydata$r=="pds" & mydata$g==0]
  s1_pds = mydata$r2[mydata$r=="pds" & mydata$g==1]
  
  s10_pss = mydata$r1[mydata$r=="pss" & mydata$g==0]
  s11_pss = mydata$r1[mydata$r=="pss" & mydata$g==1]
  s20_pss = mydata$r2[mydata$r=="pss" & mydata$g==0]
  s21_pss = mydata$r2[mydata$r=="pss" & mydata$g==1]
  
  t=c(m0_pms,d0_pds,s10_pss,m1_pms,d1_pds,s11_pss,
      s0_pms,s0_pds,s20_pss,s1_pms,s1_pds,s21_pss,
      p0,p1)
  
  return(list(t,m0_pms,d0_pds,s10_pss,m1_pms,d1_pds,s11_pss,
              s0_pms,s0_pds,s20_pss,s1_pms,s1_pds,s21_pss,
              p0,p1,c(n0,n1)))
}

#### semi-parametric approach
obj.semiparametric.pairs = function(par,t.list,n){
  #t.list = t.values[[1]]; par = ini.val.ab
  t = t.list[[1]]
  m0_pms = t.list[[2]]
  d0_pds = t.list[[3]]
  s10_pss = t.list[[4]]
  m1_pms  = t.list[[5]]
  d1_pds  = t.list[[6]]
  s11_pss = t.list[[7]]
  s0_pms = t.list[[8]]
  s0_pds = t.list[[9]]
  s20_pss = t.list[[10]]
  s1_pms  = t.list[[11]]
  s1_pds  = t.list[[12]]
  s21_pss = t.list[[13]]
  p0  = t.list[[14]]
  p1 = t.list[[15]]
  n0 = t.list[[16]][1]
  n1 = t.list[[16]][2]
  
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
  
  #### daughter part
  #lambda for Gp=1
  lambda1_1d.par = 1-(a.par*(1-a.par)^2/(a.par*(2-a.par)))
  lambda1_0d.par = 1-lambda1_1d.par
  cut1d.par = c(lambda1_0d.par,lambda1_1d.par)
  
  #lambda for Gp=0
  lambda0_1d.par = (1-a.par)
  lambda0_0d.par = 1-lambda0_1d.par
  cut0d.par = c(lambda0_0d.par,lambda0_1d.par)
  
  #proband-mother-sister
  #log likelihood part for Gp=0
  t0ms=(cut0m.par[1]+cut0m.par[2]*exp(alp+m0_pms*bet))*(cut0s.par[1]+cut0s.par[2]*exp(alp+s0_pms*bet))
  #log likelihood part for Gp=1
  t1ms=(cut1m.par[1]+cut1m.par[2]*exp(alp+m1_pms*bet))*(cut1s.par[1]+cut1s.par[2]*exp(alp+s1_pms*bet))
  
  #proband-daughter-sister
  #log likelihood part for Gp=0
  t0ds=(cut0d.par[1]+cut0d.par[2]*exp(alp+d0_pds*bet))*(cut0s.par[1]+cut0s.par[2]*exp(alp+s0_pds*bet))
  #log likelihood part for Gp=1
  t1ds=(cut1d.par[1]+cut1d.par[2]*exp(alp+d1_pds*bet))*(cut1s.par[1]+cut1s.par[2]*exp(alp+s1_pds*bet))
  
  #proband-sister1-sister2
  #log likelihood part for Gp=0
  t0ss=(cut0s.par[1]+cut0s.par[2]*exp(alp+s10_pss*bet))*(cut0s.par[1]+cut0s.par[2]*exp(alp+s20_pss*bet))
  #log likelihood part for Gp=1
  t1ss=(cut1s.par[1]+cut1s.par[2]*exp(alp+s11_pss*bet))*(cut1s.par[1]+cut1s.par[2]*exp(alp+s21_pss*bet))
  
  tp=alp+p1*bet
  
  value= n1*log(theta.par)+n0*log(1-theta.par)+sum(log(t1ms))+sum(log(t0ms))+sum(log(t1ds))+sum(log(t0ds))+
    sum(log(t1ss))+sum(log(t0ss))+sum(tp)-sum(log(1+tt*te))
  
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
  t.values[[j]] = data.generation(param.val)
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
#write.csv(est.semiparametric.pairs,"aUnknown_est.semiparametric.pairs_proband_genotype_n1000_beta-1_a0.1.csv")
