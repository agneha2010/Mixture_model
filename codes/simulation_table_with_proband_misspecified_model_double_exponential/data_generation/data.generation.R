library(dplyr)
library(invgamma)
library(ExtDist)

param.val = list(n = 1000,
                 a = 0.1, ### allele frequency
                 mu0 = 0,
                 mu1 = 1,
                 b = 1,
                 true.beta = 1
)

param.val$true.alpha = 0

param.val$alpha.ini.lb = param.val$true.alpha-2
param.val$alpha.ini.ub = param.val$true.alpha+2
param.val$beta.ini.lb = param.val$true.beta-2
param.val$beta.ini.ub = param.val$true.beta+2
param.val$mu0.ini.lb = param.val$true.mu0-2
param.val$mu0.ini.ub = param.val$true.mu0+2
param.val$mu1.ini.lb = param.val$true.mu1-2
param.val$mu1.ini.ub = param.val$true.mu1+2


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
  
  #k = 3*(seed-1)+1
  
  #set.seed(k)
  Gp = rbinom(n,1,theta)
  n1 = length(Gp[Gp==1])
  n0 = length(Gp[Gp==0])
  
  ### data for proband 
  #z0 = rLaplace(n=n0,param.val$mu0,param.val$b)
  
  z0 = rnorm(n=n0,param.val$mu0,param.val$b)
  z = rnorm(n=10000,param.val$mu0,param.val$b)
  
  #z = rLaplace(n=10000,param.val$mu0,param.val$b)
  p = exp(param.val$true.beta*z)/sum(exp(param.val$true.beta*z))
  D = rmultinom(n1,1,p)
  id = apply(D,2,function(x) which(x==1))
  z1 = z[id]
  
  
  library(dplyr)
  library(ggplot2)
  df = data.frame(data=c(z0,z1),cat=as.factor(c(rep(0,length(z0)),rep(1,length(z1)))))
  ggplot(df, aes(data, fill = cat)) + geom_histogram(alpha = 0.2)
  