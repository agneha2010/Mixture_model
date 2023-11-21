library(chemometrics)

est = read.csv("aUnknown_est.semiparametric_proband_genotype_1000c_1000nc_beta0.3_a0.1.csv")
est.parametric.trio = read.csv("aUnknown_est.parametric.trio_proband_genotype_1000c_1000nc_beta0.3_a0.1.csv")
est.parametric.pairs = read.csv("aUnknown_est.parametric.pairs_proband_genotype_1000c_1000nc_beta0.3_a0.1.csv")
est.parametric.just.moth = read.csv("aUnknown_est.parametric.just.moth_proband_genotype_1000c_1000nc_beta0.3_a0.1.csv")

est = read.csv("aUnknown_est.semiparametric_proband_genotype_1000c_1000nc_beta1_a0.1.csv")
est.parametric.trio = read.csv("aUnknown_est.parametric.trio_proband_genotype_1000c_1000nc_beta1_a0.1.csv")
est.parametric.pairs = read.csv("aUnknown_est.parametric.pairs_proband_genotype_1000c_1000nc_beta1_a0.1.csv")
est.parametric.just.moth = read.csv("aUnknown_est.parametric.just.moth_proband_genotype_1000c_1000nc_beta1_a0.1.csv")

est = read.csv("aUnknown_est.semiparametric_proband_genotype_n2000_beta0.3_a0.1.csv")
est.parametric.trio = read.csv("aUnknown_est.parametric.trio_proband_genotype_n2000_beta0.3_a0.1.csv")
est.parametric.pairs = read.csv("aUnknown_est.parametric.pairs_proband_genotype_n2000_beta0.3_a0.1.csv")
est.parametric.just.moth = read.csv("aUnknown_est.parametric.just.moth_proband_genotype_n2000_beta0.3_a0.1.csv")

est = read.csv("aUnknown_est.semiparametric_proband_genotype_n2000_beta1_a0.1.csv")
est.parametric.trio = read.csv("aUnknown_est.parametric.trio_proband_genotype_n2000_beta1_a0.1.csv")
est.parametric.pairs = read.csv("aUnknown_est.parametric.pairs_proband_genotype_n2000_beta1_a0.1.csv")
est.parametric.just.moth = read.csv("aUnknown_est.parametric.just.moth_proband_genotype_n2000_beta1_a0.1.csv")

est = est[,-1]
est.parametric.trio = est.parametric.trio[,-1]
est.parametric.pairs = est.parametric.pairs[,-1]
est.parametric.just.moth = est.parametric.just.moth[,-1]
est$a.ini = est$a.ini^2/(1+est$a.ini^2)
est.parametric.trio$a.ini = est.parametric.trio$a.ini^2/(1+est.parametric.trio$a.ini^2)
est.parametric.pairs$a.ini = est.parametric.pairs$a.ini^2/(1+est.parametric.pairs$a.ini^2)
est.parametric.just.moth$a.ini = est.parametric.just.moth$a.ini^2/(1+est.parametric.just.moth$a.ini^2)

mu2beta.f = function(par){
  mu0 = par[1]
  mu1 = par[2]
  sigma = par[3]
  a = par[4]
  beta = (mu1-mu0)/sigma^2
  alpha = (-mu1^2+mu0^2)/(2*sigma^2)
  return(c(alpha,beta,a))
}
#mu2beta.f(c(2,1.7,1,0.3))

#### find non-convergent from both and remove
ind.est = which(est[,9] !=0)
ind.parametric.est.trio = which(est.parametric.trio[,11] !=0)
ind.parametric.est.pairs = which(est.parametric.pairs[,11] !=0)
ind.parametric.est.just.moth = which(est.parametric.just.moth[,11] !=0)
ind = unique(c(ind.est,ind.parametric.est.trio,ind.parametric.est.pairs,ind.parametric.est.just.moth))
if(length(ind)==0){
  est1 = est
  est.parametric1.trio = est.parametric.trio
  est.parametric1.pairs = est.parametric.pairs
  est.parametric1.just.moth = est.parametric.just.moth
}else{
  est1 = est[-ind,]
  est.parametric1.trio = est.parametric.trio[-ind,]
  est.parametric1.pairs = est.parametric.pairs[-ind,]
  est.parametric1.just.moth = est.parametric.just.moth[-ind,]
}
est.converted.trio = t(apply(est.parametric1.trio[,5:8],1,mu2beta.f))
est.converted.pairs = t(apply(est.parametric1.pairs[,5:8],1,mu2beta.f))
est.converted.just.moth = t(apply(est.parametric1.just.moth[,5:8],1,mu2beta.f))

m1 = round(c(colMeans(est1,na.rm = T)[5],colMeans(est.converted.trio,na.rm = T)[2],
             colMeans(est.converted.pairs,na.rm = T)[2],colMeans(est.converted.just.moth,na.rm = T)[2]),2)

m2 = round(c(apply(est1,2, function(x) sd(x,na.rm=T))[5],apply(est.converted.trio,2,function(x) sd(x,na.rm=T))[2],
             apply(est.converted.pairs, 2, function(x) sd(x,na.rm=T))[2],
             apply(est.converted.just.moth,2,function(x) sd(x,na.rm=T))[2]),2)

c(paste0(paste(m1[1]),paste("("),paste(m2[1]),paste(")")),
  paste0(paste(m1[2]),paste("("),paste(m2[2]),paste(")")),
  paste0(paste(m1[3]),paste("("),paste(m2[3]),paste(")")),
  paste0(paste(m1[4]),paste("("),paste(m2[4]),paste(")")))

m1 = round(c(colMeans(est1,na.rm = T)[6],colMeans(est.converted.trio,na.rm = T)[3],
             colMeans(est.converted.pairs,na.rm = T)[3],colMeans(est.converted.just.moth,na.rm = T)[3]),2)

m2 = round(c(apply(est1,2, function(x) sd(x,na.rm=T))[6],apply(est.converted.trio,2,function(x) sd(x,na.rm=T))[3],
             apply(est.converted.pairs, 2, function(x) sd(x,na.rm=T))[3],
             apply(est.converted.just.moth,2,function(x) sd(x,na.rm=T))[3]),3)

c(paste0(paste(m1[1]),paste("("),paste(m2[1]),paste(")")),
  paste0(paste(m1[2]),paste("("),paste(m2[2]),paste(")")),
  paste0(paste(m1[3]),paste("("),paste(m2[3]),paste(")")),
  paste0(paste(m1[4]),paste("("),paste(m2[4]),paste(")")))



### semi-parametric
est1 = est[-ind,]
#est1 = est
colMeans(est1,na.rm = T)
apply(est1,2, function(x) sd(x,na.rm=T))
cor = cor(est1[,4:6])
hist(est1[,4],main = "Histogram of alpha, true alpha = 2",xlab="alpha")
hist(est1[,5],main = "Histogram of beta, true beta = -2",xlab="beta")
hist(est1[,6],main = "Histogram of a, true a = 0.3",xlab="a")

#apply(est.converted.trio,2, function(x) mean(x,trim=0.5,na.rm=T))

### parametric both relatives 
### mu,sigma
est.parametric1.trio = est.parametric.trio[-ind,]
#est.parametric1.trio = est.parametric.trio
colMeans(est.parametric1.trio[,5:8],na.rm = T)
apply(est.parametric1.trio[,5:8], 2, function(x) sd(x,na.rm=T))
cor.parametric = cor(est.parametric1.trio[,5:8])
hist(est.parametric1.trio[,5],main = "Histogram of mu0, true mu0 = 2",xlab="alpha")
hist(est.parametric1.trio[,6],main = "Histogram of mu1, true mu1 = 0",xlab="beta")
hist(est.parametric1.trio[,7],main = "Histogram of sigma, true sigma = 1",xlab="beta")
hist(est.parametric1.trio[,8],main = "Histogram of a, true a = 0.3",xlab="a")

### alpha,beta
est.converted.trio = t(apply(est.parametric1.trio[,5:8],1,mu2beta.f))
colnames(est.converted.trio) = c("alpha.est","beta.est","a.est")
colMeans(est.converted.trio,na.rm = T)
apply(est.converted.trio, 2, function(x) sd(x,na.rm=T))
cor.converted.trio = cor(est.converted.trio)
hist(est.converted.trio[,1],main = "Histogram of alpha, true alpha = 2",xlab="alpha")
hist(est.converted.trio[,2],main = "Histogram of beta, true beta = -2",xlab="beta")
hist(est.converted.trio[,3],main = "Histogram of a, true a = 0.3",xlab="a")

### parametric only pairs 
### mu,sigma
est.parametric1.pairs = est.parametric.pairs[-ind,]
#est.parametric1.pairs = est.parametric.pairs
colMeans(est.parametric1.pairs[,5:8],na.rm = T)
apply(est.parametric1.pairs[,5:8], 2, function(x) sd(x,na.rm=T))
cor.parametric.pairs = cor(est.parametric1.pairs[,5:8])
hist(est.parametric1.pairs[,5],main = "Histogram of mu0, true mu0 = 2",xlab="alpha")
hist(est.parametric1.pairs[,6],main = "Histogram of mu1, true mu1 = 0",xlab="beta")
hist(est.parametric1.pairs[,7],main = "Histogram of sigma, true sigma = 1",xlab="beta")
hist(est.parametric1.pairs[,8],main = "Histogram of a, true a = 0.3",xlab="a")

### alpha,beta
est.converted.pairs = t(apply(est.parametric1.pairs[,5:8],1,mu2beta.f))
colnames(est.converted.pairs) = c("alpha.est","beta.est","a.est")
colMeans(est.converted.pairs,na.rm = T)
apply(est.converted.pairs, 2, function(x) sd(x,na.rm=T))
cor.converted.pairs = cor(est.converted.pairs)
hist(est.converted.pairs[,1],main = "Histogram of alpha, true alpha = 2",xlab="alpha")
hist(est.converted.pairs[,2],main = "Histogram of beta, true beta = -2",xlab="beta")
hist(est.converted.pairs[,3],main = "Histogram of a, true a = 0.3",xlab="a")

### parametric only pairs 
### mu,sigma
est.parametric1.just.moth = est.parametric.just.moth[-ind,]
#est.parametric1.just.moth = est.parametric.just.moth
colMeans(est.parametric1.just.moth[,5:8],na.rm = T)
apply(est.parametric1.just.moth[,5:8], 2, function(x) sd(x,na.rm=T))
cor.parametric.just.moth = cor(est.parametric1.just.moth[,5:8])
hist(est.parametric1.just.moth[,5],main = "Histogram of mu0, true mu0 = 2",xlab="alpha")
hist(est.parametric1.just.moth[,6],main = "Histogram of mu1, true mu1 = 0",xlab="beta")
hist(est.parametric1.just.moth[,7],main = "Histogram of sigma, true sigma = 1",xlab="beta")
hist(est.parametric1.just.moth[,8],main = "Histogram of a, true a = 0.3",xlab="a")

### alpha,beta
est.converted.just.moth = t(apply(est.parametric1.just.moth[,5:8],1,mu2beta.f))
colnames(est.converted.just.moth) = c("alpha.est","beta.est","a.est")
colMeans(est.converted.just.moth,na.rm = T)
apply(est.converted.just.moth, 2, function(x) sd(x,na.rm=T))
cor.converted.just.moth = cor(est.converted.just.moth)
hist(est.converted.just.moth[,1],main = "Histogram of alpha, true alpha = 2",xlab="alpha")
hist(est.converted.just.moth[,2],main = "Histogram of beta, true beta = -2",xlab="beta")
hist(est.converted.just.moth[,3],main = "Histogram of a, true a = 0.3",xlab="a")

round(c(colMeans(est.converted.trio,na.rm = T),colMeans(est.converted.pairs,na.rm = T),colMeans(est.converted.just.moth,na.rm = T),
        colMeans(est1,na.rm = T)[4:6]),3)

round(c(apply(est.converted.trio,2,function(x) sd(x,na.rm=T)),apply(est.converted.pairs, 2, function(x) sd(x,na.rm=T)),
        apply(est.converted.just.moth,2,function(x) sd(x,na.rm=T)),apply(est1,2, function(x) sd(x,na.rm=T))[4:6]),3)

round(c(colMeans(est.parametric1.trio[,5:8],na.rm = T),colMeans(est.parametric1.pairs[,5:8],na.rm = T),
        colMeans(est.parametric1.just.moth[,5:8],na.rm = T)),3)

round(c(apply(est.parametric1.trio[,5:8],2,function(x) sd(x,na.rm=T)),apply(est.parametric1.pairs[,5:8], 2, function(x) sd(x,na.rm=T)),
        apply(est.parametric1.just.moth[,5:8],2,function(x) sd(x,na.rm=T))),3)

apply(est1[,4:6],2,function(x) mean(x,trim=0.8))
apply(est1[,4:6],2,function(x) sd_trim(x,trim=0.8,const=F))

### beta 
m1 = round(c(colMeans(est1,na.rm = T)[5],colMeans(est.converted.trio,na.rm = T)[2],
             colMeans(est.converted.pairs,na.rm = T)[2],colMeans(est.converted.just.moth,na.rm = T)[2]),2)

m2 = round(c(apply(est1,2, function(x) sd(x,na.rm=T))[5],apply(est.converted.trio,2,function(x) sd(x,na.rm=T))[2],
             apply(est.converted.pairs, 2, function(x) sd(x,na.rm=T))[2],
             apply(est.converted.just.moth,2,function(x) sd(x,na.rm=T))[2]),2)

c(paste0(paste(m1[1]),paste("("),paste(m2[1]),paste(")")),
  paste0(paste(m1[2]),paste("("),paste(m2[2]),paste(")")),
  paste0(paste(m1[3]),paste("("),paste(m2[3]),paste(")")),
  paste0(paste(m1[4]),paste("("),paste(m2[4]),paste(")")))

#### a

m1 = round(c(colMeans(est1,na.rm = T)[6],colMeans(est.converted.trio,na.rm = T)[3],
             colMeans(est.converted.pairs,na.rm = T)[3],colMeans(est.converted.just.moth,na.rm = T)[3]),2)

m2 = round(c(apply(est1,2, function(x) sd(x,na.rm=T))[6],apply(est.converted.trio,2,function(x) sd(x,na.rm=T))[3],
             apply(est.converted.pairs, 2, function(x) sd(x,na.rm=T))[3],
             apply(est.converted.just.moth,2,function(x) sd(x,na.rm=T))[3]),3)

c(paste0(paste(m1[1]),paste("("),paste(m2[1]),paste(")")),
  paste0(paste(m1[2]),paste("("),paste(m2[2]),paste(")")),
  paste0(paste(m1[3]),paste("("),paste(m2[3]),paste(")")),
  paste0(paste(m1[4]),paste("("),paste(m2[4]),paste(")")))



