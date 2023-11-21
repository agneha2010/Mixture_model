est = read.csv("aUnknown_est.semiparametric_proband_genotype_n1000_beta-1_a0.1.csv")
est.parametric.trio = read.csv("aUnknown_est.parametric.trio_proband_genotype_n1000_beta-1_a0.1.csv")
est.parametric.pairs = read.csv("aUnknown_est.parametric.pairs_proband_genotype_n1000_beta-1_a0.1.csv")
est.parametric.just.moth = read.csv("aUnknown_est.parametric.just.moth_proband_genotype_n1000_beta-1_a0.1.csv")

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


