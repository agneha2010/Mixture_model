mm = function(l){
   # l = c("aUnknown_n5000_beta0.3_a0.1/aUnknown_est.semiparametric_theta0_genotype_n5000_beta0.3_a0.1.csv",
   #                        "aUnknown_n5000_beta0.3_a0.1/aUnknown_est.parametric.trio_theta0_genotype_n5000_beta0.3_a0.1.csv",
   #                        "aUnknown_n5000_beta0.3_a0.1/aUnknown_est.parametric.pairs_theta0_genotype_n5000_beta0.3_a0.1.csv",
   #                        "aUnknown_n5000_beta0.3_a0.1/aUnknown_est.parametric.just.moth_theta0_genotype_n5000_beta0.3_a0.1.csv")
  est = read.csv(l[1])
  est.parametric.trio = read.csv(l[2])
  est.parametric.pairs = read.csv(l[3])
  est.parametric.just.moth = read.csv(l[4])
  
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
  
  ### find non-convergent from both and remove
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
  
  beta.mean = round(c(colMeans(est1,na.rm = T)[5],colMeans(est.converted.trio,na.rm = T)[2],
                      colMeans(est.converted.pairs,na.rm = T)[2],colMeans(est.converted.just.moth,na.rm = T)[2]),2)
  
  beta.sd = round(c(apply(est1,2, function(x) sd(x,na.rm=T))[5],apply(est.converted.trio,2,function(x) sd(x,na.rm=T))[2],
                    apply(est.converted.pairs, 2, function(x) sd(x,na.rm=T))[2],
                    apply(est.converted.just.moth,2,function(x) sd(x,na.rm=T))[2]),2)
  
  
  a.mean = round(c(colMeans(est1,na.rm = T)[6],colMeans(est.converted.trio,na.rm = T)[3],
                   colMeans(est.converted.pairs,na.rm = T)[3],colMeans(est.converted.just.moth,na.rm = T)[3]),2)
  
  a.sd = round(c(apply(est1,2, function(x) sd(x,na.rm=T))[6],apply(est.converted.trio,2,function(x) sd(x,na.rm=T))[3],
                 apply(est.converted.pairs, 2, function(x) sd(x,na.rm=T))[3],
                 apply(est.converted.just.moth,2,function(x) sd(x,na.rm=T))[3]),3)
  
  return(c(paste0(paste(beta.mean[1]),paste("("),paste(beta.sd[1]),paste(")")),
           paste0(paste(beta.mean[2]),paste("("),paste(beta.sd[2]),paste(")")),
           paste0(paste(beta.mean[3]),paste("("),paste(beta.sd[3]),paste(")")),
           paste0(paste(beta.mean[4]),paste("("),paste(beta.sd[4]),paste(")")),
           paste0(paste(a.mean[1]),paste("("),paste(a.sd[1]),paste(")")),
           paste0(paste(a.mean[2]),paste("("),paste(a.sd[2]),paste(")")),
           paste0(paste(a.mean[3]),paste("("),paste(a.sd[3]),paste(")")),
           paste0(paste(a.mean[4]),paste("("),paste(a.sd[4]),paste(")"))))
}


#### creating table for non-carriers increasing sample size
df1 = noquote(mm(c("aUnknown_n1000_beta0.3_a0.1_51_initial_values/aUnknown_est.semiparametric_theta0_genotype_n1000_beta0.3_a0.1.csv",
                   "aUnknown_n1000_beta0.3_a0.1_51_initial_values/aUnknown_est.parametric.trio_theta0_genotype_n1000_beta0.3_a0.1.csv",
                   "aUnknown_n1000_beta0.3_a0.1_51_initial_values/aUnknown_est.parametric.pairs_theta0_genotype_n1000_beta0.3_a0.1.csv",
                   "aUnknown_n1000_beta0.3_a0.1_51_initial_values/aUnknown_est.parametric.just.moth_theta0_genotype_n1000_beta0.3_a0.1.csv")))

df2 = noquote(mm(c("aUnknown_n5000_beta0.3_a0.1/aUnknown_est.semiparametric_theta0_genotype_n5000_beta0.3_a0.1.csv",
                   "aUnknown_n5000_beta0.3_a0.1/aUnknown_est.parametric.trio_theta0_genotype_n5000_beta0.3_a0.1.csv",
                   "aUnknown_n5000_beta0.3_a0.1/aUnknown_est.parametric.pairs_theta0_genotype_n5000_beta0.3_a0.1.csv",
                   "aUnknown_n5000_beta0.3_a0.1/aUnknown_est.parametric.just.moth_theta0_genotype_n5000_beta0.3_a0.1.csv")))

df3 = noquote(mm(c("aUnknown_n5000_beta0.3_a0.1_51_initial_values/aUnknown_est.semiparametric_theta0_genotype_n5000_beta0.3_a0.1.csv",
                   "aUnknown_n5000_beta0.3_a0.1_51_initial_values/aUnknown_est.parametric.trio_theta0_genotype_n5000_beta0.3_a0.1.csv",
                   "aUnknown_n5000_beta0.3_a0.1_51_initial_values/aUnknown_est.parametric.pairs_theta0_genotype_n5000_beta0.3_a0.1.csv",
                   "aUnknown_n5000_beta0.3_a0.1_51_initial_values/aUnknown_est.parametric.just.moth_theta0_genotype_n5000_beta0.3_a0.1.csv")))

df4 = noquote(mm(c("aUnknown_n10000_beta0.3_a0.1/aUnknown_est.semiparametric_theta0_genotype_n10000_beta0.3_a0.1.csv",
                   "aUnknown_n10000_beta0.3_a0.1/aUnknown_est.parametric.trio_theta0_genotype_n10000_beta0.3_a0.1.csv",
                   "aUnknown_n10000_beta0.3_a0.1/aUnknown_est.parametric.pairs_theta0_genotype_n10000_beta0.3_a0.1.csv",
                   "aUnknown_n10000_beta0.3_a0.1/aUnknown_est.parametric.just.moth_theta0_genotype_n10000_beta0.3_a0.1.csv")))

df5 = noquote(mm(c("aUnknown_n10000_beta0.3_a0.1_51_initial_values/aUnknown_est.semiparametric_theta0_genotype_n10000_beta0.3_a0.1.csv",
                   "aUnknown_n10000_beta0.3_a0.1_51_initial_values/aUnknown_est.parametric.trio_theta0_genotype_n10000_beta0.3_a0.1.csv",
                   "aUnknown_n10000_beta0.3_a0.1_51_initial_values/aUnknown_est.parametric.pairs_theta0_genotype_n10000_beta0.3_a0.1.csv",
                   "aUnknown_n10000_beta0.3_a0.1_51_initial_values/aUnknown_est.parametric.just.moth_theta0_genotype_n10000_beta0.3_a0.1.csv")))

dt = as.matrix(noquote(cbind(rep("Non-carriers",5),rbind(df1,df2,df3,df4,df5))))
colnames(dt) = c("Sample Type","Trio","Trio","Pairs","Sister","Trio","Trio","Pairs","Sister")
rownames(dt) = c(1,2,3,4,5)
dt


options(kableExtra.latex.load_packages = FALSE) 
library(kableExtra)
kbl(dt,"latex")
kbl(dt, booktabs = T, linesep = "","latex",align="c") %>% 
  kable_styling(latex_options = c("striped")) %>% 
  add_header_above(c(" "=1,"semi parametric" = 1, "parametric" = 3, "semi parametric" = 1, "parametric" = 3)) %>% 
  add_header_above(c(" "=1,"$\\beta$" = 4, "a" = 4)) %>% 
  kable_styling(position = "center") %>% 
  column_spec(c(2,6), border_left = T) %>%
  pack_rows("n = 1000; k = 50", 1,1) %>%
  pack_rows("n = 5000; k = 20", 2,2) %>%
  pack_rows("n = 5000; k = 50", 3,3) %>%
  pack_rows("n = 10000; k = 20", 4,4) %>%
  pack_rows("n = 10000; k = 50", 5,5) %>%
  column_spec(1, width = "4em")

kbl(dt, booktabs = T, linesep = "",align="c") %>% 
  kable_styling(latex_options = c("striped")) %>% 
  add_header_above(c(" "=1,"semi parametric" = 1, "parametric" = 3, "semi parametric" = 1, "parametric" = 3)) %>% 
  add_header_above(c(" "=1,"beta" = 4, "a" = 4)) %>% 
  kable_styling(position = "center") %>%
  column_spec(2, border_left = T) %>%
  column_spec(6, border_left = T) %>%
  pack_rows("n = 1000; k = 50", 1,1) %>%
  pack_rows("n = 5000; k = 20", 2,2) %>%
  pack_rows("n = 5000; k = 50", 3,3) %>%
  pack_rows("n = 10000; k = 20", 4,4) %>%
  pack_rows("n = 10000; k = 50", 5,5) 
