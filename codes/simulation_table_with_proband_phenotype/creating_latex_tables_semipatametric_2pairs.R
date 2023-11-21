mm = function(l){
  # l = c("proband_genotype_prospective/aUnknown_n1333_beta0.3_a0.1/aUnknown_est.semiparametric_proband_genotype_n1333_beta0.3_a0.1.csv",
  #   "proband_genotype_prospective/aUnknown_n1333_beta0.3_a0.1/aUnknown_est.semiparametric.pairs_proband_genotype_n1333_beta0.3_a0.1.csv",
  #   "proband_genotype_prospective/aUnknown_n1333_beta0.3_a0.1/aUnknown_est.parametric.trio_proband_genotype_n1333_beta0.3_a0.1.csv",
  #   "proband_genotype_prospective/aUnknown_n1333_beta0.3_a0.1/aUnknown_est.parametric.pairs_proband_genotype_n1333_beta0.3_a0.1.csv")
  est = read.csv(l[1])
  est.semiparametric.pairs = read.csv(l[2])
  est.parametric.trio = read.csv(l[3])
  est.parametric.pairs = read.csv(l[4])
  
  est = est[,-1]
  est.semiparametric.pairs = est.semiparametric.pairs[,-1]
  est.parametric.trio = est.parametric.trio[,-1]
  est.parametric.pairs = est.parametric.pairs[,-1]
  est$a.ini = est$a.ini^2/(1+est$a.ini^2)
  est.semiparametric.pairs$a.ini = est.semiparametric.pairs$a.ini^2/(1+est.semiparametric.pairs$a.ini^2)
  est.parametric.trio$a.ini = est.parametric.trio$a.ini^2/(1+est.parametric.trio$a.ini^2)
  est.parametric.pairs$a.ini = est.parametric.pairs$a.ini^2/(1+est.parametric.pairs$a.ini^2)
    
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
  ind.semiparametric.est.pairs = which(est.semiparametric.pairs[,9] !=0)
  ind.parametric.est.trio = which(est.parametric.trio[,11] !=0)
  ind.parametric.est.pairs = which(est.parametric.pairs[,11] !=0)
  ind = unique(c(ind.est,ind.semiparametric.est.pairs,ind.parametric.est.trio,ind.parametric.est.pairs))
  
  if(length(ind)==0){
    est1 = est
    est.semiparametric1.pairs = est.semiparametric.pairs
    est.parametric1.trio = est.parametric.trio
    est.parametric1.pairs = est.parametric.pairs
  }else{
    est1 = est[-ind,]
    est.parametric1.trio = est.parametric.trio[-ind,]
    est.parametric1.pairs = est.parametric.pairs[-ind,]
    est.semiparametric1.pairs = est.semiparametric.pairs[-ind,]
  }
  est.converted.trio = t(apply(est.parametric1.trio[,5:8],1,mu2beta.f))
  est.converted.pairs = t(apply(est.parametric1.pairs[,5:8],1,mu2beta.f))

  beta.mean = round(c(colMeans(est1,na.rm = T)[5],
                      colMeans(est.semiparametric1.pairs,na.rm = T)[5],
                      colMeans(est.converted.trio,na.rm = T)[2],
                      colMeans(est.converted.pairs,na.rm = T)[2]),2)
  
  beta.sd = round(c(apply(est1,2, function(x) sd(x,na.rm=T))[5],
                    apply(est.semiparametric1.pairs,2,function(x) sd(x,na.rm=T))[5],
                    apply(est.converted.trio,2,function(x) sd(x,na.rm=T))[2],
                    apply(est.converted.pairs, 2, function(x) sd(x,na.rm=T))[2]),2)
  
  a.mean = round(c(colMeans(est1,na.rm = T)[6],
                   colMeans(est.semiparametric1.pairs,na.rm = T)[6],
                   colMeans(est.converted.trio,na.rm = T)[3],
                   colMeans(est.converted.pairs,na.rm = T)[3]),2)
  
  a.sd = round(c(apply(est1,2, function(x) sd(x,na.rm=T))[6],
                 apply(est.semiparametric1.pairs,2,function(x) sd(x,na.rm=T))[6],
                 apply(est.converted.trio,2,function(x) sd(x,na.rm=T))[3],
                 apply(est.converted.pairs, 2, function(x) sd(x,na.rm=T))[3]),3)
  
  return(c(paste0(paste(beta.mean[1]),paste("("),paste(beta.sd[1]),paste(")")),
           paste0(paste(beta.mean[2]),paste("("),paste(beta.sd[2]),paste(")")),
           paste0(paste(beta.mean[3]),paste("("),paste(beta.sd[3]),paste(")")),
           paste0(paste(beta.mean[4]),paste("("),paste(beta.sd[4]),paste(")")),
           paste0(paste(a.mean[1]),paste("("),paste(a.sd[1]),paste(")")),
           paste0(paste(a.mean[2]),paste("("),paste(a.sd[2]),paste(")")),
           paste0(paste(a.mean[3]),paste("("),paste(a.sd[3]),paste(")")),
           paste0(paste(a.mean[4]),paste("("),paste(a.sd[4]),paste(")"))))
}
df11 = noquote(mm(c("theta1/aUnknown_n50_beta-1_a0.1/aUnknown_est.semiparametric_theta1_genotype_n50_beta-1_a0.1.csv",
                    "theta1/aUnknown_n50_beta-1_a0.1/aUnknown_est.semiparametric.pairs_theta1_genotype_n50_beta-1_a0.1.csv",
                    "theta1/aUnknown_n50_beta-1_a0.1/aUnknown_est.parametric.trio_theta1_genotype_n50_beta-1_a0.1.csv",
                    "theta1/aUnknown_n50_beta-1_a0.1/aUnknown_est.parametric.pairs_theta1_genotype_n50_beta-1_a0.1.csv")))

df12 = noquote(mm(c("theta0/aUnknown_n50_beta-1_a0.1/aUnknown_est.semiparametric_theta0_genotype_n50_beta-1_a0.1.csv",
                    "theta0/aUnknown_n50_beta-1_a0.1/aUnknown_est.semiparametric.pairs_theta0_genotype_n50_beta-1_a0.1.csv",
                    "theta0/aUnknown_n50_beta-1_a0.1/aUnknown_est.parametric.trio_theta0_genotype_n50_beta-1_a0.1.csv",
                    "theta0/aUnknown_n50_beta-1_a0.1/aUnknown_est.parametric.pairs_theta0_genotype_n50_beta-1_a0.1.csv")))

df13 = noquote(mm(c("proband_genotype_prospective/aUnknown_n50_beta-1_a0.1/aUnknown_est.semiparametric_proband_genotype_n50_beta-1_a0.1.csv",
                    "proband_genotype_prospective/aUnknown_n50_beta-1_a0.1/aUnknown_est.semiparametric.pairs_proband_genotype_n50_beta-1_a0.1.csv",
                    "proband_genotype_prospective/aUnknown_n50_beta-1_a0.1/aUnknown_est.parametric.trio_proband_genotype_n50_beta-1_a0.1.csv",
                    "proband_genotype_prospective/aUnknown_n50_beta-1_a0.1/aUnknown_est.parametric.pairs_proband_genotype_n50_beta-1_a0.1.csv")))

df21 = noquote(mm(c("theta1/aUnknown_n1000_beta-1_a0.1/aUnknown_est.semiparametric_theta1_genotype_n1000_beta-1_a0.1.csv",
                    "theta1/aUnknown_n1000_beta-1_a0.1/aUnknown_est.semiparametric.pairs_theta1_genotype_n1000_beta-1_a0.1.csv",
                    "theta1/aUnknown_n1000_beta-1_a0.1/aUnknown_est.parametric.trio_theta1_genotype_n1000_beta-1_a0.1.csv",
                    "theta1/aUnknown_n1000_beta-1_a0.1/aUnknown_est.parametric.pairs_theta1_genotype_n1000_beta-1_a0.1.csv")))

df22 = noquote(mm(c("theta0/aUnknown_n1000_beta-1_a0.1/aUnknown_est.semiparametric_theta0_genotype_n1000_beta-1_a0.1.csv",
                    "theta0/aUnknown_n1000_beta-1_a0.1/aUnknown_est.semiparametric.pairs_theta0_genotype_n1000_beta-1_a0.1.csv",
                    "theta0/aUnknown_n1000_beta-1_a0.1/aUnknown_est.parametric.trio_theta0_genotype_n1000_beta-1_a0.1.csv",
                    "theta0/aUnknown_n1000_beta-1_a0.1/aUnknown_est.parametric.pairs_theta0_genotype_n1000_beta-1_a0.1.csv")))

df23 = noquote(mm(c("proband_genotype_prospective/aUnknown_n1000_beta-1_a0.1/aUnknown_est.semiparametric_proband_genotype_n1000_beta-1_a0.1.csv",
                    "proband_genotype_prospective/aUnknown_n1000_beta-1_a0.1/aUnknown_est.semiparametric.pairs_proband_genotype_n1000_beta-1_a0.1.csv",
                    "proband_genotype_prospective/aUnknown_n1000_beta-1_a0.1/aUnknown_est.parametric.trio_proband_genotype_n1000_beta-1_a0.1.csv",
                    "proband_genotype_prospective/aUnknown_n1000_beta-1_a0.1/aUnknown_est.parametric.pairs_proband_genotype_n1000_beta-1_a0.1.csv")))

df31 = noquote(mm(c("theta1/aUnknown_n50_beta1_a0.2/aUnknown_est.semiparametric_theta1_genotype_n50_beta1_a0.2.csv",
                    "theta1/aUnknown_n50_beta1_a0.2/aUnknown_est.semiparametric.pairs_theta1_genotype_n50_beta1_a0.2.csv",
                    "theta1/aUnknown_n50_beta1_a0.2/aUnknown_est.parametric.trio_theta1_genotype_n50_beta1_a0.2.csv",
                    "theta1/aUnknown_n50_beta1_a0.2/aUnknown_est.parametric.pairs_theta1_genotype_n50_beta1_a0.2.csv")))

df32 = noquote(mm(c("theta0/aUnknown_n50_beta1_a0.2/aUnknown_est.semiparametric_theta0_genotype_n50_beta1_a0.2.csv",
                    "theta0/aUnknown_n50_beta1_a0.2/aUnknown_est.semiparametric.pairs_theta0_genotype_n50_beta1_a0.2.csv",
                    "theta0/aUnknown_n50_beta1_a0.2/aUnknown_est.parametric.trio_theta0_genotype_n50_beta1_a0.2.csv",
                    "theta0/aUnknown_n50_beta1_a0.2/aUnknown_est.parametric.pairs_theta0_genotype_n50_beta1_a0.2.csv")))

df33 = noquote(mm(c("proband_genotype_prospective/aUnknown_n50_beta1_a0.2/aUnknown_est.semiparametric_proband_genotype_n50_beta1_a0.2.csv",
                    "proband_genotype_prospective/aUnknown_n50_beta1_a0.2/aUnknown_est.semiparametric.pairs_proband_genotype_n50_beta1_a0.2.csv",
                    "proband_genotype_prospective/aUnknown_n50_beta1_a0.2/aUnknown_est.parametric.trio_proband_genotype_n50_beta1_a0.2.csv",
                    "proband_genotype_prospective/aUnknown_n50_beta1_a0.2/aUnknown_est.parametric.pairs_proband_genotype_n50_beta1_a0.2.csv")))

df41 = noquote(mm(c("theta1/aUnknown_n1000_beta1_a0.2/aUnknown_est.semiparametric_theta1_genotype_n1000_beta1_a0.2.csv",
                    "theta1/aUnknown_n1000_beta1_a0.2/aUnknown_est.semiparametric.pairs_theta1_genotype_n1000_beta1_a0.2.csv",
                    "theta1/aUnknown_n1000_beta1_a0.2/aUnknown_est.parametric.trio_theta1_genotype_n1000_beta1_a0.2.csv",
                    "theta1/aUnknown_n1000_beta1_a0.2/aUnknown_est.parametric.pairs_theta1_genotype_n1000_beta1_a0.2.csv")))

df42 = noquote(mm(c("theta0/aUnknown_n1000_beta1_a0.2/aUnknown_est.semiparametric_theta0_genotype_n1000_beta1_a0.2.csv",
                    "theta0/aUnknown_n1000_beta1_a0.2/aUnknown_est.semiparametric.pairs_theta0_genotype_n1000_beta1_a0.2.csv",
                    "theta0/aUnknown_n1000_beta1_a0.2/aUnknown_est.parametric.trio_theta0_genotype_n1000_beta1_a0.2.csv",
                    "theta0/aUnknown_n1000_beta1_a0.2/aUnknown_est.parametric.pairs_theta0_genotype_n1000_beta1_a0.2.csv")))

df43 = noquote(mm(c("proband_genotype_prospective/aUnknown_n1000_beta1_a0.2/aUnknown_est.semiparametric_proband_genotype_n1000_beta1_a0.2.csv",
                    "proband_genotype_prospective/aUnknown_n1000_beta1_a0.2/aUnknown_est.semiparametric.pairs_proband_genotype_n1000_beta1_a0.2.csv",
                    "proband_genotype_prospective/aUnknown_n1000_beta1_a0.2/aUnknown_est.parametric.trio_proband_genotype_n1000_beta1_a0.2.csv",
                    "proband_genotype_prospective/aUnknown_n1000_beta1_a0.2/aUnknown_est.parametric.pairs_proband_genotype_n1000_beta1_a0.2.csv")))

df51 = noquote(mm(c("theta1/aUnknown_n50_beta0.3_a0.1/aUnknown_est.semiparametric_theta1_genotype_n50_beta0.3_a0.1.csv",
                    "theta1/aUnknown_n50_beta0.3_a0.1/aUnknown_est.semiparametric.pairs_theta1_genotype_n50_beta0.3_a0.1.csv",
                    "theta1/aUnknown_n50_beta0.3_a0.1/aUnknown_est.parametric.trio_theta1_genotype_n50_beta0.3_a0.1.csv",
                    "theta1/aUnknown_n50_beta0.3_a0.1/aUnknown_est.parametric.pairs_theta1_genotype_n50_beta0.3_a0.1.csv")))

df52 = noquote(mm(c("theta0/aUnknown_n50_beta0.3_a0.1/aUnknown_est.semiparametric_theta0_genotype_n50_beta0.3_a0.1.csv",
                    "theta0/aUnknown_n50_beta0.3_a0.1/aUnknown_est.semiparametric.pairs_theta0_genotype_n50_beta0.3_a0.1.csv",
                    "theta0/aUnknown_n50_beta0.3_a0.1/aUnknown_est.parametric.trio_theta0_genotype_n50_beta0.3_a0.1.csv",
                    "theta0/aUnknown_n50_beta0.3_a0.1/aUnknown_est.parametric.pairs_theta0_genotype_n50_beta0.3_a0.1.csv")))

df53 = noquote(mm(c("proband_genotype_prospective/aUnknown_n50_beta0.3_a0.1/aUnknown_est.semiparametric_proband_genotype_n50_beta0.3_a0.1.csv",
                    "proband_genotype_prospective/aUnknown_n50_beta0.3_a0.1/aUnknown_est.semiparametric.pairs_proband_genotype_n50_beta0.3_a0.1.csv",
                    "proband_genotype_prospective/aUnknown_n50_beta0.3_a0.1/aUnknown_est.parametric.trio_proband_genotype_n50_beta0.3_a0.1.csv",
                    "proband_genotype_prospective/aUnknown_n50_beta0.3_a0.1/aUnknown_est.parametric.pairs_proband_genotype_n50_beta0.3_a0.1.csv")))

df61 = noquote(mm(c("theta1/aUnknown_n1333_beta0.3_a0.1/aUnknown_est.semiparametric_theta1_genotype_n1333_beta0.3_a0.1.csv",
                    "theta1/aUnknown_n1333_beta0.3_a0.1/aUnknown_est.semiparametric.pairs_theta1_genotype_n1333_beta0.3_a0.1.csv",
                    "theta1/aUnknown_n1333_beta0.3_a0.1/aUnknown_est.parametric.trio_theta1_genotype_n1333_beta0.3_a0.1.csv",
                    "theta1/aUnknown_n1333_beta0.3_a0.1/aUnknown_est.parametric.pairs_theta1_genotype_n1333_beta0.3_a0.1.csv")))

df62 = noquote(mm(c("theta0/aUnknown_n1333_beta0.3_a0.1/aUnknown_est.semiparametric_theta0_genotype_n1333_beta0.3_a0.1.csv",
                    "theta0/aUnknown_n1333_beta0.3_a0.1/aUnknown_est.semiparametric.pairs_theta0_genotype_n1333_beta0.3_a0.1.csv",
                    "theta0/aUnknown_n1333_beta0.3_a0.1/aUnknown_est.parametric.trio_theta0_genotype_n1333_beta0.3_a0.1.csv",
                    "theta0/aUnknown_n1333_beta0.3_a0.1/aUnknown_est.parametric.pairs_theta0_genotype_n1333_beta0.3_a0.1.csv")))

df63 = noquote(mm(c("proband_genotype_prospective/aUnknown_n1333_beta0.3_a0.1/aUnknown_est.semiparametric_proband_genotype_n1333_beta0.3_a0.1.csv",
                    "proband_genotype_prospective/aUnknown_n1333_beta0.3_a0.1/aUnknown_est.semiparametric.pairs_proband_genotype_n1333_beta0.3_a0.1.csv",
                    "proband_genotype_prospective/aUnknown_n1333_beta0.3_a0.1/aUnknown_est.parametric.trio_proband_genotype_n1333_beta0.3_a0.1.csv",
                    "proband_genotype_prospective/aUnknown_n1333_beta0.3_a0.1/aUnknown_est.parametric.pairs_proband_genotype_n1333_beta0.3_a0.1.csv")))

dt = as.matrix(noquote(cbind(rep(c("Carriers","Non-carriers","Both"),6),
                             rbind(df11,df12,df13,df21,df22,df23,
                             df31,df32,df33,df41,df42,df43,
                             df51,df52,df53,df61,df62,df63))))
colnames(dt) = c("Sample Type","Trio","Pairs","Trio","Pairs","Trio","Pairs","Trio","Pairs")
rownames(dt) = 1:18
dt


options(kableExtra.latex.load_packages = FALSE) 
library(kableExtra)
kbl(dt,"latex")
kbl(dt, booktabs = T, "latex") %>% 
  kable_styling(latex_options = c("striped", "hover")) %>% 
  add_header_above(c(" "=1,"semi parametric" = 2, "parametric" = 2, "semi parametric" = 2, "parametric" = 2)) %>% 
  add_header_above(c(" "=1,"beta" = 4, "a" = 4)) %>% 
  kable_styling(position = "center") %>% 
  column_spec(c(2,6), border_left = T) %>%
  pack_rows("n = 50, beta = -1, a = 0.1", 1,3, latex_gap_space = "1em",label_row_css = "background-color: #999; color: #fff;") %>%
  pack_rows("n = 1000, beta = -1, a = 0.1", 4,6, latex_gap_space = "1em",label_row_css = "background-color: #999; color: #fff;") %>%
  pack_rows("n = 50, beta = 1, a = 0.2", 7,9, latex_gap_space = "1em",label_row_css = "background-color: #999; color: #fff;") %>%
  pack_rows("n = 1000, beta = 1, a = 0.2", 10,12, latex_gap_space = "1em",label_row_css = "background-color: #999; color: #fff;") %>%
  pack_rows("n = 50, beta = 0.3, a = 0.1", 13,15, latex_gap_space = "1em",label_row_css = "background-color: #999; color: #fff;") %>%
  pack_rows("n = 1000, beta = 0.3, a = 0.1", 16,18, latex_gap_space = "1em",label_row_css = "background-color: #999; color: #fff;")%>%
  column_spec(1, width = "4em") %>%
  row_spec(c(3,6,9,12,15), hline_after = T) %>%
  row_spec(c(1,4,7,10,13), hline_after = F)

kbl(dt, booktabs = T, linesep = "") %>% 
  kable_styling(latex_options = c("striped"), font_size = 14) %>% 
  add_header_above(c(" "=1,"semi parametric" = 2, "parametric" = 2, "semi parametric" = 2, "parametric" = 2)) %>% 
  add_header_above(c(" "=1,"beta" = 4, "a" = 4)) %>% 
  kable_styling(position = "center") %>% 
  column_spec(c(2,6), border_left = T) %>%
  pack_rows("n = 50, beta = -1, a = 0.1", 1,3, latex_gap_space = "2em",label_row_css = "background-color: #999; color: #fff;") %>%
  pack_rows("n = 1000, beta = -1, a = 0.1", 4,6, latex_gap_space = "2em") %>%
  pack_rows("n = 50, beta = 1, a = 0.2", 7,9, latex_gap_space = "2em") %>%
  pack_rows("n = 1000, beta = 1, a = 0.2", 10,12, latex_gap_space = "2em") %>%
  pack_rows("n = 50, beta = 0.3, a = 0.1", 13,15, latex_gap_space = "2em") %>%
  pack_rows("n = 1000, beta = 0.3, a = 0.1", 16,18, latex_gap_space = "2em") %>%
  column_spec(1, width = "4em")  %>%
  row_spec(c(3,6,9,12,15), hline_after = T) %>%
  row_spec(c(1,4,7,10,13), hline_after = F)
