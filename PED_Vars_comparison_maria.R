library(abind)

library(readr)
library(doParallel)
library(foreach)
library(dplyr)
library(tidyr)
library(glmnet)
library(ggplot2)
library(purrr)

#Add your path to the SIC_Paper_Function Library
source('PATH TO LIBRARY /SIC_Paper_Function_Library.R')


library(gridExtra)

#Add your path to the folder where designs are kept
design_path = "PATH TO DESIGN FOLDER/n14_k24_designs"


beta = 3



d1<- as.matrix(read_in_design(paste0(design_path, "/d1 (1).txt")))
Var_s_80 <- as.matrix(read_in_design(paste0(design_path, "/n14_k24_vars_pos_Eff80_100.txt")))
Var_s_40 <- as.matrix(read_in_design(paste0(design_path, "/n14_k24_vars_pos_Eff40_100.txt")))
Var_s_20 <- as.matrix(read_in_design(paste0(design_path, "/n14_k24_vars_pos_Eff20_100.txt")))

n=14
k=24
k_star=5
sign_all_pos = matrix(rep(1,k_star), nrow=1, ncol=k_star, byrow=TRUE)

log_lam = seq(from = -4.5, to= 2, by = 0.1)
submodel_samples = stack_NBIBD(p=k, k=k_star, s_1=32, num_stacks = 14)

A_val = submodel_samples$A_final

log_lam = seq(-4.5, 2, by = 0.05)
# 
designs = abind( d1, Var_s_80, Var_s_40, Var_s_20, along =3)
start_time <- Sys.time()
cl <- makeCluster(4)
registerDoParallel(cl)
# 
design_names = c("PED", "Var(s+)_80", "Var(s+)_40", "Var(s+)_20")
# 
prob_results_n14_p24 = foreach(start=1:dim(designs)[3], .combine = rbind, .packages =c('MASS',"mvtnorm","dplyr"),.export = ls(globalenv())) %dopar%{
  if (start ==1){
    design_name = "PED"
  }
  if (start ==2){
    design_name = "Var(s+)_80"
  }
  if (start ==3){
    design_name = "Var(s+)_40"
  }
  else{ design_name = "Var(s+)_20"}
  F_0 = designs[,,start]
  V_0_half= diag(lapply(as.data.frame(F_0),fn))
  F_0_cs = as.data.frame(lapply(as.data.frame(F_0), cent_scale))
  joint_results= data.frame(design=character(0), lambda=numeric(0), P_joint=numeric(0))
  joint_results$design = factor(joint_results$design, levels=c("PED", "Var(s+)_80", "Var(s+)_40", "Var(s+)_20"))
  if(any(is.na(F_0_cs))){
    prob=0}
  else{
    for (i in c(1:length(log_lam))){
      log_lambda = log_lam[i]

      prob_joint <- measure_list(Joint_prob_general(F_cs =F_0_cs,V_0_half,
                                                    B_mag=rep(3,k_star),k=k_star,
                                                    sign_vects=sign_all_pos,
                                                    log_lambda = log_lambda,
                                                    submodels = A_val,
                                                    log_lambda_strategy="fixed",
                                                    sigma = 1)$results)$mean

      results_joint <- c(design=design_names[start], lambda = exp(log_lambda), prob_joint)
      


      joint_results[i,]= results_joint
    }

  }
  joint_results


}

stopCluster(cl)
end_time = Sys.time()
print(end_time - start_time)
# 
write.csv(prob_results_n14_p24, file ="./PED_Vars_comparison_n14_p24_k5_SN3.csv", row.names =FALSE, col.names = FALSE)


df_n14_p24_SN3 <- read_csv("PED_Vars_comparison_n14_p24_k5_SN3_maria.csv")
# 
ggjoint_n14_p24 = ggplot(data = df_n14_p24_SN3, aes(x=log(as.numeric(lambda)), y= as.numeric(P_joint), color= design, linetype = design))+
  geom_line(size=1.2)+
  xlab(latex2exp::TeX("$log(\\lambda)"))+
  ylab(latex2exp::TeX("$\\phi_{\\lambda}"))+
  scale_color_manual(values =c("blue","darkgreen", "purple", "red"))+
  theme_bw()+
  scale_x_continuous(limits = c(-2,2))+
  #scale_linetype_manual(values =c("solid", "dashed"))+
  theme( legend.title=element_blank(), text=element_text(size=14))



  
ggsave("PED_vars_comp_sect4.pdf",ggjoint_n14_p24)

  