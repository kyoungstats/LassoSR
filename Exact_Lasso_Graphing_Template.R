### This is a template R Script for the evaluation of 4 generic designs main effect designs
### using the exact lasso sign recovery probabilities
### Output is a graphical comparison ggplot object






library(abind)
library(readr)
library(doParallel)
library(foreach)
library(dplyr)
library(tidyr)
library(glmnet)
library(ggplot2)
library(purrr)
library(gridExtra)


#Add your path to the Function Library
source('PATH TO LIBRARY /Design_Comparison_Paper_Function_Library.R')



#Add your path to the folder where designs are kept
design_path = "PATH TO DESIGN FOLDER/n14_k24_designs"

# Set effect size
beta = 3

# Read in designs as matrices without an intercept
# Can read in more or less but will need to adjust code later in script to account for it

d1<- as.matrix(read_in_design(paste0(design_path, "/d1.txt")))
d2 <- as.matrix(read_in_design(paste0(design_path, "/d2.txt")))
d3 <- as.matrix(read_in_design(paste0(design_path, "/d3.txt")))
d4 <- as.matrix(read_in_design(paste0(design_path, "/d4.txt")))

# Set number of runs
n=14
# Set number of main effects
k=24
# Set number of active main effects
k_star=5

# Generate sign vector in this case it is all positive 
# if you want to average over all signs, change the sign_vects argument to NULL in the Joint_prob_general function below 
sign_all_pos = matrix(rep(1,k_star), nrow=1, ncol=k_star, byrow=TRUE)

# Sample the submodels to average over
# If you want to average over all possible k choose k_star submodels, change the submodels argument to NULL in the Joint_prob_general function below

submodel_samples = stack_NBIBD(p=k, k=k_star, s_1=32, num_stacks = 14)

A_val = submodel_samples$A_final

# Set grid of tuning parameter values to graph over
log_lam = seq(-4.5, 2, by = 0.05)

# Create list of designs
designs = abind( d1, d2, d3, d4, along =3)

# start the evaluation this is done in parallel

start_time <- Sys.time()
cl <- makeCluster(4)
registerDoParallel(cl)
# 
design_names = c("d1", "d2", "d3", "d4")
# 
df_results = foreach(start=1:dim(designs)[3], .combine = rbind, .packages =c('MASS',"mvtnorm","dplyr"),.export = ls(globalenv())) %dopar%{
  if (start ==1){
    design_name = "d1"
  }
  if (start ==2){
    design_name = "d2"
  }
  if (start ==3){
    design_name = "d3"
  }
  else{ design_name = "d4"}
  
  F_0 = designs[,,start]
  # Center and scale designs
  V_0_half= diag(lapply(as.data.frame(F_0),fn))
  F_0_cs = as.data.frame(lapply(as.data.frame(F_0), cent_scale))
  # Create empty df to fill for graphing
  joint_results= data.frame(design=character(0), lambda=numeric(0), P_joint=numeric(0))
  joint_results$design = factor(joint_results$design, levels=c("d1", "d2", "d3", "d4"))
  
  # If the centered and scaled design has a na in it, that means there is a constant column in the design
  # assign it a probability of 0
  if(any(is.na(F_0_cs))){
    prob=0}
  else{
    for (i in c(1:length(log_lam))){
      log_lambda = log_lam[i]
      
      
      prob_joint <- measure_list(Joint_prob_general(F_cs =F_0_cs,V_0_half,
                                                    B_mag=rep(beta,k_star),
                                                    k=k_star,
                                                    sign_vects=sign_all_pos,#change to NULL to average over all sign vects
                                                    log_lambda = log_lambda,
                                                    submodels = A_val,#change to NULL to average over all possible submodels
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

# 
ggjoint = ggplot(data = df_results, aes(x=log(as.numeric(lambda)), y= as.numeric(P_joint), color= design, linetype = design))+
  geom_line(size=1.2)+
  xlab(latex2exp::TeX("$log(\\lambda)"))+
  ylab(latex2exp::TeX("$\\phi_{\\lambda}"))+
  scale_color_manual(values =c("blue","darkgreen", "purple", "red"))+
  theme_bw()+
  scale_x_continuous(limits = c(-2,2))+
  #scale_linetype_manual(values =c("solid", "dashed"))+
  theme( legend.title=element_blank(), text=element_text(size=14))






