# This template shows how the lasso sign recovery probabilites can be simulated for a 
#screening design with main effects and 2 factor interactions

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

# First we add a few functions


#' Function to generate active main effects and active 2FI under weak heredity
#'
#' @param k Integer: number of main effects in design
#' @param me Integer: number of active main effects to sample
#' @param twoFI Integer: Number of active 2FI under weak heredity
#'
#' @return
#' @export
#'
#' @examples
active_int <- function(k, me, twoFI)
{
  arr <- matrix(2, k, 1)   #2 effect not active
  mes=sample(1:k,me)
  arr[mes]<-1
  
  inactive_mes=which(arr %in% c(2))
  
  #find which interactions are active depending on activity of main effects
  #4 both me's inactive, 2 one of the me's active, 1 both me's active
  ttt <- arr%*%t(arr)
  res <- c((diag(ttt)))
  
  
  for (i in 1:(k-1))
  {
    for (j in (i+1):k)
    {
      res <- c(res, ttt[i,j])
    }
  }
  
  
  full=c(res)
  
  activeYesorNo <- numeric(length(full)) #finds location of active me
  for (i in 1:k)
  {
    if (res[i]==1) {activeYesorNo[i] <- 1}
  }
  
  
  active=which(activeYesorNo %in% c(1))
  inactive= which(activeYesorNo %in% c(0))
  
  for (i in (k+1):length(full))  #finds active interactions--strong heredity
  {
    if (full[i]==1) {activeYesorNo[i] <-1}
  }
  
  active_strong=which(activeYesorNo %in% c(1))
  remove<-c(mes)
  strong<-setdiff(active_strong, remove)
  
  for (i in (k+1):length(full))  #finds active interactions--weak heredity
  {
    if (full[i]==2) {activeYesorNo[i] <-1}
  }
  
  
  
  bob=which(activeYesorNo %in% c(1))
  remove<-c(mes,strong)
  weak<-setdiff(bob, remove)
  
  weak_strong = c(weak, strong)
  sample_twoFI = sample(1:length(weak_strong), twoFI)
  activeweak_strong= weak_strong[sample_twoFI]
  
  remove =c(inactive_mes, activeweak_strong)
  inactive_2FI = setdiff(inactive, remove)
  
  
  
  
  
  #ppp=c(mes,activestrong,activeweak,activenone)
  ppp=c(mes,activeweak_strong)
  effects=numeric(length(activeYesorNo))
  effects[ppp]<-1
  output<-list(effects,mes,activeweak_strong,inactive_mes, inactive_2FI)
  
  return(output)
}


# Function to select Beta, this can be done in any way you want
# This follows the strategy of Mee et. al. 2017
# X is model matrix with MEs and 2FIs (no intercept)
# k is number of total MEs
# Output is a vector of effect sizes

select_beta = function(X, k){
  
  effect_size_me = c(2,2.5,3,3.5,-2,-2.5,-3, -3.5)
  #effect_size_me = c(3.5,-3.5)
  effect_size_2FI = c(0.5,1,1.5,2,-0.5,-1,-1.5,-2)
  
  beta_me=sample(effect_size_me,k, replace = TRUE)
  beta_2FI=sample(effect_size_2FI,dim(X)[2]-k, replace = TRUE)
  beta_all=c(beta_me, beta_2FI)
  return(beta_all)
}





# Function to generate a random response 

true_y<-function(X,k, active_factorsYesNo,active_me, active_2FI, inactive_me, inactive_2FI, beta_all, error=1){
  effects=1:dim(X)[2]
  n = dim(X)[1]
  
 
  active = c(active_me, active_2FI)
  
  true_beta = active_factorsYesNo*beta_all
  #true=effects%in%active
  inactive=c(inactive_me, inactive_2FI)
  #browser
  Y=X%*%true_beta+rnorm(n,mean=0, sd=error)
  
  return(list(Y=Y, active=active, inactive=inactive, active_me = active_me, 
              active_2FI=active_2FI, inactive_me = inactive_me, inactive_2FI=inactive_2FI))
  
}











### Fit lasso for one value of lambda
### X=model matrix with 2FIs and no intercept column
### Y=true response vector
### lambda=value of lambda for the fit

lasso.fit<-function(X, Y){
  require(broom)
  require(glmnet)
  n = dim(X)[1]
  Yc=Y-mean(Y)
  Xc<-apply(X, 2, function(X) X - mean(X))
  vlength<-sqrt((1/n)*(diag(t(Xc)%*%Xc)))
  Xs<-t(t(Xc)*(1/(vlength)))
  log_lam = seq(from = -4.5, to= 1.5, by = 0.1)
  lambda = exp(log_lam)
  fit<-tidy(glmnet(Xs, Yc, standardize = FALSE, intercept = FALSE, 
                   lambda=lambda))
  
  return(fit)
}







### simulation loop

## set working directory
setwd("path to folder with designs")

##
# Getting the names of desings
paths<-dir(pattern=".csv")
names(paths)<-basename(paths)
names<-names(paths)
names <- names
print(names)

#number active ME
me_set=c(2, 4)
#number active 2FI
twoFI_set = seq(1,7,1)

# Creating simulation grid

grid<-expand.grid(me_set,twoFI_set)

# Creating lists to house sim results

Des.perfect<-list()
Des.perfect_me<-list()
Des.perfect_2FI<-list()

# Set number of runs
n = 20
# Number of total main effect factors
k = 7



for (j in 1:length(names)){
  # read in the main effect designs
  X_me=as.data.frame(read_csv(names[j], col_names=FALSE))
  #compute the number of total effects in the model, main effects and 2FIs
  num_tot_eff = k + choose(k,2)
  
  Y_dummy= rnorm(nrow(X_me))
  form <- Y_dummy~.^2
  
  #Model matrix with all 2FI excluding intercept
  X = model.matrix(form, X_me)[,2:(num_tot_eff+1)]
  colnames(X)=sprintf("X%d", seq(1:num_tot_eff))
  
  k=ncol(X_me)
  n=nrow(X_me)
  print(names[j])
  
  # Generating a matrix to hold perfect sign recovery results for each lambda considered. Each row will be a different data generation
  # 61 columns for 61 lambdas considered in the lasso.fit function
  
  mat<-matrix(ncol=61, nrow=0)
  
  # Seperate results into overall sign recovery, and sign recovery for MEs and 2FI, respectively
  scenario.perfect<-data.frame(mat)
  scenario.perfect_me<-data.frame(mat)
  scenario.perfect_2FI<-data.frame(mat)
  
  for (h in 1:nrow(grid)){
    
    
    ### This is in parallel Double check how many clusters you want to use
    cl=makeCluster(4)
    registerDoParallel(cl)
    # 500 samples of betas and active sumbodels
    sim<-foreach(i=1:500, .combine = rbind) %dopar% {
      
      perfect_df =data.frame(mat)
      
      
      perfect_df_me =data.frame(mat)
      
      
      perfect_df_2FI =data.frame(mat)
      
      
      active_fun_output=active_int(k, me=grid$Var1[h], twoFI= grid$Var2)
      beta_all = select_beta(as.matrix(X), k)
      # Loop over 15 generations of Y
      for( l in 1:15){
        resp<-true_y(as.matrix(X),k, active_fun_output[[1]], active_fun_output[[2]], 
                     active_fun_output[[3]], active_fun_output[[4]], 
                     active_fun_output[[5]], beta_all)
        fit<-lasso.fit(X, Y=resp$Y)
        power<-c()
        type1<-c()
        power_me<-c()
        type1_me<-c()
        power_2FI<-c()
        type1_2FI<-c()
        
        
        # Evaluate Lasso sign recovery results
        for (d in 1:max(fit$step)) {
          sel<-as.numeric(gsub("X", "", subset(fit, step==d)$term))
          power[d]<-length(which(sel%in%resp$active)) #correct effects
          type1[d]<-length(which(sel%in%resp$inactive)) #incorrect effects
          power_me[d]<-length(which(sel%in%resp$active_me)) #correct effects
          type1_me[d]<-length(which(sel%in%resp$inactive_me)) #incorrect effects
          power_2FI[d]<-length(which(sel%in%resp$active_2FI)) #correct effects
          type1_2FI[d]<-length(which(sel%in%resp$inactive_2FI)) #incorrect effects
        }
        
        
        
        perfect<-length(which(power==length(resp$active))) #models with the correct support
        perfect.type1<-type1[which(power==length(resp$active))] #type1 errors for models with correct support
        perfect_vec = ifelse(power==length(resp$active) & type1==0, 1,0)
        
        errors<-length(which(perfect.type1>0)) #count of models with correct support and type1 errors
        perfect=perfect-errors #adjust for models with type1 errors
        
        
        
        perfect_me<-length(which(power_me==length(resp$active_me))) #models with the correct main effects
        perfect.type1_me<-type1_me[which(power_me==length(resp$active_me))] #type1 errors for models with correct main effects
        perfect_vec_me = ifelse(power_me==length(resp$active_me) & type1_me==0, 1,0)
        
        errors_me<-length(which(perfect.type1_me>0)) #count of models with correct support and type1 errors
        perfect_me=perfect_me-errors_me #adjust for models with type1 errors
        
        
        perfect_2FI<-length(which(power_2FI==length(resp$active_2FI))) #models with the correct 2FI
        perfect.type1_2FI<-type1_2FI[which(power_2FI==length(resp$active_2FI))] #type1 errors for models with correct main effects
        perfect_vec_2FI = ifelse(power_2FI==length(resp$active_2FI) & type1_2FI==0, 1,0)
        
        errors_2FI<-length(which(perfect.type1_2FI>0)) #count of models with correct support and type1 errors
        perfect_2FI=perfect_2FI-errors_2FI #adjust for models with type1 errors
        
        
        
        perfect_df[l,]=perfect_vec
       
        
        perfect_df_me[l,]=perfect_vec_me
       
        
        perfect_df_2FI[l,]=perfect_vec_2FI
        
        
      }
      return(list(perfect= colMeans(perfect_df), 
                  perfect_me = colMeans(perfect_df_me), 
                  perfect_2FI = colMeans(perfect_df_2FI)))
      
    }
    
    stopCluster(cl)
    perfect_sum<-sim[1,1][[1]]
    perfect_me_sum<-sim[1,2][[1]]
    perfect_2FI_sum<-sim[1,3][[1]]
 
    for (g in 2:500) {
      
      perfect_sum =perfect_sum+ sim[g,1][[1]]
      perfect_me_sum =perfect_me_sum+ sim[g,2][[1]]
      perfect_2FI_sum =perfect_2FI_sum+ sim[g,3][[1]]
    }
    
    scenario.perfect[h,]<-perfect_sum/500
    scenario.perfect_me[h,]<-perfect_me_sum/500
    scenario.perfect_2FI[h,]<-perfect_2FI_sum/500
   
    
  }#close h loop
  
  Des.perfect[[j]]<-scenario.perfect
  Des.perfect_me[[j]]<-scenario.perfect_me
  Des.perfect_2FI[[j]]<-scenario.perfect_2FI
  
  
  print(j)
}#close j loop




# Match log lam sequence in lasso.fit
log_lam = seq(from = -4.5, to= 1.5, by = 0.1)
lambda = exp(log_lam)

designs = names


df_results = merge( designs,grid, by=NULL)

df_results = merge(lambda, df_results, by=NULL)

colnames(df_results)= c("lambda", "design", "ME", "two_FI")

colnames(grid)= c("ME", "two_FI")

perfect_results = cbind(grid, Des.perfect[[1]])
perfect_results$design = designs[1]

perfect_results_me = cbind(grid, Des.perfect_me[[1]])
perfect_results_me$design = designs[1]

perfect_results_2FI = cbind(grid, Des.perfect_2FI[[1]])
perfect_results_2FI$design = designs[1]

for(j in c(2:6)){
  temp = cbind(grid, Des.perfect[[j]])
  temp$design=designs[j]
  perfect_results = rbind(perfect_results, temp)
  
  temp = cbind(grid, Des.perfect_me[[j]])
  temp$design=designs[j]
  perfect_results_me= rbind(perfect_results_me, temp)
  
  temp = cbind(grid, Des.perfect_2FI[[j]])
  temp$design=designs[j]
  perfect_results_2FI = rbind(perfect_results_2FI, temp)
  
  
}

perfect_results_long = perfect_results %>% pivot_longer(!c("ME","two_FI","design"), names_to="temp", values_to="perfect")

perfect_results_long_me = perfect_results_me %>% pivot_longer(!c("ME","two_FI","design"), names_to="temp", values_to="perfect_me")


perfect_results_long_2FI = perfect_results_2FI %>% pivot_longer(!c("ME","two_FI","design"), names_to="temp", values_to="perfect_2FI")



df_results = cbind(df_results, perfect_results_long$perfect, perfect_results_long_me$perfect_me,perfect_results_long_2FI$perfect_2FI)
colnames(df_results)[5]="perfect_prob"
colnames(df_results)[6]="perfect_me"
colnames(df_results)[7]="perfect_2FI"





q_perfect= ggplot(data = df_results, aes(x=log(as.numeric(lambda)), y= as.numeric(perfect_prob), color= design, linetype = design))+
  geom_line()+
  facet_grid(vars(ME), vars(two_FI))+
  xlab(latex2exp::TeX("$log(\\lambda)"))+
  ylab(latex2exp::TeX(" $\\phi_{\\lambda}"))+
  #scale_color_manual(values =c("blue","red", "green"))+
  theme_bw()+
  #scale_x_continuous(limits = c(-2,2))+
  # scale_linetype_manual(values =c("solid", "dashed", "solid"))+
  theme( legend.title=element_blank(), text=element_text(size=16))






q_me= ggplot(data = df_results, aes(x=log(as.numeric(lambda)), y= as.numeric(perfect_me), color= design))+
  geom_line()+
  facet_grid(vars(ME), vars(two_FI))+
  xlab(latex2exp::TeX("$log(\\lambda)"))+
  ylab(latex2exp::TeX("ME $\\phi_{\\lambda}"))+
  #scale_color_manual(values =c("blue","red", "green"))+
  theme_bw()+
  #scale_x_continuous(limits = c(-2,2))+
  # scale_linetype_manual(values =c("solid", "dashed", "solid"))+
  theme( legend.title=element_blank(), text=element_text(size=16))










q_2FI= ggplot(data = df_results, aes(x=log(as.numeric(lambda)), y= as.numeric(perfect_2FI), color= design, linetype = design))+
  geom_line()+
  facet_grid(vars(ME), vars(two_FI))+
  xlab(latex2exp::TeX("$log(\\lambda)"))+
  ylab(latex2exp::TeX("2FI $\\phi_{\\lambda}"))+
  #scale_color_manual(values =c("blue","red", "green"))+
  theme_bw()+
  #scale_x_continuous(limits = c(-2,2))+
  # scale_linetype_manual(values =c("solid", "dashed", "solid"))+
  theme( legend.title=element_blank(), text=element_text(size=16))


