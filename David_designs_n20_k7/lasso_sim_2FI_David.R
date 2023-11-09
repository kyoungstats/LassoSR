# Lasso Exact sign recovery simulation with 2FI from David's paper

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
source('/Users/hkyoung/Desktop/SIC_Paper_Code/SIC_Paper_Function_Library.R')




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




select_beta = function(X, k){
  
  effect_size_me = c(2,2.5,3,3.5,-2,-2.5,-3, -3.5)
  #effect_size_me = c(3.5,-3.5)
  effect_size_2FI = c(0.5,1,1.5,2,-0.5,-1,-1.5,-2)
  
  beta_me=sample(effect_size_me,k, replace = TRUE)
  beta_2FI=sample(effect_size_2FI,dim(X)[2]-k, replace = TRUE)
  beta_all=c(beta_me, beta_2FI)
  return(beta_all)
}







true_y<-function(X,k, active_factorsYesNo,active_me, active_2FI, inactive_me, inactive_2FI, beta_all, error=1){
  effects=1:dim(X)[2]
  n = dim(X)[1]
  
  #beta_vec = rep(0,dim(X)[2])
  #positive or negative effects from same size set
  effect_size_me = c(2,2.5,3,3.5,-2,-2.5,-3, -3.5)
  effect_size_2FI = c(0.5,1,1.5,2,-0.5,-1,-1.5,-2)
  
  # beta_me=sample(effect_size_me,k, replace = TRUE)
  # beta_2FI=sample(effect_size_2FI,dim(X)[2]-k, replace = TRUE)
  # beta_all=c(beta_me, beta_2FI)
  #browser()
  #active_fun_output = active_int(k, me, twoFI)
  #active_factorsYesNo = active_fun_output[[1]]
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
### X=model matrix with 2FIs
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

## The line of code below needs to change
setwd("./David_designs_n20_k7")

##

paths<-dir(pattern=".csv")
names(paths)<-basename(paths)
names<-names(paths)
names <- names
print(names)
#p=c(0, 0.5) #all positive set 0, mixed set 0.5

me_set=c(2, 4)#number active ME

twoFI_set = seq(1,7,1)



grid<-expand.grid(me_set,twoFI_set)

Des.perfect<-list()
Des.perfect_me<-list()
Des.perfect_2FI<-list()
Des.oneFP<-list()
Des.twoFP<-list()

n = 20
k = 7



for (j in 1:length(names)){
  
  X_me=as.data.frame(read_csv(names[j], col_names=FALSE))
  
  Y_dummy= rnorm(nrow(X_me))
  form <- Y_dummy~.^2
  #matrix with all 2FI
  X = model.matrix(form, X_me)[,2:29]
  colnames(X)=sprintf("X%d", seq(1:28))
  k=ncol(X_me)
  n=nrow(X_me)
  print(names[j])
  
  # Generating a matrix to hold perfect, one FP and two FP results for each lambda considered. Each row will be a different data generation
  mat<-matrix(ncol=61, nrow=0)
  
  scenario.perfect<-data.frame(mat)
  scenario.perfect_me<-data.frame(mat)
  scenario.perfect_2FI<-data.frame(mat)
  scenario.oneFP<-c()
  scenario.twoFP<-c()
  for (h in 1:nrow(grid)){
    # if(grid$Var3[h]=="med"){
    #   A_val = A_val_med
    # } else if(grid$Var3[h]=="hi"){
    #   A_val = A_val_hi
    # }else{
    #   A_val = A_val_low
    # }
    
    ### Double check how many clusters you want to use
    cl=makeCluster(2)
    registerDoParallel(cl)
    
    sim<-foreach(i=1:500, .combine = rbind) %dopar% {
      perfect_df =data.frame(mat)
      oneFP_df = data.frame(mat)
      twoFP_df = data.frame(mat)
      
      perfect_df_me =data.frame(mat)
      oneFP_df_me = data.frame(mat)
      twoFP_df_me = data.frame(mat)
      
      perfect_df_2FI =data.frame(mat)
      oneFP_df_2FI = data.frame(mat)
      twoFP_df_2FI = data.frame(mat)
      
      #active_fun_output=active_int(k, me=grid$Var1[h], twoFI= grid$Var2[h])
      active_fun_output=active_int(k, me=grid$Var1[h], twoFI= grid$Var2)
      beta_all = select_beta(as.matrix(X), k)
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
        oneFP_vec =ifelse((power==length(resp$active) & type1<=1), 1,0)
        twoFP_vec =ifelse((power==length(resp$active) & type1<=2), 1,0)
        errors<-length(which(perfect.type1>0)) #count of models with correct support and type1 errors
        perfect=perfect-errors #adjust for models with type1 errors
        oneFP<-length(which(perfect.type1==1)) #models with correct support and one type1 error
        twoFP<-length(which(perfect.type1==2)) #models with correct support and two type1 error
        
        
        perfect_me<-length(which(power_me==length(resp$active_me))) #models with the correct main effects
        perfect.type1_me<-type1_me[which(power_me==length(resp$active_me))] #type1 errors for models with correct main effects
        perfect_vec_me = ifelse(power_me==length(resp$active_me) & type1_me==0, 1,0)
        oneFP_vec_me =ifelse((power_me==length(resp$active_me) & type1_me<=1), 1,0)
        twoFP_vec_me =ifelse((power_me==length(resp$active_me) & type1_me<=2), 1,0)
        errors_me<-length(which(perfect.type1_me>0)) #count of models with correct support and type1 errors
        perfect_me=perfect_me-errors_me #adjust for models with type1 errors
        oneFP_me<-length(which(perfect.type1_me==1)) #models with correct support and one type1 error
        twoFP_me<-length(which(perfect.type1_me==2)) #models with correct support and two type1 error
        
        perfect_2FI<-length(which(power_2FI==length(resp$active_2FI))) #models with the correct 2FI
        perfect.type1_2FI<-type1_2FI[which(power_2FI==length(resp$active_2FI))] #type1 errors for models with correct main effects
        perfect_vec_2FI = ifelse(power_2FI==length(resp$active_2FI) & type1_2FI==0, 1,0)
        oneFP_vec_2FI =ifelse((power_2FI==length(resp$active_2FI) & type1_2FI<=1), 1,0)
        twoFP_vec_2FI =ifelse((power_2FI==length(resp$active_2FI) & type1_2FI<=2), 1,0)
        errors_2FI<-length(which(perfect.type1_2FI>0)) #count of models with correct support and type1 errors
        perfect_2FI=perfect_2FI-errors_2FI #adjust for models with type1 errors
        oneFP_2FI<-length(which(perfect.type1_2FI==1)) #models with correct support and one type1 error
        twoFP_2FI<-length(which(perfect.type1_2FI==2)) #models with correct support and two type1 error
        
        
        perfect_df[l,]=perfect_vec
        oneFP_df[l,]=oneFP_vec
        twoFP_df[l,]=twoFP_vec
        
        perfect_df_me[l,]=perfect_vec_me
        oneFP_df_me[l,]=oneFP_vec_me
        twoFP_df_me[l,]=twoFP_vec_me
        
        perfect_df_2FI[l,]=perfect_vec_2FI
        oneFP_df_2FI[l,]=oneFP_vec_2FI
        twoFP_df_2FI[l,]=twoFP_vec_2FI
        
      }
      return(list(perfect= colMeans(perfect_df), oneFP= colMeans(oneFP_df), twoFP=colMeans(twoFP_df),
                  perfect_me = colMeans(perfect_df_me), oneFP_me = colMeans(oneFP_df_me), twoFP_me=colMeans(twoFP_df_me),
                  perfect_2FI = colMeans(perfect_df_2FI), oneFP_2FI = colMeans(oneFP_df_2FI), twoFP_2FI=colMeans(twoFP_df_2FI)))
      #return(list(perfect= (sum(colMeans(perfect_df)/61)), oneFP= (sum(colMeans(oneFP_df))/61), twoFP=(sum(colMeans(twoFP_df))/61)))
    }
    
    stopCluster(cl)
    perfect_sum<-sim[1,1][[1]]
    perfect_me_sum<-sim[1,4][[1]]
    perfect_2FI_sum<-sim[1,7][[1]]
    # oneFP<-c()
    # twoFP<-c()
    for (g in 2:500) {
      # perfect[g]<-sim[g,]$perfect
      # oneFP[g]<-sim[g,]$oneFP
      # twoFP[g]<-sim[g,]$twoFP
      perfect_sum =perfect_sum+ sim[g,1][[1]]
      perfect_me_sum =perfect_me_sum+ sim[g,4][[1]]
      perfect_2FI_sum =perfect_2FI_sum+ sim[g,7][[1]]
    }
    
    scenario.perfect[h,]<-perfect_sum/500
    scenario.perfect_me[h,]<-perfect_me_sum/500
    scenario.perfect_2FI[h,]<-perfect_2FI_sum/500
    # scenario.oneFP[h]<-mean(oneFP)
    # scenario.twoFP[h]<-mean(twoFP)
    
  }#close h loop
  
  Des.perfect[[j]]<-scenario.perfect
  Des.perfect_me[[j]]<-scenario.perfect_me
  Des.perfect_2FI[[j]]<-scenario.perfect_2FI
  
  # Des.oneFP[[j]]<-scenario.oneFP
  # Des.twoFP[[j]]<-scenario.twoFP
  print(j)
}#close j loop





log_lam = seq(from = -4.5, to= 1.5, by = 0.1)
lambda = exp(log_lam)

designs = c("BayesD", "MEPI", "OA18", "OA1a", "OA1b", "PEC")


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

write.csv(df_results, "david_paper_lasso_sim_sampled_ME.csv", col.names=TRUE, row.names=FALSE)
#write.csv(df_results, "david_paper_lasso_sim_large_ME.csv", col.names=TRUE, row.names=FALSE)

df_results

q= ggplot(data = df_results, aes(x=log(as.numeric(lambda)), y= as.numeric(perfect_prob), color= design, linetype = design))+
  geom_line()+
  facet_grid(vars(ME), vars(two_FI))+
  xlab(latex2exp::TeX("$log(\\lambda)"))+
  ylab(latex2exp::TeX(" $\\phi_{\\lambda}"))+
  #scale_color_manual(values =c("blue","red", "green"))+
  theme_bw()+
  #scale_x_continuous(limits = c(-2,2))+
 # scale_linetype_manual(values =c("solid", "dashed", "solid"))+
  theme( legend.title=element_blank(), text=element_text(size=16))



q_trimmed= ggplot(data = df_results[(df_results$two_FI%in% c(1,4,7)),], aes(x=log(as.numeric(lambda)), y= as.numeric(perfect_prob), color= design, linetype=design))+
  geom_line(size=1.2)+
  facet_grid(vars(ME), vars(two_FI))+
  xlab(latex2exp::TeX("$log(\\lambda)"))+
  ylab(latex2exp::TeX("$\\phi_{\\lambda}"))+
  #scale_color_manual(values =c("blue","red", "green"))+
  theme_bw()+
  #scale_x_continuous(limits = c(-2,2))+
  # scale_linetype_manual(values =c("solid", "dashed", "solid"))+
  theme( legend.title=element_blank(), text=element_text(size=22))

ggsave("david_lasso_sim_overall_trimmed.pdf",q_trimmed)
  

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
  
ggsave("david_lasso_sim_me.pdf",q_me)


q_me_trimmed= ggplot(data = df_results[(df_results$two_FI%in% c(1,4,7)),], aes(x=log(as.numeric(lambda)), y= as.numeric(perfect_me), color= design, linetype=design))+
  geom_line(size=1.2)+
  facet_grid(vars(ME), vars(two_FI))+
  xlab(latex2exp::TeX("$log(\\lambda)"))+
  ylab(latex2exp::TeX("ME $\\phi_{\\lambda}"))+
  #scale_color_manual(values =c("blue","red", "green"))+
  theme_bw()+
  #scale_x_continuous(limits = c(-2,2))+
  # scale_linetype_manual(values =c("solid", "dashed", "solid"))+
  theme( legend.title=element_blank(), text=element_text(size=22))

ggsave("david_lasso_sim_me_trimmed.pdf",q_me_trimmed)






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


q_2FI_trimmed= ggplot(data = df_results[(df_results$two_FI%in% c(1,4,7)),], aes(x=log(as.numeric(lambda)), y= as.numeric(perfect_2FI), color= design, linetype=design))+
  geom_line(size=1.2)+
  facet_grid(vars(ME), vars(two_FI))+
  xlab(latex2exp::TeX("$log(\\lambda)"))+
  ylab(latex2exp::TeX("2FI $\\phi_{\\lambda}"))+
  #scale_color_manual(values =c("blue","red", "green"))+
  theme_bw()+
  #scale_x_continuous(limits = c(-2,2))+
  # scale_linetype_manual(values =c("solid", "dashed", "solid"))+
  theme( legend.title=element_blank(), text=element_text(size=22))

ggsave("david_lasso_sim_2FI_trimmed.pdf",q_2FI_trimmed)

q_2_1 =ggplot(data = df_results[df_results$ME==2 & df_results$two_FI==1,], aes(x=log(as.numeric(lambda)), y= as.numeric(perfect_prob), color= design, linetype = design))+
  geom_line()+
xlab(latex2exp::TeX("$log(\\lambda)"))+
  ylab(latex2exp::TeX("$\\phi_{\\lambda}"))+
  #scale_color_manual(values =c("blue","red", "green"))+
  theme_bw()+
  #scale_x_continuous(limits = c(-2,2))+
  # scale_linetype_manual(values =c("solid", "dashed", "solid"))+
  theme( legend.title=element_blank(), text=element_text(size=16), legend.position = c(0.1,.75))
  
  
q_2_1_me =ggplot(data = df_results[df_results$ME==2 & df_results$two_FI==1,], aes(x=log(as.numeric(lambda)), y= as.numeric(perfect_me), color= design, linetype = design))+
  geom_line()+
  xlab(latex2exp::TeX("$log(\\lambda)"))+
  ylab(latex2exp::TeX("$\\phi_{\\lambda}"))+
  #scale_color_manual(values =c("blue","red", "green"))+
  theme_bw()+
  #scale_x_continuous(limits = c(-2,2))+
  # scale_linetype_manual(values =c("solid", "dashed", "solid"))+
  theme( legend.title=element_blank(), text=element_text(size=16), legend.position = c(0.1,.75))
  
  df_results_early_example = df_results[df_results$ME==2 & df_results$two_FI==1,]
    
  df_results_early_example = df_results_early_example[df_results_early_example$design=="BayesD"| df_results_early_example$design=="PEC",]
    
  df_results_early_example = df_results_early_example %>% mutate(design = recode(design, "BayesD"="Design A", "PEC"="Design B"))  
    
    
  
  q_early =ggplot(data = df_results_early_example, aes(x=log(as.numeric(lambda)), y= as.numeric(perfect_me), color= design))+
    geom_line(size=1.2)+
    xlab(latex2exp::TeX("$log(\\lambda)"))+
    ylab(latex2exp::TeX("$\\phi_{\\lambda}"))+
    scale_color_manual(values =c("blue","red"))+
    theme_bw()+
    #scale_x_continuous(limits = c(-2,2))+
    # scale_linetype_manual(values =c("solid", "dashed", "solid"))+
    theme( legend.title=element_blank(), text=element_text(size=14))
  
  
  ggsave("lasso_prob_early_example.pdf",q_early)
  
    
    