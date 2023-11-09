
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


candidate_designs = abind(d1, Var_s_20, Var_s_40, Var_s_80, along=3)





dantzig=function(X,Y,delta,gamma, vhalf){
  k=dim(X)[2]
  c1=matrix(1,k,1)
  c2=matrix(0,k,1)
  c=rbind(c1,c2)
  A11=t(X)%*%X
  A12=-t(X)%*%X
  A1=cbind(A11,A12)
  A21=-t(X)%*%X
  A22=t(X)%*%X
  A2=cbind(A21,A22)
  A31=2*diag(k)
  A32=-diag(k)
  A3=cbind(A31,A32)
  A=rbind(A1,A2,A3)
  b1=-t(X)%*%Y-delta*matrix(1,k,1)
  b2=t(X)%*%Y-delta*matrix(1,k,1)
  b3=matrix(0,k,1)
  b=rbind(b1,b2,b3)
  dir=matrix(">=",dim(b)[1],1)
  
  vars=lp("min", c, A, dir, b)$solution
  vars=as.vector(vars)
  u=vars[1:k]
  beta=(vars[(k+1):length(vars)]-u)/vhalf
  terms=(abs(beta)>gamma)
  answer=list(sum(terms),terms, beta, delta)
  return(answer)
  
  
}

dantzig.selection <- function(X,Y, sc, l.grid, thresh)  {
  require(lpSolve)
  Yc=Y-mean(Y)
  Xc<-apply(X, 2, function(X) X - mean(X))
  vlength<-sqrt(diag(t(Xc)%*%Xc))
  Xs<-t(t(Xc)*(1/vlength))
  n=length(Y)
  i=0
  #d0=round(max(abs((t(Xs)%*%Yc))))
  d0=16
  
  
  if (l.grid=="fine"){
    
    if (thresh=="data"){
      out=dantzig(Xs,Yc,delta=0,gamma=0, vlength)
      gamma<-0.1*max(abs(out[[3]]))
    } else {
      gamma<-1
    }
    
    output=matrix(0,nrow=100,ncol=3)
    est=matrix(0, nrow=100, ncol=ncol(X))
    for (delta in seq(0,d0,length.out=100)){
      i=i+1
      out=dantzig(Xs,Yc,delta,gamma, vlength)
      beta=out[[3]]
      beta=as.vector(beta)
      p=out[[1]]+1
      if (sum(out[[2]])==0) { Xmodel=rep(1,dim(X)[1])
      } else { Xmodel=X[,out[[2]]] }
      lin=lm(Y~Xmodel)
      SSE=deviance(lin)
      if (sc==1){
        crit=n*log(SSE/n)+(2*p)+((2*p*(p+1))/(n-p-1)) #AICc
      } else {
        crit<-n*log(SSE/n)+p*log(n) #BIC
      }
      est[i,]=beta
      output[i,1]=delta
      output[i,2]=p
      output[i,3]=crit
    }
    
  } else {
    rakhi<-seq(0, d0, length.out=12)
    rakhi<-rakhi[c(-1, -12)]
    output=matrix(0,length(rakhi),ncol=3)
    est=matrix(0, length(rakhi), ncol=ncol(X))
    
    if (thresh=="data"){
      i=0
      for (delta in rakhi){
        i=i+1
        out=dantzig(Xs,Yc,delta,gamma=0, vlength)
        est[i,]<-out[[3]]
      }
      gamma<-0.1*max(abs(est))
    } else {
      gamma<-1
    }
    
    output=matrix(0,length(rakhi),ncol=3)
    est=matrix(0, length(rakhi), ncol=ncol(X))
    i=0
    for (delta in rakhi){
      i=i+1
      out=dantzig(Xs,Yc,delta,gamma, vlength)
      beta=out[[3]]
      beta=as.vector(beta)
      p=out[[1]]+1
      if (sum(out[[2]])==0) { Xmodel=rep(1,dim(X)[1])
      } else { Xmodel=X[,out[[2]]] }
      lin=lm(Y~Xmodel)
      SSE=deviance(lin)
      if (sc==1){
        crit=n*log(SSE/n)+(2*p)+((2*p*(p+1))/(n-p-1)) #AICc
      } else {
        crit<-n*log(SSE/n)+p*log(n) #BIC
      }
      est[i,]=beta
      output[i,1]=delta
      output[i,2]=p
      output[i,3]=crit
      
    }
  }
  crit.best=min(output[,3])
  delta.best=output[,1][which(output[,3]==min(output[,3]))[1]]
  model.terms=dantzig(Xs,Yc,delta.best,gamma, vlength)
  answer=list(terms=model.terms, delta=delta.best, estimates=est)#, plot=p, gamma=gamma)
  return(answer)
  
}










library(doParallel)
library(foreach)
library(glmnet)
library(ggplot2)

## Function to return true coefficient vector
### a=number active
### k=number of linear effect columns
### p=probability of sign flipping
### active= index of active columns
### aceoff= acitve coef vector
### inacoeff=inactive effect coeff vector
### X=model matrix, linear effects
### error=sigma for norm error

true_y<-function(X, k, size, p, sn, error,active){
  effects=1:k
  n = dim(X)[1]
  a = length(active)
  #if(size=="med") { a=floor(n/3) } else if (size=="hi") {a=3} else {a=floor(0.75*n)}
  # active=sample(k, a, replace = FALSE, prob = NULL)
  true=effects%in%active
  inactive=which(true %in% c(0))
  R=(rbinom(a,1,p)*-2)+1 #creates a vector of 1 and -1 of signs
  #acoeff=(rexp(a, rate=1)+sn)*R
  acoeff=(rep(sn,a))*R
  inacoeff=rep(0,(k-a))
  Y=X[,active]%*%acoeff+X[,inactive]%*%inacoeff+rnorm(n,mean=0, sd=error)
  
  return(list(Y=Y, active=active, inactive=inactive))
  
}

### Fit lasso for one value of lambda
### X=design matrix, linear effects
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

gen_all_submodels<-function(p,k){
  return(t(combn(p,k)))
}




### simulation loop

## The line of code below needs to change
# setwd(".your_directory_here")

##

#paths<-dir(pattern=".csv")
# names(paths)<-basename(paths)
# names<-names(paths)
# names <- names
names<- c("PED", "Var_s_20", "Var_s_40", "Var_s_80")
print(names)
p=c(0) #all positive set 0, mixed set 0.5

sn=c(1, 3)#coefficient magnitude

size=c("med","hi", "low")  #number active

grid<-expand.grid(p,sn,size)
grid$Var3<-as.character(grid$Var3)
Des.perfect<-list()
Des.oneFP<-list()
Des.twoFP<-list()

n = 14
p = 24

### Generate the same sample of submodels for each sparsity setting so that we compare the same submodels for each scenario.

A_val_low = gen_all_submodels(p=p, k=ceiling(0.75*n))[sample(1:choose(p,ceiling(0.75*n)), 2000),]
A_val_med = gen_all_submodels(p=p, k=7)[sample(1:choose(p,7), 1500),]
A_val_hi =gen_all_submodels(p=p, k=ceiling(n/4))[sample(1:choose(p,ceiling(n/4)), 1000),]


for (j in 1:length(names)){
  
  X=as.matrix(candidate_designs[,,j])
  k=ncol(X)
  n=nrow(X)
  print(names[j])
  
  Xc<-apply(X, 2, function(X) X - mean(X))
  vlength<-sqrt(diag(t(Xc)%*%Xc))
  Xs<-t(t(Xc)*(1/vlength))
  
  # Generating a matrix to hold perfect, one FP and two FP results for each lambda considered. Each row will be a different data generation
  mat<-matrix(ncol=100, nrow=0)
  
  scenario.perfect<-data.frame(mat)
  scenario.oneFP<-data.frame(mat)
  scenario.twoFP<-data.frame(mat)
  for (h in 1:nrow(grid)){
    if(grid$Var3[h]=="med"){
      A_val = A_val_med
    } else if(grid$Var3[h]=="hi"){
      A_val = A_val_hi
    }else{
      A_val = A_val_low
    }
    
    ### Double check how many clusters you want to use
    cl=makeCluster(16)
    registerDoParallel(cl)
    
    sim<-foreach(i=1:nrow(A_val), .combine = rbind) %dopar% {
      perfect_df =data.frame(mat)
      oneFP_df=data.frame(mat)
      twoFP_df =  data.frame(mat)
      for( l in 1:25){
        resp<-true_y(X, k, size=grid$Var3[h], p=grid$Var1[h], sn=grid$Var2[h], error=1, active =A_val[i,] )
        Y=resp$Y
        Yc=Y-mean(Y)
        
        fit<- dantzig.selection(Xc, Yc, sc=1, l.grid="fine", thresh =0)
        #fit_lasso<-lasso.fit(X, Y=resp$Y)
        power<-c()
        type1<-c()
        
        for (d in 1:100) {
          sel = which(fit$estimates[d,]!=0)
          #sel<-as.numeric(gsub("X", "", subset(fit_lasso, step==d)$term))
          power[d]<-length(which(sel%in%resp$active)) #correct effects
          type1[d]<-length(which(sel%in%resp$inactive)) #incorrect effects
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
        
        
        perfect_df[l,]=perfect_vec
        oneFP_df[l,]=oneFP_vec
        twoFP_df[l,]=twoFP_vec
        
      }
      
      return(list(perfect= (colMeans(perfect_df)), oneFP= (colMeans(oneFP_df)), twoFP=(colMeans(twoFP_df))))
    }
    
    stopCluster(cl)
    perfect_sum<-sim[1,1][[1]]
    oneFP_sum<-sim[1,2][[1]]
    twoFP_sum<-sim[1,3][[1]]
    # oneFP<-c()
    # twoFP<-c()
    for (g in 2:nrow(A_val)) {
      # perfect[g]<-sim[g,]$perfect
      # oneFP[g]<-sim[g,]$oneFP
      # twoFP[g]<-sim[g,]$twoFP
      perfect_sum =perfect_sum+ sim[g,1][[1]]
      oneFP_sum =oneFP_sum+ sim[g,2][[1]]
      twoFP_sum =twoFP_sum+ sim[g,3][[1]]
    }
    
    scenario.perfect[h,]<-perfect_sum/nrow(A_val)
    scenario.oneFP[h,]<-oneFP_sum/nrow(A_val)
    scenario.twoFP[h,]<-twoFP_sum/nrow(A_val)
    
    
    
    
    
    # perfect<-c()
    # oneFP<-c()
    # twoFP<-c()
    # for (g in 1:nrow(A_val)) {
    #   perfect[g]<-sim[g,]$perfect
    #   oneFP[g]<-sim[g,]$oneFP
    #   twoFP[g]<-sim[g,]$twoFP
    # }
    # 
    # scenario.perfect[h]<-mean(perfect)
    # scenario.oneFP[h]<-mean(oneFP)
    # scenario.twoFP[h]<-mean(twoFP)
    
  }#close h loop
  
  Des.perfect[[j]]<-scenario.perfect
  
  Des.oneFP[[j]]<-scenario.oneFP
  Des.twoFP[[j]]<-scenario.twoFP
  print(j)
}#close j loop






# log_lam = seq(from = -4.5, to= 1.5, by = 0.1)
# lambda = exp(log_lam)

delta = seq(0,16,length.out=100)

designs =names


df_results = merge( designs,grid, by=NULL)

df_results = merge(delta, df_results, by=NULL)

colnames(df_results)= c("delta", "design", "Signs", "Beta", "Sparsity")

colnames(grid)= c("Signs", "Beta", "Sparsity")

perfect_results = cbind(grid, Des.perfect[[1]])
perfect_results$design = designs[1]

oneFP_results = cbind(grid, Des.oneFP[[1]])
oneFP_results$design = designs[1]

twoFP_results = cbind(grid, Des.twoFP[[1]])
twoFP_results$design = designs[1]

for(j in c(2:4)){
  temp = cbind(grid, Des.perfect[[j]])
  temp$design=designs[j]
  perfect_results = rbind(perfect_results, temp)
  
  temp = cbind(grid, Des.oneFP[[j]])
  temp$design=designs[j]
  oneFP_results= rbind(oneFP_results, temp)
  
  temp = cbind(grid, Des.twoFP[[j]])
  temp$design=designs[j]
  twoFP_results = rbind(twoFP_results, temp)
  
  
}

perfect_results_long = perfect_results %>% pivot_longer(!c("Signs", "Beta", "Sparsity","design"), names_to="temp", values_to="perfect")

oneFP_results_long = oneFP_results %>% pivot_longer(!c("Signs", "Beta", "Sparsity","design"), names_to="temp", values_to="oneFP")


twoFP_results_long = twoFP_results %>% pivot_longer(!c("Signs", "Beta", "Sparsity","design"), names_to="temp", values_to="twoFP")



df_results = cbind(df_results, perfect_results_long$perfect, oneFP_results_long$oneFP,twoFP_results_long$twoFP)
colnames(df_results)[6]="perfect_prob"
colnames(df_results)[7]="oneFP_prob"
colnames(df_results)[8]="twoFP_prob"






write.csv(df_results, "dantzig_sim_perfect_SR_n14_k24.csv", col.names=TRUE, row.names=FALSE)
#write.csv(df_results, "david_paper_lasso_sim_large_ME.csv", col.names=TRUE, row.names=FALSE)

#df_results
# 
# q= ggplot(data = df_results, aes(x=log(as.numeric(lambda)), y= as.numeric(perfect_prob), color= design, linetype = design))+
#   geom_line()+
#   facet_grid(vars(ME), vars(two_FI))+
#   xlab(latex2exp::TeX("$log(\\lambda)"))+
#   ylab(latex2exp::TeX(" $\\phi_{\\lambda}"))+
#   #scale_color_manual(values =c("blue","red", "green"))+
#   theme_bw()+
#   #scale_x_continuous(limits = c(-2,2))+
#   # scale_linetype_manual(values =c("solid", "dashed", "solid"))+
#   theme( legend.title=element_blank(), text=element_text(size=16))
# 
# 
# 
# q_trimmed= ggplot(data = df_results[(df_results$two_FI%in% c(1,4,7)),], aes(x=log(as.numeric(lambda)), y= as.numeric(perfect_prob), color= design, linetype=design))+
#   geom_line(size=1.2)+
#   facet_grid(vars(ME), vars(two_FI))+
#   xlab(latex2exp::TeX("$log(\\lambda)"))+
#   ylab(latex2exp::TeX("$\\phi_{\\lambda}"))+
#   #scale_color_manual(values =c("blue","red", "green"))+
#   theme_bw()+
#   #scale_x_continuous(limits = c(-2,2))+
#   # scale_linetype_manual(values =c("solid", "dashed", "solid"))+
#   theme( legend.title=element_blank(), text=element_text(size=22))
# 
# ggsave("david_lasso_sim_overall_trimmed.pdf",q_trimmed)
# 
# 
# q_me= ggplot(data = df_results, aes(x=log(as.numeric(lambda)), y= as.numeric(perfect_me), color= design))+
#   geom_line()+
#   facet_grid(vars(ME), vars(two_FI))+
#   xlab(latex2exp::TeX("$log(\\lambda)"))+
#   ylab(latex2exp::TeX("ME $\\phi_{\\lambda}"))+
#   #scale_color_manual(values =c("blue","red", "green"))+
#   theme_bw()+
#   #scale_x_continuous(limits = c(-2,2))+
#   # scale_linetype_manual(values =c("solid", "dashed", "solid"))+
#   theme( legend.title=element_blank(), text=element_text(size=16))
# 
# ggsave("david_lasso_sim_me.pdf",q_me)
# 
# 
# q_me_trimmed= ggplot(data = df_results[(df_results$two_FI%in% c(1,4,7)),], aes(x=log(as.numeric(lambda)), y= as.numeric(perfect_me), color= design, linetype=design))+
#   geom_line(size=1.2)+
#   facet_grid(vars(ME), vars(two_FI))+
#   xlab(latex2exp::TeX("$log(\\lambda)"))+
#   ylab(latex2exp::TeX("ME $\\phi_{\\lambda}"))+
#   #scale_color_manual(values =c("blue","red", "green"))+
#   theme_bw()+
#   #scale_x_continuous(limits = c(-2,2))+
#   # scale_linetype_manual(values =c("solid", "dashed", "solid"))+
#   theme( legend.title=element_blank(), text=element_text(size=22))
# 
# ggsave("david_lasso_sim_me_trimmed.pdf",q_me_trimmed)
# 
# 
# 
# 
# 
# 
# q_2FI= ggplot(data = df_results, aes(x=log(as.numeric(lambda)), y= as.numeric(perfect_2FI), color= design, linetype = design))+
#   geom_line()+
#   facet_grid(vars(ME), vars(two_FI))+
#   xlab(latex2exp::TeX("$log(\\lambda)"))+
#   ylab(latex2exp::TeX("2FI $\\phi_{\\lambda}"))+
#   #scale_color_manual(values =c("blue","red", "green"))+
#   theme_bw()+
#   #scale_x_continuous(limits = c(-2,2))+
#   # scale_linetype_manual(values =c("solid", "dashed", "solid"))+
#   theme( legend.title=element_blank(), text=element_text(size=16))
# 
# 
# q_2FI_trimmed= ggplot(data = df_results[(df_results$two_FI%in% c(1,4,7)),], aes(x=log(as.numeric(lambda)), y= as.numeric(perfect_2FI), color= design, linetype=design))+
#   geom_line(size=1.2)+
#   facet_grid(vars(ME), vars(two_FI))+
#   xlab(latex2exp::TeX("$log(\\lambda)"))+
#   ylab(latex2exp::TeX("2FI $\\phi_{\\lambda}"))+
#   #scale_color_manual(values =c("blue","red", "green"))+
#   theme_bw()+
#   #scale_x_continuous(limits = c(-2,2))+
#   # scale_linetype_manual(values =c("solid", "dashed", "solid"))+
#   theme( legend.title=element_blank(), text=element_text(size=22))
# 
# ggsave("david_lasso_sim_2FI_trimmed.pdf",q_2FI_trimmed)
# 
# q_2_1 =ggplot(data = df_results[df_results$ME==2 & df_results$two_FI==1,], aes(x=log(as.numeric(lambda)), y= as.numeric(perfect_prob), color= design, linetype = design))+
#   geom_line()+
#   xlab(latex2exp::TeX("$log(\\lambda)"))+
#   ylab(latex2exp::TeX("$\\phi_{\\lambda}"))+
#   #scale_color_manual(values =c("blue","red", "green"))+
#   theme_bw()+
#   #scale_x_continuous(limits = c(-2,2))+
#   # scale_linetype_manual(values =c("solid", "dashed", "solid"))+
#   theme( legend.title=element_blank(), text=element_text(size=16), legend.position = c(0.1,.75))
# 
# 
# q_2_1_me =ggplot(data = df_results[df_results$ME==2 & df_results$two_FI==1,], aes(x=log(as.numeric(lambda)), y= as.numeric(perfect_me), color= design, linetype = design))+
#   geom_line()+
#   xlab(latex2exp::TeX("$log(\\lambda)"))+
#   ylab(latex2exp::TeX("$\\phi_{\\lambda}"))+
#   #scale_color_manual(values =c("blue","red", "green"))+
#   theme_bw()+
#   #scale_x_continuous(limits = c(-2,2))+
#   # scale_linetype_manual(values =c("solid", "dashed", "solid"))+
#   theme( legend.title=element_blank(), text=element_text(size=16), legend.position = c(0.1,.75))
# 
# df_results_early_example = df_results[df_results$ME==2 & df_results$two_FI==1,]
# 
# df_results_early_example = df_results_early_example[df_results_early_example$design=="BayesD"| df_results_early_example$design=="PEC",]
# 
# df_results_early_example = df_results_early_example %>% mutate(design = recode(design, "BayesD"="Design A", "PEC"="Design B"))  
# 
# 
# 
# q_early =ggplot(data = df_results_early_example, aes(x=log(as.numeric(lambda)), y= as.numeric(perfect_me), color= design))+
#   geom_line(size=1.2)+
#   xlab(latex2exp::TeX("$log(\\lambda)"))+
#   ylab(latex2exp::TeX("$\\phi_{\\lambda}"))+
#   scale_color_manual(values =c("blue","red"))+
#   theme_bw()+
#   #scale_x_continuous(limits = c(-2,2))+
#   # scale_linetype_manual(values =c("solid", "dashed", "solid"))+
#   theme( legend.title=element_blank(), text=element_text(size=14))
# 
# 
# ggsave("lasso_prob_early_example.pdf",q_early)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
