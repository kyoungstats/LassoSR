### This script houses all functions used for the SIC paper

library(mvtnorm)
library(dplyr)
library(tidyr)
library(MASS)
library(ibd)
library(doParallel)
library(foreach)
library(HadamardR)
library(readr)
library(ggpubr)
#library(MixMatrix)





read_in_design<-function(path){
  df <- read_table2(path,
                    col_names = FALSE)
  
  df <- df[, 1:ncol(df)-1]
  return(df)
}
#library(doSNOW)

#' Title
#'
#' @param F_A_cs  A centered and scaled submatrix of the active effects
#' @param V_half  Scaling matrix as definded in the paper
#' @param lambda  Scalar; LASSO Tuning Param
#' @param B_mag Vector of magnitude of active effects
#' @param Z_A Sign vector corresponding to Bvec
#' @param sigma Error variance set to 1
#'
#' @return P_S, probability of S event
#' @export
#'
#' @examples
Power_local <-function(F_A_cs, V_half, lambda, B_mag, Z_A, sigma = 1){
  # Consider adding option for lambda to be on log scale
  Bvec = B_mag*Z_A
  n = dim(F_A_cs)[1]
  k = dim(F_A_cs)[2]
  # The next line checks if  F_A is full column rank if not, we give it a 0 probability
  F_A_cs= as.matrix(F_A_cs)
  if (!is.null(dim(tryCatch(solve(t(F_A_cs)%*%F_A_cs), error = function(u)0)))){
    C_AA_inv = n*solve(t(F_A_cs)%*%F_A_cs)
    Cov_S = (sigma^2)* diag(Z_A)%*%C_AA_inv%*%diag(Z_A)
    upper_S = sqrt(n)*diag(Z_A)%*%V_half%*%Bvec -
      lambda*sqrt(n)* diag(C_AA_inv%*%Z_A%*%t(Z_A))
    if(any(is.na(upper_S))){
      P_S =0
    }
    else{
      P_S= mvtnorm:: pmvnorm(lower= -Inf, upper = as.vector(upper_S), mean= rep(0,k), sigma = Cov_S)[1]
    }
  }
  else {
    P_S = 0
  }
  return(P_S)
}


#' Power_local_log_lam
#'
#' @param F_A_cs  A centered and scaled submatrix of the active effects
#' @param V_half  Scaling matrix as definded in the paper
#' @param log_lambda  Scalar; log of LASSO Tuning Param
#' @param B_mag Vector of magnitude of active effects
#' @param Z_A Sign vector corresponding to Bvec
#' @param sigma Error variance set to 1
#'
#' @return P_S, probability of S event
#' @export
#'
#' @examples
Power_local_log_lam <-function(F_A_cs, V_half, log_lambda, B_mag, Z_A, sigma = 1){
  # THis is the same as Power_local, but with lambda on the log scale
  Bvec = B_mag*Z_A
  n = dim(F_A_cs)[1]
  k = dim(F_A_cs)[2]
  lambda = exp(log_lambda)
  # The next line checks if  F_A is full column rank if not, we give it a 0 probability
  F_A_cs= as.matrix(F_A_cs)
  if (!is.null(dim(tryCatch(solve(t(F_A_cs)%*%F_A_cs), error = function(u)0)))){
    C_AA_inv = n*solve(t(F_A_cs)%*%F_A_cs)
    Cov_S = (sigma^2)* diag(Z_A)%*%C_AA_inv%*%diag(Z_A)
    upper_S = sqrt(n)*diag(Z_A)%*%V_half%*%Bvec -
      lambda*sqrt(n)* diag(C_AA_inv%*%Z_A%*%t(Z_A))
    P_S= mvtnorm:: pmvnorm(lower= -Inf, upper = as.vector(upper_S), mean= rep(0,k), sigma = Cov_S)[1]
  }
  else {
    P_S = 0
  }
  return(P_S)
}
### Function for creating the diagonals of V_half, use lapply to find V
fn<-function(x){
  sqrt(sum((x-mean(x))^2)/length(x))
}
### Function for centering and scaling, use lapply to center and scale F
cent_scale <- function(x){
  (x - mean(x))/ (sqrt((sum((x-mean(x))^2))/length(x)))
}

get_cs_and_V <- function(mat){
  F_cs = as.data.frame(lapply(mat, cent_scale))
  V_half=diag(lapply(mat,fn))
  return(list("F_cs"= F_cs, "V_half"=V_half))
}
#' Title
#'
#' @param F_A_start Initial uncentered, unscaled design matrix to start coord exchange from
#' @param lambda scalar, Lasso tuning param
#' @param Bvec Vector of magnitude of active effects
#' @param sigma scalar; error variance, default to 1
#' @param Z_A sign vector of active effects, if set to NULL, then average over all sign vectors
#'
#' @return list of best design and the probability of S associated with that design
#' @export
#'
#' @examples
Local_coord_ex <- function(F_A_start, lambda, B_mag, sigma =1,Z_A =NULL){
  F_0 = F_A_start
  n = dim(F_0)[1]
  k = dim(F_0)[2]
  V_0 = diag(lapply(as.data.frame(F_0),fn))
  F_0_cs = as.data.frame(lapply(as.data.frame(F_0), cent_scale))
  #browser()
  if(is.null(Z_A)){
    print("Using Avg Prob over all sign vects as metric to construct design. ")
    P_S_0 = mean(search_sign_vects(F_0_cs,V_0, B_mag, log(lambda), sigma)$results)
    iterations = 0
    count_between_flips=0
    stop = FALSE
    while ( iterations <= 5 ){
      #print(paste0("Iteration: ", iterations))
      for (i in c(1:n)){
        for (j in c(1:k)){
          #print(paste0("Count Between Flips: ", count_between_flips))
          F_1 = F_0
          #browser()
          F_1[i,j]= -F_0[i,j]
          #V_1 <- update_V(V_0, i, j, F_0, n)
          V_1 = V_0
          V_1[j,j]<- fn(F_1[,j])
          F_1_cs = F_0_cs
          F_1_cs[,j]= cent_scale(F_1[,j])
          #browser()
          P_S_1 = mean(search_sign_vects(F_1_cs,V_1, B_mag, log(lambda), sigma)$results)
          if( P_S_1 > P_S_0) {
            F_0 = F_1
            V_0 =V_1
            F_0_cs=F_1_cs
            P_S_0 = P_S_1
            count_between_flips = 0
          }
          else{
            count_between_flips = count_between_flips+1
          }
          if(count_between_flips> n*k){
            #print("break")
            stop = TRUE
            break
          }
          #Breaking out of outter loops

        }
        if(stop){break}

      }
      if(stop){break}
      iterations=iterations+1
    }
  }
    else{
      P_S_0 = Power_local(F_0_cs, V_half= V_0, lambda, B_mag, Z_A, sigma)
      iterations = 0
      count_between_flips=0
      stop = FALSE
      while ( iterations <= 10 ){
        #print(paste0("Iteration: ", iterations))
        for (i in c(1:n)){
          for (j in c(1:k)){
            F_1 = F_0
            #browser()
            F_1[i,j]= -F_0[i,j]
            #V_1 <- update_V(V_0, i, j, F_0, n)
            V_1 = V_0
            V_1[j,j]<- fn(F_1[,j])
            F_1_cs = F_0_cs
            F_1_cs[,j]= cent_scale(F_1[,j])
            #browser()
            P_S_1 = Power_local(F_1_cs, V_half= V_1, lambda, B_mag, Z_A, sigma)
            if( P_S_1 > P_S_0) {
              #print("FLIPPED")
              F_0 = F_1
              V_0 =V_1
              F_0_cs=F_1_cs
              P_S_0 = P_S_1
              count_between_flips = 0
            }
            else{
              count_between_flips = count_between_flips+1
            }
            if(count_between_flips> n*k){
              #print("break")
              stop = TRUE
              break
            }
            # Breaking from the next two outter loops

          }
          if(stop){break}

        }
        if(stop){break}
        iterations=iterations+1
      }
    }


  return(list("design"= F_0, "prob"= P_S_0))
}


#' Title
#'
#' @param F_A_cs A centered and scaled submatrix of the active effects
#' @param V_half  Scaling matrix as definded in the paper
#' @param B_mag Vector of magnitude of active effects
#' @param log_lambda Scalar; LASSO Tuning Param
#' @param sigma Scalar; Error variance set to 1
#'
#' @return dataframe with a first k columns denoting the sign vector and results column giving cooresponding prob
#' @export
#'
#' @examples
search_sign_vects <- function( F_A_cs,V_half, B_mag, log_lambda, sigma = 1){

  k = dim(F_A_cs)[2]
  l=rep(list(c(-1,1)), as.integer(k))
  # expand.grid will give a df where each row is a unique sign vector. The first half of the rows are the negartives of the second half
  # a is the index of the half way point of the rows
  #b is the largest row index
  a = ((2^(k-1))+1)
  b = 2^k
  sign_vec_grid = expand.grid(l)[a:b,]
  #sign_vec_grid = expand.grid(l)
  colnames(sign_vec_grid) <- paste("s", 1:k,sep="")
  results <- apply(sign_vec_grid, 1, Power_local_log_lam, F_A_cs = F_A_cs, V_half = V_half,
                   log_lambda= log_lambda, B_mag = B_mag, sigma=sigma )
  output = cbind(sign_vec_grid, results)
  return(output)
}

#' integral_over_lam
#'
#' @param F_A_cs A centered and scaled submatrix of the active effects
#' @param V_half Scaling matrix as definded in the paper
#' @param lambda Scalar; LASSO Tuning Param. If log_scale is TRUE, this is the log of the LASSO tuning param
#' @param Bvec Vector of active effects, sign included
#' @param Z_A Sign vector corresponding to Bvec
#' @param sigma Scalar; Error variance set to 1
#' @param lam_min Scalar;Minimum of lambda range to integrate over, or log_lambda if log scale is TRUE
#' @param lam_max Scalar; Maximum of lambda range to integrate over, or log_lambda if log scale is TRUE
#' @param step_size Scalar; stepsize between points on reimman integration grid
#' @param log_scale Boolean; if TRUE, all lambdas are in the log scale

#'
#' @return Scalar; intergal of P_S over lambda or log_lambda grid.
#' @export
#'
#' @examples
integral_over_lam<- function(F_A_cs, V_half, lambda, Bvec, Z_A, sigma = 1, lam_min=0, lam_max = 4,
                             step_size = 0.01, log_scale = FALSE){
  #This function is the Reimman integral approximation of the power of a local optimal design.
  if (log_scale){
    log_lambda_grid <- seq(lam_min, lam_max, by=step_size)

    return(sum(step_size*as.data.frame(lapply(lambda_grid, Power_local_log_lam, F_A_cs=F_A_cs,V_half=V_half, Bvec=rep(3,7), Z_A= rep(1,7), sigma=1 ))))

  }
  else{
    lambda_grid<- seq(lam_min, lam_max, by=step_size)
    return(sum(step_size*as.data.frame(lapply(lambda_grid, Power_local, F_A_cs=F_A_cs,V_half=V_half, Bvec=rep(3,7), Z_A= rep(1,7), sigma=1 ))))

  }
}




##### Non-local Coordinate Exchange Functions

check_inverse <- function(M){
  class(try(solve(M), silent = T))=="matrix"
}

fix_neg_diag <- function(M){
  #browser()
  if (any(diag(M)<=0)){
    # When we have negative diags in covariance matrix for I, they are very close to zero.
    # So, using the absolute value will flip them to slightly positive
    diag(M)=abs(diag(M))+ 1e-5
  }
  return(M)
}


#'Joint_prob_all_fixed
#'
#' @param A Vector of column indexes in the active set
#' @param F_cs Full centered and scaled design matrix
#' @param V_half Full sqrt of the scale matrix
#' @param B_mag Vector of magnitude of active effects, should be all positive
#' @param Z_A Sign vector of active effects
#' @param sigma Scalar, variance parameter, default is 1
#' @param log_lambda scalar, log of tuning parameter
#'
#' @return P_joint, scalar probability of sign recovery
#' @export
#'
#' @examples
Joint_prob_all_fixed <-function(A, F_cs, V_half, B_mag, Z_A, log_lambda, sigma = 1){
  # Subsetting design into active and inactive columns
  print(A)
  F_A=as.matrix(F_cs[,A])
  F_I = as.matrix(F_cs[,-A])
  V_half=diag(diag(V_half)[A])

  lambda = exp(log_lambda)

  n = dim(F_cs)[1]
  k = length(Z_A)
  one_k = rep(1,k)
  one_q = rep(1, dim(F_I)[2])
  B_min_vec = B_mag*Z_A
  C_IA = (1/n)*t(F_I)%*%F_A
  C_II = (1/n)*t(F_I)%*%F_I
  C_AA = (1/n)*t(F_A)%*%F_A



  # Check if C_AA has an inverse, if not, use moore-penrose. (Discuss with group)
    C_AA_inv = n*ginv(t(F_A)%*%F_A)
    Mean_I = lambda*sqrt(n)*C_IA%*%C_AA_inv%*%Z_A
    Cov_I = (sigma^2)*(C_II - (C_IA%*%C_AA_inv%*%t(C_IA)))

    # Fix the diagonals of Cov_I so that they are slightly positive, if they are slightly negative
    Cov_I = fix_neg_diag(Cov_I)
    P_I = mvtnorm::pmvnorm(lower = -sqrt(n)*lambda*one_q, upper = sqrt(n)*lambda*one_q,
                           mean = as.vector(Mean_I), sigma = as.matrix(Cov_I), algorithm = GenzBretz() )[1]
    Cov_S = (sigma^2)* diag(Z_A)%*%C_AA_inv%*%diag(Z_A)
    upper_S = sqrt(n)*diag(Z_A)%*%V_half%*%B_min_vec -
      lambda*sqrt(n)* diag(C_AA_inv%*%Z_A%*%t(Z_A))
    #browser()
    P_S= mvtnorm:: pmvnorm(lower= -Inf, upper = as.vector(upper_S), mean= rep(0,k), sigma = Cov_S)[1]
    P_joint = P_S*P_I

 return(P_joint)

}


Joint_prob_all_fixed_deconstructed<-function(lower_I_WO_lam, upper_I_WO_lam, mean_I_WO_lam, Cov_I, Cov_S, upper_S,mean_S_WO_lam, log_lambda){

  lambda = exp(log_lambda)
  k = dim(Cov_S)[1]


  Mean_I = mean_I_WO_lam*lambda
  

  P_I = mvtnorm::pmvnorm(lower = lambda*as.vector(lower_I_WO_lam), upper = lambda * as.vector(upper_I_WO_lam),
                         mean = as.vector(Mean_I), sigma = as.matrix(Cov_I), algorithm = GenzBretz() )[1]
  #browser()

  P_S= mvtnorm:: pmvnorm(lower= -Inf, upper = as.vector(upper_S)-lambda*as.vector(mean_S_WO_lam), mean= rep(0,k), sigma = Cov_S)[1]


  P_joint = P_S*P_I
  return(P_joint)

}
Joint_prob_all_fixed_deconstructed_more_output<-function(lower_I_WO_lam, upper_I_WO_lam, mean_I_WO_lam, Cov_I, Cov_S, upper_S,mean_S_WO_lam, log_lambda, output_option="joint"){
  
  lambda = exp(log_lambda)
  k = dim(Cov_S)[1]
  
  
  Mean_I = mean_I_WO_lam*lambda
  
  P_I = mvtnorm::pmvnorm(lower = lambda*as.vector(lower_I_WO_lam), upper = lambda * as.vector(upper_I_WO_lam),
                         mean = as.vector(Mean_I), sigma = as.matrix(Cov_I), algorithm = GenzBretz() )[1]
  #browser()
  
  P_S= mvtnorm:: pmvnorm(lower= -Inf, upper = as.vector(upper_S)-lambda*as.vector(mean_S_WO_lam), mean= rep(0,k), sigma = Cov_S)[1]
  
  
  P_joint = P_S*P_I
  if (output_option =="joint"){ return(P_joint)}
  if (output_option =="I"){return(P_I)}
  else{ return(P_S)}
  
}

### Strategies to handle lambda


#' Title
#'
#' @param log_lam_min Scalar, minimum for grid of log lambdas for a warm start
#' @param log_lam_max Scalar, maximum for grid of log lambdas for a warm start
#' @param A Vector of column indexes in the active set
#' @param F_cs Full centered and scaled design matrix
#' @param V_half Full sqrt of the scale matrix
#' @param B_mag Vector of magnitude of active effects, should be all positive
#' @param Z_A Sign vector of active effects
#' @param sigma Scalar, variance parameter, default is 1
#' @param stepsize Scalar, stepsize between grid points for warm start, default is 0.05
#'
#' @return list with one item denoting the optimal log lambda value, and another giving the P_joint value at the optimal log lambda
#' @export
#'
#' @examples
opt_log_lambda <- function( log_lam_min, log_lam_max, F_cs, V_half, B_mag, k, submodels, sign_vects, sigma =1, stepsize= 0.05){
  # Get warm start for lamnda values
  #browser()
  log_lambda_grid = seq(log_lam_min,log_lam_max,stepsize)
  grid_probs <-sapply(log_lambda_grid, Joint_prob_general_fixed_lambda, F_cs=F_cs, V_half = V_half, B_mag=B_mag, k=k, submodels=submodels, sign_vects=sign_vects,
                      log_lambda_strategy="fixed", sigma=sigma, output_option="joint")
  warm_start <- cbind(log_lambda_grid,grid_probs)
  #print(max(warm_start[,2]))
  #browser()

  if (is.na(max(warm_start[,2]))|| max(warm_start[,2])<0.0001){
    start_lam <- runif(1, 0, 3)
  }
  else{
    start_lam <- warm_start[which.max(warm_start[,2]),1]
  }
  # maybe consider doing box-constaints, but since log-lam can take any real number, I thought other solvers could be faster
  #browser()
  opt_lam <-optim(par=start_lam, fn=Joint_prob_general_fixed_lambda, F_cs=F_cs, V_half = V_half, B_mag=B_mag, k=k, submodels=submodels, sign_vects=sign_vects,
                           log_lambda_strategy="fixed", sigma=sigma, output_option="joint",
                           method="Brent",lower=-5, upper =5,
                           control = list("fnscale"=-1,warn.1d.NelderMead=FALSE))
  # opt_lam <-tryCatch(optim(par=start_lam, fn=Joint_prob_all_fixed_deconstructed,lower_I_WO_lam=lower_I_WO_lam,
  #                 upper_I_WO_lam=upper_I_WO_lam, mean_I_WO_lam=mean_I_WO_lam,
  #                 Cov_I=Cov_I, Cov_S = Cov_S, upper_S= upper_S, mean_S_WO_lam= mean_S_WO_lam,
  #                 method="L-BFGS-B", upper = 10, lower = -Inf,
  #                 control = list("fnscale"=-1,warn.1d.NelderMead=FALSE)), error = function(c){return(list("par"=start_lam, "value"=0))})
  return(list("opt_log_lam"= opt_lam$par, "opt_val" = opt_lam$value))


}



#' This gives the area under the joint probability curve as it depends on log lambda
#'
#' @param A Vector of column indexes in the active set
#' @param F_cs Full centered and scaled design matrix
#' @param V_half Full sqrt of the scale matrix
#' @param B_mag Vector of magnitude of active effects, should be all positive
#' @param Z_A Sign vector of active effects
#' @param sigma Scalar, variance parameter, default is 1
#' @param int_lower_bound Scalar, lower bound on region of integration, default is infinity
#' @param int_upper_bound Scalar, upper bound on region of integration, default is 3
#' @param method String,  Quadrature or Sum, Sum indicates Riemann sum
#'
#'
#' @return resulting value of integral
#' @export
#'
#' @examples
integrate_log_lambda <- function(lower_I_WO_lam, upper_I_WO_lam,
                                 mean_I_WO_lam, Cov_I, Cov_S, upper_S, mean_S_WO_lam, int_lower_bound=-5, int_upper_bound = 1.5, method = "Quadrature", step_size = 0.02){
 

  
  int_methods = c("Quadrature", "Sum")
  if(!method%in% int_methods){
    stop("The input integration method must be either  'Quadrature', or 'Sum'")
  }
  if (method == "Quadrature"){
    Joint_prob_all_fixed_vec = Vectorize(Joint_prob_all_fixed_deconstructed, "log_lambda")
    #browser()
    #print("Integrating via quadrature")
    
    integrate_results = integrate(Joint_prob_all_fixed_vec, lower= int_lower_bound, upper= int_upper_bound,
                                  lower_I_WO_lam=lower_I_WO_lam,
                                  upper_I_WO_lam=upper_I_WO_lam,
                                  mean_I_WO_lam=mean_I_WO_lam,
                                  Cov_I=Cov_I, Cov_S = Cov_S,
                                  upper_S= upper_S, mean_S_WO_lam= mean_S_WO_lam)
    return(integrate_results$value)
  }
  else{
    #print("Integrating via Riemann Sum")
    log_lambda_grid <- seq(int_lower_bound, int_upper_bound, by=step_size)
  
    
    above_0_yet = FALSE
    sum_prob = Joint_prob_all_fixed_deconstructed(lower_I_WO_lam, upper_I_WO_lam, mean_I_WO_lam, Cov_I, Cov_S, upper_S, mean_S_WO_lam, log_lambda = int_lower_bound)
    for( j in c(1:1000)){
      log_lambda = int_lower_bound + step_size*j
      #print(log_lambda)
      prob = tryCatch({Joint_prob_all_fixed_deconstructed(lower_I_WO_lam, upper_I_WO_lam, mean_I_WO_lam, Cov_I, Cov_S, upper_S, mean_S_WO_lam, log_lambda)}, error =function(c){return(0)})
      sum_prob = sum_prob + prob
      if (is.na(prob)){
        #something with MVN algorithm givens NaNs sometimes
        prob =0
      }
     
      #print(paste0("prob: ", prob))
      if( prob > 0.0001 && above_0_yet == FALSE){ 
        #print("prob is above 0 for first time")
        above_0_yet= TRUE}
      if( prob < 0.0001 && above_0_yet == TRUE){ 
        #print("breaking loop")
        break}
      
    }
    return(step_size* sum_prob)
  }
  


}


#'This is a function for a fixed submodel and sign vector, with an option for
#'
#' @param A Vector of column indexes in the active set
#' @param F_cs Full centered and scaled design matrix
#' @param V_half Full sqrt of the scale matrix
#' @param B_mag Vector of magnitude of active effects, should be all positive
#' @param Z_A Sign vector of active effects
#' @param sigma Scalar, variance parameter, default is 1
#' @param log_lambda Scalar, fixed value of log_lambda, this is set to NULL, if the log_lambda strategy is fixed, this must be defined
#' @param log_lambda_strategy string, strategy to handle log_lambdas, takes in either "fixed", "optimal", "integrate"
#' @param int_lower_bound
#' @param int_upper_bound
#'
#' @return P_joint, scalar probability of sign recovery
#' @export
#'
#' @examples
Joint_prob_submodel_sign_fixed <-function(A, F_cs, V_half, B_mag, Z_A, log_lambda=NULL, log_lambda_strategy="integrate",
                                          sigma = 1, int_lower_bound=NULL, int_upper_bound=NULL, int_method = "Quadrature", output_option = "joint"){
  strategies = c("fixed", "optimal", "integrate")
  if(!log_lambda_strategy%in% strategies){
    stop("The input log_lam strategy must be either 'fixed, 'optimal', or 'integrate'")
  }
  # Subsetting design into active and inactive columns
  F_A=as.matrix(F_cs[,A])
  F_I = as.matrix(F_cs[,-A])
  V_half=diag(diag(V_half)[A])

  #lambda = exp(log_lambda)

  n = dim(F_cs)[1]
  k = length(Z_A)
  one_k = rep(1,k)
  one_q = rep(1, dim(F_I)[2])
  B_min_vec = B_mag*Z_A
  C_IA = (1/n)*t(F_I)%*%F_A
  C_II = (1/n)*t(F_I)%*%F_I
  C_AA = (1/n)*t(F_A)%*%F_A



  # Check if C_AA has an inverse, if not, set the joint prob to zero. (Discuss with group)
  
  C_AA_inv = n*ginv(t(F_A)%*%F_A)
  
  mean_I_WO_lam = sqrt(n)*C_IA%*%C_AA_inv%*%Z_A
  Cov_I = (sigma^2)*(C_II - (C_IA%*%C_AA_inv%*%t(C_IA)))
  upper_I_WO_lam = sqrt(n)*one_q
  lower_I_WO_lam = -upper_I_WO_lam
  
  # Fix the diagonals of Cov_I so that they are slightly positive, if they are slightly negative
  
  Cov_I = fix_neg_diag(Cov_I)
  Cov_S = (sigma^2)* diag(Z_A)%*%C_AA_inv%*%diag(Z_A)
  #browser()
  upper_S = sqrt(n)*diag(Z_A)%*%V_half%*%B_min_vec
  mean_S_WO_lam = sqrt(n)* diag(C_AA_inv%*%Z_A%*%t(Z_A))
  if( log_lambda_strategy=="fixed"){
    if(is.null(log_lambda)){
      stop("For a fixed lambda strategy, you must specify what the fixed log_lambda is. It is currently NULL ")
    }
    else{
      return(Joint_prob_all_fixed_deconstructed_more_output(lower_I_WO_lam, upper_I_WO_lam, mean_I_WO_lam, Cov_I, Cov_S, upper_S, mean_S_WO_lam, log_lambda, output_option))
    }
  }
  if( log_lambda_strategy=="integrate"){
    if(any(is.null(c(int_lower_bound, int_lower_bound)))){
      warning(" Integral upper and/or lower bounds are not specified, using defualt of lower bound at -5.")
      return(integrate_log_lambda(lower_I_WO_lam, upper_I_WO_lam,
                                  mean_I_WO_lam, Cov_I, Cov_S,
                                  upper_S= upper_S, mean_S_WO_lam= mean_S_WO_lam, method =int_method ))
    }
    else{ 
      return(integrate_log_lambda(lower_I_WO_lam, upper_I_WO_lam,
                                  mean_I_WO_lam, Cov_I, Cov_S,
                                  upper_S= upper_S, mean_S_WO_lam= mean_S_WO_lam,
                                  int_lower_bound=int_lower_bound ,
                                  int_upper_bound = int_upper_bound, method= int_method))
    }
    
    }

    if(log_lambda_strategy=="optimal"){
      #print(A)
      #browser()
      return(opt_log_lambda( log_lam_min=-2, log_lam_max=3, lower_I_WO_lam, upper_I_WO_lam,
                                         mean_I_WO_lam, Cov_I, Cov_S, upper_S= upper_S, mean_S_WO_lam= mean_S_WO_lam, stepsize= 0.05)$opt_val)
    }

  

}

### Now we move to a fixed submodel, but the ability to loop over sign vectors


gen_all_sign_vects<-function(A){
  k = length(A)
  l=rep(list(c(-1,1)), as.integer(k))
  # expand.grid will give a df where each row is a unique sign vector. The first half of the rows are the negartives of the second half
  # a is the index of the half way point of the rows
  #b is the largest row index
  a = ((2^(k-1))+1)
  b = 2^k
  sign_vec_grid = expand.grid(l)[a:b,]
  return(sign_vec_grid)
}

#'This is a function for a fixed submodel with the option to specify sign vectors and lambda strategy
#'
#' @param A Vector of column indexes in the active set
#' @param F_cs Full centered and scaled design matrix
#' @param V_half Full sqrt of the scale matrix
#' @param B_mag Vector of magnitude of active effects, should be all positive
#' @param sign_vects Matrix, with |A| columns where each row is a sign vector, if NULL then it loops over all sign vectors. Defualt is NULL.
#' @param sigma Scalar, variance parameter, default is 1
#' @param log_lambda Scalar, fixed value of log_lambda, this is set to NULL, if the log_lambda strategy is fixed, this must be defined
#' @param log_lambda_strategy string, strategy to handle log_lambdas, takes in either "fixed", "optimal", "integrate"
#' @param int_lower_bound
#' @param int_upper_bound
#'
#' @return Dataframe where the first |A| columns represent the sign vector and the last column is the joint probablity
#' @export
#'
#' @examples
Joint_prob_submodel_fixed <-function(A, F_cs, V_half, B_mag, sign_vects=NULL, log_lambda=NULL, log_lambda_strategy="integrate",
                                          sigma = 1, int_lower_bound=NULL, int_upper_bound=NULL, int_method="Quadrature", output_option = "joint"){
  #print(A)
  #browser()
  if(is.null(sign_vects)){
    sign_vects= gen_all_sign_vects(A)
  }
  colnames(sign_vects) <- paste("s", 1:length(A),sep="")
  sign_vects = as.data.frame(sign_vects)
  #browser()
  results <- apply(as.data.frame(sign_vects), 1, Joint_prob_submodel_sign_fixed, A=A, F_cs=F_cs, V_half=V_half, B_mag=B_mag,
                   log_lambda=log_lambda, log_lambda_strategy=log_lambda_strategy,
                   sigma = sigma, int_lower_bound=int_lower_bound, int_upper_bound=int_upper_bound, int_method=int_method, output_option=output_option)
  #browser()
  output = cbind(sign_vects, results)
  return(output)
}


### Now we move to the completely general case, looping over all (or a given subset of) submodels and sign vectors



gen_all_submodels<-function(p,k){
  return(t(combn(p,k)))
}

sample_submodels <-function(p, k, s_1){
  
  # this Uses NBIBD approach
  return(ibd(v=p, b= s_1, k=k)$design)
}

#' Title
#'
#' @param p number of main effect factors
#' @param k number of factors in active set
#' @param s_1 The number for one NBIBD design, recommend no more than 64
#' @param num_stacks Intieger, number of times the NBIBD design of size s_1 will be permuted and stacked to give the final sample submodel set
#'
#' @return
#' @export
#'
#' @examples
stack_NBIBD <-function(p,k,s_1,num_stacks){
  A_NBIBD = sample_submodels(p=p,k=k,s_1=s_1)
  A_final = A_NBIBD
  for(s in 1:(num_stacks-1)){
    A_perm = A_NBIBD
    mapping = sample(1:p,size =p)
    for (i in 1:s_1){
      for (j in 1:k){
        A_perm[i,j]= mapping[A_NBIBD[i,j]]
        
      }
    }
    A_final = rbind(A_final, A_perm)
  }
  
  return(list("A_final"= A_final, "A_NBIBD"=A_NBIBD))
}

Joint_prob_general_fixed_lambda <-function( F_cs, V_half, B_mag,k=NULL, submodels=NULL, sign_vects=NULL, log_lambda=1, log_lambda_strategy="fixed",
                               sigma = 1,  output_option = "joint"){
  p=dim(F_cs)[2]
  n = dim(F_cs)[1]
  if(any(is.na(F_cs))){
    return(c(0,0,0)) 
  }
  if(is.null(submodels)){
    if(is.null(k)){
      k = floor(n/3)
    }
    submodels= gen_all_submodels(p,k)
  }
  else{ k = dim(submodels)[2]}
  if(is.null(sign_vects)){
    num_sign_vects = 2^(k-1)
  }
  else{
    num_sign_vects = dim(sign_vects)[1]
  }
  #browser()
  colnames(submodels) <- paste("A", 1:k,sep="")
  submodels.expanded = submodels[rep(seq_len(nrow(submodels)), rep(num_sign_vects,nrow(submodels))), 1:k]
  #browser()
  results <- do.call(rbind,apply(submodels,1, Joint_prob_submodel_fixed,sign_vects=sign_vects, F_cs=F_cs, V_half=V_half, B_mag=B_mag,
                                 log_lambda=log_lambda, log_lambda_strategy=log_lambda_strategy,
                                 sigma = sigma, output_option=output_option))
  #browser()
  output = cbind(submodels.expanded, results)
  return(mean(output$results))
} 


Joint_prob_general_optimal_lambda <-function( F_cs, V_half, B_mag,k=NULL, submodels=NULL, sign_vects=NULL,
                                            sigma = 1,  output_option = "joint", log_lam_min=-2, log_lam_max= 2, stepsize=0.1){
  p=dim(F_cs)[2]
  n = dim(F_cs)[1]
  if(any(is.na(F_cs))){
    return(c(0,0,0)) 
  }
  if(is.null(submodels)){
    if(is.null(k)){
      k = floor(n/3)
    }
    submodels= gen_all_submodels(p,k)
  }
  else{ k = dim(submodels)[2]}
  if(is.null(sign_vects)){
    num_sign_vects = 2^(k-1)
    sign_vects = gen_all_sign_vects(seq(1,k, by =1))
  }
  else{
    num_sign_vects = dim(sign_vects)[1]
  }
  # #browser()
  # colnames(submodels) <- paste("A", 1:k,sep="")
  # submodels.expanded = submodels[rep(seq_len(nrow(submodels)), rep(num_sign_vects,nrow(submodels))), 1:k]
  # #browser()
  # results <- do.call(rbind,apply(submodels,1, Joint_prob_submodel_fixed,sign_vects=sign_vects, F_cs=F_cs, V_half=V_half, B_mag=B_mag,
  #                                log_lambda=log_lambda, log_lambda_strategy=log_lambda_strategy,
  #                                sigma = sigma, output_option=output_option))
  # #browser()
  # output = cbind(submodels.expanded, results)
  results = opt_log_lambda(log_lam_min = log_lam_min, log_lam_max = log_lam_max, F_cs, V_half, k, B_mag=B_mag, submodels, sign_vects,sigma,stepsize)
  return(results$opt_val)
}


#'This is a function that gives the joint probability, or the intergral across log lambda, for all submodels and sign vectors supplied by the user. 
#'
#'
#' @param F_cs Full centered and scaled design matrix
#' @param V_half Full sqrt of the scale matrix
#' @param B_mag Vector of magnitude of active effects, should be all positive
#' @param k Integer, size of submodels, only used if submodels=NULL, if k=NULL, it will be set to floor of n/3.
#' @param submodels Matrix, with k columns where each row represents a particular submodel column index, if NULL then it loops over all p choose k submodels. Defualt is NULL.
#' @param sign_vects Matrix, with k columns where each row is a sign vector, if NULL then it loops over all sign vectors. Defualt is NULL.
#' @param sigma Scalar, variance parameter, default is 1
#' @param log_lambda Scalar, fixed value of log_lambda, this is set to NULL, if the log_lambda strategy is fixed, this must be defined
#' @param log_lambda_strategy string, strategy to handle log_lambdas, takes in either "fixed", "optimal", "integrate"
#' @param int_lower_bound Scalar, log lambda lower bound of integration
#' @param int_upper_bound Scalar, log lambda upper bound of integration
#' @param int_method string, "Sum" or "Quadrature", Quadrature takes a long time
#' @param output_option string, must be one of "I","S",or "joint", specifies which event, or joint event the output is in terms of. 
#'
#' @return Dataframe where the first |A| columns represent the submodel index, the next |A| column represent the sign vector and the last column is the output_option probablity
#' @export
#'
#' @examples
Joint_prob_general <-function( F_cs, V_half, B_mag,k=NULL, submodels=NULL, sign_vects=NULL, log_lambda=NULL, log_lambda_strategy="integrate",
                                     sigma = 1, int_lower_bound=NULL, int_upper_bound=NULL, int_method="Sum", output_option = "joint"){
  p=dim(F_cs)[2]
  n = dim(F_cs)[1]
  if(any(is.na(F_cs))){
    return(c(0,0,0)) 
  }
  if(is.null(submodels)){
    if(is.null(k)){
      k = floor(n/3)
    }
    submodels= gen_all_submodels(p,k)
  }
  else{ k = dim(submodels)[2]}
  if(is.null(sign_vects)){
    num_sign_vects = 2^(k-1)
  }
  else{
    num_sign_vects = dim(sign_vects)[1]
  }
  #browser()
  colnames(submodels) <- paste("A", 1:k,sep="")
  submodels.expanded = submodels[rep(seq_len(nrow(submodels)), rep(num_sign_vects,nrow(submodels))), 1:k]
  #browser()
  results <- do.call(rbind,apply(submodels,1, Joint_prob_submodel_fixed,sign_vects=sign_vects, F_cs=F_cs, V_half=V_half, B_mag=B_mag,
                   log_lambda=log_lambda, log_lambda_strategy=log_lambda_strategy,
                   sigma = sigma, int_lower_bound=int_lower_bound, int_upper_bound=int_upper_bound, int_method=int_method, output_option=output_option))
  #browser()
  output = cbind(submodels.expanded, results)
  return(output)
}

 



### Measures of output of Joint_prob_general

#' Title
#'
#' @param results vector of the results for the set of submodels and sign vectors. This is the last column of the output from Joint_prob_general.
#'
#' @return List of several summary measures.
#' @export
#'
#' @examples
measure_list <-function(results){
  percent_0 = sum(results ==0)/length(results)
  percent_LT_0.1 = sum(results <0.1)/length(results)
  percent_GT_0.8 = sum(results >0.8)/length(results)
  summary_results = summary(results)
  percentile_25=as.numeric(summary_results[2])
  percentile_75=as.numeric(summary_results[5])
  mean_results = as.numeric(summary_results[4])
  median_results = as.numeric(summary_results[3])
  return(list("percent_0"= percent_0, "percent_LT_0.1"= percent_LT_0.1,
              "percent_GT_0.8"=percent_GT_0.8, "percentile_25"=percentile_25, 
              "percentile_75"=percentile_75, "mean"=mean_results,
              "median"=median_results))
}












Joint_prob_all_fixed_thresh <- function(A,I1, F_cs, V_half, B_mag, Z_A, Z_I1, log_lambda, sigma = 1, thresh = Inf){
  if(any(I1%in%A)){
    return(0)
  }
  else{
    A_hat = union(A, I1)
    
    #print(A_hat)
    # F_A=as.matrix(F_cs[,A])
    # F_I1 = as.matrix(F_cs[,I1])
    F_A_hat = as.matrix(F_cs[,A_hat])
    F_I2 = as.matrix(F_cs[,-A_hat])
    # V_half_A=diag(diag(V_half)[A])
    # V_half_I1=diag(diag(V_half)[I1])
    V_half_A_hat=diag(diag(V_half)[A_hat])
    
    lambda = exp(log_lambda)
    
    n = dim(F_cs)[1]
    k = length(Z_A)
    q1 = length(Z_I1)
    one_k = rep(1,k)
    one_q2 = rep(1, dim(F_I2)[2])
    B_min_vec_A = B_mag*Z_A
    B_min_vec_I1 = rep(0,q1)
    B_min_vec = c(B_min_vec_A,B_min_vec_I1)
    C_I2A_hat = (1/n)*t(F_I2)%*%F_A_hat
    C_I2I2 = (1/n)*t(F_I2)%*%F_I2
    C_AA_hat = (1/n)*t(F_A_hat)%*%F_A_hat
    z_A_hat = as.matrix(c(Z_A, as.matrix(Z_I1)))
    #browser()
    
   # print(z_A_hat)
    #browser()
    # Check if C_AA has an inverse, if not, use moore-penrose. (Discuss with group)
    C_AA_inv = n*ginv(t(F_A_hat)%*%F_A_hat)
    Mean_I = lambda*sqrt(n)*C_I2A_hat%*%C_AA_inv%*%as.matrix(z_A_hat)
    Cov_I = (sigma^2)*(C_I2I2 - (C_I2A_hat%*%C_AA_inv%*%t(C_I2A_hat)))
    
    # Fix the diagonals of Cov_I so that they are slightly positive, if they are slightly negative
    Cov_I = fix_neg_diag(Cov_I)
    P_I = mvtnorm::pmvnorm(lower = -sqrt(n)*lambda*one_q2, upper = sqrt(n)*lambda*one_q2,
                           mean = as.vector(Mean_I), sigma = as.matrix(Cov_I), algorithm = GenzBretz() )[1]
    #browser()
    Cov_S = (sigma^2)* diag(as.vector(z_A_hat))%*%C_AA_inv%*%diag(as.vector(z_A_hat))
    mean_S = sqrt(n) * diag(as.vector(z_A_hat)) %*%( (V_half_A_hat%*%B_min_vec) -(lambda* C_AA_inv%*%z_A_hat))
    # upper_S = sqrt(n)*diag(Z_A)%*%V_half%*%B_min_vec -
    #   lambda*sqrt(n)* diag(C_AA_inv%*%Z_A%*%t(Z_A))
    upper_S = c(rep(Inf,k), rep(thresh,q1))
    lower_S = rep(0, k+q1)
    #browser()
    P_S= mvtnorm:: pmvnorm(lower= as.vector(lower_S), upper = as.vector(upper_S), mean= as.vector(mean_S), sigma = Cov_S)[1]
    P_joint = P_S*P_I
    if(is.na(P_joint)){
      return(0)
    }
    else{return(P_joint)}
    }
}

# This function for a fixed I1 gets all possible A sets
get_possible_A<-function(I1,k, p){
  possible_A_cols = setdiff(seq(1,p,1), I1)
  submodels_A = t(combn(possible_A_cols,k))
  return(submodels_A)
}



# Function to generate all possible sign vects for I1, and A. We can't just use gen_all_sign_vects, because that assumes sign symmetry, which we don't have

gen_all_sign_vects_thresh<-function(k){
  #k = length(A)
  l=rep(list(c(-1,1)), as.integer(k))
  # expand.grid will give a df where each row is a unique sign vector. The first half of the rows are the negartives of the second half
  sign_vec_grid = expand.grid(l)
  return(sign_vec_grid)
}

# Now we fix I1 and Z_I1 and average over all possible A and Z_A

Joint_prob_avg_fixed_I1 <-function(I1,Z_I1, F_cs, V_half, B_mag, submodels_A, sign_vects_A, log_lambda, sigma = 1, thresh = Inf){
  
  k= dim(submodels_A)[2]
  
  sum = 0
  count = 0
  
  for(i in seq(1, dim(submodels_A)[1])){
    A = submodels_A[i,]
    #print(A)
    for(j in seq(1, dim(sign_vects_A)[1])){
      Z_A = sign_vects_A[j,]
      P_joint = Joint_prob_all_fixed_thresh(A,I1, F_cs, V_half, B_mag, Z_A, Z_I1, log_lambda, sigma, thresh)
      #print(P_joint)
      sum = sum + P_joint
      count=count+1
    }
  }
  #print(count)
  return(sum/count)
  
  
}




# Now we set the numbner of false positives to include and get the probability of getting EXACTLY q1 false Positives

Joint_prob_q1_FP <-function(k,q1, F_cs, V_half, B_mag, sign_vects_A, log_lambda, sigma = 1, thresh = Inf){
  
  p= dim(F_cs)[2]
  #get all possible I1 sets
  I1_set= gen_all_submodels(p,q1)
  submodels_A = gen_all_submodels(p,k)
  sign_vects_I1 = gen_all_sign_vects_thresh(q1)
  sum = 0
  
  for(i in seq(1, dim(I1_set)[1])){
    I1 = I1_set[i,]
    #print(I1)
    #submodels_A = get_possible_A(I1, k=k, p=p)
    
    for(j in seq(1, dim(sign_vects_I1)[1])){
      Z_I1 = sign_vects_I1[j,]
      P_joint = Joint_prob_avg_fixed_I1(I1,Z_I1, F_cs, V_half, B_mag, submodels_A, sign_vects_A, log_lambda, sigma, thresh)
      #print(P_joint)
      sum = sum + P_joint
    }
  }
  
  return(sum)
  
  
}
  
  
  
# Now we set a max number of FPs to include and sum them up

Joint_prob_gen_thresh <-function(k,q1_max, F_cs, V_half, B_mag, sign_vects_A, log_lambda, sigma = 1, thresh = Inf){
  # first, the q1=0 case
  result_df = data.frame(q1=numeric(0), joint=numeric(0))
  sum = measure_list(Joint_prob_general(F_cs, V_half, B_mag, k, sign_vects = sign_vects_A, log_lambda, submodels = NULL, log_lambda_strategy = "fixed", sigma=sigma)$results)$mean
  result_df[1,]<-c(0, sum)
  if(q1_max >0){
    for(q1 in (1:q1_max)){
      #print(q1)
      P_joint_q1 = Joint_prob_q1_FP(k,q1, F_cs, V_half, B_mag, sign_vects_A, log_lambda, sigma, thresh)
      sum = sum + P_joint_q1
      print(P_joint_q1)
      result_df[q1+1,]<-c(q1, sum)
    }
  }
  
  return(result_df)
  
  
}


# Fix A and Z_A and loop over all I1 and Z_I1

Joint_prob_fixed_A_thresh <-function(A,I1_set, F_cs, V_half, B_mag, Z_A, Z_I1_set, log_lambda, sigma = 1, thresh = Inf){
  
  
  sum = 0
  
  
  for(i in seq(1, dim(I1_set)[1])){
    I1 = I1_set[i,]
    #print(I1)
    for(j in seq(1, dim(Z_I1_set)[1])){
      Z_I1 = Z_I1_set[j,]
      P_joint = Joint_prob_all_fixed_thresh(A,I1, F_cs, V_half, B_mag, Z_A, Z_I1, log_lambda, sigma, thresh)
      #print(P_joint)
      sum = sum + P_joint
      #count=count+1
    }
  }
  
  return(sum)
}












#' 
#' #' Title
#' #'
#' #' @param all_model_results dataframe output of the Joint_prob_general function for all submodels and a set of sign vectors
#' #' @param step_size How many submodels to skip between sampling. If overall design has 330 submodels, and a step size of 10 is chosen, this will sample 33
#' #' @param k_star Integer: size of support
#' #'
#' #' @return Matrix where each row is a submodel to sample
#' #' @export
#' #'
#' #' @examples
#' dist_sample <-function(all_model_results,step_size = 10, k_star=4){
#'   
#'   all_model_results_ave = all_model_results %>% 
#'     group_by(across(c(1:k_star))) %>% 
#'     summarise(mean_results=mean(results), .groups ="keep")
#'   
#'   #browser()
#'   sorted_results = sort(all_model_results_ave$mean_results)
#'   sample_result = sorted_results[seq(1, length(sorted_results), step_size)]
#'   
#'   
#'   zero_index = 0
#'   index_vec = c()
#'   for (val in sample_result){
#'     if (val == 0){
#'       zero_index = zero_index + 1
#'       index = which(val == all_model_results_ave$mean_results)
#'       index_vec = c(index[zero_index], index_vec)
#'       #browser()
#'     }
#'     else{
#'       index = which(val == all_model_results_ave$mean_results)
#'       index_vec = c(index_vec, index[1])
#'       #browser()
#'     }
#'     
#'     
#'   }
#'   
#'   #print(index_vec)
#'   return(as.matrix(all_model_results_ave[index_vec,c(1:k_star)]))
#'   
#'   
#' }
### Hueristics and Hueristic design construction

UES <- function(F_0){
  p = dim(F_0)[2]
  n = dim(F_0)[1]
  F_0_w_int =as.data.frame(F_0)
  F_0_w_int$int = 1
  
  FTF= t(as.matrix(F_0_w_int))%*%as.matrix(F_0_w_int)
  FTF_upper <- FTF[upper.tri(FTF, diag = FALSE)]
  UES <- (2/(p*(p+1))) * sum(FTF_upper)
  return(UES)
}

UES_squared <- function(F_0){
  p = dim(F_0)[2]
  n = dim(F_0)[1]
  F_0_w_int =as.data.frame(F_0)
  F_0_w_int$int = 1
  
  FTF= t(as.matrix(F_0_w_int))%*%as.matrix(F_0_w_int)
  FTF_upper <- FTF[upper.tri(FTF, diag = FALSE)]
  UES_sq = (2/(p*(p+1))) * sum(FTF_upper^2)
  return(UES_sq)
}

Var_S_Plus <- function(F_0, UES_sq_star, eff = 0.8, UES_const = 0){
  p = dim(F_0)[2]
  n = dim(F_0)[1]
  F_0_w_int =as.data.frame(F_0)
  F_0_w_int$int = 1
  
  FTF= t(as.matrix(F_0_w_int))%*%as.matrix(F_0_w_int)
  FTF_upper <- FTF[upper.tri(FTF, diag = FALSE)]
  UES_sq = (2/(p*(p+1))) * sum(FTF_upper^2)
  UES <- (2/(p*(p+1))) * sum(FTF_upper)
  #browser()
  if(UES >UES_const & (UES_sq_star/UES_sq) > eff){
    valid_result <- TRUE
  }
  else{ valid_result = FALSE}
  Var_s = UES_sq - (UES^2)
  return(list("Var_s"=Var_s, "valid"= valid_result, "UES"=UES))
}


UES_sq_construction<-function(n,main_effects){
  p = main_effects+1
  if(p%%4 ==0){
    H_p = Normcol(Hadamard_Matrix(p))
    X_0 =  H_p[sample(p, size = n, replace = FALSE),]
  }
  if (p%%4 ==1){
    H = Normcol(Hadamard_Matrix(p-1))
    V = H[sample(p-1, size = n, replace = FALSE),]
    phi = sample(c(-1,1), n, replace=TRUE)
    X_0 = as.matrix(cbind(V, phi))
    
  }
  if (p%%4 ==2){
    if(n%%2==0){
      m = n/2
      H = Normcol(Hadamard_Matrix(p-2))
      X_star = H[sample(p-2, size = n, replace = FALSE),]
      G = matrix(c(1,1,-1,-1,-1,1,1,-1), nrow = 4, ncol=2, byrow = TRUE)
      U_11 = G[sample(2,m, replace = TRUE),]
      U_12 = G[sample(3:4,m, replace = TRUE),]
      U_1 = rbind(U_11, U_12)
      X_0 = as.matrix(cbind(X_star, U_1))
    }
    if(n%%2!=0){
      m = (n-1)/2
      H = Normcol(Hadamard_Matrix(p-2))
      X_star = H[sample(p-2, size = n, replace = FALSE),]
      G = matrix(c(1,1,-1,-1,-1,1,1,-1), nrow = 4, ncol=2, byrow = TRUE)
      U_21 = G[sample(2,m, replace = TRUE),]
      U_22 = G[sample(3:4,m+1, replace = TRUE),]
      U_2 = rbind(U_21, U_22)
      X_0 = as.matrix(cbind(X_star, U_2))
       }
    
  }
  if (p%%4 ==3){
    H = t(Hadamard_Matrix(p+1))
    X_star = H[sample(p+1, size = n, replace = FALSE),]
    X_0 = X_star[,1:p]
  }
  return(X_0)
}

min_UES_sq_val <- function(n, main_effects){
  p = main_effects+1
  if(p%%2 ==1){
    min_UES_sq = ((n*(n-1)) + (n*p*(p-n)))/(p*(p-1))
  }
  if (p%%4 ==0){
    min_UES_sq = (n*(p-n))/(p-1)
    
  }
  if (p%%4 ==2){
    if(n%%2==0){
      min_UES_sq = ((2*n*(n-2)) + (n*p*(p-n)))/(p*(p-1))
    }
    if(n%%2!=0){
      min_UES_sq = ((2*((n-1)**2)) + (n*p*(p-n)))/(p*(p-1))
    }
    
  }
  return(min_UES_sq)
}

coord_exchange_UES_sq<-function( F_init, max_iter=10000){
  
  F_0 = F_init
  n = dim(F_0)[1]
  p = dim(F_0)[2]
  
  
  P_0 =  as.numeric(UES_squared(F_0))
  
  
  iterations = 0
  # Intitialize the count between finding a design with a higher measure.
  count_between_flips=0
  # exchange_data = data.frame( exchange_ind = numeric(0), sum_correl=numeric(0), max_correl=numeric(0), 
  #                             min_correl= numeric(0), improvement=numeric(0))
  stop = FALSE
  while ( iterations <= max_iter ){
    #print(paste0("Iteration: ", iterations))
    for (i in c(1:n)){
      for (j in c(1:p)){
        # print(i)
        # print(j)
        #print(paste0("Count Between Flips: ", count_between_flips))
        F_1 = F_0
        # Change Coordinate and update scaling
        F_1[i,j]= -F_0[i,j]
        
        P_1 = as.numeric(UES_squared(F_1))
        if( P_1 < P_0) {
          #Capture the sum and max absolute column wise correlations of the factor that is flipped
          #browser()
          # abs_correl_mat = (abs(t(as.matrix(F_0_cs))%*%as.matrix(F_0_cs))/n)- diag(1, nrow=p, ncol=p)
          # sum_correl = rowSums(abs_correl_mat)[j]
          # max_correl = max(abs_correl_mat[j,])
          # min_correl = min(abs_correl_mat[j,])
          # imp= P_1-P_0
          # exchange_data[nrow(exchange_data)+1,] <- c(1, sum_correl,max_correl , min_correl, imp)
          #browser()
          F_0 = F_1
          P_0 = P_1
          # We just flipped designs so now we reset count
          count_between_flips = 0
          
        }
        else{
          #Capture the sum and max absolute column wise correlations of the factor that is NOT flipped
          count_between_flips = count_between_flips+1
        }
        if(count_between_flips> n*p){
          # If we have looped through the whole matrix without changing designs break. 
          stop = TRUE
          break
        }
        #Breaking out of outter loops
        
      }
      if(stop){break}
      
    }
    if(stop){break}
    iterations=iterations+1
  }
  
  return(list("design"= F_0, "P_0"=P_0))
}

coord_exchange_Var_s<-function( n, p, eff = 0.8,UES_const = 0, max_iter=10000){
  UES_sq_standard = coord_exchange_UES_sq(matrix(sample(c(-1,1), n*p, TRUE), nrow = n, ncol = p), max_iter = 1000)
  UES_sq_star = UES_sq_standard$P_0
  #print(UES_sq_star)
  
  
  invalid = TRUE
  # First exchange to get a valid var(s +) design
  while(invalid){
    # Indicator for negative average off diagonals
    neg_off <- TRUE
    #browser()
    # Start with random matrices until we get one with UES >0
    while( neg_off){
      
      F_init = matrix(sample(c(-1,1), n*p, TRUE), nrow = n, ncol = p)
      
      if (UES(F_init)>UES_const){
        neg_off = FALSE
        #print("found positive matrix")
      }
      else{ neg_off = TRUE}
    }
    F_0 = F_init
    # Now we exchange to get up to the efficiency limit
    eff_0 = UES_sq_star / UES_squared(F_0)
    row = 1
    col = 1
    count = 0
    overall_count = 1000
    #browser()
    while( eff_0 < eff){
      F_1= F_0
      #exchange coordinates
      F_1[row,col] = -F_1[row,col]
      eff_1 = UES_sq_star / UES_squared(F_1)
      if(eff_1>eff_0){
        if( UES(F_1)>UES_const){
          F_0 = F_1
          eff_0 = eff_1
          count = 0
        }
      }
      row = (row + 1)%%n
      col = (col+ 1)%%p
      count = count +1
      if(count >n*p){
        #print("Not efficient enough, trying again.")
        #browser()
        break}
    }
    if(eff_0 >= eff){
      invalid = FALSE
    }
    
  }
  
  P_0 = as.numeric(Var_S_Plus(F_0, UES_sq_star, eff,UES_const = UES_const)$Var_s)
  
  #browser()
  
  
  iterations = 0
  # Intitialize the count between finding a design with a higher measure.
  count_between_flips=0
  # exchange_data = data.frame( exchange_ind = numeric(0), sum_correl=numeric(0), max_correl=numeric(0), 
  #                             min_correl= numeric(0), improvement=numeric(0))
  stop = FALSE
  while ( iterations <= max_iter ){
    #print(paste0("Iteration: ", iterations))
    for (i in c(1:n)){
      for (j in c(1:p)){
        # print(i)
        # print(j)
        #print(paste0("Count Between Flips: ", count_between_flips))
        F_1 = F_0
        # Change Coordinate and update scaling
        F_1[i,j]= -F_0[i,j]
        
        Var_s_results = Var_S_Plus(F_1, UES_sq_star, eff)
        if(Var_s_results$valid){
          P_1 = Var_s_results$Var_s
        }
        else{P_1 = P_0 +0.5}
        if( P_1 <= P_0) {
          #Capture the sum and max absolute column wise correlations of the factor that is flipped
          #browser()
          # abs_correl_mat = (abs(t(as.matrix(F_0_cs))%*%as.matrix(F_0_cs))/n)- diag(1, nrow=p, ncol=p)
          # sum_correl = rowSums(abs_correl_mat)[j]
          # max_correl = max(abs_correl_mat[j,])
          # min_correl = min(abs_correl_mat[j,])
          # imp= P_1-P_0
          # exchange_data[nrow(exchange_data)+1,] <- c(1, sum_correl,max_correl , min_correl, imp)
          #browser()
          F_0 = F_1
          P_0 = P_1
          # We just flipped designs so now we reset count
          count_between_flips = 0
          
        }
        else{
          #Capture the sum and max absolute column wise correlations of the factor that is NOT flipped
          count_between_flips = count_between_flips+1
        }
        if(count_between_flips> n*p){
          # If we have looped through the whole matrix without changing designs break. 
          stop = TRUE
          break
        }
        #Breaking out of outter loops
        
      }
      if(stop){break}
      
    }
    if(stop){break}
    iterations=iterations+1
  }
  #print(F_init - F_0)
  return(list("design"= F_0, "P_0"=P_0, "UES"= UES(F_0)))
}


MSE_c <- function(F_0, c_opt){
  
  n=dim(F_0)[1]
  k=dim(F_0)[2]
  
  F_0_cs = as.data.frame(lapply(as.data.frame(F_0), cent_scale))
  C = (1/n)*t(as.matrix(F_0_cs))%*%as.matrix(F_0_cs)
  
  MSE =mean((C[upper.tri(C)]-c_opt)^2)
  return(MSE)
  
}
#'This is a function  for the coordinate exchange for one initial matrix. This will be wrapped in parallel across a number of random starts in a subsequent function
#'
#'
#' @param F_cs_init Full centered and scaled design matrix
#' @param V_half Full sqrt of the scale matrix
#' @param B_mag Vector of magnitude of active effects, should be all positive
#' @param k Integer, size of submodels, only used if submodels=NULL, if k=NULL, it will be set to floor of n/3.
#' @param submodels Matrix, with k columns where each row represents a particular submodel column index, if NULL then it loops over all p choose k submodels. Defualt is NULL.
#' @param sign_vects Matrix, with k columns where each row is a sign vector, if NULL then it loops over all sign vectors. Defualt is NULL.
#' @param sigma Scalar, variance parameter, default is 1
#' @param log_lambda Scalar, fixed value of log_lambda, this is set to NULL, if the log_lambda strategy is fixed, this must be defined
#' @param log_lambda_strategy string, strategy to handle log_lambdas, takes in either "fixed", "optimal", "integrate"
#' @param int_lower_bound
#' @param int_upper_bound
#'
#' @return Dataframe where the first |A| columns represent the sign vector and the last column is the joint probablity
#' @export
#'
#' @examples
coord_exchange_1_run <-function( F_init, B_mag,k=NULL, submodels=NULL, sign_vects=NULL, log_lambda=NULL, log_lambda_strategy="integrate",
                               sigma = 1, int_lower_bound=NULL, int_upper_bound=NULL, int_method="Sum", measure_decsision= "mean", max_iter=10){
  
  F_0 = F_init
  n = dim(F_0)[1]
  p = dim(F_0)[2]
  #center and scale
  V_0_half= diag(lapply(as.data.frame(F_0),fn))
  F_0_cs = as.data.frame(lapply(as.data.frame(F_0), cent_scale))
  #browser()
  # This gives a list of several design measures in a named list
  measure_list_0 = measure_list(Joint_prob_general(F_cs = F_0_cs, V_half=V_0_half, 
                                                   B_mag=B_mag,
                                                   sign_vects=sign_vects,
                                                   log_lambda = log_lambda,
                                                   submodels = submodels,
                                                   log_lambda_strategy=log_lambda_strategy,
                                                   sigma = sigma, int_lower_bound =int_lower_bound, int_upper_bound = int_upper_bound)$results)
  # Choose which of these measures to make a decision on. 
  P_0 = as.numeric(measure_list_0[measure_decsision])
  
  
  iterations = 0
  # Intitialize the count between finding a design with a higher measure.
  count_between_flips=0
  exchange_data = data.frame( exchange_ind = numeric(0), sum_correl=numeric(0), max_correl=numeric(0), 
                              min_correl= numeric(0), improvement=numeric(0))
  stop = FALSE
  while ( iterations <= max_iter ){
    #print(paste0("Iteration: ", iterations))
    for (i in c(1:n)){
      for (j in c(1:p)){
        print(i)
        print(j)
        #print(paste0("Count Between Flips: ", count_between_flips))
        F_1 = F_0
        # Change Coordinate and update scaling
        F_1[i,j]= -F_0[i,j]
        V_1_half = V_0_half
        V_1_half[j,j]<- fn(F_1[,j])
        F_1_cs = F_0_cs
        F_1_cs[,j]= cent_scale(F_1[,j])
        measure_list_1 = measure_list(Joint_prob_general(F_cs = F_1_cs, V_half=V_1_half, 
                                                         B_mag=B_mag,
                                                         sign_vects=sign_vects,
                                                         log_lambda = log_lambda,
                                                         submodels = submodels,
                                                         log_lambda_strategy=log_lambda_strategy,
                                                         sigma = sigma, int_lower_bound =int_lower_bound, int_upper_bound = int_upper_bound)$results)
        # Choose which of these measures to make a decision on. 
        #browser()
        P_1 = as.numeric(measure_list_1[measure_decsision])
        if( P_1 > P_0) {
          #Capture the sum and max absolute column wise correlations of the factor that is flipped
          #browser()
          abs_correl_mat = (abs(t(as.matrix(F_0_cs))%*%as.matrix(F_0_cs))/n)- diag(1, nrow=p, ncol=p)
          sum_correl = rowSums(abs_correl_mat)[j]
          max_correl = max(abs_correl_mat[j,])
          min_correl = min(abs_correl_mat[j,])
          imp= P_1-P_0
          exchange_data[nrow(exchange_data)+1,] <- c(1, sum_correl,max_correl , min_correl, imp)
          #browser()
          F_0 = F_1
          V_0_half =V_1_half
          F_0_cs=F_1_cs
          measure_list_0=measure_list_1
          P_0 = P_1
          # We just flipped designs so now we reset count
          count_between_flips = 0
          
        }
        else{
          #Capture the sum and max absolute column wise correlations of the factor that is NOT flipped
          abs_correl_mat = (abs(t(as.matrix(F_0_cs))%*%as.matrix(F_0_cs))/n)- diag(1, nrow=p, ncol=p)
          sum_correl = rowSums(abs_correl_mat)[j]
          max_correl = max(abs_correl_mat[j,])
          min_correl = min(abs_correl_mat[j,])
          
          exchange_data[nrow(exchange_data)+1,] <- c(0, sum_correl,max_correl , min_correl, 0)
          #browser()
          count_between_flips = count_between_flips+1
        }
        if(count_between_flips> n*p){
          # If we have looped through the whole matrix without changing designs break. 
          stop = TRUE
          break
        }
        #Breaking out of outter loops
        
      }
      if(stop){break}
      
    }
    if(stop){break}
    iterations=iterations+1
  }
  
  return(list("design"= F_0, "measure_list"= measure_list_0, "P_0"=P_0, "exchange_df"=exchange_data))
}


# Now we need to make this run in parallel across different random starts.





### Column Exchange Algorithm 
get_max_corr<-function(F_cs, taboo_index){
  FTF_abs= abs(t(as.matrix(F_cs))%*%as.matrix(F_cs))
  # Get the upper triangle of the information matrix in absolute value
  FTF_abs[lower.tri(FTF_abs, diag=TRUE)]<-0
  #max_corr_cols <- which(FTF_abs ==max(FTF_abs), arr.ind=TRUE)
  corr_order = unique(sort(FTF_abs, decreasing=TRUE))
  
  for (corr in corr_order){
    #browser()
    corr_indices = which(FTF_abs==corr, arr.ind=TRUE)
    for (index in seq(1:nrow(corr_indices))){
      #browser()
      if(Position(function(x) identical(x, corr_indices[index,]), taboo_index, nomatch=0) > 0){
        print("Found a taboo index")
      }
      else{
        return(corr_indices[index,])
      }
    }
  }

}

update_tab_list <- function(corr_cols, tab_list){
  tab_list$tab1 <- tab_list$tab2
  tab_list$tab2 <- tab_list$tab3
  tab_list$tab3 <- corr_cols
  
  storage.mode(tab_list$tab1)="integer"
  storage.mode(tab_list$tab2)="integer"
  storage.mode(tab_list$tab3)="integer"
  return(tab_list)
}


#' Title
#'
#' @param F_0 
#' @param B_mag 
#' @param k 
#' @param submodels 
#' @param sign_vects 
#' @param log_lambda 
#' @param log_lambda_strategy 
#' @param sigma 
#' @param int_lower_bound 
#' @param int_upper_bound 
#' @param int_method 
#' @param measure_decsision 
#' @param max_iter 
#' @param tol 
#'
#' @return
#' @export
#'
#' @examples
pairwise_exchange_1_run <-function( F_0, B_mag,k=NULL, submodels=NULL, sign_vects=NULL, log_lambda=NULL, log_lambda_strategy="integrate",
                                    sigma = 1, int_lower_bound=NULL, int_upper_bound=NULL, int_method="Sum", measure_decsision= "mean", 
                                    max_iter=10, tol=0.0005, max_lookback=1, early_stopping = 0.01){
  
  # F_0 = F_init
  n = dim(F_0)[1]
  p = dim(F_0)[2]
  #center and scale
  V_0_half= diag(lapply(as.data.frame(F_0),fn))
  F_0_cs = as.data.frame(lapply(as.data.frame(F_0), cent_scale))
  #browser()
  # This gives a list of several design measures in a named list
  measure_list_0 = measure_list(Joint_prob_general(F_cs = F_0_cs, V_half=V_0_half, 
                                                   B_mag=B_mag,k=k,
                                                   sign_vects=sign_vects,
                                                   log_lambda = log_lambda,
                                                   submodels = submodels,
                                                   log_lambda_strategy=log_lambda_strategy,
                                                   sigma = sigma, int_lower_bound =int_lower_bound,int_method=int_method, int_upper_bound = int_upper_bound)$results)
  # Choose which of these measures to make a decision on. 
  P_0 = as.numeric(measure_list_0[measure_decsision])
  
  
  iterations = 0
  lookback = 0
  #tracking_data = data.frame(col_ind_1= numeric(0), col_index_2 = numeric(0), P_prev=numeric(0), P_new=numeric(0))
  #stop = FALSE
  #taboo_index = NULL
  tab_list = list("tab1"= NULL, "tab2"=NULL, "tab3"=NULL)
  #print(P_0)
  P_prev=P_0
  P_new=P_0
  P_early_stopping = P_0
  #browser()
  while ( iterations <= max_iter ){
    if( P_new - P_prev < tol){
      lookback = lookback +1
    }
    if( lookback > max_lookback){
      break
    }
    else{
      P_prev=P_new
      lookback = 0 
    }
    if(iterations !=0 &&iterations %% 5 ==0){
      if(P_prev - P_early_stopping < early_stopping){
          break
      }
    }
    
    
    max_corr_cols = get_max_corr(F_0_cs, taboo_index= tab_list)
    #print(max_corr_cols)
    tab_list = update_tab_list(max_corr_cols, tab_list) 
    #browser()
    #print(tab_list)
    #Subset the design to the most correlated columns
    for( i in c(1:n)){
      #print(i)
      # Generate three new designs based on the different pairwise combinations of one row of 2 cols of -1 and 1\
      # This might be a very sloppy way of doing this
      F_1=F_0
      F_2=F_0
      F_3=F_0
      #browser()
      F_1[i,max_corr_cols]<- F_0[i,max_corr_cols]*c(-1,1)
      F_2[i,max_corr_cols]<- F_0[i,max_corr_cols]*c(1,-1)
      F_3[i,max_corr_cols]<- F_0[i,max_corr_cols]*c(-1,-1)
      
      V_1_half= diag(lapply(as.data.frame(F_1),fn))
      F_1_cs = as.data.frame(lapply(as.data.frame(F_1), cent_scale))
      V_2_half= diag(lapply(as.data.frame(F_2),fn))
      F_2_cs = as.data.frame(lapply(as.data.frame(F_2), cent_scale))
      V_3_half= diag(lapply(as.data.frame(F_3),fn))
      F_3_cs = as.data.frame(lapply(as.data.frame(F_3), cent_scale))
      
      # Calculate the measure list for all designs
      measure_list_1 = measure_list(Joint_prob_general(F_cs = F_1_cs, V_half=V_1_half, 
                                                       B_mag=B_mag,k=k,
                                                       sign_vects=sign_vects,
                                                       log_lambda = log_lambda,
                                                       submodels = submodels,
                                                       log_lambda_strategy=log_lambda_strategy,
                                                       sigma = sigma, int_lower_bound =int_lower_bound, 
                                                       int_method=int_method, int_upper_bound = int_upper_bound)$results)
      measure_list_2 = measure_list(Joint_prob_general(F_cs = F_2_cs, V_half=V_2_half, 
                                                       B_mag=B_mag,k=k,
                                                       sign_vects=sign_vects,
                                                       log_lambda = log_lambda,
                                                       submodels = submodels,
                                                       log_lambda_strategy=log_lambda_strategy,
                                                       sigma = sigma, int_lower_bound =int_lower_bound, int_method=int_method,int_upper_bound = int_upper_bound)$results)
      measure_list_3 = measure_list(Joint_prob_general(F_cs = F_3_cs, V_half=V_3_half, 
                                                       B_mag=B_mag,k=k,
                                                       sign_vects=sign_vects,
                                                       log_lambda = log_lambda,
                                                       submodels = submodels,
                                                       log_lambda_strategy=log_lambda_strategy,
                                                       sigma = sigma, int_lower_bound =int_lower_bound,int_method=int_method, int_upper_bound = int_upper_bound)$results)
      # Choose which of these measures to make a decision on. 
      #browser()
      P_1 = as.numeric(measure_list_1[measure_decsision])
      P_2 = as.numeric(measure_list_2[measure_decsision])
      P_3 = as.numeric(measure_list_3[measure_decsision])
      result_table <-tribble(~measure, ~design,~design_cs,~V_half,
                             P_0, F_0,F_0_cs, V_0_half,
                             P_1, F_1, F_1_cs, V_1_half,
                             P_2, F_2, F_2_cs, V_2_half,
                             P_3, F_3, F_3_cs, V_3_half)
      max_index = which(result_table$measure==max(result_table$measure), arr.ind=TRUE)
      #Update the optimal design
      #browser()
      F_0 = matrix(unlist(result_table[max_index,]$design), nrow = n, ncol = p)
      F_0_cs = matrix(unlist(result_table[max_index,]$design_cs), nrow = n, ncol = p)
      V_0_half= matrix(unlist(result_table[max_index,]$V_half), nrow = p, ncol = p)
      P_new = max(result_table$measure)
      #print(P_new)
      P_0=P_new
      #tracking_data[nrow(tracking_data)+1,]<-c(max_corr_cols[1], max_corr_cols[2], P_prev, P_new)
      #browser()
    }
    #tracking_data[nrow(tracking_data)+1,]<-c(max_corr_cols[1], max_corr_cols[2], P_prev, P_new)
    iterations =iterations +1
  }
  
  
  return(list("design"= F_0,"P_0"=P_new, "final_iteration"= iterations))
}


pairwise_exchange_1_run_more_output <-function( F_0, B_mag,k=NULL, submodels=NULL, sign_vects=NULL, log_lambda=NULL, log_lambda_strategy="integrate",
                                    sigma = 1, int_lower_bound=NULL, int_upper_bound=NULL, int_method="Sum", measure_decsision= "mean", max_iter=10, tol=0.0005){
  
  # F_0 = F_init
  n = dim(F_0)[1]
  p = dim(F_0)[2]
  #center and scale
  V_0_half= diag(lapply(as.data.frame(F_0),fn))
  F_0_cs = as.data.frame(lapply(as.data.frame(F_0), cent_scale))
  #browser()
  # This gives a list of several design measures in a named list
  measure_list_0 = measure_list(Joint_prob_general(F_cs = F_0_cs, V_half=V_0_half, 
                                                   B_mag=B_mag,k=k,
                                                   sign_vects=sign_vects,
                                                   log_lambda = log_lambda,
                                                   submodels = submodels,
                                                   log_lambda_strategy=log_lambda_strategy,
                                                   sigma = sigma, int_lower_bound =int_lower_bound, int_upper_bound = int_upper_bound)$results)
  # Choose which of these measures to make a decision on. 
  P_0 = as.numeric(measure_list_0[measure_decsision])
  
  
  iterations = 0
  
  #tracking_data = data.frame(col_ind_1= numeric(0), col_index_2 = numeric(0), P_prev=numeric(0), P_new=numeric(0))
  #stop = FALSE
  tab_list = list("tab1"= NULL, "tab2"=NULL, "tab3"=NULL)
  P_prev=P_0
  P_new=P_0
  #browser()
  while ( iterations <= max_iter ){
    print(iterations)
    print(P_new-P_prev)
    if( iterations != 0 && P_new-P_prev <tol){
      break
    }
    else{
      P_prev=P_new
    }
    
    
    max_corr_cols = get_max_corr(F_0_cs, taboo_index=tab_list)
    #print(max_corr_cols)
    taboo_index = max_corr_cols
    #Subset the design to the most correlated columns
    for( i in c(1:n)){
      print(i)
      # Generate three new designs based on the different pairwise combinations of one row of 2 cols of -1 and 1\
      # This might be a very sloppy way of doing this
      F_1=F_0
      F_2=F_0
      F_3=F_0
      #browser()
      F_1[i,max_corr_cols]<- F_0[i,max_corr_cols]*c(-1,1)
      F_2[i,max_corr_cols]<- F_0[i,max_corr_cols]*c(1,-1)
      F_3[i,max_corr_cols]<- F_0[i,max_corr_cols]*c(-1,-1)
      
      V_1_half= diag(lapply(as.data.frame(F_1),fn))
      F_1_cs = as.data.frame(lapply(as.data.frame(F_1), cent_scale))
      print(F_1)
      print(F_1_cs)
      V_2_half= diag(lapply(as.data.frame(F_2),fn))
      F_2_cs = as.data.frame(lapply(as.data.frame(F_2), cent_scale))
      print(F_2)
      print(F_2_cs)
      V_3_half= diag(lapply(as.data.frame(F_3),fn))
      F_3_cs = as.data.frame(lapply(as.data.frame(F_3), cent_scale))
      print(F_3)
      print(F_3_cs)
      
      # Calculate the measure list for all designs
      measure_list_1 = measure_list(Joint_prob_general(F_cs = F_1_cs, V_half=V_1_half, 
                                                       B_mag=B_mag,k=k,
                                                       sign_vects=sign_vects,
                                                       log_lambda = log_lambda,
                                                       submodels = submodels,
                                                       log_lambda_strategy=log_lambda_strategy,
                                                       sigma = sigma, int_lower_bound =int_lower_bound, int_upper_bound = int_upper_bound, int_method = int_method)$results)
      measure_list_2 = measure_list(Joint_prob_general(F_cs = F_2_cs, V_half=V_2_half, 
                                                       B_mag=B_mag,k=k,
                                                       sign_vects=sign_vects,
                                                       log_lambda = log_lambda,
                                                       submodels = submodels,
                                                       log_lambda_strategy=log_lambda_strategy,
                                                       sigma = sigma, int_lower_bound =int_lower_bound, int_upper_bound = int_upper_bound)$results)
      measure_list_3 = measure_list(Joint_prob_general(F_cs = F_3_cs, V_half=V_3_half, 
                                                       B_mag=B_mag,k=k,
                                                       sign_vects=sign_vects,
                                                       log_lambda = log_lambda,
                                                       submodels = submodels,
                                                       log_lambda_strategy=log_lambda_strategy,
                                                       sigma = sigma, int_lower_bound =int_lower_bound, int_upper_bound = int_upper_bound)$results)
      # Choose which of these measures to make a decision on. 
      #browser()
      P_1 = as.numeric(measure_list_1[measure_decsision])
      P_2 = as.numeric(measure_list_2[measure_decsision])
      P_3 = as.numeric(measure_list_3[measure_decsision])
      result_table <-tribble(~measure, ~design,~design_cs,~V_half,
                             P_0, F_0,F_0_cs, V_0_half,
                             P_1, F_1, F_1_cs, V_1_half,
                             P_2, F_2, F_2_cs, V_2_half,
                             P_3, F_3, F_3_cs, V_3_half)
      max_index = which(result_table$measure==max(result_table$measure), arr.ind=TRUE)
      #Update the optimal design
      #browser()
      F_0 = matrix(unlist(result_table[max_index,]$design), nrow = n, ncol = p)
      F_0_cs = matrix(unlist(result_table[max_index,]$design_cs), nrow = n, ncol = p)
      V_0_half= matrix(unlist(result_table[max_index,]$V_half), nrow = p, ncol = p)
      P_new = max(result_table$measure)
      P_0=P_new
      #tracking_data[nrow(tracking_data)+1,]<-c(max_corr_cols[1], max_corr_cols[2], P_prev, P_new)
      #browser()
    }
    #tracking_data[nrow(tracking_data)+1,]<-c(max_corr_cols[1], max_corr_cols[2], P_prev, P_new)
    iterations =iterations +1
  }
  
  
  return(list("design"= F_0,"P_0"=P_new))
}


pairwise_exchange_parallel <-function(start_design_list, cores, B_mag,k=NULL, submodels=NULL, sign_vects=NULL, log_lambda=NULL, log_lambda_strategy="integrate",
                                    sigma = 1, int_lower_bound=NULL, int_upper_bound=NULL, int_method="Sum", measure_decsision= "mean", max_iter=100, tol=0.0005, early_stopping = 0.0051){
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  print("Starting Parallel For Loop")
  
  final_results = foreach(start=1:dim(start_design_list)[3], .combine = rbind, .packages =c('MASS',"mvtnorm","dplyr"), .export = ls(globalenv())) %dopar%{
    pairwise_1_step <- pairwise_exchange_1_run( F_0=start_design_list[,,start], B_mag,k, 
                                              submodels, sign_vects, log_lambda=1, 
                                              log_lambda_strategy=log_lambda_strategy,
                                              sigma = 1, int_lower_bound=int_lower_bound, int_upper_bound=int_upper_bound , 
                                              int_method=int_method, measure_decsision= measure_decsision, max_iter=max_iter, tol = tol, early_stopping = early_stopping)
    pairwise_1_step
    
  }
  
  stopCluster(cl)
  
  final_results_df = as.data.frame(final_results, row.names=FALSE)
  
  final_design = final_results_df$design[which(unlist(final_results_df$P_0) == max(unlist(final_results_df$P_0)))]
  final_measure = max(unlist(final_results_df$P_0))
  
  return(list("final_design"= final_design[[1]], "final_measure"= final_measure))
  
}

pairwise_exchange_parallel_val <-function(start_design_list, cores, B_mag,k=NULL, submodels=NULL,val_submodels = NULL, sign_vects=NULL, log_lambda=NULL, log_lambda_strategy="integrate",
                                      sigma = 1, int_lower_bound=NULL, int_upper_bound=NULL, int_method="Sum", measure_decsision= "mean", max_iter=100, tol=0.0005, early_stopping = 0.0051){
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  print("Starting Parallel For Loop")
  
  final_results = foreach(start=1:dim(start_design_list)[3], .combine = rbind, .packages =c('MASS',"mvtnorm","dplyr"), .export = ls(globalenv())) %dopar%{
    pairwise_1_step <- pairwise_exchange_1_run( F_0=start_design_list[,,start], B_mag,k, 
                                                submodels, sign_vects, log_lambda=1, 
                                                log_lambda_strategy=log_lambda_strategy,
                                                sigma = 1, int_lower_bound=int_lower_bound, int_upper_bound=int_upper_bound , 
                                                int_method=int_method, measure_decsision= measure_decsision, max_iter=max_iter, tol = tol, early_stopping = early_stopping)
    pairwise_1_step
    
  }
  
  stopCluster(cl)
  
  final_results_df = as.data.frame(final_results, row.names=FALSE)
  
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  prob_results = foreach(start=1:4, .combine = rbind, .packages =c('MASS',"mvtnorm","dplyr"),.export = ls(globalenv())) %dopar%{
    F_0 = final_results_df$design[start][[1]]
    V_0_half= diag(lapply(as.data.frame(F_0),fn))
    F_0_cs = as.data.frame(lapply(as.data.frame(F_0), cent_scale))
    
    prob <- as.numeric(measure_list(Joint_prob_general(F_cs =F_0_cs,V_0_half, 
                                            B_mag=B_mag,k=k,
                                            sign_vects=sign_vects,
                                            log_lambda = log_lambda,
                                            submodels = val_submodels,
                                            log_lambda_strategy= log_lambda_strategy,
                                            int_method=int_method,
                                            int_lower_bound=int_lower_bound, int_upper_bound=int_upper_bound,
                                            sigma = 1)$results)[measure_decsision])
    prob
    
  }
  
  stopCluster(cl)
  print("final design validation measures")
  print(sort(prob_results,decreasing = TRUE, index.return = TRUE)$x[1:4])
  max_index = sort(prob_results,decreasing = TRUE, index.return = TRUE)$ix[1]
  
  print(max_index)
  
  final_design = final_results_df$design[max_index]
  final_measure = sort(prob_results,decreasing = TRUE, index.return = TRUE)$x[1]
  
  return(list("final_design"= final_design[[1]], "final_measure"= final_measure))
  
}
  

pairwise_exchange_parallel_dist_sampling <-function(start_design_list, cores, B_mag,k=NULL, sign_vects=NULL, log_lambda=NULL, log_lambda_strategy="integrate",
                                      sigma = 1, int_lower_bound=NULL, int_upper_bound=NULL, int_method="Sum", measure_decsision= "mean", max_iter=100, tol=0.0005, step_size=10){
  #browser()
  print("Starting Parallel For Loop")
  # V_0_half= diag(lapply(as.data.frame(start_design_list[,,1]),fn))
  # F_0_cs = as.data.frame(lapply(as.data.frame(start_design_list[,,1]), cent_scale))
  # JointProbResults = Joint_prob_general(F_cs = F_0_cs, V_half=V_0_half,
  #                                       B_mag=B_mag,k=k,
  #                                       sign_vects=sign_vects,
  #                                       log_lambda = log_lambda,
  #                                       submodels = NULL,
  #                                       log_lambda_strategy=log_lambda_strategy,
  #                                       sigma = sigma, int_lower_bound =int_lower_bound, int_upper_bound = int_upper_bound, int_method = int_method)
  # A_dist_list = dist_sample(JointProbResults, step_size=step_size, k_star = length(B_mag))
  # for(i in c(2:dim(start_design_list)[3])){
  #   V_0_half= diag(lapply(as.data.frame(start_design_list[,,i]),fn))
  #   F_0_cs = as.data.frame(lapply(as.data.frame(start_design_list[,,i]), cent_scale))
  #   JointProbResults = Joint_prob_general(F_cs = F_0_cs, V_half=V_0_half,
  #                                         B_mag=B_mag,k=k,
  #                                         sign_vects=sign_vects,
  #                                         log_lambda = log_lambda,
  #                                         submodels = NULL,
  #                                         log_lambda_strategy=log_lambda_strategy,
  #                                         sigma = sigma, int_lower_bound =int_lower_bound, int_upper_bound = int_upper_bound, int_method = int_method)
  #   A_i =dist_sample(JointProbResults, step_size=step_size, k_star = length(B_mag))
  #   abind(A_dist_list,A_i)
  # }
  # 
  # browser()
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  final_results = foreach(start=1:dim(start_design_list)[3], .combine = rbind, .packages =c('MASS',"mvtnorm","dplyr"), .export = ls(globalenv())) %dopar%{
    n = dim(start_design_list[,,start])[1]
    p = dim(start_design_list[,,start])[2]
    # center and scale
    V_0_half= diag(lapply(as.data.frame(start_design_list[,,start]),fn))
    F_0_cs = as.data.frame(lapply(as.data.frame(start_design_list[,,start]), cent_scale))
    #browser()
    # This gives a list of several design measures in a named list
    # ( F_cs, V_half, B_mag,k=NULL, submodels=NULL, sign_vects=NULL, log_lambda=NULL, log_lambda_strategy="integrate",
    #   sigma = 1, int_lower_bound=NULL, int_upper_bound=NULL, int_method="Sum"){
    JointProbResults = Joint_prob_general(F_cs = F_0_cs, V_half=V_0_half,
                                                     B_mag=B_mag,k=k,
                                                     sign_vects=sign_vects,
                                                     log_lambda = log_lambda,
                                                     submodels = NULL,
                                                     log_lambda_strategy=log_lambda_strategy,
                                                     sigma = sigma, int_lower_bound =int_lower_bound, int_upper_bound = int_upper_bound, int_method = int_method)
    A_dist = dist_sample(JointProbResults, step_size=step_size, k_star = length(B_mag))
    pairwise_1_step <- pairwise_exchange_1_run( F_0=start_design_list[,,start], k,B_mag,
                                                submodels=A_dist, sign_vects, log_lambda=1,
                                                log_lambda_strategy=log_lambda_strategy,
                                                sigma = 1, int_lower_bound=int_lower_bound, int_upper_bound=int_upper_bound ,
                                                int_method=int_method, measure_decsision= measure_decsision, max_iter=max_iter, tol = tol)

    pairwise_1_step
    # JointProbResults = Joint_prob_general(F_cs=F_0_cs, V_half = V_0_half, B_mag=B_mag,k=k,submodels =NULL,sign_vects=sign_vects,
    #                                       log_lambda = 1, log_lambda_strategy = log_lambda_strategy,sigma=sigma,int_lower_bound = int_lower_bound, int_upper_bound = int_upper_bound, int_method = int_method)
    
  }
  
  close(pb)
  stopCluster(cl)
  #browser()
  final_results_df = as.data.frame(final_results, row.names=FALSE)
  
  final_design = final_results_df$design[which(unlist(final_results_df$P_0) == max(unlist(final_results_df$P_0)))]
  final_measure = max(unlist(final_results_df$P_0))
  
  return(list("final_design"= final_design[[1]], "final_measure"= final_measure, "final_iters"= final_results_df$final_iteration))
  
}





pairwise_exchange_1_step <-function( F_0, B_mag,k=NULL, submodels=NULL, sign_vects=NULL, log_lambda=NULL, log_lambda_strategy="integrate",
                                    sigma = 1, int_lower_bound=NULL, int_upper_bound=NULL, int_method="Sum", measure_decsision= "mean", 
                                    max_iter=10, tol=0.0005, max_lookback=1, early_stopping = 0.01){
  
  # F_0 = F_init
  n = dim(F_0)[1]
  p = dim(F_0)[2]
  #center and scale
  V_0_half= diag(lapply(as.data.frame(F_0),fn))
  F_0_cs = as.data.frame(lapply(as.data.frame(F_0), cent_scale))
  #browser()
  # This gives a list of several design measures in a named list
  measure_list_0 = measure_list(Joint_prob_general(F_cs = F_0_cs, V_half=V_0_half, 
                                                   B_mag=B_mag,k=k,
                                                   sign_vects=sign_vects,
                                                   log_lambda = log_lambda,
                                                   submodels = submodels,
                                                   log_lambda_strategy=log_lambda_strategy,
                                                   sigma = sigma, int_lower_bound =int_lower_bound,int_method=int_method, int_upper_bound = int_upper_bound)$results)
  # Choose which of these measures to make a decision on. 
  P_0 = as.numeric(measure_list_0[measure_decsision])
  
  
  iterations = 0
  lookback = 0
  change = 0
  #tracking_data = data.frame(col_ind_1= numeric(0), col_index_2 = numeric(0), P_prev=numeric(0), P_new=numeric(0))
  #stop = FALSE
  #taboo_index = NULL
  tab_list = list("tab1"= NULL, "tab2"=NULL, "tab3"=NULL)
  #print(P_0)
  P_prev=P_0
  P_new=P_0
  P_early_stopping = P_0
  #browser()
  while ( iterations <= max_iter & change == 0){
    if( P_new - P_prev < tol){
      lookback = lookback +1
    }
    if( lookback > max_lookback){
      break
    }
    else{
      P_prev=P_new
      lookback = 0 
    }
    if(iterations !=0 &&iterations %% 5 ==0){
      if(P_prev - P_early_stopping < early_stopping){
        break
      }
    }
    
    
    max_corr_cols = get_max_corr(F_0_cs, taboo_index= tab_list)
    #print(max_corr_cols)
    tab_list = update_tab_list(max_corr_cols, tab_list) 
    #browser()
    #print(tab_list)
    #Subset the design to the most correlated columns
    for( i in c(1:n)){
      #print(i)
      # Generate three new designs based on the different pairwise combinations of one row of 2 cols of -1 and 1\
      # This might be a very sloppy way of doing this
      F_1=F_0
      F_2=F_0
      F_3=F_0
      #browser()
      F_1[i,max_corr_cols]<- F_0[i,max_corr_cols]*c(-1,1)
      F_2[i,max_corr_cols]<- F_0[i,max_corr_cols]*c(1,-1)
      F_3[i,max_corr_cols]<- F_0[i,max_corr_cols]*c(-1,-1)
      
      V_1_half= diag(lapply(as.data.frame(F_1),fn))
      F_1_cs = as.data.frame(lapply(as.data.frame(F_1), cent_scale))
      V_2_half= diag(lapply(as.data.frame(F_2),fn))
      F_2_cs = as.data.frame(lapply(as.data.frame(F_2), cent_scale))
      V_3_half= diag(lapply(as.data.frame(F_3),fn))
      F_3_cs = as.data.frame(lapply(as.data.frame(F_3), cent_scale))
      
      # Calculate the measure list for all designs
      measure_list_1 = measure_list(Joint_prob_general(F_cs = F_1_cs, V_half=V_1_half, 
                                                       B_mag=B_mag,k=k,
                                                       sign_vects=sign_vects,
                                                       log_lambda = log_lambda,
                                                       submodels = submodels,
                                                       log_lambda_strategy=log_lambda_strategy,
                                                       sigma = sigma, int_lower_bound =int_lower_bound, 
                                                       int_method=int_method, int_upper_bound = int_upper_bound)$results)
      measure_list_2 = measure_list(Joint_prob_general(F_cs = F_2_cs, V_half=V_2_half, 
                                                       B_mag=B_mag,k=k,
                                                       sign_vects=sign_vects,
                                                       log_lambda = log_lambda,
                                                       submodels = submodels,
                                                       log_lambda_strategy=log_lambda_strategy,
                                                       sigma = sigma, int_lower_bound =int_lower_bound, int_method=int_method,int_upper_bound = int_upper_bound)$results)
      measure_list_3 = measure_list(Joint_prob_general(F_cs = F_3_cs, V_half=V_3_half, 
                                                       B_mag=B_mag,k=k,
                                                       sign_vects=sign_vects,
                                                       log_lambda = log_lambda,
                                                       submodels = submodels,
                                                       log_lambda_strategy=log_lambda_strategy,
                                                       sigma = sigma, int_lower_bound =int_lower_bound,int_method=int_method, int_upper_bound = int_upper_bound)$results)
      # Choose which of these measures to make a decision on. 
      #browser()
      P_1 = as.numeric(measure_list_1[measure_decsision])
      P_2 = as.numeric(measure_list_2[measure_decsision])
      P_3 = as.numeric(measure_list_3[measure_decsision])
      result_table <-tribble(~measure, ~design,~design_cs,~V_half,
                             P_0, F_0,F_0_cs, V_0_half,
                             P_1, F_1, F_1_cs, V_1_half,
                             P_2, F_2, F_2_cs, V_2_half,
                             P_3, F_3, F_3_cs, V_3_half)
      max_index = which(result_table$measure==max(result_table$measure), arr.ind=TRUE)
      #Update the optimal design
      if (max_index!=1){
        change = 1
      }
      #browser()
      F_0 = matrix(unlist(result_table[max_index,]$design), nrow = n, ncol = p)
      F_0_cs = matrix(unlist(result_table[max_index,]$design_cs), nrow = n, ncol = p)
      V_0_half= matrix(unlist(result_table[max_index,]$V_half), nrow = p, ncol = p)
      P_new = max(result_table$measure)
      #print(P_new)
      P_0=P_new
      #tracking_data[nrow(tracking_data)+1,]<-c(max_corr_cols[1], max_corr_cols[2], P_prev, P_new)
      #browser()
    }
    #tracking_data[nrow(tracking_data)+1,]<-c(max_corr_cols[1], max_corr_cols[2], P_prev, P_new)
    iterations =iterations +1
  }
  
  
  return(list("design"= F_0,"P_0"=P_new, "final_iteration"= iterations))
}



#### Analog to the code above just for S only. This was used for section 2. For future work, build in flexibility on what event or joint even to optimize over and delete the folowing code. 

S_prob_all_fixed <-function(A, F_cs, V_half, B_mag, Z_A, log_lambda, sigma = 1){
  # Subsetting design into active and inactive columns
  print(A)
  F_A=as.matrix(F_cs[,A])
  #F_I = as.matrix(F_cs[,-A])
  V_half=diag(diag(V_half)[A])
  
  lambda = exp(log_lambda)
  
  n = dim(F_cs)[1]
  k = length(Z_A)
  one_k = rep(1,k)
  B_min_vec = B_mag*Z_A
  C_AA = (1/n)*t(F_A)%*%F_A
  
  
  
  # Check if C_AA has an inverse, if not, use moore-penrose. (Discuss with group)
  C_AA_inv = n*ginv(t(F_A)%*%F_A)
  
  
  Cov_S = (sigma^2)* diag(Z_A)%*%C_AA_inv%*%diag(Z_A)
  upper_S = sqrt(n)*diag(Z_A)%*%V_half%*%B_min_vec -
    lambda*sqrt(n)* diag(C_AA_inv%*%Z_A%*%t(Z_A))
  #browser()
  P_S= mvtnorm:: pmvnorm(lower= -Inf, upper = as.vector(upper_S), mean= rep(0,k), sigma = Cov_S)[1]
  
  return(P_S)
  
}

S_prob_all_fixed_deconstructed<-function( Cov_S, upper_S,mean_S_WO_lam, log_lambda){
  
  lambda = exp(log_lambda)
  k = dim(Cov_S)[1]
  
  
  
  P_S= mvtnorm:: pmvnorm(lower= -Inf, upper = as.vector(upper_S)-lambda*as.vector(mean_S_WO_lam), mean= rep(0,k), sigma = Cov_S)[1]
  
  
  return(P_S)
  
}

opt_log_lambda_S <- function( log_lam_min, log_lam_max, Cov_S, upper_S, mean_S_WO_lam, stepsize= 0.05){
  # Get warm start for lamnda values
  #browser()
  log_lambda_grid = seq(log_lam_min,log_lam_max,stepsize)
  grid_probs <-sapply(log_lambda_grid, S_prob_all_fixed_deconstructed,
                      Cov_S = Cov_S, upper_S= upper_S, mean_S_WO_lam= mean_S_WO_lam)
  warm_start <- cbind(log_lambda_grid,grid_probs)
  #print(max(warm_start[,2]))
  
  if (is.na(max(warm_start[,2]))|| max(warm_start[,2])<0.0001){
    start_lam <- runif(1, 0, 3)
  }
  else{
    start_lam <- warm_start[which.max(warm_start[,2]),1]
  }
  # maybe consider doing box-constaints, but since log-lam can take any real number, I thought other solvers could be faster
  #browser()
  opt_lam <-tryCatch(optim(par=start_lam, fn=S_prob_all_fixed_deconstructed, Cov_S = Cov_S, upper_S= upper_S, mean_S_WO_lam= mean_S_WO_lam,
                           method="L-BFGS-B", upper = 10, lower = -Inf,
                           control = list("fnscale"=-1,warn.1d.NelderMead=FALSE)), error = function(c){return(list("par"=start_lam, "value"=0))})
  return(list("opt_log_lam"= opt_lam$par, "opt_val" = opt_lam$value))
  
  
}


integrate_log_lambda_S <- function( Cov_S, upper_S, mean_S_WO_lam, int_lower_bound=-4.5, int_upper_bound = 1.5, method = "Quadrature", step_size = 0.1){

  int_methods = c("Quadrature", "Sum")
  ROI = int_upper_bound-int_lower_bound
  if(!method%in% int_methods){
    stop("The input integration method must be either  'Quadrature', or 'Sum'")
  }
  if (method == "Quadrature"){
    S_prob_all_fixed_vec = Vectorize(S_prob_all_fixed_deconstructed, "log_lambda")
    #browser()
    #print("Integrating via quadrature")
    
    integrate_results = integrate(S_prob_all_fixed_vec, lower= int_lower_bound, upper= int_upper_bound,
                                   Cov_S = Cov_S,
                                  upper_S= upper_S, mean_S_WO_lam= mean_S_WO_lam)
    return(integrate_results$value/ROI)
  }
  else{
    #print("Integrating via Riemann Sum")
    log_lambda_grid <- seq(int_lower_bound, int_upper_bound, by=step_size)
    #browser()
    return((sum(step_size*as.data.frame(lapply(log_lambda_grid, S_prob_all_fixed_deconstructed, Cov_S = Cov_S,
                                               upper_S= upper_S, mean_S_WO_lam= mean_S_WO_lam  ))))/ROI)
    
    
  }
  
  
  
}



S_prob_submodel_sign_fixed <-function(A, F_cs, V_half, B_mag, Z_A, log_lambda=NULL, log_lambda_strategy="integrate",
                                          sigma = 1, int_lower_bound=NULL, int_upper_bound=NULL, int_method = "Quadrature"){
  strategies = c("fixed", "optimal", "integrate")
  if(!log_lambda_strategy%in% strategies){
    stop("The input log_lam strategy must be either 'fixed, 'optimal', or 'integrate'")
  }
  # Subsetting design into active
  F_A=as.matrix(F_cs[,A])
  V_half=diag(diag(V_half)[A])
  
  #lambda = exp(log_lambda)
  
  n = dim(F_cs)[1]
  k = length(Z_A)
  one_k = rep(1,k)

  B_min_vec = B_mag*Z_A
  
  C_AA = (1/n)*t(F_A)%*%F_A
  
  
  
  # Check if C_AA has an inverse, if not, set the joint prob to zero. (Discuss with group)
  
  C_AA_inv = n*ginv(t(F_A)%*%F_A)
  
 
  
  # Fix the diagonals of Cov_I so that they are slightly positive, if they are slightly negative
  
  
  Cov_S = (sigma^2)* diag(Z_A)%*%C_AA_inv%*%diag(Z_A)
  #browser()
  upper_S = sqrt(n)*diag(Z_A)%*%V_half%*%B_min_vec
  mean_S_WO_lam = sqrt(n)* diag(C_AA_inv%*%Z_A%*%t(Z_A))
  if( log_lambda_strategy=="fixed"){
    if(is.null(log_lambda)){
      stop("For a fixed lambda strategy, you must specify what the fixed log_lambda is. It is currently NULL ")
    }
    else{
      return(S_prob_all_fixed_deconstructed( Cov_S, upper_S, mean_S_WO_lam, log_lambda))
    }
  }
  if( log_lambda_strategy=="integrate"){
    if(any(is.null(c(int_lower_bound, int_lower_bound)))){
      warning(" Integral upper and/or lower bounds are not specified, using defualt of lower bound at -4.5, and upper at 1.5.")
      return(integrate_log_lambda_S( Cov_S,
                                  upper_S= upper_S, mean_S_WO_lam= mean_S_WO_lam, method =int_method ))
    }
    else{ 
      return(integrate_log_lambda_S( Cov_S,
                                  upper_S= upper_S, mean_S_WO_lam= mean_S_WO_lam,
                                  int_lower_bound=int_lower_bound ,
                                  int_upper_bound = int_upper_bound, method= int_method))
    }
    
  }
  
  if(log_lambda_strategy=="optimal"){
    #print(A)
    #browser()
    return(opt_log_lambda_S( log_lam_min=-2, log_lam_max=3, Cov_S, upper_S= upper_S, mean_S_WO_lam= mean_S_WO_lam, stepsize= 0.05)$opt_val)
  }
  
  
  
}

### Now we move to a fixed submodel, but the ability to loop over sign vectors

S_prob_submodel_fixed <-function(A, F_cs, V_half, B_mag, sign_vects=NULL, log_lambda=NULL, log_lambda_strategy="integrate",
                                     sigma = 1, int_lower_bound=NULL, int_upper_bound=NULL, int_method="Quadrature"){
  
  if(is.null(sign_vects)){
    sign_vects= gen_all_sign_vects(A)
  }
  colnames(sign_vects) <- paste("s", 1:length(A),sep="")
  sign_vects = as.data.frame(sign_vects)
  #browser()
  results <- apply(as.data.frame(sign_vects), 1, S_prob_submodel_sign_fixed, A=A, F_cs=F_cs, V_half=V_half, B_mag=B_mag,
                   log_lambda=log_lambda, log_lambda_strategy=log_lambda_strategy,
                   sigma = sigma, int_lower_bound=int_lower_bound, int_upper_bound=int_upper_bound, int_method=int_method)
  #browser()
  output = cbind(sign_vects, results)
  return(output)
}
S_prob_general <-function( F_cs, V_half, B_mag,k=NULL, submodels=NULL, sign_vects=NULL, log_lambda=NULL, log_lambda_strategy="integrate",
                               sigma = 1, int_lower_bound=NULL, int_upper_bound=NULL, int_method="Sum"){
  p=dim(F_cs)[2]
  n = dim(F_cs)[1]
  if(any(is.na(F_cs))){
    return(c(0,0,0)) 
  }
  if(is.null(submodels)){
    if(is.null(k)){
      k = floor(n/3)
    }
    submodels= gen_all_submodels(p,k)
  }
  else{ k = dim(submodels)[2]}
  if(is.null(sign_vects)){
    num_sign_vects = 2^(k-1)
  }
  else{
    num_sign_vects = dim(sign_vects)[1]
  }
  #browser()
  colnames(submodels) <- paste("A", 1:k,sep="")
  submodels.expanded = submodels[rep(seq_len(nrow(submodels)), rep(num_sign_vects,nrow(submodels))), 1:k]
  #browser()
  results <- do.call(rbind,apply(submodels,1, S_prob_submodel_fixed,sign_vects=sign_vects, F_cs=F_cs, V_half=V_half, B_mag=B_mag,
                                 log_lambda=log_lambda, log_lambda_strategy=log_lambda_strategy,
                                 sigma = sigma, int_lower_bound=int_lower_bound, int_upper_bound=int_upper_bound, int_method=int_method))
  #browser()
  output = cbind(submodels.expanded, results)
  return(output)
}

pairwise_exchange_1_run_S <-function( F_0, B_mag,k=NULL, submodels=NULL, sign_vects=NULL, log_lambda=NULL, log_lambda_strategy="integrate",
                                    sigma = 1, int_lower_bound=NULL, int_upper_bound=NULL, int_method="Sum", measure_decsision= "mean", 
                                    max_iter=10, tol=0.0005, max_lookback=1, early_stopping = 0.01){
  
  # F_0 = F_init
  n = dim(F_0)[1]
  p = dim(F_0)[2]
  #center and scale
  V_0_half= diag(lapply(as.data.frame(F_0),fn))
  F_0_cs = as.data.frame(lapply(as.data.frame(F_0), cent_scale))
  #browser()
  # This gives a list of several design measures in a named list
  measure_list_0 = measure_list(S_prob_general(F_cs = F_0_cs, V_half=V_0_half, 
                                                   B_mag=B_mag,k=k,
                                                   sign_vects=sign_vects,
                                                   log_lambda = log_lambda,
                                                   submodels = submodels,
                                                   log_lambda_strategy=log_lambda_strategy,
                                                   sigma = sigma, int_lower_bound =int_lower_bound,int_method=int_method, int_upper_bound = int_upper_bound)$results)
  # Choose which of these measures to make a decision on. 
  P_0 = as.numeric(measure_list_0[measure_decsision])
  
  
  iterations = 0
  lookback = 0
  #tracking_data = data.frame(col_ind_1= numeric(0), col_index_2 = numeric(0), P_prev=numeric(0), P_new=numeric(0))
  #stop = FALSE
  #taboo_index = NULL
  tab_list = list("tab1"= NULL, "tab2"=NULL, "tab3"=NULL)
  #print(P_0)
  P_prev=P_0
  P_new=P_0
  P_early_stopping = P_0
  #browser()
  while ( iterations <= max_iter ){
    if( P_new - P_prev < tol){
      lookback = lookback +1
    }
    if( lookback > max_lookback){
      break
    }
    else{
      P_prev=P_new
      lookback = 0 
    }
    if(iterations !=0 &&iterations %% 5 ==0){
      if(P_prev - P_early_stopping < early_stopping){
        break
      }
    }
    
    
    max_corr_cols = get_max_corr(F_0_cs, taboo_index= tab_list)
    #print(max_corr_cols)
    tab_list = update_tab_list(max_corr_cols, tab_list) 
    #browser()
    #print(tab_list)
    #Subset the design to the most correlated columns
    for( i in c(1:n)){
      #print(i)
      # Generate three new designs based on the different pairwise combinations of one row of 2 cols of -1 and 1\
      # This might be a very sloppy way of doing this
      F_1=F_0
      F_2=F_0
      F_3=F_0
      #browser()
      F_1[i,max_corr_cols]<- F_0[i,max_corr_cols]*c(-1,1)
      F_2[i,max_corr_cols]<- F_0[i,max_corr_cols]*c(1,-1)
      F_3[i,max_corr_cols]<- F_0[i,max_corr_cols]*c(-1,-1)
      
      V_1_half= diag(lapply(as.data.frame(F_1),fn))
      F_1_cs = as.data.frame(lapply(as.data.frame(F_1), cent_scale))
      V_2_half= diag(lapply(as.data.frame(F_2),fn))
      F_2_cs = as.data.frame(lapply(as.data.frame(F_2), cent_scale))
      V_3_half= diag(lapply(as.data.frame(F_3),fn))
      F_3_cs = as.data.frame(lapply(as.data.frame(F_3), cent_scale))
      
      # Calculate the measure list for all designs
      measure_list_1 = measure_list(S_prob_general(F_cs = F_1_cs, V_half=V_1_half, 
                                                       B_mag=B_mag,k=k,
                                                       sign_vects=sign_vects,
                                                       log_lambda = log_lambda,
                                                       submodels = submodels,
                                                       log_lambda_strategy=log_lambda_strategy,
                                                       sigma = sigma, int_lower_bound =int_lower_bound, 
                                                       int_method=int_method, int_upper_bound = int_upper_bound)$results)
      measure_list_2 = measure_list(S_prob_general(F_cs = F_2_cs, V_half=V_2_half, 
                                                       B_mag=B_mag,k=k,
                                                       sign_vects=sign_vects,
                                                       log_lambda = log_lambda,
                                                       submodels = submodels,
                                                       log_lambda_strategy=log_lambda_strategy,
                                                       sigma = sigma, int_lower_bound =int_lower_bound, int_method=int_method,int_upper_bound = int_upper_bound)$results)
      measure_list_3 = measure_list(S_prob_general(F_cs = F_3_cs, V_half=V_3_half, 
                                                       B_mag=B_mag,k=k,
                                                       sign_vects=sign_vects,
                                                       log_lambda = log_lambda,
                                                       submodels = submodels,
                                                       log_lambda_strategy=log_lambda_strategy,
                                                       sigma = sigma, int_lower_bound =int_lower_bound,int_method=int_method, int_upper_bound = int_upper_bound)$results)
      # Choose which of these measures to make a decision on. 
      #browser()
      P_1 = as.numeric(measure_list_1[measure_decsision])
      P_2 = as.numeric(measure_list_2[measure_decsision])
      P_3 = as.numeric(measure_list_3[measure_decsision])
      result_table <-tribble(~measure, ~design,~design_cs,~V_half,
                             P_0, F_0,F_0_cs, V_0_half,
                             P_1, F_1, F_1_cs, V_1_half,
                             P_2, F_2, F_2_cs, V_2_half,
                             P_3, F_3, F_3_cs, V_3_half)
      max_index = which(result_table$measure==max(result_table$measure), arr.ind=TRUE)
      #Update the optimal design
      #browser()
      F_0 = matrix(unlist(result_table[max_index,]$design), nrow = n, ncol = p)
      F_0_cs = matrix(unlist(result_table[max_index,]$design_cs), nrow = n, ncol = p)
      V_0_half= matrix(unlist(result_table[max_index,]$V_half), nrow = p, ncol = p)
      P_new = max(result_table$measure)
      #print(P_new)
      P_0=P_new
      #tracking_data[nrow(tracking_data)+1,]<-c(max_corr_cols[1], max_corr_cols[2], P_prev, P_new)
      #browser()
    }
    #tracking_data[nrow(tracking_data)+1,]<-c(max_corr_cols[1], max_corr_cols[2], P_prev, P_new)
    iterations =iterations +1
  }
  
  
  return(list("design"= F_0,"P_0"=P_new, "final_iteration"= iterations))
}
pairwise_exchange_parallel_S <-function(start_design_list, cores, B_mag,k=NULL, submodels=NULL, sign_vects=NULL, log_lambda=NULL, log_lambda_strategy="integrate",
                                      sigma = 1, int_lower_bound=NULL, int_upper_bound=NULL, int_method="Sum", measure_decsision= "mean", max_iter=100, tol=0.0005, early_stopping = 0.0051){
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  print("Starting Parallel For Loop")
  
  final_results = foreach(start=1:dim(start_design_list)[3], .combine = rbind, .packages =c('MASS',"mvtnorm","dplyr"), .export = ls(globalenv())) %dopar%{
    pairwise_1_step <- pairwise_exchange_1_run_S( F_0=start_design_list[,,start], B_mag,k, 
                                                submodels, sign_vects, log_lambda=1, 
                                                log_lambda_strategy=log_lambda_strategy,
                                                sigma = 1, int_lower_bound=int_lower_bound, int_upper_bound=int_upper_bound , 
                                                int_method=int_method, measure_decsision= measure_decsision, max_iter=max_iter, tol = tol, early_stopping = early_stopping)
    pairwise_1_step
    
  }
  
  stopCluster(cl)
  
  final_results_df = as.data.frame(final_results, row.names=FALSE)
  
  final_design = final_results_df$design[which(unlist(final_results_df$P_0) == max(unlist(final_results_df$P_0)))]
  final_measure = max(unlist(final_results_df$P_0))
  
  return(list("final_design"= final_design[[1]], "final_measure"= final_measure))
  
}







#Compound Symmetric Code


make_cs_matrix <-function(p,c){
  # This generates a completely symmetric matrix with off diagonals c
  mat = matrix(rep(c,p*p), nrow=p, ncol=p)
  diag(mat)= rep(1,p)
  return(mat)
  #return(CSgenerate(p,c))
}

prob_conditon_CS_joint <- function(c, sigma_sq, n, kstar,q,lambda,B1, sign_vec){
  # This generates the probabilites of inclusion and exclusion for a CS X1TX1 and X2TX2
  
  
  X1TX1= n*make_cs_matrix(p=kstar,c=c)
  X2TX2 =n* make_cs_matrix(p=q,c=c)
  
  X2TX1 = n*c*matrix(rep(1,kstar*q), nrow=q, ncol=kstar)
  
  X1TX2=t(X2TX1)
  
  Sig1<-diag(1/(1-c),kstar)-(c/(1-c))*(1/(1+c*(kstar-1)))*matrix(1,nrow=kstar,ncol=kstar)
  
  Var_inc<-sigma_sq*diag(sign_vec)%*%Sig1%*%diag(sign_vec)
  #Var_inc = sigma_sq * n *(diag(sign_vec)%*% solve(X1TX1)%*% diag(sign_vec))
  
  
  one_vec_inc = rep(1,kstar)
  
  z = t(one_vec_inc)%*%sign_vec
  #browser()
  
  
  mean_inc = sqrt(n)*lambda*( ( (1/(1-c))*one_vec_inc ) - (z*(c/((1-c)*(1+c*(kstar-1))))*sign_vec) )
  
  
  #rhs_inc = n*lambda* solve(X1TX1)%*%one_vec_inc-B1
  
  #mean_vec = rhs_inc
  
  #browser()
  rhs_vec = sqrt(n) * one_vec_inc*B1
  prob_inc =pmvnorm(upper = c(rhs_vec), mean = c(mean_inc), sigma=Var_inc)
  
  #browser()
  
  I_q = diag(rep(1,q))
  
  J_q = matrix(rep(1,q*q), nrow=q, ncol=q)
  
  Var_ex = sigma_sq* ( ((1-c)*I_q) + ((c*(1-c)/(1+c*(kstar-1)))*J_q) )
  one_vec_ex =rep(1,q)
  
  mean_ex = sqrt(n)*lambda*z* (c/(1+c*(kstar-1)))*one_vec_ex
  
  #bound = n*lambda*(one_vec_ex-abs( ((c*pstar)/(1+c*(pstar-1))) *one_vec_ex))
  bound = sqrt(n)*lambda*one_vec_ex
  if (any(bound < 0)){
    prob_ex =0
  }
  else{prob_ex=pmvnorm(lower=-bound, upper=bound,mean = c(mean_ex), sigma=Var_ex)}
  # return(list('prob_inclusion'=prob_inc[1], 'prob_exclusion' = prob_ex[1], "joint"= prob_inc[1]*prob_ex[1]))
  joint = prob_inc[1]*prob_ex[1]
  return(joint)
}


prob_CS_joint_sum_lambda <-function(c, sigma_sq, n, kstar,q,lambda_list,B1, sign_vec){
  return(sum(sapply(lambda_list, prob_conditon_CS_joint, c=c, sigma_sq=sigma_sq,n=n, kstar=kstar, q=q,B1=B1, sign_vec = sign_vec)))
}

prob_CS_joint_max_lambda <-function(c, sigma_sq, n, kstar,q,lambda_list,B1, sign_vec){
  
  return(max(sapply(lambda_list, prob_conditon_CS_joint, c=c, sigma_sq=sigma_sq,n=n, kstar=kstar, q=q,B1=B1, sign_vec = sign_vec)))
}

opt_CS_c <- function(c_min, c_max, sigma_sq, n, kstar,q,lambda_list,B1, sign_vec, int=TRUE){
  if (int){
    c_opt = optim(par=0, prob_CS_joint_sum_lambda, sigma_sq=sigma_sq, n=n, kstar=kstar, q=q, 
                  lambda_list=lambda_list,B1=B1, sign_vec=sign_vec, 
                  lower=c_min, upper=c_max, control= list(fnscale=-1))
  }
  else{
    c_opt = optim(par=0, prob_CS_joint_max_lambda, sigma_sq=sigma_sq, n=n, kstar=kstar, q=q, 
                  lambda_list=lambda_list,B1=B1, sign_vec=sign_vec, 
                  lower=c_min, upper=c_max, control= list(fnscale=-1))
  }
  return(c_opt)
}


eval_exact_design_with_CS<-function(X, c_opt){
  n= nrow(X)
  p =ncol(X)
  V= (diag(lapply(as.data.frame(X),fn)))^2
  m_V = sum(diag(V))/p
  
  C_opt = make_cs_matrix(p,c_opt)
  X_cs = as.data.frame(lapply(as.data.frame(X), cent_scale))
  C= t(as.matrix(X_cs))%*%as.matrix(X_cs)/n
  frob_norm = norm(abs(C-C_opt), type="F")
  return(list("m_V"= m_V, "frob_norm"=frob_norm))
}


# 
# 
# 
# ##########################################################
# # Evalfunc: a function for calculating the criteria values for EPCEA
# ##########################################################
# 
evalfunc = function(des, c_opt) {
  
  nrun = nrow(des)
  #colnames(des) = c("x1","x2","x3")
  #desfr = as.data.frame(des)
  #M = crossprod(model.matrix(~ (x1+x2+x3)^2+I(x1^2)+I(x2^2)+I(x3^2), desfr))
  
  results = eval_exact_design_with_CS(des,c_opt)
  
  if(is.na(results$frob_norm)){ return(c(0,-1000))}
  else{ return(c(results$m_V, -results$frob_norm))}
}

compare = function(newpt, newdes, curpf, curpfdes, nfactor)
{
  #browser()
  g1 = round(newpt[1L], 8L) > round(curpf[,1L], 8L)
  g2 = round(newpt[2L], 8L) > round(curpf[,2L], 8L)
  
  ge1 = round(newpt[1L], 8L) >= round(curpf[,1L], 8L)
  ge2 = round(newpt[2L], 8L) >= round(curpf[,2L], 8L)
  
  l1 = round(newpt[1L], 8L) < round(curpf[,1L], 8L)
  l2 = round(newpt[2L], 8L) < round(curpf[,2L], 8L)
  
  le1 = round(newpt[1L], 8L) <= round(curpf[,1L], 8L)
  le2 = round(newpt[2L], 8L) <= round(curpf[,2L], 8L)
  
  eq1 = round(newpt[1L], 8L) == round(curpf[,1L], 8L)
  eq2 = round(newpt[2L], 8L) == round(curpf[,2L], 8L)
  
  cond1 = g1*ge2 + g2*ge1 == 0L
  cond2 = sum(l1*le2 + l2*le1 + eq1*eq2)
  cond3 = seq(1L, nrow(curpf))[cond1]
  cond4 = rep(0L, nfactor*length(cond3))
  
  if(length(cond3) > 0L)
  {
    #browser()
    for(i in 1L:length(cond3))
    {
      cond4[i*nfactor-nfactor + 1L:nfactor]=cond3[i]*nfactor-nfactor + 1L:nfactor
    }
  }
  
  newpf = curpf[cond1,]
  newpfdes = curpfdes[,cond4]
  if(cond2 == 0L)
  {
    newpf = rbind(newpf, newpt)
    newpfdes = cbind(newpfdes, newdes)
  }
  
  list(matrix(newpf, ncol = 2L), newpfdes)
}


cnt = function(A, B){
  ind = rep(TRUE, nrow(A))
  for (i in 1L:nrow(A)){
    for (j in 1L:nrow(B)){
      if(all(round(A[i,], 4L) == round(B[j,], 4L))){
        ind[i] = FALSE
      }
    }
  }
  #browser()
  ind
}

###########################################
# Pareto based coordinate exchange operator
###########################################

coexch1 = function(des, c_opt, nfactor) {
  
  cdesmat = des
  cbestvals = evalfunc(des, c_opt)
  cpfvals = matrix(cbestvals, ncol = 2L)
  cchange = TRUE
  cpfdes = des
  
  while(cchange)
  {
    # cdesmat = des
    # cbestvals = evalfunc(des, c_opt)
    # cpfvals = matrix(cbestvals, ncol = 2L)
    # cpfdes = des
    #
    cchange = FALSE
    
    for (i in 1L:length(cdesmat))
    {
      restart=FALSE
      for (j in level)
      {
        cnewdesmat = cdesmat
        #browser()
        if(cnewdesmat[i] == j) next
        cnewdesmat[i] = j
        cnewbestval = evalfunc(cnewdesmat, c_opt)
        
        # If new design dominates old discard the old and start over?
        if (any(round(cnewbestval, 8L) > round(cbestvals, 8L)) & all(round(cnewbestval, 8L) >= round(cbestvals, 8L)))
        {
          ctemppf = compare(t(cnewbestval), cnewdesmat, cpfvals, cpfdes, nfactor)
          cpfvals = ctemppf[[1L]]
          cpfdes = ctemppf[[2L]]
          cbestvals = cnewbestval
          cchange = TRUE
          cdesmat = cnewdesmat
          #browser()
          restart = TRUE
          break
          
        }
        #browser()
        ctemppf = compare(t(cnewbestval), cnewdesmat, cpfvals, cpfdes, nfactor)
        cpfvals = ctemppf[[1L]]
        cpfdes = ctemppf[[2L]]
        
      }
      if(restart){break}
    }
  }
  
  ctemppf
}

################################################
# EPCEA for a single random start
################################################

coexch2 = function(level, nrun, nfactor, c_opt) {
  
  des = matrix(sample(level, nrun*nfactor, replace = TRUE), ncol = nfactor)
  # colnames(des) = c("x1","x2","x3")
  # desfr = as.data.frame(des)
  # dcheck = qr(crossprod(model.matrix(~ (x1+x2+x3)^2+I(x1^2)+I(x2^2)+I(x3^2), desfr)))$rank
  #
  # while (dcheck < 10L){
  #   des = matrix(sample(level, nrun*nfactor, replace = TRUE), ncol = nfactor)
  #   colnames(des) = c("x1","x2","x3")
  #   desfr = as.data.frame(des)
  #   dcheck = qr(crossprod(model.matrix(~ (x1+x2+x3)^2+I(x1^2)+I(x2^2)+I(x3^2), desfr)))$rank}
  #
  temp = coexch1(des, c_opt, nfactor)
  psold = psnew = temp[[2L]]
  pfold = pfnew = temp[[1L]]
  
  for (i in 1L:nrow(pfold))
  {
    des = psold[,(i-1L)*nfactor+1L:nfactor]
    temp = coexch1(des, c_opt, nfactor)
    for(j in 1L:nrow(temp[[1L]]))
    {
      temppf = compare(matrix(temp[[1L]][j,], ncol = 2L), temp[[2L]][,j*nfactor-nfactor+1L:nfactor], pfnew, psnew, nfactor)
      pfnew = temppf[[1L]]
      psnew = temppf[[2L]]
      
    }
    
  }
  #browser()
  ind = cnt(pfnew, pfold)
  
  while(sum(ind) > 0)
  {
    dind = rep(ind, each = nfactor)
    pfdiff = matrix(pfnew[ind,], ncol = 2L)
    #browser()
    psdiff = psnew[,dind]
    
    pfold = pfnew
    psold = psnew
    
    for (i in 1L:nrow(pfdiff))
    {
      desmat = psdiff[,(i-1L)*nfactor+1L:nfactor]
      temp = coexch1(desmat, c_opt, nfactor)
      for(j in 1L:nrow(temp[[1L]]))
      {
        temppf = compare(matrix(temp[[1L]][j,], ncol = 2L), temp[[2L]][,j*nfactor-nfactor+1L:nfactor], pfnew, psnew, nfactor)
        pfnew = temppf[[1L]]
        psnew = temppf[[2L]]
        
      }
      
    }
    
    ind = cnt(pfnew, pfold)
    
  }
  
  list(pfnew, psnew)
}

##################################################
# EPCEA with multiple random starts
##################################################

coexchpf = function(level, nrun, nfactor, ns, c_opt)
{
  
  temp = coexch2(level, nrun, nfactor, c_opt)
  paretodes = temp[[2]]
  paretodescrit = temp[[1]]
  
  for (icount in 2L:ns)
  {
    
    temp = coexch2(level, nrun, nfactor, c_opt)
    for(i in 1L:nrow(temp[[1L]]))
    {
      temppf = compare(matrix(temp[[1L]][i,], ncol = 2L),temp[[2L]][,i*nfactor-nfactor+1L:nfactor], paretodescrit, paretodes, nfactor)
      paretodescrit=temppf[[1L]]
      paretodes=temppf[[2L]]
      
    }
    
  }
  
  list(paretodescrit, paretodes)
}




