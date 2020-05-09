#' K-fold Cross-Validation
#' @description Divide the whole dataset into K folds. Take K-1 folds as training dataset 
#' and take the rest one as the testing dataset. Return mean square error (MSE).
#' @param data precessed data
#' @param lambda the tuning parameter
#' @param K number of folds
#' @param ka number of pcs for the first response variable
#' @param kb number of pcs for the second response variable
#'
#' @return the sum of mean square error of K folds
#' @export
#'
#' @examples
#' rawdata = gen_data(n=50)
#' visit = seq(0,100,20)
#' data = predata(nobs = rawdata$nobs, 
#'                time = rawdata$time, 
#'                y = rawdata$y, 
#'                z = rawdata$z, 
#'                knots = visit, 
#'                order = 3)
#' MSE = matrix(nrow=5,ncol=5)
#' for(i in 1:5)
#'   for(j in 1:5)
#'   {
#'      lambda = c(200*i, 200*i, 4000*j, 4000*j)
#'      MSE[i,j] = Kf_CV(data, lambda, K=10, ka=1,kb=2)
#'   }
#'
Kf_CV<-function(data,lambda,K,ka,kb)
{
  ## Split the data set into K folders
  ## Use likelihood as K-folder cross-validation value
  ## source('cirteria.R')
  n    = data$n
  fold = round(n/K)
  ind  = c(0, cumsum(data$obs_times))
  sumn = sum(data$obs_times)
  
  ## Each iteration, 1 folder as testing data, K-1 folder as training data
  value= 0
  for(i in 1:K)
  {
    ## Index of training and testing dataset
    if(i!=K){
      ind_remove = ((i-1)*fold+1):(i*fold)
      ind_remain = setdiff(1:n,ind_remove)
      dat_remove = (ind[(i-1)*fold+1]+1):ind[i*fold+1]
      dat_remain = setdiff(1:sumn,dat_remove)
    }else{
      ind_remove = ((K-1)*fold+1):n
      ind_remain = setdiff(1:n,ind_remove)
      dat_remove = (ind[(K-1)*fold+1]+1):ind[length(ind)]
      dat_remain = setdiff(1:sumn,dat_remove)
    }
    
    
    ## Define training set
    n_train        = length(ind_remain)
    dataset_train  = data$dataset[dat_remain,]
    obs_times_train= data$obs_times[ind_remain]
    bas            = orthbasis(dataset_train[,2],knots = data$knots)
    B              = bas$B
    Omega          = bas$Omega
    train_data     = list( n         = n_train,
                           dataset   = dataset_train,
                           obs_times = obs_times_train,
                           B         = B,
                           Omega     = Omega)
    rm("n_train","dataset_train","obs_times_train","bas","B","Omega","ind_remain","dat_remain")
    para = minEM(train_data, lambda, ka,kb, tol = 1e-4, maxiter = 50)
    
    ## Testing Data Set
    n_test         = length(ind_remove)
    dataset_test   = data$dataset[dat_remove,]
    obs_times_test = data$obs_times[ind_remove]
    bas            = orthbasis(dataset_test[,2],knots = data$knots)
    B              = bas$B
    Omega          = bas$Omega
    test_data      = list( n         = n_test,
                           dataset   = dataset_test,
                           obs_times = obs_times_test,
                           B         = B,
                           Omega     = Omega)
    rm("n_test","dataset_test","obs_times_test","bas","B","Omega","ind_remove","dat_remove")
    
    ## Compute the Criteria Value
    value = value + MSE(test_data,para)
  }
  return(value)
}
