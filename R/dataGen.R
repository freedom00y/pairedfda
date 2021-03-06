#' Generate Data
#' @description Generate the simulation data described in the paper.
#' @param n number of subjects
#'
#' @return a list of 4 components including nobs, time, y and z. The meaning of the 4 parameters can be found from \code{\link{predata}}.
#' @importFrom MASS mvrnorm stats
#' @export
#'
#' @examples
#' rawdata = gen_data(50)
gen_data <- function(n){
  Da = 36
  Db = diag(c(36, 16))
  C  = c(-28.8, -10.8)
  Sigma = rbind( cbind(Da,t(C)), cbind(C,Db) )
  score = mvrnorm(n,rep(0,3),Sigma)
  
  nobs=rep(0,n)
  time=c()
  y=c()
  z=c()
  for(i in 1:n){
    tt = simu_t()
    nobs[i] = length(tt)
    time = c(time,tt)
    yi = mu_t(tt) + score[i,1]*fy_t(tt)
    y = c(y,yi)
    zi = nu_t(tt) + score[i,2]*fz1_t(tt) + score[i,3]*fz2_t(tt)
    z = c(z,zi)
  }
  y = y + rnorm(length(y),0,0.5)
  z = z + rnorm(length(z),0,0.5)
  
  return(list("nobs"=nobs, "time"=time, "y"=y, "z"=z))
}

#' Simulation function mu(t)
#' mean function of the first functional variable
#' @param t covariate
#' @return mean function mu(t)
mu_t <- function(t){
  len = length(t)
  return( rep(1,len) + t/100 + exp(-(t-60)**2/500) ) 
}

#' Simulation function nu(t)
#' mean function of the second functional variable
#' @param t covariate
#' @return mean function nu(t)
nu_t <- function(t){
  len = length(t)
  return( rep(1,len) - t/100 - exp(-(t-30)**2/500) ) 
}

#' Simulation function 1st pc of Y
#' 1st pc of Y
#' @param t covariate
#' @return f_1(t)
fy_t <- function(t){
  sin(2*pi*t/100)/sqrt(50)
}

#' Simulation function 1st pc of Z
#' 1st pc of Z
#' @param t covariate
#' @return g_1(t)
fz1_t <- function(t){
  fy_t(t)
}

#' Simulation function 2nd pc of Z
#' 2nd pc of Z
#' @param t covariate
#' @return g_2(t)
fz2_t <- function(t){
  cos(2*pi*t/100)/sqrt(50)
}

#' One simulation attempt
simu_t <- function(){
  res = c(0)
  k=1
  while(k<=4){
    cur = res[k]+rnorm(1,mean=30,sd=10)
    if(cur<=100){
      res=c(res,cur)
      k = k+1
    }else{
      break
    }
  }
  return(res)
}

