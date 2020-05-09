#' Preprocessing data
#' 
#' @description Change the raw data to a sturcture that can be used in this package. Create the design matrix and penalty matrix.
#'
#' @param nobs A vector whose component is the number of the observations of each subject
#' @param time A vector records the observation time of all subjects
#' @param y A vector records first response variable of all subjects
#' @param z A vector records second response variable of all subjects
#' @param knots The full set of knots used to define the basis functions. It can be a number or a vector. If it is a number, the knots = seq(0,1,1/knots); if it is a vector, knots is this vector.
#' @param order Order of bsplines
#'
#' @return 
#' A list including number of subjects (n), obs times (nobs), paired data (dataset), the full set of knots (knots), design matrix (B) and penalty matrix (Omega) 
#' 
#' @details Suppose the data include the information of n subjects. 
#' @details For each subject, there are m_i (i=1,...n) observations which are observed at time \eqn{t_i =  (t_{i1},...t_{i m_i})}.
#' @details Suppose the data are generated from potential curves Y_i and Z_i. Then 
#' @details nobs = (m_1,...,m_n), time = (t_1,...t_n), y = (Y_1(t_1),...,Y_n(t_n)), z = (Z_1(t_1),...,Z_n(t_n)).
#' 
#' @export
#'
#' @examples
#' n = 5
#' nobs = 1+rbinom(n,5,0.9)
#' time=c()
#' for(i in 1:n)
#' {
#'   time=c(time,0,runif(nobs[i]-1))
#' }
#' sumobs = sum(nobs)
#' y = rnorm(sumobs,0,1.5)
#' z = rnorm(sumobs,1,1)
#' data = predata(nobs,time,y,z,knots = 8,order=3)
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

predata <- function(nobs,time,y,z,knots = 8, order =3)
{
  nsum = sum(nobs)
  if(length(y)!=nsum | length(z)!=nsum | length(time)!=nsum )
  {
    print("The obs times and the length of data don't match!")
    break
  }
  n = length(nobs)
  ind = rep(1:n-1,nobs)
  dataset = cbind(ind,time,y,z)
  colnames(dataset)  = c("subject","obs_time","y_t","z_t")
  bas        = orthbasis(time,knots,order)
  
  data = list( n = n,
               obs_times = nobs,
               dataset   = dataset,
               knots     = bas$knots,
               B         = bas$B,
               Omega     = bas$Omega)

  return(data)
}