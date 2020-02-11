#' Plot the raw data
#' @description display two plots of the paired raw data
#' @param data Processed data. Use "predata" to preprocess the raw data first.
#'
#'
#' @examples
#' plt.data(data)
plt.data <- function(data)
{
  n = data$n
  nobs = data$obs_times
  ind = c(0,cumsum(nobs))
  color = rainbow(n)
  par(mfrow=c(1,2))
  t = data$dataset[,2]
  y = data$dataset[,3]
  z = data$dataset[,4]
  
  indrange = (ind[1]+1):ind[2]
  plot(t[indrange],y[indrange],'l',col=color[1],ylim=range(y), xlab='Time', ylab='1st Response Variable')
  for(i in 2:n)
  {
    indrange = (ind[i]+1):ind[i+1]
    lines(t[indrange],y[indrange],'l',col=color[i])
  }
  
  indrange = (ind[1]+1):ind[2]
  plot(t[indrange],z[indrange],'l',col=color[1],ylim=range(z), xlab='Time', ylab='2nd Response Variable')
  for(i in 2:n)
  {
    indrange = (ind[i]+1):ind[i+1]
    lines(t[indrange],z[indrange],'l',col=color[i])
  }
  par(mfrow=c(1,1))
}