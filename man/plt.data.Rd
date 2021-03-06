% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plt.data.R
\name{plt.data}
\alias{plt.data}
\title{Plot the raw data}
\usage{
plt.data(data)
}
\arguments{
\item{data}{Processed data. Use "predata" to preprocess the raw data first.}
}
\description{
display two plots of the paired raw data
}
\examples{
rawdata = gen_data(n=50)
visit = seq(0,100,20)
data = predata(nobs = rawdata$nobs, 
               time = rawdata$time,  
               y = rawdata$y, 
               z = rawdata$z, 
               knots = visit, 
               order = 3)
plt.data(data)
}
