\name{pairedfda-package}
\alias{pairedfda-package}
\alias{pairedfda}
\docType{package}
\title{
\packageTitle{pairedfda}
}
\description{
\packageDescription{pairedfda}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{pairedfda}
\packageIndices{pairedfda}
~~ An overview of how to use the package, including the most important ~~
~~ functions ~~
}
\author{
\packageAuthor{pairedfda}

Maintainer: \packageMaintainer{pairedfda}
}
\references{
Zhou, L., Huang, J. Z., & Carroll, R. J. (2008). Joint modelling of paired sparse functional data using principal components. Biometrika, 95(3), 601-619.
}

\keyword{ package }
\seealso{
\code{\link{predata}, \link{plt.data}, \link{Kf_CV}, \link{minEM}}
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
## without penalty
lambda = c(0,0,0,0)
pt_nopen = minEM(data, lambda, ka=1, kb=2, tol = 1e-4, maxiter = 100)
## with penalty
lambda_pen = c(6000,6000,16000,16000,16000)
pt_pen = minEM(data, lambda_pen, ka=1, kb=2, tol = 1e-4, maxiter = 100)
}
