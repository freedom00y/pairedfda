# pairedfda
Joint Modeling for Paired Functional Data

## Install packages:
* Check the operating system version, clang (for MacOS) or gcc (for Windows) version and gfortran version. They should match with each other so that you can compile rpp in R.
* Install devtools in R
You need the R package "devtools" to install "pairedfda". Use the following command to check if you have installed "devtools"
```r
require("devtools")
```
If you do not have it, install using the following command
```r
library(devtools)
install.packages("devtools")
```
* install rbpdfda in R:
```r
devtools::install_github("freedom00y/pairedfda")
```
## Use pairedfda package:
The following gives the description of each function in this package. The details can be found in the help documentation.

* predata: convert your raw paired data to the form this package can be used
* plt.data: plot the paired raw data
* Kf_CV: evaluate the mean square error using K-fold cross-validation
* minEM: get the estimation of the model parameters

