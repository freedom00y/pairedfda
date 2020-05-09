#include <RcppArmadillo.h>
#include "EMinner.hpp"
#include "loglike.hpp"

using namespace Rcpp;
using namespace arma;
//[[Rcpp::depends(RcppArmadillo)]]
//' Estimations
//' 
//' Use this function to estimate parameters set 
//' 
//' @param data Processed data. Use "predata" to preprocess the raw data first.
//' @param lambda Tuning parameter, a vector with 4 components. The tuning parameters shows as the following order, the first mean curve \eqn{(\lambda_\mu)}, the second mean curve \eqn{(\lambda_\nu)}, the first pcs \eqn{(\lambda_f)}, the second pcs \eqn{(\lambda_g)}.
//' @param tol Tolerance of the EM algorithm
//' @param ka Number of pcs for the first reponse variable
//' @param kb Number of pcs for the second reponse variable
//' @param maxiter Maximum iteration time, the default is 100 times
//' 
//' @return the estimation of the parameters
//' 
//' @details  We suppose the model is 
//' \deqn{Y_i = B_i \theta_\mu + B_i f \alpha_i + \epsilon_i, Z_i = B_i \theta_\nu + B_i g \beta_i + \xi_i,}
//' where \eqn{(\alpha_i, \beta_i)} and residuals follow normal distribution. We denote that 
//' \deqn{(\alpha_i, \beta_i) \sim N(0,\Sigma_{\alpha\beta}), \epsilon_i\sim N(0,\sigma_\epsilon^2), \xi_i\sim N(0,\sigma_\xi^2),}
//' and \eqn{\Sigma_{\alpha\beta} = (D_a & C\\C^T & Db)}.
//' This function returns the estimation of \eqn{\sigma^2_\epsilon,\sigma^2_\xi,\theta_\mu,\theta_\nu,\theta_f,\theta_g,Da,Db,C}.
//' 
//' @examples 
//' rawdata = gen_data(n=50)
//' visit = seq(0,100,20)
//' data = predata(nobs = rawdata$nobs, 
//'                time = rawdata$time,  
//'                y = rawdata$y, 
//'                z = rawdata$z, 
//'                knots = visit, 
//'                order = 3)
//' ## without penalty
//' lambda = c(0,0,0,0)
//' pt_nopen = minEM(data, lambda, ka=1, kb=2, tol = 1e-4, maxiter = 100)
//' ## with penalty
//' lambda_pen = c(6000,6000,16000,16000,16000)
//' pt_pen = minEM(data, lambda_pen, ka=1, kb=2, tol = 1e-4, maxiter = 100)
//' 
//' @export
//[[Rcpp::export]]
const List minEM(const List data, const arma::vec lambda, const int ka, const int kb, const double tol=1e-4, int maxiter = 100)
{
  // Data
  mat B       = data["B"];
  int ncolB   = B.n_cols;
  
  // Initialize Parameters
  vec theta_mu0 = ones<vec>(ncolB);
  vec theta_nu0 = ones<vec>(ncolB);
  mat X(ncolB,ncolB,fill::eye);
  mat theta_f0  = X.head_cols(ka);
  mat theta_g0  = X.head_cols(kb);
  double sig_eps0 = 1;
  double sig_xi0 = 1;
  mat Da0(ka,ka,fill::eye);
  mat Db0(kb,kb,fill::eye);
  mat C0(ka,kb,fill::zeros);
  vec llv = zeros<vec>(maxiter);
  
  List para0 = List::create(Named("sig_eps")  = sig_eps0,
                            Named("sig_xi")   = sig_xi0,
                            Named("theta_mu") = theta_mu0,
                            Named("theta_nu") = theta_nu0,
                            Named("theta_f")  = theta_f0,
                            Named("theta_g")  = theta_g0,
                            Named("Da")       = Da0,
                            Named("Db")       = Db0,
                            Named("C")        = C0);
  
  // EM steps
  List para1    = EM(para0,data,lambda);
  double value0 = loglike(data,para0);
  double value1 = loglike(data,para1);
  double differ = 1.0-value1/value0;
  llv(0) = value0;
  
  int iter = 1;
  llv(iter) = value1;
  while( iter < maxiter-1 & fabs(differ)>tol)
  {
    para0  = para1;
    value0 = value1;
    para1  = EM(para0,data,lambda);
    value1 = loglike(data,para1);
    differ = (value0-value1)/fabs(value0);
    if(differ<-0.05)
    {
      para1 = para0;
      value1 = value0;
      break;
    }
    iter   = iter+1;
    llv(iter) = value1;
  }
  para1["-2ll"] = llv.head(iter+1);
  para1["iter"] = iter+1;
  
  return para1;
}
