#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::export]]
const List orth_algo(arma::mat Th, arma::mat V)
{
  vec eigval1;
  mat eigvec1;
  eig_sym(eigval1, eigvec1, V);
  eigval1 = reverse(eigval1);
  eigvec1 = fliplr(eigvec1);
  mat L = eigvec1*diagmat(sqrt(eigval1));
  mat A = Th*L;
  mat AtA = A.t()*A;
  vec eigval;
  mat eigvec;
  eig_sym(eigval, eigvec, AtA);
  eigval = reverse(eigval);
  eigvec = fliplr(eigvec);
  mat D = diagmat(eigval);
  mat temp = 1/sqrt(eigval);
  mat Q = A*eigvec*diagmat(1/sqrt(eigval));
  return List::create(Named("D")=D,
                      Named("Q")=Q);
}
