#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
//[[Rcpp::depends(RcppArmadillo)]]
//' Evaluate the loglikelihood value
//' @param data postprocessed data
//' @param para parameter set
//' @return -2*log-likelihood value
//' @keywords internal
//[[Rcpp::export]]
double loglike(const List data, const List para)
{
  // Data
  int n       = data["n"];
  mat dataset = data["dataset"];
  vec y       = dataset.col(2);
  vec z       = dataset.col(3);
  vec obs     = data["obs_times"];
  mat B       = data["B"];
  
  // Parameter
  vec    mu   = para["theta_mu"];
  vec    nu   = para["theta_nu"];
  mat    f    = para["theta_f"];
  mat    g    = para["theta_g"];
  mat    Da   = para["Da"];
  mat    Db   = para["Db"];
  mat    C    = para["C"];
  double eps  = para["sig_eps"];
  double xi   = para["sig_xi"];
  
  int ka = Da.n_rows;
  int kb = Db.n_rows;
  int kab= ka+kb;
  double total = 0;
  vec temp0 = "0";
  vec ind = join_vert(temp0,cumsum(obs));
  mat Bf = B*f;
  mat Bg = B*g;
  vec res_y = y - B*mu;
  vec res_z = z - B*nu;
  
  mat sig_ab = zeros(kab,kab);
  sig_ab(span(0,ka-1),span(0,ka-1))     = Da;
  sig_ab(span(0,ka-1),span(ka,kab-1))   = C;
  sig_ab(span(ka,kab-1),span(0,ka-1))   = C.t();
  sig_ab(span(ka,kab-1),span(ka,kab-1)) = Db;
  mat inv_sig_ab = inv(sig_ab); 
  
  
  for(int i=0;i<n;i++)
  {
    int ni  = obs(i);
    int ind1 = ind(i);
    int ind2 = ind(i+1)-1;
    mat Bfi  = Bf.rows(ind1,ind2);
    mat Bgi  = Bg.rows(ind1,ind2);
    mat Bth  = zeros<mat>(2*ni,ka+kb);
    Bth(span(0,ni-1),span(0,ka-1))        = Bfi;
    Bth(span(ni,2*ni-1),span(ka,ka+kb-1)) = Bgi;
    vec resid  = join_cols(res_y.subvec(ind1,ind2),res_z.subvec(ind1,ind2));
    vec sigvec = join_cols(eps*ones<vec>(ni),xi*ones<vec>(ni) );
    
    mat sig_i  = zeros<mat>(2*ni,2*ni);
    sig_i(span(0,ni-1),span(0,ni-1))      = Bfi*Da*Bfi.t() + eps*eye(ni,ni);
    sig_i(span(0,ni-1),span(ni,2*ni-1))   = Bfi*C*Bgi.t();
    sig_i(span(ni,2*ni-1),span(ni,2*ni-1))= Bgi*Db*Bgi.t() + xi*eye(ni,ni);
    sig_i(span(ni,2*ni-1),span(0,ni-1))   = Bgi*C.t()*Bfi.t();
    
    mat inv_Ssig     = diagmat(1/sigvec);
    mat BthSs        = Bth.t()*inv_Ssig;
    mat inner        = BthSs*Bth+inv_sig_ab;
    mat inv_inner    = inv(inner);
    inv_inner.elem( find(abs(inv_inner)<1e-10) ).zeros();
    mat inv_Sigi     = inv_Ssig-BthSs.t()*inv_inner*BthSs;
    mat sig_inv      = inv_Sigi*resid;
    double delta2    = sum(resid%sig_inv);
    total += log(det(sig_i))+delta2+2*ni*log(2*datum::pi);
  }
  return total;
}
