#include <RcppArmadillo.h>
#include "orth_algo.h"
using namespace Rcpp;
using namespace arma;
//[[Rcpp::depends(RcppArmadillo)]]
//' One EM iterate
//' @param oldpara parameter set from the last step
//' @param data postprocessed data
//' @param lambda tuning parameter
//' @return new parameter set
//' @keywords internal
//[[Rcpp::export]]
Rcpp::List EM(const Rcpp::List oldpara, const Rcpp::List data, const arma::vec lambda){
  // Parameters
  double eps0 = oldpara["sig_eps"];
  double xi0  = oldpara["sig_xi"];
  vec thmu0   = oldpara["theta_mu"];
  vec thnu0   = oldpara["theta_nu"];
  mat thf0    = oldpara["theta_f"];
  mat thg0    = oldpara["theta_g"];
  mat Da0     = oldpara["Da"];
  mat Db0     = oldpara["Db"];
  mat C0      = oldpara["C"];
  
  // Data
  int n = data["n"];
  mat dataset = data["dataset"];
  vec y=dataset.col(2);
  vec z=dataset.col(3);
  vec obs_times=data["obs_times"];
  mat B=data["B"];
  mat Omega=data["Omega"];
  int ncolB=B.n_cols;
  
  // E Step
  mat Bf = B*thf0;
  mat Bg = B*thg0;
  vec res_y = y - B*thmu0;
  vec res_z = z - B*thnu0;
  int ka = Da0.n_rows;
  int kb = Db0.n_rows;
  int kab= ka+kb;
  
  //Initialize sig_ab, declare (u, logu, delta2), (alpha, beta), bar_sig
  mat sig_ab = zeros(kab,kab);
  sig_ab(span(0,ka-1),span(0,ka-1))     = Da0;
  sig_ab(span(0,ka-1),span(ka,kab-1))   = C0;
  sig_ab(span(ka,kab-1),span(0,ka-1))   = C0.t();
  sig_ab(span(ka,kab-1),span(ka,kab-1)) = Db0;
  mat inv_sig_ab = inv(sig_ab); 
  
  vec  delta2 = zeros<vec>(n);
  mat  alpha  = zeros<mat>(n,ka);
  mat  beta   = zeros<mat>(n,kb);
  cube bar_sig(kab,kab,n);
  vec  temp0  = "0";
  vec  ind    = join_vert(temp0,cumsum(obs_times));

  
  // For each subject, compute their latent variables u,alpha,beta
  for(int i=0;i<n;i++)
  {
    int ni   = obs_times(i);
    int ind1 = ind(i);
    int ind2 = ind(i+1)-1;
    mat Bfi  = Bf.rows(ind1,ind2);
    mat Bgi  = Bg.rows(ind1,ind2);
    mat Bth  = zeros<mat>(2*ni,ka+kb);
    Bth(span(0,ni-1),span(0,ka-1))        = Bfi;
    Bth(span(ni,2*ni-1),span(ka,ka+kb-1)) = Bgi;
    vec resid  = join_cols(res_y.subvec(ind1,ind2),res_z.subvec(ind1,ind2));
    vec sigvec = join_cols(eps0*ones<vec>(ni),xi0*ones<vec>(ni) );
    mat sig_i  = zeros<mat>(2*ni,2*ni);
    sig_i(span(0,ni-1),span(0,ni-1))      = Bfi*Da0*Bfi.t() + eps0*eye(ni,ni);
    sig_i(span(0,ni-1),span(ni,2*ni-1))   = Bfi*C0*Bgi.t();
    sig_i(span(ni,2*ni-1),span(ni,2*ni-1))= Bgi*Db0*Bgi.t() + xi0*eye(ni,ni);
    sig_i(span(ni,2*ni-1),span(0,ni-1))   = Bgi*C0.t()*Bfi.t();
    mat sig_abi = sig_ab*Bth.t();
    
    // Compute inverse of Sigma_i
    // Sigma_i = Bth*Sigma_ab*Bth.t()+Sigma_sigma, where Sigma_sigma is a diagonal matrix
    // inverse(Sigma_i) = inv(Sigma_sigma) - inv(Sigma_sigma)*Bth* 
    // ( inv(Sigma_ab + Bth.t()*inv(Sigma_sigma)*Bth) *Bth.t()*inv(Sigma_sigma))
    mat inv_Ssig     = diagmat(1/sigvec);
    mat BthSs        = Bth.t()*inv_Ssig;
    mat inner        = BthSs*Bth+inv_sig_ab;
    mat inv_inner    = inv(inner);//good
    inv_inner.elem( find(abs(inv_inner)<1e-10) ).zeros();
    bar_sig.slice(i) = inv_inner;
    mat inv_Sigi     = inv_Ssig-BthSs.t()*bar_sig.slice(i)*BthSs;
    mat sig_inv      = inv_Sigi*resid;
    mat bar_ab       = sig_abi*sig_inv;
    alpha.row(i)     = bar_ab.rows(0,ka-1).t();
    beta.row(i)      = bar_ab.rows(ka,kab-1).t();
    delta2(i)        = sum(resid%sig_inv);
   
  }

  mat ab = join_rows(alpha,beta);
  
  // Center a and b
  rowvec a_mean = mean(alpha,0);
  rowvec b_mean = mean(beta,0);
  alpha = alpha.each_row()-a_mean;
  beta  = beta.each_row()-b_mean;
  
  // M step
  // Compute residual variance
  double sum_eps  = 0;
  double sum_xi   = 0;
  mat    sum_ubtb = zeros<mat>(ncolB,ncolB);
  vec    sum_uby  = zeros<vec>(ncolB);
  vec    sum_ubz  = zeros<vec>(ncolB);
  vec sum_btbfua  = zeros<vec>(ncolB);
  vec sum_btbgub  = zeros<vec>(ncolB);
  for(int i=0;i<n;i++)
  {
    int ind1 = ind(i);
    int ind2 = ind(i+1)-1;
    mat Bi   = B.rows(ind1,ind2);
    mat btb  = Bi.t()*Bi;
    mat Bfi  = Bf.rows(ind1,ind2);
    mat Bgi  = Bg.rows(ind1,ind2);
    vec yi   = y.subvec(ind1,ind2);
    vec zi   = z.subvec(ind1,ind2);
    vec res1 = yi-Bi*thmu0-Bfi*alpha.row(i).t();
    vec res2 = zi-Bi*thnu0-Bgi*beta.row(i).t();
    mat sigaa= bar_sig(0,0,i,size(ka,ka,1));
    mat sigbb= bar_sig(ka,ka,i,size(kb,kb,1));
    mat A1=Bfi*sigaa*Bfi.t();
    double val1 = sum( A1.diag() ) + sum(res1%res1);
    mat A2=Bgi*sigbb*Bgi.t();
    double val2 = sum( A2.diag() ) + sum(res2%res2);
    sum_eps += val1;
    sum_xi  += val2;
    sum_ubtb+= btb;
    sum_uby += Bi.t()*yi;
    sum_ubz += Bi.t()*zi;
    sum_btbfua += btb*thf0*alpha.row(i).t();
    sum_btbgub += btb*thg0*beta.row(i).t();
  }
  double eps1 = sum_eps/sum(obs_times);
  double xi1  = sum_xi/sum(obs_times);
  
  // Compute theta_mu, theta_nu
  vec thmu1   = solve(sum_ubtb/n+lambda[0]*eps1*Omega,sum_uby/n-sum_btbfua/n);
  vec thnu1   = solve(sum_ubtb/n+lambda[1]*xi1 *Omega,sum_ubz/n-sum_btbgub/n);
  
  // Compute theta_f, theta_g
  cube sum_f1 = zeros<cube>(ncolB,ncolB,ka);
  cube sum_f2 = zeros<cube>(ncolB,1,ka);
  cube sum_g1 = zeros<cube>(ncolB,ncolB,kb);
  cube sum_g2 = zeros<cube>(ncolB,1,kb);
  
  for(int i=0;i<n;i++)
  {
    int ind1 = ind(i);
    int ind2 = ind(i+1)-1;
    mat Bi   = B.rows(ind1,ind2);
    mat btb  = Bi.t()*Bi;
    vec yi   = y.subvec(ind1,ind2);
    vec zi   = z.subvec(ind1,ind2);
    
    mat uab = ab.row(i).t()*ab.row(i) + bar_sig.slice(i);
    mat uaa = uab(span(0,ka-1),span(0,ka-1));
    mat ubb = uab(span(ka,kab-1),span(ka,kab-1));
    for(int j=0;j<ka;j++)
    {
      vec tfua = thf0*uaa.col(j);
      sum_f1.slice(j) += uaa(j,j)*btb;
      sum_f2.slice(j) += alpha(i,j)*Bi.t()*(yi-Bi*thmu1)+uaa(j,j)*btb*thf0.col(j)-btb*tfua;
    }
    
    for(int k=0;k<kb;k++)
    {
      vec tgub = thg0*ubb.col(k);
      sum_g1.slice(k) += ubb(k,k)*btb;
      sum_g2.slice(k) += beta(i,k)*Bi.t()*(zi-Bi*thnu1)+ubb(k,k)*btb*thg0.col(k)-btb*tgub;
    }
  }
    
  mat thf1(ncolB,ka);
  for(int j=0;j<ka;j++)
  {
    mat temp = sum_f1.slice(j)/n+eps1*lambda[2]*Omega;
    vec la;
    mat V;
    eig_sym(la, V, temp);
    thf1.col(j) = V*diagmat(1/la)*V.t()*sum_f2.slice(j)/n;
  }
  
  mat thg1(ncolB,kb);
  for(int k=0;k<kb;k++)
  {
    mat temp = sum_g1.slice(k)/n+xi1*lambda[3]*Omega;
    vec la;
    mat V;
    eig_sym(la, V, temp);
    thg1.col(k) = V*diagmat(1/la)*V.t()*sum_g2.slice(k)/n;
  }
  
  // Orthogonal
  mat sum_bar = sum(bar_sig,2);
  mat hat_sigma = (ab.t()*ab+sum_bar)/n;
  mat Vaa = hat_sigma(span(0,ka-1),span(0,ka-1));
  mat Vbb = hat_sigma(span(ka,kab-1),span(ka,kab-1));
  mat Vab = hat_sigma(span(0,ka-1),span(ka,kab-1));
  List temp1 = orth_algo(thf1,Vaa);
  mat qf = temp1["Q"];
  mat Da1 = temp1["D"];
  List temp2 = orth_algo(thg1,Vbb);
  mat qg = temp2["Q"];
  mat Db1 = temp2["D"];
  mat C1  = qf.t()*thf1*Vab*thg1.t()*qg;
  mat after = zeros(kab,kab);
  after(span(0,ka-1),span(0,ka-1))     = Da1;
  after(span(0,ka-1),span(ka,kab-1))   = C1;
  after(span(ka,kab-1),span(0,ka-1))   = C1.t();
  after(span(ka,kab-1),span(ka,kab-1)) = Db1;
  
  // Change sign in order to promise the second element in each column is positive
  thf1=qf*diagmat(sign(qf.row(1)));
  thg1=qg*diagmat(sign(qg.row(1)));

  return List::create(Named("sig_eps") = eps1,
                      Named("sig_xi")  = xi1,
                      Named("theta_mu")= thmu1,
                      Named("theta_nu")= thnu1,
                      Named("theta_f") = thf1,
                      Named("theta_g") = thg1,
                      Named("Da")      = Da1,
                      Named("Db")      = Db1,
                      Named("C")       = C1,
                      Named("alpha")   = alpha,
                      Named("beta")    = beta
  );
}
