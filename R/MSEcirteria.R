MSE <- function(data,para)
{
  ## Load dataset
  n = data$n
  dataset=data$dataset
  times = data$obs_times
  y = dataset[,"y_t"]
  z = dataset[,"z_t"]
  B = data$B
  
  
  ## Load parameters
  mu = para$theta_mu
  nu = para$theta_nu
  f  = para$theta_f
  g  = para$theta_g
  ka = ncol(f)
  kb = ncol(g)
  Da = para$Da
  Db = para$Db
  C  = para$C
  eps = para$sig_eps
  xi  = para$sig_xi
  n_total = nrow(B)
  
  kab=ka+kb
  ind_ka = as.vector(1:ka)
  ind_kb = as.vector((ka+1):kab)
  sig_ab = matrix(nrow = kab, ncol = kab)
  sig_ab[ind_ka, ind_ka] = Da
  sig_ab[ind_ka, ind_kb] = C
  sig_ab[ind_kb, ind_ka] = t(C)
  sig_ab[ind_kb, ind_kb] = Db
  eig = eigen(sig_ab)
  inv_sig_ab = eig$vectors%*%diag(1/eig$values)%*%t(eig$vectors)
  
  
  ## Estimate alpha & beta for test curves
  Bf    = B%*%f
  Bg    = B%*%g
  res_y = y - B%*%mu
  res_z = z - B%*%nu
  alpha = matrix(nrow = n, ncol = ka)
  beta  = matrix(nrow = n, ncol = kb)
  
  ind = c(0, cumsum(times))                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   
  for (i in 1 : n )
  {
    ni   = times[i]
    ind1 = ind[i] + 1
    ind2 = ind[i+1]
    Bfi  = as.matrix(Bf[ind1:ind2,])
    Bgi  = as.matrix(Bg[ind1:ind2,])
    inv_sig = diag(rep(c(1/eps,1/xi),each=ni))
    ind3    = 1:ni;
    ind4    = (ni+1):(2*ni);
    Bth              = matrix(0,ncol = (ka+kb), nrow = 2*ni)
    Bth[ind3,ind_ka] = Bfi
    Bth[ind4,ind_kb] = Bgi
    resid     = c(res_y[ind1:ind2],res_z[ind1:ind2])
    
    inner     = inv_sig_ab + t(Bth)%*%inv_sig%*%Bth
    eig=eigen(inner)
    inv_inner = eig$vectors%*%diag(1/eig$values)%*%t(eig$vectors)
    inv_sigi  = inv_sig-inv_sig%*%Bth%*%inv_inner%*%t(Bth)%*%inv_sig
    bar_ab    = sig_ab%*% t(Bth) %*%inv_sigi%*%resid
    alpha[i,] = bar_ab[ind_ka]
    beta[i,]  = bar_ab[ind_kb]
  }
  
  
  ## Create Augmented matrix Alpha and Beta
  Alpha = matrix(ncol = ka, nrow = n_total)
  for(i in 1:ka)
  {
    Alpha[,i]=rep(alpha[,i],times)
  }
  
  Beta = matrix(ncol = kb, nrow = n_total)
  for(i in 1:kb)
  {
    Beta[,i]=rep(beta[,i],times)
  }
  
  Bmu  = B%*%mu
  aBf  = apply(Alpha*(B%*%f),1,sum)
  yres = Bmu+aBf-y 
  sy2  = sum(yres^2)
  
  Bnu  = B%*%nu
  bBg  = apply(Beta*(B%*%g),1,sum)
  zres = Bnu+bBg-z 
  sz2  = sum(zres^2)
  
  return(sy2+sz2)
}