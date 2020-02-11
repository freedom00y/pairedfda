#' Create bsplines
#'
#' @param x the points you are interested in
#' @param knots if it is a number, the knots = seq(0,1,1/nknots); if it is a vector, knots is this vector
#' @param order default is 3, it means cubic spline
#'
#' @return a list including two elements. B(x) and covariance matrix of bsplines 
#' @import orthogonalsplinebasis 
#'
orthbasis <- function(x, knots = 8, order = 3) {
  #################################################
  # This code is to get the orthnormal basis matrix for given knots
  # knots for basis functions are equally spaced on [0,1]
  # input:
  #     x: a given vector between [0,1]
  # ouput:
  #     B: a basis matrix at the points x;
  
  len = length(knots)
  if(len == 1){
    innknots = seq(0,1,1/knots)
  }else{
    innknots = knots
  }
  
  knots = expand.knots(innknots,order+1) #set knots
  obase = OBasis(knots,order + 1)              ## get orthogonal basis

  B     = evaluate(obase, x)
  Omega = OuterProdSecondDerivative(obase)

  result = list( knots = innknots,
                 B = B,
                 Omega = Omega)
  return(result)
}
