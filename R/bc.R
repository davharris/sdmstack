cor_from_upper = function(x){
  # length(x) == (dim * (dim - 1)) / 2
  dim = (sqrt(8 * length(x) + 1) + 1) / 2
  
  out = matrix(0, dim, dim)
  
  out[upper.tri(out)] = x    # Fill in upper triangle
  out = out + t(out)         # Fill in lower triangle
  diag(out) = 1              # Diagonal of correlation matrix is 1
  
  out
}


#' @importFrom BayesComm BC
#' @importFrom mvtnorm rmvnorm
#' @importFrom poibin dpoibin
#' @importFrom progress progress_bar
#' @export
bc_stack = function(p_train, y_train, p_test, its = 1000, thin = 10){
  covlist = lapply(1:n_spp, identity)
  bc = BC(
    Y = y_train, 
    X = qlogis(p_train), 
    model = "full", 
    covlist = covlist, 
    its = its, 
    thin = thin
  )
  
  thinned_its = nrow(bc$trace$R)
  
  correlations = lapply(
    1:nrow(bc$trace$R),
    function(i){
      cor_from_upper(bc$trace$R[i, ])
    }
  )
  
  test_mu_list = lapply(
    1:ncol(p_train),
    function(i) {
      cbind(1, qlogis(p_test)[ , i]) %*% t(bc$trace$B[[i]])
    }
  )
  test_mu = aperm(
    do.call(abind, c(test_mu_list, along=3)),
    c(1, 3, 2)
  )

  out = array(NA, dim(test_mu))
  for(i in 1:thinned_its){
    out[ , , i] = pnorm(
      test_mu[ , , i] + rmvnorm(nrow(test_mu), sigma = correlations[[i]]),
      0,
      1
    )
  }
  
    
  out
}
