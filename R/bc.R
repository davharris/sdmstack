probit = qnorm
inv_probit = pnorm

cor_from_upper = function(x){
  # Turn a numeric vector with upper triangle values and
  # complete the correlation matrix with this information

  # Dimension of the correlation matrix
  # dim defined by solving length(x) == (dim * (dim - 1)) / 2
  dim = (sqrt(8 * length(x) + 1) + 1) / 2

  out = matrix(0, dim, dim)  # Initialize
  out[upper.tri(out)] = x    # Fill in upper triangle
  out = out + t(out)         # Fill in lower triangle
  diag(out) = 1              # Diagonal of correlation matrix is 1

  out
}


#' @importFrom BayesComm BC
#' @importFrom mvtnorm rmvnorm
#' @importFrom poibin dpoibin
#' @importFrom progress progress_bar
#' @importFrom abind abind
#' @export
bc_stack = function(p_train, y_train, p_test, its = 1100, burn = 100, thin = 10, verbose = 1, ...){

  # Each species gets its covariate from the corresponding column
  # of probit(p_train)
  covlist = lapply(1:n_spp, identity)


  # 2 / sqrt(2) probably needed because probit model
  # doubles the variance?
  bc = BC(
    Y = y_train,
    X = sqrt(2) / 2 * probit(p_train),
    model = "full",
    covlist = covlist,
    its = its,
    thin = thin,
    verbose = verbose,
    ...
  )

  thinned_its = nrow(bc$trace$R)

  residual_cor_list = lapply(
    1:thinned_its,
    function(i){
      cor_from_upper(bc$trace$R[i, ])
    }
  )

  # Predict mean probit values for each species at test locations
  test_mu_list = lapply(
    1:ncol(p_train),
    function(i) {
      # First matrix has ones for intercept,
      cbind(1, probit(p_test)[ , i]) %*% t(bc$trace$B[[i]])
    }
  )
  test_mu = aperm(
    do.call(abind, c(test_mu_list, along = 3)),
    c(1, 3, 2)
  )

  predicted_array = array(NA, dim(test_mu))
  for (i in 1:thinned_its) {
    predicted_array[ , , i] = pnorm(
      test_mu[ , , i] + rmvnorm(nrow(test_mu), sigma = residual_cor_list[[i]]),
      0,
      1
    )
  }

  list(
    residual_cor_list = residual_cor_list,
    predicted_array = predicted_array,
    coefficients = bc$trace$B
  )
}
