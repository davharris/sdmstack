richness_pdf = function(p){
  dpoibin(kk = 0:length(p), pp = p)
}

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


#' @importFrom mvtnorm rmvnorm
#' @importFrom poibin dpoibin
#' @importFrom progress progress_bar
#' @importFrom abind abind
#' @export
bc_stack = function(p_train, y_train, p_test, iters = 1010, burn = 10, thin = 2, verbose = 1, ...){

  # Initial values
  z = y_train - predicted_p_train
  R = cor(z)

  # Multiply by sqrt(2)??
  train_mu = probit(p_train) * sqrt(2)
  test_mu  = probit(p_test) * sqrt(2)

  bc = BayesComm::BCfit(
    y = y_train,
    X = NULL,
    covlist = NULL,
    R = R,
    z = z,
    mu = train_mu,
    iters = iters,
    burn = burn,
    thin = thin,
    verbose = verbose,
    updateR = TRUE,
    updateMu = FALSE
  )


  thinned_its = nrow(bc$R)

  predicted_array = array(NA, c(dim(p_test), thinned_its))
  rho_array = array(NA, c(dim(R), thinned_its))

  pb = progress_bar$new(
    format = "Species-level predictions [:bar] :percent",
    total = thinned_its, clear = FALSE, width = 70
  )

  for (i in 1:thinned_its) {

    pb$tick()

    # Sample a correlation matrix
    rho = cor_from_upper(bc$R[i, ])


    # Samples from latent Gaussian, using rho as our covariance
    z = test_mu + rmvnorm(nrow(test_mu), sigma = rho)

    # Save rho
    rho_array[ , , i] = rho

    # Save conditional occurrence probabilities
    predicted_array[ , , i] = pnorm(z, 0, 1)
  }


  pb = progress_bar$new(
    format = "Richness distributions [:bar] :percent",
    total = nrow(p_test), clear = FALSE, width = 70
  )
  rich_p = t(
    apply(
      predicted_array,
      1,
      function(site_predictions){
        pb$tick()

        # Loop through the predictions for the site and average their
        # probability distribution functions for richness
        rowMeans(
          apply(site_predictions, 2, richness_pdf)
        )
      }
    )
  )
  dimnames(rich_p)[1:2] = dimnames(predicted_p_test)



  list(
    rho_array = rho_array,
    predicted_array = predicted_array,
    rich_p = rich_p
  )
}
