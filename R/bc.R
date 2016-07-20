richness_pdf = function(p){
  dpoibin(kk = 0:length(p), pp = p)
}

probit = qnorm
inv_probit = pnorm


#' @importFrom mvtnorm rmvnorm
#' @importFrom poibin dpoibin
#' @importFrom progress progress_bar
#' @export
bc_stack = function(p_train, y_train, p_test, iters = 1010, burn = 10, thin = 2,
                    verbose = 1, ...){

  # Initial values
  z = y_train - p_train
  R = cor(z)

  train_mu = probit(p_train)
  test_mu  = probit(p_test)

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
    updateMu = FALSE,
    priW = c(nrow(z) + ncol(z) + 1, ncol(z) + 1)
  )


  rich_sims = sapply(
    1:nrow(bc$R),
    function(i){
      # For the current correlation matrix, make a richness prediction at
      # each site
      R <- BayesComm:::upper2cor(bc$R[i, ])
      z <- test_mu + mvtnorm::rmvnorm(nrow(test_mu), sigma = R)
      rowSums(z > 0)
    }
  )

  list(
    bc_model = bc,
    sims = rich_sims
  )
}
