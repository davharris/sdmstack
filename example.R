library(poibin)
devtools::load_all()

set.seed(2)

n_spp = 10
n_loc = 10000
n_obs_x = 5

true_intercepts = rnorm(n_spp, -1)

full_x_both = replicate(
  2,
  sapply(
    1:n_spp,
    function(i){
      rnorm(n_loc, sd = sqrt(0.25 / i))
    }
  ),
  simplify = FALSE
)
full_x_train = full_x_both[[1]]
full_x_test = full_x_both[[2]]


pca = prcomp(full_x_train, scale = FALSE)
observed_x_train = predict(pca, full_x_train)[ , seq(2, n_obs_x + 1)]
observed_x_test = predict(pca, full_x_test)[ , seq(2, n_obs_x + 1)]


true_beta = matrix(rnorm(n_spp^2, mean = 1), nrow = n_spp, ncol = n_spp)
true_intercept = -1

true_p_train = inv_probit(full_x_train %*% true_beta + true_intercept)
true_p_test = inv_probit(full_x_test %*% true_beta + true_intercept)

y_train = structure(
  rbinom(length(true_p_train), size = 1, prob = true_p_train),
  dim = dim(true_p_train)
)
y_test = structure(
  rbinom(length(true_p_test), size = 1, prob = true_p_test),
  dim = dim(true_p_test)
)


glms = lapply(
  1:n_spp,
  function(i){
    glm(y_train[,i] ~ .,
        data = as.data.frame(observed_x_train),
        family = binomial)
  }
)

predicted_p_train = sapply(glms, predict, type = "response")
predicted_p_test = sapply(glms, predict, newdata = as.data.frame(observed_x_test), type = "response")


bc = bc_stack(
  p_train = predicted_p_train,
  y_train = y_train,
  p_test = predicted_p_test,
  5010, burn = 10, thin = 10
)

stack = bc$predicted_array




marginal_p = marginalize_out(stack, 3)

plot(
  marginal_p,
  predicted_p_test,
  pch = ".",
  main = "BC stacking introduces no biases",
  xlab = "Expected occurrence probability (post-stacking)",
  ylab = "Expected occurrence probability (pre-stacking)"
)
abline(0,1)

par(mfrow = c(2, 1))
plot(probit(bc$predicted_array[1, 2, ]), type = "l", main = "Great mixing for predictions")
plot(bc$rho_array[1, 2, ], type = "l", main = "Good mixing for correlations")
par(mfrow = c(1, 1))

rich_p_naive = t(apply(predicted_p_test, 1, function(x){dpoibin(0:n_spp, x)}))
rich_p = bc$rich_p


mean(log(rich_p[cbind(1:nrow(y_test), rowSums(y_test) + 1)])) # Full
mean(log(rich_p_naive[cbind(1:nrow(y_test), rowSums(y_test) + 1)])) # Naive
log(1 / (n_spp + 1)) # Very naive


par(mfrow = c(1, 2))
plot(seq(0, 1, length = nrow(y_test)), sort(sapply(1:nrow(y_test), function(i){sum(rich_p_naive[i, 1:sum(y_test[i, ])])})),
     main = "Naive stacking yields many observations in the tails of the distribution")
abline(0,1)
plot(seq(0, 1, length = nrow(y_test)), sort(sapply(1:nrow(y_test), function(i){sum(rich_p[i, 1:sum(y_test[i, ])])})),
     main = "BC stacking is less extreme, but not optimal")
abline(0,1)
par(mfrow = c(1, 1))



i = sample.int(n_spp, 1)
plot(-.075 + seq(0, n_spp), rich_p_naive[i, ], type = "h", bty = "l", yaxs = "i", lwd = 2)
lines(.055 + seq(0, n_spp), rich_p[i, ], type = "h", col = 2, lwd = 2)
points(sum(y_test[i, ]), 0, pch = 16, cex = 2)



# analyses ----------------------------------------------------------------

true_rho = cov2cor(cov(probit(true_p_test) - probit(predicted_p_test)))
est_rho = marginalize_out(bc$rho_array, 3)
plot(true_rho ~ est_rho)
abline(0,1)
