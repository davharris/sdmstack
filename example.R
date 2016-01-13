library(mvtnorm)
library(poibin)
devtools::load_all()

set.seed(2)

n_spp = 50
n_loc = 1000
n_obs_x = 5

true_rho = matrix(.7, n_spp, n_spp)
diag(true_rho) = 1

true_intercepts = rnorm(n_spp, -2)

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
  its = 500
)
stack = bc$predicted_array

str(stack)

marginal_p = marginalize_out(stack, 3)
plot(marginal_p, predicted_p_test)

rich_p = t(
  apply(
    stack,
    1,
    function(x){
      rowMeans(
        apply(
          x,
          2,
          function(p){
            dpoibin(0:n_spp, p)
          }
        )
      )
    }
  )
)
dimnames(rich_p)[1:2] = dimnames(predicted_p_test)

rich_p_naive = t(apply(predicted_p_test, 1, function(x){dpoibin(0:n_spp, x)}))



sum(log(rich_p[cbind(1:nrow(y_test), rowSums(y_test) + 1)])) # Full
sum(log(rich_p_naive[cbind(1:nrow(y_test), rowSums(y_test) + 1)])) # Naive

i = sample.int(n_spp, 1)
plot(-.05 + seq(0, n_spp), rich_p_naive[i, ], type = "h", bty = "l", yaxs = "i", lwd = 2)
lines(.05 + seq(0, n_spp), rich_p[i, ], type = "h", col = 2, lwd = 2)
points(sum(y_test[i, ]), 0, pch = 16, cex = 2)


