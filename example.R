library(poibin)
library(arm)

devtools::load_all()

set.seed(2)

n_spp = 300
n_loc = 500

true_intercepts = rnorm(n_spp, -1)

full_x_both = replicate(
  2,
  sapply(
    1:n_spp,
    function(i){
      rt(n_loc, df = 8) / (2 + i)
    }
  ),
  simplify = FALSE
)
full_x_train = full_x_both[[1]]
full_x_test = full_x_both[[2]]


pca = prcomp(full_x_train, scale = FALSE)
variances = pca$sdev^2 / sum(pca$sdev^2)
included = 2:which(cumsum(variances[-1]) > 0.7)[1]
observed_x_train = predict(pca, full_x_train)[ ,included]
observed_x_test = predict(pca, full_x_test)[ , included]


true_beta = matrix(1 + rt(n_spp^2, df = 7), nrow = n_spp, ncol = n_spp)
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
    bayesglm(y_train[,i] ~ .,
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
  1E3,
  burn = 10,
  thin = 1
)

stack = bc$sims


plot(
  rowMeans(stack),
  rowSums(predicted_p_test),
  main = "BC stacking introduces no biases",
  xlab = "Expected richness (post-stacking)",
  ylab = "Expected richness (pre-stacking)"
)
abline(0,1)


# Confidence intervals ----------------------------------------------------
quantiles = sapply(1:nrow(stack), function(i){mean(stack[i, ] > sum(y_test[i, ]))})
naive_quantiles = sapply(1:nrow(stack), function(i){ppoibin(sum(y_test[i, ]), predicted_p_test[i, ])})

bad_naive = naive_quantiles < .025 | naive_quantiles > .975
bad = quantiles < .025 | quantiles > .975

pdf("poster-CI.pdf", width = 7, height = 4.1)
par(mfrow = c(1, 2), las = 1)
plot(rowSums(predicted_p_test), rowSums(y_test), ylab = "observed richness", asp = 1, xlab = "predicted richness", xlim = c(0, max(rowSums(y_test))),
     ylim = 1.02 * c(0, max(rowSums(y_test))), main = "with independent model errors", bty = "l", xaxs = "i", yaxs = "i",
     col = ifelse(bad_naive, "red", "black"), pch = ifelse(bad_naive, 4, 1), cex = 1/2, lwd = 1/2)
abline(0,1)
mean(bad_naive)

plot(rowMeans(stack), rowSums(y_test), ylab = "observed richness", asp = 1, xlab = "predicted richness", xlim = c(0, max(rowSums(y_test))),
     ylim = 1.02 * c(0, max(rowSums(y_test))), main = "with residual correlations", bty = "l", xaxs = "i", yaxs = "i",
     col = ifelse(bad, "red", "black"), pch = ifelse(bad, 4, 1), cex = 1/2, lwd = 1/2)
abline(0,1)
mean(bad)
dev.off()

rich_p_naive = t(apply(predicted_p_test, 1, function(x){dpoibin(0:n_spp, x)}))
rich_p = t(apply(stack,
                 1,
                 function(x){sapply(0:n_spp, function(n){mean(x==n)})}
))
smoothed_rich_p = (rich_p + 1/ncol(stack)) / (1 + ncol(rich_p)/ncol(stack))

i = sample.int(n_spp, 1)
plot(-.075 + seq(0, n_spp), rich_p_naive[i, ], type = "h", bty = "l", yaxs = "i", lwd = 2)
lines(.055 + seq(0, n_spp), smoothed_rich_p[i, ], type = "h", col = "#FF000099", lwd = 2)
points(sum(y_test[i, ]), 0, pch = 16, cex = 2)

# coverage ----------------------------------------------------------------

my_quantiles = sapply(
  1:n_loc,
  function(i) {
    sum(rich_p[i, sum(y_test[i, ]):n_spp])
  }
)
my_smoothed_quantiles = sapply(
  1:n_loc,
  function(i){
    sum(smoothed_rich_p[i, sum(y_test[i, ]):n_spp])
  }
)
naive_quantiles = sapply(
  1:n_loc,
  function(i) {
    sum(rich_p_naive[i, sum(y_test[i, ]):n_spp])
  }
)

plot(density(naive_quantiles, bw = 0.02), yaxs = "i", xaxs = "i", ylim = c(0, 9))
lines(density(my_smoothed_quantiles, bw = 0.02), col = 2)

# Coverage of 95% interval
mean(naive_quantiles > .025 & naive_quantiles < .975)
mean(my_smoothed_quantiles > .025 & my_quantiles < .975)

