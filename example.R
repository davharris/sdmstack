library(mvtnorm)
library(poibin)

set.seed(1)

n_spp = 20
n_loc = 500

true_rho = matrix(.7, n_spp, n_spp)
diag(true_rho) = 1

true_intercepts = rnorm(n_spp, -2)

expected_p = matrix(
  plogis(rnorm(n_spp * n_loc, mean = true_intercepts)),
  nrow = n_loc, 
  ncol = n_spp,
  byrow = TRUE
)

true_z = qnorm(expected_p) + rmvnorm(n_loc, rep(0, n_spp), true_rho)
true_z2 = qnorm(expected_p) + rmvnorm(n_loc, rep(0, n_spp), true_rho)

y = ifelse(true_z > 0, 1, 0)
y2 = ifelse(true_z2 > 0, 1, 0)

stack = bc_stack(p_train = expected_p, y_train = y, p_test = expected_p)

str(stack)

marginal_p = marginalize_out(stack, 3)


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
colnames(rich_p) = 0:n_spp

rich_p_naive = t(apply(expected_p, 1, function(x){dpoibin(0:n_spp, x)}))

sum(log(rich_p[cbind(1:nrow(y2), rowSums(y2) + 1)]))
sum(log(rich_p_naive[cbind(1:nrow(y2), rowSums(y2) + 1)]))


i = i+1
plot((0:n_spp) + .1, rich_p_naive[i, ], type = "h", yaxs = "i", bty = "l")
lines((0:n_spp) - .1, rich_p[i, ], type = "h", col = 2)
points(sum(y2[i,]), 0, pch = 16)
