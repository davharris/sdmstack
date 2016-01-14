# rstan -------------------------------------------------------------------

library(rstan)
model = stan_model("extras/stack.stan")

z = y_train - predicted_p_train
z_pos = z[y_train == 1]
z_neg = z[y_train == 0]

L_rho = t(chol(cor(z)))

fit = sampling(
  model,
  data = list(
    D = ncol(y_train),
    N = nrow(y_train),
    y = y_train,
    N_pos = sum(y_train),
    N_neg = sum(1 - y_train),
    n_pos = row(y_train)[y_train == 1],
    n_neg = row(y_train)[y_train == 0],
    d_pos = col(y_train)[y_train == 1],
    d_neg = col(y_train)[y_train == 0],
    probit_p = probit(predicted_p_train),
    lkj_eta = 2
  ),
  init = list(
    list(
      L_rho = L_rho,
      z_pos = z_pos,
      z_neg = z_neg
    )
  ),
  pars = c("rho", "z"),
  chains = 1,
  iter = 100
)

e = extract(fit)
