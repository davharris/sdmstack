n_spp = 100
i = i + 1
x = bc$sims[i, ]


f = function(rho){
  -sum(VGAM::dbetabinom(x, size = n_spp, prob = mean(x) / n_spp, rho = rho, log = TRUE))
}

rho = optimize(f, c(0, 1))$min


plot(table(x), col = "gray", xlim = c(0, n_spp))
p = VGAM::dbetabinom(0:n_spp, size = n_spp, prob = mean(x) / n_spp, rho = rho, log = FALSE)
log_p = VGAM::dbetabinom(0:n_spp, size = n_spp, prob = mean(x) / n_spp, rho = rho, log = TRUE)
lines(
  0:n_spp,
  length(x) * p,
  col = 2,
  lwd = 2
)

# Sample from the beta-binomial
fake_x = VGAM::rbetabinom(1E5, size = n_spp, prob = mean(x) / n_spp, rho = rho)


# Entropy of the beta-binomial
-sum(p * log_p)

# Log-likelihood of beta-binomial given the samples
-mean(VGAM::dbetabinom(x, size = n_spp, prob = mean(x) / n_spp, rho = rho, log = TRUE))

# Log-likelihood of beta-binomial given  samples from the beta-binomial
-mean(VGAM::dbetabinom(fake_x, size = n_spp, prob = mean(x) / n_spp, rho = rho, log = TRUE))


# Calculate highest-density intervals
q = .95
include = cumsum(sort(p, decreasing = TRUE)) < q
sorted_indices = (0:n_spp)[order(p, decreasing = TRUE)]
abline(v = range(sorted_indices[include]), col = 2, lwd = 2)

