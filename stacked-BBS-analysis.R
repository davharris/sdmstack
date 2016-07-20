devtools::load_all()

library("GRaF")
library("parallel")
library("progress")
library("poibin")
library("glmnet")
library("MASS")

load("birds.Rdata")
x = as.data.frame(scale(x[ , 1:8]))
storage.mode(route.presence.absence) = "integer"

if(TRUE){
  p = readRDS("p.RDS")
  bc = readRDS("birds_bc.RDS")
}else{

  pb <- progress_bar$new(total = ncol(route.presence.absence))

  p = lapply(
    1:ncol(route.presence.absence),
    function(i) {
      pb$tick()
      m = graf(y = route.presence.absence[in.train, i], x = x[in.train, ])
      predict(m, x)[,1]
    }
  )
  p = do.call(cbind, p)


  bc = bc_stack(
    p_train = p[in.train, ],
    y_train = route.presence.absence[in.train, ],
    p_test = p[in.test, ],
    1E3,
    burn = 10,
    thin = 1
  )
}


# name things -------------------------------------------------------------

predicted_p_test = p[in.test, ]
y_test = route.presence.absence[in.test, ]


# compare means -----------------------------------------------------------

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


# correlations ------------------------------------------------------------

mean_cor = BayesComm:::upper2cor(colMeans(bc$bc_model$R))
eig = eigen(mean_cor)
plot(eig$values, type = "h")


# aab ---------------------------------------------------------------------

aab = read.csv("aab_traits.csv", stringsAsFactors = FALSE)

hab = aab[match(colnames(route.presence.absence), aab$common_name), ]$habitat

d = cbind(hab, as.data.frame(eig$vectors[,1:5]))
d$hab = as.character(d$hab)
d$hab[d$hab %in% c("Marsh", "Lake/Pond", "River/Stream")] = "wet"
d$hab[d$hab %in% c("Open Woodland", "")] = ""
d$hab[d$hab %in% c("Mountains", "Shore-line", "Deserts", "Town")] = ""
d = d[d$hab != "" & !is.na(d$hab), ]

da = lda(hab ~ ., data = d)
levels(factor(d$hab))
colors = c("#009E73", "#E69F00", "Darkgray", "#56B4E9")
shapes = c(17, 15, 1, 16)


# Wetland* means "Marsh" or "Lake/Pond" or "River/Stream"
pdf("LDA.pdf", width = 5, height = 5)
par(mar=c(5.1,5.1,1.1,2.1)-1)
plot(predict(da)$x[,1:2], col = colors[factor(d$hab)], pch = shapes[factor(d$hab)], axes = FALSE, xlab = "Linear Discriminant 1", ylab = "Linear Discriminant 2", lwd = 2)
abline(v = 0, h = 0, col = "#00000030")
dev.off()
pdf("LDA-legend.pdf", width = 5, height = 5)
par(mar=c(5.1,5.1,1.1,2.1)-1)
plot(predict(da)$x[,1:2], type = "n", axes = FALSE)
legend("topleft", legend = c("Forest", "Grassland", "Scrub", "Wetland*"), pch = shapes, col = colors)
dev.off()
