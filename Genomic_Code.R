library(ISLR)
library(pcaMethods)
library(class)
library(ALL)
library(Biobase)
library(MASS)

# Funzione di valutazione -------------------------------------------------

evaluate <- function(scores_train, scores_test, y_train, y_test, k_knn = 5) {
  pred <- knn(train = as.matrix(scores_train),
              test  = as.matrix(scores_test),
              cl    = y_train, k = k_knn)
  return(mean(pred == y_test))
}

# Analisi dataset ALL ---------------------------------------------------
data(ALL)

X <- data.frame(t(exprs(ALL))) # 128 x 12625
y <- as.factor(ALL$BT)        
y_bt <- as.factor(ifelse(grepl("^B", ALL$BT), "B", "T"))
# sottotipi B/T con varianti (B, B1, B2, B3, B4, T, T1, ...)

# Split train-test --------------------------------------------------------

set.seed(123)
n         <- nrow(X)
idx_train <- sample(1:n, floor(0.7 * n))
X_train   <- X[idx_train, ];  y_train <- y_bt[idx_train]
X_test    <- X[-idx_train, ]; y_test  <- y_bt[-idx_train]
cat("\nDistribuzione classi nel training set:\n")
print(table(y_train))
cat("\nDistribuzione classi nel test set:\n")
print(table(y_test))

# Prefiltraggio -----------------------------------------------------------

gene_var  <- apply(X_train, 2, var)
top_genes <- order(gene_var, decreasing = TRUE)[1:100]
X_train_f <- X_train[, top_genes]
X_test_f  <- X_test[, top_genes]

# PCA ---------------------------------------------------------------------
set.seed(123)
pca_fit <- prcomp(X_train_f, scale. = TRUE)
summary(pca_fit)
plot(pca_fit$sdev^2, type = "b",
     ylab = "Varianza spiegata")

for (k in c(2, 5, 10, 20)) {
  scores_train <- pca_fit$x[, 1:k]
  scores_test  <- predict(pca_fit, X_test_f)[, 1:k]
  acc          <- evaluate(scores_train, scores_test, y_train, y_test)
  cat("PCA  k =", k, "| Accuracy =", round(acc, 3),
      "| R2 cum =", round(sum(pca_fit$sdev[1:k]^2) /
                            sum(pca_fit$sdev^2), 3), "\n")
}

# PPCA --------------------------------------------------------------------
set.seed(123)
ppca_fit <- pca(X_train_f, method = "ppca", nPcs = 20, scale = "uv", seed = 123)

for (k in c(2, 5, 10, 20)) {
  scores_train  <- scores(ppca_fit)[, 1:k]
  X_test_scaled <- scale(X_test_f,
                         center = ppca_fit@center,
                         scale  = ppca_fit@scale)
  scores_test   <- X_test_scaled %*% loadings(ppca_fit)[, 1:k]
  acc           <- evaluate(scores_train, scores_test, y_train, y_test)
  cat("PPCA k =", k, "| Accuracy =", round(acc, 3),
      "| R2 cum =", round(sum(ppca_fit@R2[1:k]), 3), "\n")
}
lines(ppca_fit@sDev^2,col=2)

# BPCA --------------------------------------------------------------------
set.seed(123)
bpca_fit <- pca(X_train_f, method = "bpca", nPcs = 20, scale = "uv",
                maxSteps = 300, verbose = TRUE,seed=123)
summary(bpca_fit)

k_star <- sum(bpca_fit@R2 > 1e-10)
cat("\nBPCA k* dall'ARD =", k_star, "\n")

bpca_fit <- pca(X_train_f, method = "bpca", nPcs = k_star, scale = "uv",
                maxSteps = 300,seed=123)
summary(bpca_fit)

k_star <- sum(bpca_fit@R2 > 1e-10)
cat("\nBPCA k* dall'ARD =", k_star, "\n")

bpca_fit <- pca(X_train_f, method = "bpca", nPcs = k_star, scale = "uv",
                maxSteps = 300,seed=123)
summary(bpca_fit)

scores_train  <- scores(bpca_fit)[, 1:k_star, drop = FALSE]
X_test_scaled <- scale(X_test_f,
                       center = bpca_fit@center,
                       scale  = bpca_fit@scale)
scores_test   <- X_test_scaled %*% loadings(bpca_fit)[, 1:k_star, drop = FALSE]
acc           <- evaluate(scores_train, scores_test, y_train, y_test)
cat("BPCA k* =", k_star, "| Accuracy =", round(acc, 3),
    "| R2 cum =", round(sum(bpca_fit@R2[1:k_star]), 3), "\n")

# Analisi dataset NCI60 ---------------------------------------------------

data(NCI60)
X <- data.frame(NCI60$data )   
y <- as.factor(NCI60$labs)
levels(y)

# Split train-test --------------------------------------------------------

set.seed(123)
n         <- nrow(X)
idx_train <- sample(1:n, floor(0.7 * n))
X_train   <- X[idx_train, ];  y_train <- y[idx_train]
X_test    <- X[-idx_train, ]; y_test  <- y[-idx_train]
cat("\nDistribuzione classi nel training set:\n")
print(table(y_train))
cat("\nDistribuzione classi nel test set:\n")
print(table(y_test))

# Prefiltraggio -----------------------------------------------------------

gene_var  <- apply(X_train, 2, var)
top_genes <- order(gene_var, decreasing = TRUE)[1:100]
X_train_f <- X_train[, top_genes]
X_test_f  <- X_test[, top_genes]

# PCA ---------------------------------------------------------------------
set.seed(123)
pca_fit <- prcomp(X_train_f, scale. = TRUE)
plot(pca_fit$sdev^2, type = "b",
     ylab = "Varianza spiegata")

for (k in c(2, 5, 10, 20)) {
  scores_train <- pca_fit$x[, 1:k]
  scores_test  <- predict(pca_fit, X_test_f)[, 1:k]
  acc          <- evaluate(scores_train, scores_test, y_train, y_test)
  cat("PCA  k =", k, "| Accuracy =", round(acc, 3),
      "| R2 cum =", round(sum(pca_fit$sdev[1:k]^2) /
                            sum(pca_fit$sdev^2), 3), "\n")
}

# PPCA --------------------------------------------------------------------
set.seed(123)
ppca_fit <- pca(X_train_f, method = "ppca", nPcs = 20, scale = "uv", seed = 123)

for (k in c(2, 5, 10, 20)) {
  scores_train  <- scores(ppca_fit)[, 1:k]
  X_test_scaled <- scale(X_test_f,
                         center = ppca_fit@center,
                         scale  = ppca_fit@scale)
  scores_test   <- X_test_scaled %*% loadings(ppca_fit)[, 1:k]
  acc           <- evaluate(scores_train, scores_test, y_train, y_test)
  cat("PPCA k =", k, "| Accuracy =", round(acc, 3),
      "| R2 cum =", round(sum(ppca_fit@R2[1:k]), 3), "\n")
}
lines(ppca_fit@sDev^2,col=2)

# BPCA --------------------------------------------------------------------
set.seed(123)
bpca_fit <- pca(X_train_f, method = "bpca", nPcs = 20, scale = "uv",
                maxSteps = 300, verbose = TRUE,seed=123)
summary(bpca_fit)

k_star <- sum(bpca_fit@R2 > 1e-10)
cat("\nBPCA k* dall'ARD =", k_star, "\n")

bpca_fit <- pca(X_train_f, method = "bpca", nPcs = k_star, scale = "uv",
                maxSteps = 300,seed=123)
summary(bpca_fit)

k_star <- sum(bpca_fit@R2 > 1e-10)
cat("\nBPCA k* dall'ARD =", k_star, "\n")

bpca_fit <- pca(X_train_f, method = "bpca", nPcs = k_star, scale = "uv",
                maxSteps = 300,seed=123)
summary(bpca_fit)

scores_train  <- scores(bpca_fit)[, 1:k_star, drop = FALSE]
X_test_scaled <- scale(X_test_f,
                       center = bpca_fit@center,
                       scale  = bpca_fit@scale)
scores_test   <- X_test_scaled %*% loadings(bpca_fit)[, 1:k_star, drop = FALSE]
acc           <- evaluate(scores_train, scores_test, y_train, y_test)
cat("BPCA k* =", k_star, "| Accuracy =", round(acc, 3),
    "| R2 cum =", round(sum(bpca_fit@R2[1:k_star]), 3), "\n")

# Analisi dataset KHAN ---------------------------------------------------

data(Khan)

# Split train-test --------------------------------------------------------

X_train   <- Khan$xtrain;  y_train <- Khan$ytrain
X_test    <- Khan$xtest; y_test  <- Khan$ytest
cat("\nDistribuzione classi nel training set:\n")
print(table(y_train))
cat("\nDistribuzione classi nel test set:\n")
print(table(y_test))


# PCA ---------------------------------------------------------------------
set.seed(123)
pca_fit <- prcomp(X_train, scale. = TRUE)
summary(pca_fit)
plot(pca_fit$sdev^2, type = "b",
     ylab = "Varianza spiegata")

for (k in c(2, 5, 10, 20)) {
  scores_train <- pca_fit$x[, 1:k]
  scores_test  <- predict(pca_fit, X_test)[, 1:k]
  acc          <- evaluate(scores_train, scores_test, y_train, y_test)
  cat("PCA  k =", k, "| Accuracy =", round(acc, 3),
      "| R2 cum =", round(sum(pca_fit$sdev[1:k]^2) /
                            sum(pca_fit$sdev^2), 3), "\n")
}

# PPCA --------------------------------------------------------------------
set.seed(123)
ppca_fit <- pca(X_train, method = "ppca", nPcs = 20, scale = "uv", seed = 123)

for (k in c(2, 5, 10, 20)) {
  scores_train  <- scores(ppca_fit)[, 1:k]
  X_test_scaled <- scale(X_test,
                         center = ppca_fit@center,
                         scale  = ppca_fit@scale)
  scores_test   <- X_test_scaled %*% loadings(ppca_fit)[, 1:k]
  acc           <- evaluate(scores_train, scores_test, y_train, y_test)
  cat("PPCA k =", k, "| Accuracy =", round(acc, 3),
      "| R2 cum =", round(sum(ppca_fit@R2[1:k]), 3), "\n")
}
lines(ppca_fit@sDev^2,col=2)

# BPCA --------------------------------------------------------------------
set.seed(123)
bpca_fit <- pca(X_train, method = "bpca", nPcs = 20, scale = "uv",
                maxSteps = 300, verbose = TRUE,seed=123)
summary(bpca_fit)

k_star <- sum(bpca_fit@R2 > 1e-10)
cat("\nBPCA k* dall'ARD =", k_star, "\n")

bpca_fit <- pca(X_train, method = "bpca", nPcs = k_star, scale = "uv",
                maxSteps = 300,seed=123)
summary(bpca_fit)

k_star <- sum(bpca_fit@R2 > 1e-10)
cat("\nBPCA k* dall'ARD =", k_star, "\n")

bpca_fit <- pca(X_train, method = "bpca", nPcs = k_star, scale = "uv",
                maxSteps = 300,seed=123)
summary(bpca_fit)

scores_train  <- scores(bpca_fit)[, 1:k_star, drop = FALSE]
X_test_scaled <- scale(X_test,
                       center = bpca_fit@center,
                       scale  = bpca_fit@scale)
scores_test   <- X_test_scaled %*% loadings(bpca_fit)[, 1:k_star, drop = FALSE]
acc           <- evaluate(scores_train, scores_test, y_train, y_test)
cat("BPCA k* =", k_star, "| Accuracy =", round(acc, 3),
    "| R2 cum =", round(sum(bpca_fit@R2[1:k_star]), 3), "\n")
