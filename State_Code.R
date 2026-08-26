library(ggplot2)
library(tidyverse)
library(pcaMethods)
library(naniar)
library(patchwork)
library(GGally)

# Analisi dataset state.x77 con NA artificiali ----------------------------
data("state")
dati_completi <- as.data.frame(state.x77)
boxplot(dati_completi)
ggpairs(dati_completi)

mu    <- colMeans(dati_completi)
sigma <- apply(dati_completi, 2, sd)
dati_completi_scalati <- scale(dati_completi, center = mu, scale = sigma)
dati_na <- dati_completi
set.seed(123)
n_celle <- prod(dim(dati_na))
n_na <- floor(.15 * n_celle)
idx <- sample(n_celle, n_na)
m <- as.matrix(dati_na)
m[idx] <- NA
dati_na <- as.data.frame(m)
cat("% NA globale:", round(mean(is.na(dati_na)) * 100, 2), "%\n")
print(round(colMeans(is.na(dati_na)) * 100, 2))
vis_miss(dati_na, sort_miss = F) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1, size = 12),
    plot.title  = element_text(face = "bold", size = 14)
  ) +
  scale_fill_manual(
    values = c("TRUE" = "red", "FALSE" = "grey80"),
    labels = c("TRUE" = "NA", "FALSE" = "non NA")
  )

X <- dati_na  
X_std <- scale(X, center = mu, scale = sigma)
X_std_NN <- scale(na.omit(X))
dim(X_std_NN)
boxplot(X)
ggpairs(X)
glimpse(X)
round(apply(X,2,function(x) var(na.omit(x))),1)
# PCA ---------------------------------------------------------------------

pca_class <- prcomp(X_std_NN,scale=F,center=F)
summary(pca_class)
# Screeplot
plot(pca_class,type="l")
abline(h=1,col=2)
# Varianza spiegata
var_exp <- pca_class$sdev^2 / sum(pca_class$sdev^2)
cum_var <- cumsum(var_exp)
data.frame(PC=1:4,var_spiegata_cumulata=cum_var)
# Loadings
loadings_df <- as.data.frame(pca_class$rotation)
round(loadings_df, 3)
# Biplot
scores_df <- as.data.frame(pca_class$x)
biplot_pca <- ggplot(scores_df, aes(PC1, PC2)) +
  geom_point(size = 2.5, alpha = 1,col="forestgreen") +
  # frecce loadings
  geom_segment(data = loadings_df,
               aes(x = 0, y = 0,
                   xend = PC1 * 2.5,
                   yend = PC2 * 2.5),
               arrow = arrow(length = unit(0.2,"cm")),
               color = "black", inherit.aes = FALSE) +
  geom_text(data = loadings_df,
            aes(x = PC1 * 2.8, y = PC2 * 2.8,
                label = rownames(loadings_df)),
            size = 5, color = "black", inherit.aes = FALSE) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "gray70") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "gray70") +
  labs(x = "PC1", y = "PC2") +
  theme(
    panel.background = element_rect(fill = "#F5F7FA")
  )
biplot_pca

# PPCA --------------------------------------------------------------------

ppca_res <- pca(X_std, method = "ppca", nPcs = 8, seed = 123)
summary(ppca_res)
# sigma^2 stimato
sigma2_ppca <- mean(ppca_res@sDev[4:8]^2)
cat("Rumore:", sigma2_ppca, "\n")
ppca_res <- pca(X_std, method = "ppca", nPcs = 3, seed = 123)
summary(ppca_res)

# Biplot
scores_ppca <- as.data.frame(scores(ppca_res))

biplot_ppca1 <- ggplot(scores_ppca, aes(PC1, PC2)) +
  geom_point(size = 2.5, alpha = 1,col="forestgreen") +
  geom_segment(data = ppca_res@loadings,
               aes(x = 0, y = 0,
                   xend = PC1 * 2.5,
                   yend = PC2 * 2.5),
               arrow = arrow(length = unit(0.2,"cm")),
               color = "black", inherit.aes = FALSE) +
  geom_text(data = ppca_res@loadings,
            aes(x = PC1 * 2.8, y = PC2 * 2.8,
                label = rownames(loadings_df)),
            size = 5, color = "black", inherit.aes = FALSE) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "gray70") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "gray70") +
  labs(x = "PC1", y = "PC2") +
  theme(
    panel.background = element_rect(fill = "#F5F7FA")
  )
biplot_ppca2 <- ggplot(scores_ppca, aes(PC1, PC3)) +
  geom_point(size = 2.5, alpha = 1,col="forestgreen") +
  geom_segment(data = ppca_res@loadings,
               aes(x = 0, y = 0,
                   xend = PC1 * 2.5,
                   yend = PC3 * 2.5),
               arrow = arrow(length = unit(0.2,"cm")),
               color = "black", inherit.aes = FALSE) +
  geom_text(data = ppca_res@loadings,
            aes(x = PC1 * 2.8, y = PC3 * 2.8,
                label = rownames(loadings_df)),
            size = 5, color = "black", inherit.aes = FALSE) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "gray70") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "gray70") +
  labs(x = "PC1", y = "PC3") +
  theme(
    panel.background = element_rect(fill = "#F5F7FA")
  )
biplot_ppca1+biplot_ppca2

# BPCA --------------------------------------------------------------------

bpca_res <- pca(X_std, method = "bpca", nPcs = 7, maxSteps = 200,seed=123)
summary(bpca_res) 
W_hat <- loadings(bpca_res)
d     <- ncol(X_std)
alpha <- apply(W_hat, 2, function(w) d / sum(w^2))
cat("Alpha:\n"); print(round(alpha, 3))
#anche con alpha noto come il quarto sia esploso
bpca_res <- pca(X_std, method = "bpca", nPcs = 3, maxSteps = 200,seed=123)
summary(bpca_res)
bpca_res <- pca(X_std, method = "bpca", nPcs = 2, maxSteps = 200,seed=123)
summary(bpca_res)
#ARD da k=2

# Biplot
scores_bpca <- as.data.frame(scores(bpca_res))

biplot_bpca <- ggplot(scores_bpca, aes(PC1, PC2)) +
  geom_point(size = 2.5, alpha = 1,col="forestgreen") +
  geom_segment(data = bpca_res@loadings,
               aes(x = 0, y = 0,
                   xend = PC1 * 2.5,
                   yend = PC2 * 2.5),
               arrow = arrow(length = unit(0.2,"cm")),
               color = "black", inherit.aes = FALSE) +
  geom_text(data = bpca_res@loadings,
            aes(x = PC1 * 2.8, y = PC2 * 2.8,
                label = rownames(loadings_df)),
            size = 5, color = "black", inherit.aes = FALSE) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "gray70") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "gray70") +
  labs(x = "PC1", y = "PC2") +
  theme(
    panel.background = element_rect(fill = "#F5F7FA")
  )
biplot_bpca

# Confronto R^2 -----------------------------------------------------------

r2_ppca <- ppca_res@R2cum[3]
cat("R^2 PPCA        (k=3):", round(r2_ppca, 4), "\n")

r2_bpca <- bpca_res@R2cum[2]
cat("R^2 BPCA        (k=2):", round(r2_bpca, 4), "\n")

# RMSE di ricostruzione su tutti i dati ---------------------------------------------------

recon_ppca <- completeObs(ppca_res)
rmse_ppca  <- sqrt(mean((dati_completi_scalati - recon_ppca)^2))

recon_bpca <- completeObs(bpca_res)
rmse_bpca  <- sqrt(mean((dati_completi_scalati - recon_bpca)^2))

cat("PPCA         :", round(rmse_ppca,  5), "\n")
cat("BPCA         :", round(rmse_bpca,  5), "\n")

# RMSE di ricostruzione solo NA -------------------------------------------

rmse_NA_ppca  <- sqrt(mean((dati_completi_scalati[idx] - recon_ppca[idx])^2))
rmse_NA_bpca  <- sqrt(mean((dati_completi_scalati[idx] - recon_bpca[idx])^2))

cat("PPCA:", round(rmse_NA_ppca, 5), "\n")
cat("BPCA:", round(rmse_NA_bpca, 5), "\n")

#PPCA fa leggermente meglio sia sulla ricostruzione completa che su quella degli
#NA... c'è da ricordare che ppca usa 3 componenti mentre bpca 2.

# Confronto con lo stesso k=3 ---------------------------------------------

ppca_k2 <- pca(X_std, method = "ppca", nPcs = 2, seed = 123)
bpca_k2 <- pca(X_std, method = "bpca", nPcs = 2, seed = 123)

recon_ppca <- completeObs(ppca_k2)
rmse_ppca  <- sqrt(mean((dati_completi_scalati - recon_ppca)^2))

recon_bpca <- completeObs(bpca_k2)
rmse_bpca  <- sqrt(mean((dati_completi_scalati - recon_bpca)^2))

r2_ppca <- ppca_k2@R2cum[2]
cat("R^2 PPCA        (k=2):", round(r2_ppca, 4), "\n")

r2_bpca <- bpca_k2@R2cum[2]
cat("R^2 BPCA        (k=2):", round(r2_bpca, 4), "\n")


cat("PPCA         :", round(rmse_ppca,  5), "\n")
cat("BPCA         :", round(rmse_bpca,  5), "\n")

rmse_NA_ppca  <- sqrt(mean((dati_completi_scalati[idx] - recon_ppca[idx])^2))
rmse_NA_bpca  <- sqrt(mean((dati_completi_scalati[idx] - recon_bpca[idx])^2))

cat("PPCA:", round(rmse_NA_ppca, 5), "\n")
cat("BPCA:", round(rmse_NA_bpca, 5), "\n")

