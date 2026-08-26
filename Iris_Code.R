library(ggplot2)
library(tidyverse)
library(pcaMethods)
library(corrplot)
library(GGally)
library(naniar)
library(patchwork)
library(ISLR)
library(class)
library(ALL)
library(Biobase)
library(MASS)

# Iris dataset ----------------------------------------------------
rm(list=ls())
data(iris)
X      <- iris[, 1:4]
X_std  <- scale(X)
specie <- iris$Species
table(iris$Species) #perfettamente diviso in 3
ggpairs(X,aes(colour = specie), upper = list(continuous = wrap("cor", size = 9)))
round(cor(X),2)
ggplot(iris)+
  geom_boxplot(aes(x=Sepal.Length,y=Species,fill = "red"))+
  geom_boxplot(aes(x=Sepal.Width,y=Species,fill = "blue"))+
  geom_boxplot(aes(x=Petal.Length,y=Species,fill = "green"))+
  geom_boxplot(aes(x=Petal.Width,y=Species,fill = "purple"))+
  labs(x="",y="")+
  scale_fill_manual(labels=c("Sepal.Length","Sepal.Width","Petal.Length","Petal.Width"),values=c("red","blue","green","purple"))+
  guides(fill=guide_legend("Variabile"))+
  theme(text = element_text(size = 20),
        axis.text.x = element_text(angle = 90, hjust = 1)) 
summary(iris[,1:4])
round(diag(var(iris[,1:4])),2)
#evidenzia come setosa sia ben diversa dalle altre due specie 
#inoltre le scale e le varianze sono decisamente diverse, quindi si standardizza

# PCA ---------------------------------------------------------------------
pca_class <- prcomp(X_std, center = FALSE, scale. = FALSE)
summary(pca_class)
# Varianza spiegata
var_exp <- pca_class$sdev^2 / sum(pca_class$sdev^2)
cum_var <- cumsum(var_exp)
data.frame(PC=1:4,var_spiegata_cumulata=cum_var)
# Screeplot
plot(pca_class$sdev^2,type="b",ylab="var_spiegata",xlab="PC")
abline(h=1,col=2)
# Loadings
loadings_df <- as.data.frame(pca_class$rotation)
round(loadings_df, 3)
# Biplot
scores_df <- as.data.frame(pca_class$x)
scores_df$Species <- specie
biplot_pca <- ggplot(scores_df, aes(PC1, PC2, color = Species)) +
  geom_point(size = 2.5, alpha = 1) +
  geom_segment(data = loadings_df,
               aes(x = 0, y = 0,
                   xend = PC1 * 2.5,
                   yend = PC2 * 2.5),
               arrow = arrow(length = unit(0.2,"cm")),
               color = "black", inherit.aes = FALSE) +
  geom_text(data = loadings_df,
            aes(x = PC1 * 2.8, y = PC2 * 2.8,
                label = rownames(loadings_df)),
            size = 3.5, color = "black", inherit.aes = FALSE) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "gray70") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "gray70") +
  labs(x = "PC1", y = "PC2") +
  theme(
    panel.background = element_rect(fill = "#F5F7FA")
  )
biplot_pca

# PPCA --------------------------------------------------------------------

ppca_res <- pca(X_std, method = "ppca", nPcs = 4, seed = 123)
summary(ppca_res)
# Screeplot
plot(ppca_res@sDev^2,type="b")
abline(h=1,col=2)
# Sigma^2
sigma2_ppca <- mean(ppca_res@sDev[3:4]^2)
cat("Rumore:", round(sigma2_ppca, 4), "\n")
ppca_res <- pca(X_std, method = "ppca", nPcs = 2, seed = 123)
summary(ppca_res)
# Biplot
scores_ppca <- as.data.frame(scores(ppca_res))
scores_ppca$Species <- specie
biplot_ppca <- ggplot(scores_ppca, aes(PC1, PC2, color = Species)) +
  geom_point(size = 2.5, alpha = 1) +
  geom_segment(data = ppca_res@loadings,
               aes(x = 0, y = 0,
                   xend = PC1 * 2.5,
                   yend = PC2 * 2.5),
               arrow = arrow(length = unit(0.2,"cm")),
               color = "black", inherit.aes = FALSE) +
  geom_text(data = ppca_res@loadings,
            aes(x = PC1 * 2.8, y = PC2 * 2.8,
                label = rownames(loadings_df)),
            size = 3.5, color = "black", inherit.aes = FALSE) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "gray70") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "gray70") +
  labs(x = "PC1", y = "PC2") +
  theme(
    panel.background = element_rect(fill = "#F5F7FA")
  )
biplot_ppca

# BPCA --------------------------------------------------------------------

bpca_res <- pca(X_std, method = "bpca", nPcs = 3, maxSteps = 200,seed=123)
summary(bpca_res)
# Alpha
W_hat <- loadings(bpca_res)
d     <- ncol(X_std)
alpha <- apply(W_hat, 2, function(w) d / sum(w^2))
cat("Alpha ARD :\n"); print(round(alpha, 3))

bpca_res <- pca(X_std, method = "bpca", nPcs = 2, maxSteps = 200,seed=123)
summary(bpca_res)
# Alpha
W_hat <- loadings(bpca_res)
d     <- ncol(X_std)
alpha <- apply(W_hat, 2, function(w) d / sum(w^2))
cat("Alpha ARD :\n"); print(round(alpha, 3))
#si potrebbe fare bootstrap per stimare la varianza di W

# Score plot BPCA
scores_bpca <- as.data.frame(scores(bpca_res))
scores_bpca$Species <- specie

biplot_bpca <- ggplot(scores_bpca, aes(PC1, PC2, color = Species)) +
  geom_point(size = 2.5, alpha = 1) +
  geom_segment(data = bpca_res@loadings,
               aes(x = 0, y = 0,
                   xend = PC1 * 2.5,
                   yend = PC2 * 2.5),
               arrow = arrow(length = unit(0.2,"cm")),
               color = "black", inherit.aes = FALSE) +
  geom_text(data = bpca_res@loadings,
            aes(x = PC1 * 2.8, y = PC2 * 2.8,
                label = rownames(loadings_df)),
            size = 3.5, color = "black", inherit.aes = FALSE) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "gray70") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "gray70") +
  labs(x = "PC1", y = "PC2") +
  theme(
    panel.background = element_rect(fill = "#F5F7FA")
  )
biplot_bpca


# Confronto R^2 -----------------------------------------------------------
r2_classic <- sum(var_exp[1:2])
cat("R^2 PCA classica (k=2):", round(r2_classic, 4), "\n")

r2_ppca <- ppca_res@R2cum[2]
cat("R^2 PPCA        (k=2):", round(r2_ppca, 4), "\n")

r2_bpca <- bpca_res@R2cum[2]
cat("R^2 BPCA        (k=2):", round(r2_bpca, 4), "\n")

# RMSE di ricostruzione ---------------------------------------------------

k=2
recon_class <- pca_class$x[, 1:k] %*% t(pca_class$rotation[, 1:k])
recon_ppca  <- scores(ppca_res)[,1:k] %*% t(loadings(ppca_res)[,1:k])
recon_bpca  <- scores(bpca_res)[,1:k] %*% t(loadings(bpca_res)[,1:k])

rmse <- function(orig, recon) sqrt(mean((orig - recon)^2))

cat("PCA classica :", round(rmse(X_std, recon_class), 5), "\n")
cat("PPCA         :", round(rmse(X_std, recon_ppca),  5), "\n")
cat("BPCA         :", round(rmse(X_std, recon_bpca),  5), "\n")
