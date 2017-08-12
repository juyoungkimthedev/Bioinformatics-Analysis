#batch effect detection for hthgu133a chip type 
library(readr)
library(ggbiplot)
library(ggplot2)
wd0 <- "C:/Users/samsung/Desktop/University of Cape Town/MSc 2016-2017/R codes/hthgu133a/batch effect/"
setwd(wd0)
batch_csv <- readRDS("batch_csv.rds")

#PCA PLOT BEFORE BATCH EFFECT REMOVAL
#data preparation
colnames(batch_csv)[1] <- "sample"
batch_csv <- as.data.frame(batch_csv)
#normalizeQuantiles should be used for further analysis from now on
log.batch <- normalizeQuantiles(log2(batch_csv[, -c(1,2,3)]))
batch_id_data <- batch_csv[, 3] # A/B column
#apply pca
batch.pca <- prcomp(log.batch, center = T, scale = T)
plot(batch.pca, type = "l")
#ggbiplot
color_batch <- c(A = "red", B = "black")
col_batch <- color_batch[batch_id_data]
g <- ggbiplot(batch.pca, choices = 1:2, obs.scale = 1, var.scale = 1,
              groups = batch_id_data, ellipse = F, circle = F, var.axes = F, alpha = 0.01)
g <- g + scale_color_discrete(name = "") + scale_y_continuous("PC2") + scale_x_continuous("PC1")
g <- g + theme(legend.direction = "horizontal", legend.position = "none")
g <- g + geom_point(shape = 3, col = col_batch)
g <- g + ggtitle("PCA plot before batch effect removal")
print(g)

#PCA PLOT AFTER BATCH EFFECT REMOVAL [COMBAT]
library(sva)
batchid <- batch_csv$batch_id
com <- ComBat(t(log.batch), batchid, mod = NULL, par.prior = T, prior.plots = F)
#apply pca
pca.com <- prcomp(t(com), centerg = T, scale = T)
#ggbiplot
g <- ggbiplot(pca.com, choices = 1:2, obs.scale = 1, var.scale = 1,
              groups = batch_id_data, ellipse = F, circle = F, var.axes = F, alpha = 0.01)
g <- g + scale_color_discrete(name = "") + scale_y_continuous("PC2") + scale_x_continuous("PC1")
g <- g + theme(legend.direction = "horizontal", legend.position = "none")
g <- g + geom_point(shape = 3, col = col_batch)
g <- g + ggtitle("PCA plot after batch effect removal")
print(g)

#dev.off() to turn off the graphics
#student's t-test to confirm significance
t.test(log.batch, t(com))

#################################
#PCA plots for PC1,PC2,PC3 before and after batch effect removal
#PC1vsPC2 before batch effect removal
g <- ggbiplot(batch.pca, choices = 1:2, obs.scale = 1, var.scale = 1,
              groups = batch_id_data, ellipse = F, circle = F, var.axes = F, alpha = 0.01)
g <- g + scale_color_discrete(name = "") + scale_y_continuous("PC2") + scale_x_continuous("PC1")
g <- g + theme(legend.direction = "horizontal", legend.position = "none")
g <- g + geom_point(shape = 3, col = col_batch)
g <- g + ggtitle("PCA plot before batch effect removal")
print(g)

#PC1vsPC3 before batch effect removal
g <- ggbiplot(batch.pca, choices = c(1,3), obs.scale = 1, var.scale = 1,
              groups = batch_id_data, ellipse = F, circle = F, var.axes = F, alpha = 0.01)
g <- g + scale_color_discrete(name = "") + scale_y_continuous("PC3") + scale_x_continuous("PC1")
g <- g + theme(legend.direction = "horizontal", legend.position = "none")
g <- g + geom_point(shape = 3, col = col_batch)
g <- g + ggtitle("PCA plot before batch effect removal")
print(g)

#PC2vsPC3 before batch effect removal
g <- ggbiplot(batch.pca, choices = 2:3, obs.scale = 1, var.scale = 1,
              groups = batch_id_data, ellipse = F, circle = F, var.axes = F, alpha = 0.01)
g <- g + scale_color_discrete(name = "") + scale_y_continuous("PC3") + scale_x_continuous("PC2")
g <- g + theme(legend.direction = "horizontal", legend.position = "none")
g <- g + geom_point(shape = 3, col = col_batch)
g <- g + ggtitle("PCA plot before batch effect removal")
print(g)

#################################
#PC1vsPC2 after batch effect removal
g <- ggbiplot(pca.com, choices = 1:2, obs.scale = 1, var.scale = 1,
              groups = batch_id_data, ellipse = F, circle = F, var.axes = F, alpha = 0.01)
g <- g + scale_color_discrete(name = "") + scale_y_continuous("PC2") + scale_x_continuous("PC1")
g <- g + theme(legend.direction = "horizontal", legend.position = "none")
g <- g + geom_point(shape = 3, col = col_batch)
g <- g + ggtitle("PCA plot after batch effect removal")
print(g)

#PC1vsPC3 after batch effect removal
g <- ggbiplot(pca.com, choices = c(1,3), obs.scale = 1, var.scale = 1,
              groups = batch_id_data, ellipse = F, circle = F, var.axes = F, alpha = 0.01)
g <- g + scale_color_discrete(name = "") + scale_y_continuous("PC3") + scale_x_continuous("PC1")
g <- g + theme(legend.direction = "horizontal", legend.position = "none")
g <- g + geom_point(shape = 3, col = col_batch)
g <- g + ggtitle("PCA plot after batch effect removal")
print(g)

#PC2vsPC3 after batch effect removal
g <- ggbiplot(pca.com, choices = 2:3, obs.scale = 1, var.scale = 1,
              groups = batch_id_data, ellipse = F, circle = F, var.axes = F, alpha = 0.01)
g <- g + scale_color_discrete(name = "") + scale_y_continuous("PC2") + scale_x_continuous("PC3")
g <- g + theme(legend.direction = "horizontal", legend.position = "none")
g <- g + geom_point(shape = 3, col = col_batch)
g <- g + ggtitle("PCA plot after batch effect removal")
print(g)
dev.off()