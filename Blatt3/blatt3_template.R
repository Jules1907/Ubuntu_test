# Blatt 3 – Visualisierung (Template)

library(pheatmap)

# 1) Einlesen
expr <- read.csv("hafer_expression.csv", stringsAsFactors = FALSE)
rownames(expr) <- expr$Gene; expr$Gene <- NULL
mat <- as.matrix(expr)
traits <- read.csv("hafer_traits.csv", stringsAsFactors = FALSE)

# Optional: Z-Matrix aus Modul 2 bevorzugen
if (file.exists("blatt2_mat_z.rds")) mat <- readRDS("blatt2_mat_z.rds")

# 2) Export-Device
out_dir <- "Blatt3_plots"; dir.create(out_dir, showWarnings = FALSE)
pdf(file.path(out_dir, "blatt3_plots.pdf"), width = 7, height = 6)

# 3) Boxplot + Dichte
boxplot(mat, main="Expression per Sample", las=2)

cols <- rainbow(ncol(mat))
plot(density(mat[,1]), col=cols[1], main="Dichten je Sample", xlab="log2(Expression)")
for (j in 2:ncol(mat)) lines(density(mat[,j]), col=cols[j], lwd=2)
legend("topright", legend=colnames(mat), col=cols, lty=1, lwd=2, bty="o")

# 4) PCA
pca <- prcomp(t(mat), scale. = FALSE)
var_exp <- (pca$sdev^2) / sum(pca$sdev^2)
trt <- traits$Treatment[match(colnames(mat), traits$Sample)]
col_trt <- setNames(c("steelblue","tomato"), sort(unique(trt)))
plot(pca$x[,1], pca$x[,2], pch=19, col=col_trt[trt],
     xlab=paste0("PC1 (", round(100*var_exp[1],1), "%)"),
     ylab=paste0("PC2 (", round(100*var_exp[2],1), "%)"),
     main="PCA nach Treatment")
text(pca$x[,1], pca$x[,2], labels=colnames(mat), pos=3)
legend("topright", legend=names(col_trt), col=col_trt, pch=19, title="Treatment", bty="n")

# 5) Heatmap
ann <- data.frame(Variety = traits$Variety,
                  Treatment = traits$Treatment,
                  Batch = traits$Batch)
row.names(ann) <- traits$Sample
pheatmap(mat, scale="row", annotation_col = ann)

dev.off()
