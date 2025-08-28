# Blatt 4 – Korrelation & Dendrogramm (Template)

library(pheatmap)

# 1) Daten laden
if (file.exists("blatt2_mat_z.rds")) {
  mat <- readRDS("blatt2_mat_z.rds")
} else {
  expr <- read.csv("hafer_expression.csv", stringsAsFactors = FALSE)
  rownames(expr) <- expr$Gene; expr$Gene <- NULL
  mat <- t(scale(t(as.matrix(expr))))
}
traits <- read.csv("hafer_traits.csv", stringsAsFactors = FALSE)

# 2) Korrelation & Distanz
cor_mat <- cor(t(mat), method = "pearson", use = "pairwise.complete.obs")
dist_mat <- as.dist(1 - cor_mat)

# 3) Clustering
hc <- hclust(dist_mat, method = "average")

# 4) Export
out_dir <- "Blatt4_plots"; dir.create(out_dir, showWarnings = FALSE)
pdf(file.path(out_dir, "blatt4_plots.pdf"), width = 7, height = 6)

plot(hc, main = "Hierarchisches Clustering der Samples")

ann <- data.frame(Variety = traits$Variety,
                  Treatment = traits$Treatment,
                  Batch = traits$Batch)
row.names(ann) <- traits$Sample
pheatmap(cor_mat, annotation_col = ann, annotation_row = ann,
         main = "Sample–Sample Pearson-Korrelation")

grp <- cutree(hc, k = 2)
plot(hc, main = "Dendrogramm mit Cluster (k=2)"); rect.hclust(hc, k = 2, border = "red")
dev.off()

# Tabellen zur Zuordnung
tab_trt <- table(Cluster = grp, Treatment = traits$Treatment[match(names(grp), traits$Sample)])
tab_batch <- table(Cluster = grp, Batch = traits$Batch[match(names(grp), traits$Sample)])
tab_var <- table(Cluster = grp, Variety = traits$Variety[match(names(grp), traits$Sample)])
print(tab_trt); print(tab_batch); print(tab_var)
