# Blatt 2 – Datenaufbereitung (Template)

# 1) Einlesen
expr <- read.csv("hafer_expression.csv", stringsAsFactors = FALSE)
rownames(expr) <- expr$Gene; expr$Gene <- NULL
mat <- as.matrix(expr)
stopifnot(is.numeric(mat))

# 2) Varianz filtern (oberes 50%-Quantil)
gene_var <- apply(mat, 1, var, na.rm = TRUE)
thr <- quantile(gene_var, 0.5)
keep <- gene_var >= thr
mat_filt <- mat[keep, ]

# 3) Z-Score je Gen
mat_z <- t(scale(t(mat_filt)))

# 4) Transponieren (Samples als Zeilen)
mat_t <- t(mat_z)

# 5) Speichern
saveRDS(mat_filt, "blatt2_mat_filt.rds")
saveRDS(mat_z,    "blatt2_mat_z.rds")
saveRDS(mat_t,    "blatt2_mat_t.rds")

# (Optional) Schnellcheck
# boxplot(mat_filt, las=2, main="Nach Varianzfilter")
# pca <- prcomp(mat_t, scale.=FALSE)
# plot(pca$x[,1], pca$x[,2], pch=19, main="PCA nach Aufbereitung"); text(pca$x[,1], pca$x[,2], labels=rownames(mat_t), pos=3)
