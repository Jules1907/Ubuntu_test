# Übungsblatt 3 – Visualisierung (R)

**Ziel:** Saubere, reproduzierbare Visualisierungen erstellen und speichern.
- Boxplots & Dichtekurven je Sample
- PCA mit Prozentanteilen und Trait-Färbung
- Heatmap mit Spalten-Annotation
- Export als mehrseitige PDF

**Voraussetzung:** `hafer_expression.csv`, `hafer_traits.csv` (aus Modul 1).
Optional: Ergebnisse aus Modul 2 (`blatt2_mat_filt.rds`, `blatt2_mat_z.rds`).

---

## Aufgaben

1) **Einlesen**
   - Lies `hafer_expression.csv` ein, setze `Gene` als Zeilennamen, erzeuge `mat`.
   - Falls vorhanden: nutze `blatt2_mat_z.rds` (Z-Score pro Gen) bevorzugt für Plots.

2) **Boxplot & Dichte**
   - Boxplot über alle Samples.
   - Dichtekurven je Sample, farbig nach Sample und mit Legende.

3) **PCA**
   - PCA auf `t(mat)` (oder `mat_t` aus Modul 2).
   - Achsen mit **% erklärter Varianz** beschriften.
   - Punkte nach **Treatment** färben, optional nach **Batch/Variety** formen.

4) **Heatmap**
   - `pheatmap` mit `scale="row"`.
   - Spaltenannotation aus `traits` (Variety, Treatment, Batch).

5) **Export**
   - Lege Ordner `Blatt3_plots/` an.
   - Erzeuge **eine PDF** `blatt3_plots.pdf` mit allen Grafiken (mehrseitig).

---

## Startercode (du kannst direkt in RStudio kopieren)

```r
library(pheatmap)

# --- 1) Einlesen ---
expr <- read.csv("hafer_expression.csv", stringsAsFactors = FALSE)
rownames(expr) <- expr$Gene; expr$Gene <- NULL
mat <- as.matrix(expr)
traits <- read.csv("hafer_traits.csv", stringsAsFactors = FALSE)
stopifnot(is.numeric(mat))

# Optional: Z-Score-Matrix aus Modul 2 verwenden, falls vorhanden
mat_z_path <- "blatt2_mat_z.rds"
if (file.exists(mat_z_path)) {
  mat <- readRDS(mat_z_path)
}

# --- 2) Boxplot & Dichte ---
out_dir <- "Blatt3_plots"; dir.create(out_dir, showWarnings = FALSE)
pdf(file.path(out_dir, "blatt3_plots.pdf"), width = 7, height = 6)

boxplot(mat, main="Expression per Sample", las=2)

cols <- rainbow(ncol(mat))
plot(density(mat[,1]), col=cols[1], main="Dichten je Sample", xlab="log2(Expression)")
for (j in 2:ncol(mat)) lines(density(mat[,j]), col=cols[j], lwd=2)
legend("topright", legend=colnames(mat), col=cols, lty=1, lwd=2, bty="o")

# --- 3) PCA ---
pca <- prcomp(t(mat), scale. = FALSE) # mat ist ggf. schon z-transformiert
var_exp <- (pca$sdev^2) / sum(pca$sdev^2)

# Farben nach Treatment
trt <- traits$Treatment[match(colnames(mat), traits$Sample)]
col_trt <- setNames(c("steelblue","tomato"), sort(unique(trt)))
plot(pca$x[,1], pca$x[,2], pch=19, col=col_trt[trt],
     xlab=paste0("PC1 (", round(100*var_exp[1],1), "%)"),
     ylab=paste0("PC2 (", round(100*var_exp[2],1), "%)"),
     main="PCA nach Treatment")
text(pca$x[,1], pca$x[,2], labels=colnames(mat), pos=3)
legend("topright", legend=names(col_trt), col=col_trt, pch=19, title="Treatment", bty="n")

# --- 4) Heatmap mit Annotation ---
ann <- data.frame(Variety = traits$Variety,
                  Treatment = traits$Treatment,
                  Batch = traits$Batch)
row.names(ann) <- traits$Sample

pheatmap(mat, scale="row", annotation_col = ann)

dev.off()
```

**Hinweise**
- Für publikationsreife Grafiken eignen sich **PDFs** (vektorbasiert). PNG mit `res=300` ist gut für Slides.
- Bei sehr vielen Genen: Heatmap auf Top-Variablen beschränken (Varianzfilter aus Modul 2).

Viel Erfolg!
