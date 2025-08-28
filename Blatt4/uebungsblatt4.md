# Übungsblatt 4 – Korrelation & Dendrogramm (R)

**Ziel:** Sample-zu-Sample-Ähnlichkeiten verstehen und als Cluster-Bäume darstellen.
- Korrelationsmatrix der Samples
- Distanzen aus Korrelation ableiten
- Hierarchisches Clustering & Dendrogramm
- Heatmap der Korrelationen mit Annotation
- Cluster schneiden und mit Traits vergleichen

**Voraussetzung:** Aufbereitete Matrix (idealerweise `blatt2_mat_z.rds`).

---

## Aufgaben

1) **Daten laden**
   - Nutze `blatt2_mat_z.rds` (Z-Score je Gen). Falls nicht vorhanden, nimm `hafer_expression.csv` und zentriere/skaliere zeilenweise.

2) **Korrelation & Distanz**
   - Berechne **Pearson-Korrelation** zwischen Samples: `cor(t(mat_z))` → Matrix `S×S`.
   - Wandle zu einer **Distanz**: `dist_mat <- as.dist(1 - cor_mat)`.

3) **Hierarchisches Clustering**
   - `hc <- hclust(dist_mat, method="average")` (oder `complete`/`ward.D2`).
   - Dendrogramm plotten; optional farbige Balken für Traits.

4) **Heatmap der Korrelationen**
   - `pheatmap(cor_mat, annotation_col=traits, annotation_row=traits)`.

5) **Cluster schneiden & vergleichen**
   - `cutree(hc, k=2 or h=...)` → Gruppenzuordnung.
   - Vergleiche Cluster mit `Treatment/Batch/Variety` (Tabellen, Randtests).

6) **Export**
   - Lege `Blatt4_plots/` an und speichere alle Grafiken in eine PDF `blatt4_plots.pdf`.

---

## Startercode

```r
library(pheatmap)

# 1) Laden
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

# 4) Export-PDF
out_dir <- "Blatt4_plots"; dir.create(out_dir, showWarnings = FALSE)
pdf(file.path(out_dir, "blatt4_plots.pdf"), width = 7, height = 6)

plot(hc, main = "Hierarchisches Clustering der Samples")

# 5) Heatmap der Korrelationen mit Annotation
ann <- data.frame(Variety = traits$Variety,
                  Treatment = traits$Treatment,
                  Batch = traits$Batch)
row.names(ann) <- traits$Sample
pheatmap(cor_mat, annotation_col = ann, annotation_row = ann,
         main = "Sample–Sample Pearson-Korrelation")

# 6) Cluster schneiden
grp <- cutree(hc, k = 2)
plot(hc, main = "Dendrogramm mit Cluster (k=2)"); rect.hclust(hc, k = 2, border = "red")
print(grp)

dev.off()

# 7) Vergleich mit Traits
tab_trt <- table(Cluster = grp, Treatment = traits$Treatment[match(names(grp), traits$Sample)])
tab_batch <- table(Cluster = grp, Batch = traits$Batch[match(names(grp), traits$Sample)])
tab_var <- table(Cluster = grp, Variety = traits$Variety[match(names(grp), traits$Sample)])
tab_trt; tab_batch; tab_var
```

**Hinweise**
- Für WGCNA ist die Sample-Korrelation und das Dendrogramm ein früher QC-Schritt (Ausreißer erkennen).
- `1 - cor` ist eine *Ähnlichkeit → Distanz*-Heuristik; für Ward-Clustering besser `dist(t(mat))` auf z-transformierten Daten testen.

Viel Erfolg!
