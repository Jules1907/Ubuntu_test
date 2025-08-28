# Übungsblatt 1 – Daten einlesen & erste Schritte (R)

**Dateien:** 
- `hafer_expression.csv` (50 Gene × 6 Samples, log2-Expressionswerte)
- `hafer_traits.csv` (Metadaten zu den 6 Samples; für spätere Übungen)

## Aufgaben
1) **Einlesen**
   - Lies `hafer_expression.csv` als data.frame ein.
   - Prüfe die Dimensionen (`dim()`), die Spaltennamen (`colnames()`), und die ersten Zeilen (`head()`).
2) **Aufräumen & Struktur**
   - Setze die `Gene`-Spalte als Zeilennamen (rownames).
   - Konvertiere die restlichen Spalten in eine numerische Matrix.
3) **Deskriptive Statistik**
   - Berechne **Spaltenmittelwerte** und **Spalten-Standardabweichungen** (pro Sample).
   - Finde die **Top-5 Gene** mit der höchsten *Varianz* über die 6 Samples.
4) **Visualisierung (QC)**
   - Erstelle einen **Boxplot** der 6 Samples.
   - Erstelle ein **Density-Plot** (Dichtekurven) aller Samples in einem Plot.
5) **Bonus**
   - Erstelle eine **Heatmap** der Top-20 variabelsten Gene.
   - Erstelle eine **PCA** der Samples auf Basis der Matrix.

## R-Startercode (optional)
```r
# Arbeitsverzeichnis anpassen, falls nötig:
# setwd("/Pfad/zu/deinem/Repo")

expr <- read.csv("hafer_expression.csv", stringsAsFactors = FALSE)

# 1) Überblick
dim(expr); colnames(expr); head(expr)

# 2) Gene als rownames, Matrix extrahieren
rownames(expr) <- expr$Gene
expr$Gene <- NULL
mat <- as.matrix(expr)

# 3) Deskriptive Statistik
colMeans(mat)
apply(mat, 2, sd)  # sd pro Spalte
gene_var <- apply(mat, 1, var)
head(sort(gene_var, decreasing = TRUE), 5)

# 4) Visualisierung
boxplot(mat, main = "Expression per Sample", las = 2)
plot(density(mat[,1]), main = "Dichten je Sample")
for (i in 2:ncol(mat)) lines(density(mat[,i]))

# 5) Bonus: Heatmap & PCA
# install.packages("pheatmap")  # falls nötig
# library(pheatmap)
# top20 <- names(sort(gene_var, decreasing = TRUE))[1:20]
# pheatmap(mat[top20, ], scale = "row")

pca <- prcomp(t(mat), scale. = TRUE)
plot(pca$x[,1], pca$x[,2], xlab = "PC1", ylab = "PC2", main = "PCA der Samples")
text(pca$x[,1], pca$x[,2], labels = colnames(mat), pos = 3)
```
