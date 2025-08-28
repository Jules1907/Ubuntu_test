# Übungsblatt 2 – Datenaufbereitung (R)

**Ziel:** Du übst typische Aufbereitungsschritte vor Clustering/WGCNA:
- Varianz-basiertes Filtern
- (falls nötig) Log-Transformation
- Z-Scoring (Skalierung je Gen)
- Transponieren (Samples in Zeilen)
- Zwischenergebnisse speichern

**Dateien (aus Modul 1):**
- `hafer_expression.csv` (log2-Expression, 50 Gene × 6 Samples)
- `hafer_traits.csv`

---

## Aufgaben

### 1) Einlesen & Struktur prüfen
1. Lies `hafer_expression.csv` ein und setze `Gene` als Zeilennamen.
2. Erzeuge eine **numerische** Matrix `mat`.
3. Prüfe auf NAs/nicht-numerische Spalten.

### 2) Varianz-basiertes Filtern
1. Berechne die Zeilenvarianz für alle Gene.
2. Wähle die **obersten 50%** variabelsten Gene (oder alternativ Top-20).
3. Erzeuge `mat_filt` mit nur diesen Genen.

### 3) (Optional) Log-Transformation
Unsere Übungsdaten sind bereits **log2**. Für rohe Counts: `log2(counts + 1)` oder besser `vst()`/`rlog()` (DESeq2).

### 4) Skalierung (Z-Score)
Skaliere **pro Gen** auf Mittelwert 0, SD 1: `t(scale(t(mat_filt)))` → Ergebnis `mat_z`.

### 5) Transponieren
Für Methoden mit Beobachtungen in Zeilen (z. B. PCA): `mat_t <- t(mat_z)`.

### 6) Speichern
Speichere `mat_filt`, `mat_z`, `mat_t` als `.rds`.

---

## Startercode

```r
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

# Alternative: Top-N
# topN <- 30
# idx <- order(gene_var, decreasing = TRUE)[1:topN]
# mat_filt <- mat[idx, ]

# 3) Z-Score je Gen
mat_z <- t(scale(t(mat_filt)))  # center=TRUE, scale=TRUE per default

# 4) Transponieren
mat_t <- t(mat_z)

# 5) Speichern
saveRDS(mat_filt, file = "blatt2_mat_filt.rds")
saveRDS(mat_z,    file = "blatt2_mat_z.rds")
saveRDS(mat_t,    file = "blatt2_mat_t.rds")
```

## Kontrollfragen
- Warum filtern wir nach Varianz? Was passiert bei zu streng/zu locker?
- Wieso `t(scale(t(X)))`?
- Warum transponieren wir vor PCA?
