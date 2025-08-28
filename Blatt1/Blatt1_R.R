expr <- read.csv("hafer_expression.csv", stringsAsFactors = FALSE)
traits <- read.csv("hafer_traits.csv", stringsAsFactors = FALSE)
dim(expr)
head(expr)
colnames(expr)

rownames(expr) <- expr$Gene
expr$Gene <- NULL
mat <- as.matrix(expr)

colMeans(mat) #Mittelwert pro Sample
sd <- apply(mat, 2, sd)
gene_var <- apply(mat, 1, var)
head(sort(gene_var, decreasing=TRUE), 5) # hohe Varianz = pot. informatives Gen

####QC####

#Boxplot
boxplot(mat, main="Expression per sample", las=2)

#Dichtekurve
cols <- rainbow(ncol(mat))
plot(density(mat[,1]), col=cols[1], main="Dichten je Sample", xlab="log2(Expression)")
for (i in 2:ncol(mat)) lines(density(mat[,i]), col=cols[i]) 
legend("topright", colnames(mat),col=cols,lty=1,bty="o")

#Heatmap
library(pheatmap)
top20 <- names(sort(gene_var, decreasing=TRUE))[1:20]
pheatmap(mat[top20, ], scale="row")

#PCA
pca <- prcomp(t(mat), scale. =TRUE)
plot(pca$x[,1], pca$x[,2], xlab="PC1", ylab="PC2", main="PCA der Samples")
text(pca$x[,1], pca$x[,2], labels = colnames(mat), pos=3)








