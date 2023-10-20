#import library
library(DropletTestFiles)
library(DropletUtils)
library(ggplot2)
library(Seurat)

my.theme <- theme(axis.title = element_text(size = 12), axis.text.x = element_text(angle = 0, hjust=0.5),
                  axis.text = element_text(size=12, color='black'), plot.title = element_text(size=12, face="plain"),
                  legend.position = "none")

# Import 10x data
data_dir <- '/path/data'
list.files(data_dir)
mat <- Read10X(data.dir = data_dir)

# Histogram Fig 3c
hist(colSums(mat), main="Histogram of UMIs/cell")
hist(colSums(mat>0), main="Histogram of Genes/cell", xlab='Genes per barcode', ylab='Cell frequency')

#*********************************************#
# Discard cells with more than 10,000 UMIs
mat <- mat[, colSums(mat)<=10e3]

# # Discard cells with fewer than 1,000 UMIs
# mat <- mat[, colSums(mat)>=1e3]

# Discard cells with fewer than 500 genes
mat <- mat[, colSums(mat>0)>=500]

# Keep genes detected in at least 5 cells
mat <- mat[rowSums(mat>0)>=5,]

#*********************************************#
#------- use seurat --------
pbmc <- CreateSeuratObject(mat)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
#Fig 3d
VlnPlot(pbmc, features=c('nCount_RNA','nFeature_RNA','percent.mt'), adjust=2)

#*********************************************#
#Creat knee plot based on filtered_cell_barcode
sce <- read10xCounts(data_dir)

colData(sce)[1:2,]
rowData(sce)[1:2,]

bcrank <- barcodeRanks(counts(sce))
bcrank[1:2,]

uniq <- !duplicated(bcrank$rank)
#Fig 3b
plot(bcrank$rank[uniq], bcrank$total[uniq], log = "xy", xlab = "Cell barcode rank", ylab = "Total UMI counts", 
     main='Knee plot', cex = 0.8, cex.lab = 1.2, cex.main = 1.5)
abline(h = metadata(bcrank)$inflection, col = "darkgreen", lty = 2)
abline(h = metadata(bcrank)$knee, col = "dodgerblue", lty = 2)