require(DECENT)
require(limma)
require(edgeR)
require(Seurat)
source('func_de_methods.R')

### Read data ###
# through personal contact
cluster <- read.table('NBID_Rh41_clust.info.txt', as.is = T, sep='\t')
colnames(cluster) <- cluster[1, ]
cluster <- cluster[-1, ]
cluster$barcode <- substr(cluster$barcode, 1, 16)

# UMI count matrix can be downloaded from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE113660
sc <- Seurat::Read10X('mex/') # '/mex' is the directory containing the 3 files in "GSE113660_RAW.tar"
sc <- sc[, colnames(sc) %in% cluster$barcode]
sc <- sc[rowSums(sc) > 0, ]
cut <- log10(rowSums(sc)) > 2
sc <- sc[cut, ]
data.obs <- as.matrix(sc)
cell.type <- as.numeric(cluster$`Cluster ID`)
rm(sc)

### Derive bulk reference DEG list ###
# bulk DE
bulk <- read.table('GSE113660_bulk_htseq_count.txt', row.names = 1, sep='\t', as.is = T)
colnames(bulk) <- bulk[1, ]
bulk <- bulk[-1, -c(1,2)]
bulk <- as.matrix(bulk)
mode(bulk) <- 'numeric'

cpm <- edgeR::cpm(bulk)
keep.exprs <- rowSums(cpm>1)>=3
bulk <- bulk[keep.exprs, ]
group <- factor(c(2, 1, 1, 1, 2, 2))
batch <- factor(c(1, 1, 2, 3, 2, 3))
design <- model.matrix(~group+batch)
y <- DGEList(counts = bulk[, 2:7], group = group)
v <- voom(y, design, plot=TRUE)
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, coefficients =2)
efit <- eBayes(vfit)
plotSA(efit, main="Final model: Mean−variance trend")
lv.table <- topTable(efit, number = Inf)
ref500 <- rownames(lv.table)[order(lv.table$P.Value)][1:500]
ref500 <- intersect(rownames(data.obs), ref500)

### single-cell DE analysis ###
# ZINB-WaVE, SCDE and DECENT are recommended to be run on server
er.table <- runEdgeR(data.obs, cell.type)
zer.table <- runZinbwaveEdgeR(data.obs, cell.type, J = 1e3)
mnc.table <- runMonocleDE(data.obs, cell.type)
mst.table <- runMASTDE(data.obs, cell.type)
scde.table <- runSCDE(data.obs, cell.type)
decent.table <- decent(data.obs, ~cell.type, CE.range=c(0.04,0.2))

# calculate and save curves
rank <- list()
rank$er <- rownames(er.table)[order(er.table$PValue)]
rank$zer <- rownames(zer.table)[order(zer.table$PValue)]
rank$mnc <- rownames(mnc.table)[order(mnc.table$pval)]
rank$mst <- mst.table$primerid[order(mst.table$`Pr(>Chisq)`)]
rank$scde <- rownames(scde.table)[order(abs(scde.table$Z), decreasing = T)]
rank$decent <- decent.table$gene[order(decent.table$pvalue)]

c.curves <- list()
c.curves[['fpr']] <- c.curves[['tpr']] <- c.curves[['pr']] <- list()

for (met in c('er', 'zer', 'mnc', 'mst', 'scde', 'tasc', 'decent')){
  if (met != 'tasc') {
    c.curves$fpr[[met]] <-  calcFPR(rank[[met]], ref500)
    c.curves$tpr[[met]] <-  calcTPR(rank[[met]], ref500)
    c.curves$pr[[met]] <-  calcPR(rank[[met]], ref500)
  } else {
    c.curves$fpr[[met]] <- c.curves$tpr[[met]] <- rep(0, length(rank$er))
    c.curves$pr[[met]] <- rep(1, length(rank$er))
  }
}

save(c.curves, file='Chen_benckmark.rda')
