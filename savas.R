require(DECENT)
require(edgeR)
source('func_de_methods.R')

### Read data ###
# single-cell data were clustered as described in the original paper. 
# We provide the subsetted data here.
data.obs <- readRDS('data.savas.rds')
cell.type <- readRDS('ct.savas.rds')

# bulk reference list is available as the Supp table 4 in the original paper. First export as .csv.
bulk.de <- read.csv("41591_2018_78_MOESM5_ESM.csv", skip = 1, row.names = 1, header = T)
ref <- substring(rownames(bulk.de), 2)

### single-cell DE analysis ###
# ZINB-WaVE, SCDE and DECENT are recommended to be run on server
er.table <- runEdgeR(data.obs, cell.type)
zer.table <- runZinbwaveEdgeR(data.obs, cell.type, J = 1e13)
mnc.table <- runMonocleDE(data.obs, cell.type)
mst.table <- runMASTDE(data.obs, cell.type)
scde.table <- runSCDE(data.obs, cell.type)
decent.table <- decent(data.obs, ~cell.type)

# calculate and save curves
rank <- list()
rank$er <- rownames(er.table)[order(er.table$PValue)]
rank$zer <- rownames(zer.table)[order(zer.table$PValue)]
rank$mnc <- rownames(mnc.table)[order(mnc.table$pval)]
rank$mst <- mst.table$primerid[order(mst.table$`Pr(>Chisq)`)]
rank$scde <- rownames(scde.table)[order(abs(scde.table$Z), decreasing = T)]
rank$decent <- decent.table$gene[order(decent.table$pvalue)]

sv.curves <- list()
sv.curves[['fpr']] <- sv.curves[['tpr']] <- sv.curves[['pr']] <- list()

for (met in c('er', 'zer', 'mnc', 'mst', 'scde', 'tasc', 'decent')){
  if (met != 'tasc') {
    sv.curves$fpr[[met]] <-  calcFPR(rank[[met]], ref)
    sv.curves$tpr[[met]] <-  calcTPR(rank[[met]], ref)
    sv.curves$pr[[met]] <-  calcPR(rank[[met]], ref)
  } else {
    sv.curves$fpr[[met]] <- sv.curves$tpr[[met]] <- rep(0, nrow(data.obs))
    sv.curves$pr[[met]] <- rep(1, nrow(data.obs))
  }
}

save(sv.curves, file='Savas_benckmark.rda')
