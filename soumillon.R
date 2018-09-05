require(DECENT)
require(edgeR)
require(ggpubr)

source('func_de_methods.R')

### Read data ###
# downloaded from GEO entry: GSE53638
sc <- read.delim('GSE53638_D3_UMI.dat')
bulk <- read.delim('GSE53638_D3_Bulk_UMI.dat')


### Subset single-cell data ###
time.point <- substr(colnames(sc), 1, 4)

# two group benckmark
data.obs  <- as.matrix(sc[,time.point %in% c('D3T0', 'D3T7')])
cell.type <- as.factor(time.point[time.point %in% c('D3T0', 'D3T7')])
data.obs <- data.obs[rowSums(data.obs) > 0, ]
cut <- median(log(rowSums(data.obs))) - mad(log(rowSums(data.obs)))
data.obs <- data.obs[log(rowSums(data.obs)) > cut, ]

# three group benchmark
data.obs.3g  <- as.matrix(sc[,time.point != 'X'])
cell.type.3g <- as.factor(time.point[time.point != 'X'])
data.obs.3g <- data.obs.3g[rowSums(data.obs.3g) > 0, ]
cut <- median(log(rowSums(data.obs.3g))) - mad(log(rowSums(data.obs.3g)))
data.obs.3g <- data.obs.3g[log(rowSums(data.obs.3g)) > cut, ]


### Derive bulk reference gene list ###
# two group benckmark
bulk.g <- edgeR::cpm(bulk[, c(2,4)], log=T, prior.count = 5)
f <- rowSums(bulk[, c(2, 4)] > 0) == 2
logfc <- bulk.g[f, 1] - bulk.g[f, 2]
ref500 <- as.character(bulk$X)[f][order(abs(logfc), decreasing = T)][1:500]
ref500 <- intersect(rownames(data.obs), ref500)

# three group benchmark
bulk.3g <- edgeR::cpm(bulk[, c(2, 3, 4)], log=T, prior.count = 5)
f.3g <- rowSums(bulk[, c(2, 3, 4)] > 0) == 3
var.3g <- rowVars(bulk.3g[f, ])
ref500.3g <- as.character(bulk$X)[f.3g][order(var.3g, decreasing = T)][1:500]
ref500.3g <- intersect(rownames(data.obs), ref500.3g)


### Single-cell DE analysis ###
# ZINB-WaVE, SCDE and DECENT are recommended to be run on server

# two group
er.table <- runEdgeR(data.obs, cell.type)
zer.table <- runZinbwaveEdgeR(data.obs, cell.type, J = 1e7)
mnc.table <- runMonocleDE(data.obs, cell.type)
mst.table <- runMASTDE(data.obs, cell.type)
scde.table <- runSCDE(data.obs, cell.type) #
decent.table <- decent(data.obs, ~cell.type, normalize = 'TMM')
decent.table.ml <- decent(data.obs, ~cell.type, normalize = 'ML')

# three group
er.table.3g <- runEdgeR(data.obs.3g, cell.type.3g)
zer.table.3g <- runZinbwaveEdgeR(data.obs.3g, cell.type.3g, J = 1e7)
mnc.table.3g <- runMonocleDE(data.obs.3g, cell.type.3g)
mst.table.3g <- runMASTDE3(data.obs.3g, cell.type.3g)
decent.table.3g <- decent(data.obs.3g, ~cell.type.3g, normalize = 'TMM')

# calculate curves
rank <- list()
rank$er <- rownames(er.table)[order(er.table$PValue)]
rank$zer <- rownames(zer.table)[order(zer.table$PValue)]
rank$mnc <- rownames(mnc.table)[order(mnc.table$pval)]
rank$mst <- mst.table$primerid[order(mst.table$`Pr(>Chisq)`)]
rank$scde <- rownames(scde.table)[order(abs(scde.table$Z), decreasing = T)]
rank$decent <- decent.table$gene[order(decent.table$pvalue)]
rank$decentml <- decent.table.ml$gene[order(decent.table.ml$pvalue)]

sm.curves <- list()
sm.curves[['fpr']] <- sm.curves[['tpr']] <- sm.curves[['pr']] <- list()

for (met in c('er', 'zer', 'mnc', 'mst', 'scde', 'tasc', 'decent', 'decentml')){
  if (met != 'tasc') {
    sm.curves$fpr[[met]] <-  calcFPR(rank[[met]], ref500)
    sm.curves$tpr[[met]] <-  calcTPR(rank[[met]], ref500)
    sm.curves$pr[[met]] <-  calcPR(rank[[met]], ref500)
  } else {
    sm.curves$fpr[[met]] <- sm.curves$tpr[[met]] <- rep(0, nrow(data.obs))
    sm.curves$pr[[met]] <- rep(1, nrow(data.obs))
  }
}

rank.3g <- list()
rank.3g$er <- rownames(er.table.3g)[order(er.table.3g$PValue)]
rank.3g$zer <- rownames(zer.table.3g)[order(zer.table.3g$PValue)]
rank.3g$mnc <- rownames(mnc.table.3g)[order(mnc.table.3g$pval)]
rank.3g$mst <- mst.table.3g$primerid[order(mst.table.3g$`Pr(>Chisq)`)]
rank.3g$decent <- decent.table.3g$gene[order(decent.table.3g$pvalue)]

sm.curves.3g <- list()
sm.curves.3g[['fpr']] <- sm.curves.3g[['tpr']] <- sm.curves.3g[['pr']] <- list()

for (met in c('er', 'zer', 'mnc', 'mst', 'decent')){
  sm.curves.3g$fpr[[met]] <-  calcFPR(rank.3g[[met]], ref500.3g)
  sm.curves.3g$tpr[[met]] <-  calcTPR(rank.3g[[met]], ref500.3g)
  sm.curves.3g$pr[[met]] <-  calcPR(rank.3g[[met]], ref500.3g)
}

sm.logfc.tmm <- decent.table$logfc
sm.logfc.ml <- decent.table.ml$logfc
save(sm.curves, sm.curves.3g, sm.logfc.tmm, sm.logfc.ml, file='Soumillon_benckmark.rda')


### Type I error control ###
data.base <- data.obs[, cell.type == 'D3T0'] #
data.base <- data.base[rowSums(data.base) > 10, ] #

smp_size <- 400
er.ctrl <- zer.ctrl <- decent.ctrl <- mnc.ctrl <- mst.ctrl <- scde.ctrl <- decent.ctrl <- list()

for (i in 1:20) {
  set.seed(7*i)
  ind <- sample(seq_len(ncol(data.base)), size = smp_size) #
  data.split <- data.base[, ind]
  split <- c(rep(1, 200), rep(2, 200))
  
  er.ctrl[[i]] <- runEdgeR(data.split, split)
  zer.ctrl[[i]] <- runZinbwaveEdgeR(data.split, split, J = 1e7)
  mnc.ctrl[[i]] <- runMonocleDE(data.split, split)
  mst.ctrl[[i]] <- runMASTDE(data.split, split)
  scde.ctrl[[i]] <- runSCDE(data.split, split)
  scde.ctrl[[i]]$pval <- 2*pnorm(abs(scde.ctrl[[i]]$Z),lower.tail=F)
  decent.ctrl[[i]] <- decent(data.split, X = ~as.factor(split), normalize = 'TMM')
  message(i)
}

# qqplot
qq <- function(pvector, group, title="Quantile-quantile plot of p-values") {
  o = -log10(sort(pvector,decreasing=F)); print(length(o)); print(o[1:5]);
  e = -log10( 1:length(o)/length(o) )
  data <- data.frame(o = o, e = e)
  ggscatter(data = data, x = "e", y = "o",
            ggtheme = theme_bw(), title = title, size = 1,
            xlab = '', #expression(Expected~~-log[10](italic(p))),
            ylab = '')+ #expression(Observed~~-log[10](italic(p)))) +
    geom_abline(intercept=0,slope=1, col="#7AA6DCFF", linetype='dashed', size=0.6) +
    xlim(c(-log10(0.5), max(c(o, e)))) + ylim(c(-log10(0.5), max(c(o, e)))) +
    theme(title = element_text(size = 10))
}

p1 <-
ggarrange(qq(decent.ctrl[[8]]$pval[decent.ctrl[[8]]$pval > 1e-10], title = 'DECENT'),
          qq(mst.ctrl[[8]]$`Pr(>Chisq)`[mst.ctrl[[8]]$`Pr(>Chisq)` > 1e-10], title = 'MAST'),
          qq(scde.ctrl[[8]]$pval[scde.ctrl[[8]]$pval > 1e-10], title = 'SCDE'),
          qq(mnc.ctrl[[8]]$pval[mnc.ctrl[[8]]$pval > 1e-10], title = 'Monocle2'),
          qq(er.ctrl[[8]]$PValue[er.ctrl[[8]]$PValue > 1e-10], title='Z-edgeR'),
          qq(zer.ctrl[[8]]$PValue[zer.ctrl[[8]]$PValue > 1e-10], title='edgeR'),
          ncol = 3, nrow = 2)
p1.fpr.s <-
  annotate_figure(p1, 
                  left = text_grob(expression(Observed~~-log[10](italic(p))), color = "black", size = 10, rot = 90, vjust = 1),
                  bottom = text_grob(expression(Expected~~-log[10](italic(p))), color = "black", size = 10, vjust = -0.3))

# boxplot
p2.fpr.s <- 
ggboxplot(data = data.frame(FPR=c(sapply(1:20, function(i) mean(decent.ctrl[[i]]$pval < 0.05)),
                                        sapply(1:20, function(i) mean(mst.ctrl[[i]]$`Pr(>Chisq)` < 0.05)), 
                                        sapply(1:20, function(i) mean(scde.ctrl[[i]]$pval < 0.05)), 
                                        sapply(1:20, function(i) mean(mnc.ctrl[[i]]$pval < 0.05)), 
                                        sapply(1:20, function(i) mean(zer.ctrl[[i]]$PValue < 0.05)), 
                                        sapply(1:20, function(i) mean(er.ctrl[[i]]$PValue < 0.05))),
                                  Method=c(rep('DECENT', 20), 
                                           rep('MAST', 20), 
                                           rep('SCDE', 20), 
                                           rep('Monocle2', 20), 
                                           rep('Z-edgeR', 20), 
                                           rep('edgeR', 20))
                                  ), 
                x = 'Method', y = 'FPR', width = 0.7,
                xlab = '', ylab = 'Type I error rate')+ ylim(c(0,0.17))+
  theme_pubclean() + border()+
  theme(axis.text.x = element_text(angle=40, hjust = 1, size = 11)) +
  geom_hline(aes(yintercept = 0.05), color = '#A73030FF', linetype = "dashed", size=0.7)

save(p1.fpr.s, p2.fpr.s, file = 'fpr.plots.soumillon.rda')


# ggarrange(p1 + theme(plot.margin = margin(l = 5, t = 5)), 
#           p2 + theme(plot.margin = margin(l = 5, t = 10, r= 10, b=30)), 
#           nrow=2, ncol=1, labels = 'auto', font.label = list(size=16))

# 6 9
# 539 810