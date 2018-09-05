require(DECENT)
require(edgeR)
require(ggpubr)
source('func_de_methods.R')

### Read data ###
# downloaded from https://github.com/jdblischak/singleCellSeq
sc <- read.table('molecule-filtered.txt', sep='\t', header = T, row.names = 1)
# downloaded from GSE77288
bulk <- read.table('GSE77288_reads-raw-bulk-per-sample.txt', sep='\t', header = T)
bulk <- t(bulk)
bulk <- bulk[-(1:3), ]
mode(bulk) <- 'numeric'

### Subset single-cell data ###
ind <- substr(colnames(sc), 1, 7)
data.obs <- cbind(sc[, ind == 'NA19101'], sc[, ind == 'NA19239'])
cell.type <- c(rep(sum(ind == 'NA19101')), rep(sum(ind == 'NA19239')))


### Derive bulk reference gene list ###
# limma-voom (v3.36.1)
cpm <- edgeR::cpm(bulk)
keep.exprs <- rowSums(cpm>1)>=3
bulk <- bulk[keep.exprs, ]
group <- factor(c(1, 1, 1, 2, 2, 2))
y <- DGEList(counts = bulk[, 4:9], group = group)
y <- calcNormFactors(y)
y <- estimateDisp(y, model.matrix(~group))
v <- voom(y, model.matrix(~group), plot=TRUE)
vfit <- lmFit(v, model.matrix(~group))
vfit <- contrasts.fit(vfit, coefficients =2)
efit <- eBayes(vfit)
plotSA(efit, main="Final model: Mean−variance trend")
lv.table <- topTable(efit, number = Inf)
ref500 <- rownames(lv.table)[order(lv.table$P.Value)][1:500]
ref500 <- intersect(rownames(data.obs), ref500)


### single-cell DE analysis ###
# ZINB-WaVE, SCDE and DECENT are recommended to be run on server

er.table <- runEdgeR(data.obs, cell.type)
zer.table <- runZinbwaveEdgeR(data.obs, cell.type, J = 1e5)
mnc.table <- runMonocleDE(data.obs, cell.type)
mst.table <- runMASTDE(data.obs, cell.type)
scde.table <- runSCDE(data.obs, cell.type)

# DECENT w/o spike-ins
decent.table.nsp <- decent(data.obs, ~cell.type, tau.global = F, dir = 'tung/nsp')
decent.table.nsp.0.5x <- decent(data.obs, ~cell.type, CE.range = c(0.01, 0.05), tau.global = F, dir = 'tung/nsp_0.5x')
decent.table.nsp.1.5x <- decent(data.obs, ~cell.type, CE.range = c(0.03, 0.15), tau.global = F, dir = 'tung/nsp_1.5x')
  
# DECETN w/ spike-ins, get spike-in information first
# available at https://www.thermofisher.com/order/catalog/product/4456740
sp.ref <- read.table("./cms_095046.txt", sep='\t', row.names=1, header=T) 
sp.obs <- data.obs[grepl('^ERCC', rownames(data.obs)), ]
sp.nom <- getSpikeInMolecule(sp.ref$concentration.in.Mix.1..attomoles.ul., 13, 50000)
names(sp.nom) <- sp.ref$ERCC.ID
sp.nom <- sp.nom[match(rownames(sp.obs), names(sp.nom))]
decent.table <- decent(data.obs, ~cell.type, tau.global = F, use.spikes = T, spikes = sp.obs, spike.conc = sp.nom, dir = 'tung/sp')

# TASC is not a R package. It was run separately and has no tunable parameter. We just provide the result here.
tasc.table <- read.table('tasc/tung.de', row.names = 1, header = T, sep = '\t')

# calculate and save curves
rank <- list()
rank$er <- rownames(er.table)[order(er.table$PValue)]
rank$zer <- rownames(zer.table)[order(zer.table$PValue)]
rank$mnc <- rownames(mnc.table)[order(mnc.table$pval)]
rank$mst <- mst.table$primerid[order(mst.table$`Pr(>Chisq)`)]
rank$scde <- rownames(scde.table)[order(abs(scde.table$Z), decreasing = T)]
rank$tasc <- rownames(tasc.table)[order(tasc.table$group.LRT.PVal)]
rank$decent <- decent.table$gene[order(decent.table$pvalue)]
rank$decentnsp <- decent.table.nsp$gene[order(decent.table.nsp$pvalue)]
rank$decentnspl <- decent.table.nsp.0.5x$gene[order(decent.table.nsp.0.5x$pvalue)]
rank$decentnsph <- decent.table.nsp.1.5x$gene[order(decent.table.nsp.1.5x$pvalue)]

t.curves <- list()
t.curves[['fpr']] <- t.curves[['tpr']] <- t.curves[['pr']] <- list()

for (met in c('er', 'zer', 'mnc', 'mst', 'scde', 'tasc', 'decent', 'decentnsp', 'decentnspl', 'decentnsph')){
    t.curves$fpr[[met]] <-  calcFPR(rank[[met]], ref500)
    t.curves$tpr[[met]] <-  calcTPR(rank[[met]], ref500)
    t.curves$pr[[met]] <-  calcPR(rank[[met]], ref500)
}

save(t.curves, file='Tung_benckmark.rda')


### Before and after dropout ###

obsVsInf <- function(g, data.obs, data.imp, cell.type) {
  cell.type <- as.factor(cell.type)
  levels(cell.type) <- c('NA19101', 'NA19239')
  cell.type.imp <- rep(cell.type, each=50)
  p1 <- ggviolin(data = data.frame(y = c(data.obs[g, cell.type=='NA19101'], data.obs[g, cell.type=='NA19239']),
                                   x = cell.type),
                 y="y", x="x", fill = "#00AFBB", trim = T,
                 ylab='Normalized expression', xlab = 'Observed',
                 add = c("mean_sd"), add.params = list(size=0.4)) +
    theme(axis.text = element_text(size=10), axis.title = element_text(size=11))
  p2 <- ggviolin(data = data.frame(y = c(data.imp[g, cell.type.imp=='NA19101'], data.imp[g, cell.type.imp=='NA19239']),
                                   x = cell.type.imp),
                 y="y", x="x", fill = "#E7B800", trim = T,
                 ylab='Normalized expression', xlab = 'Inferred pre-dropout',
                 add = c("mean_sd"), add.params = list(size=0.4)) +
    theme(axis.text = element_text(size=10), axis.title = element_text(size=11))
  p <- ggarrange(p1, p2, ncol=2, nrow=1, common.legend = T, legend = 'none')
  annotate_figure(p, top = text_grob(g, color = "black", size = 14, vjust = 0.3))
}

# Select genes
u.gene <- rank.neo1b[(rank.neo1b %in% intersect(rank.neo1b[1:500], ref500)) & 
                       (!rank.neo1b %in% Reduce(x = list(rank.mst[1:500], rank.er[1:500]), f = union))]
u.gene <- u.gene[c(1, 4, 5, 6, 7, 8, 9, 10, 12, 15)]

# draw samples from pre-dropout distribution
set.seed(7)

out.lrt <- readRDS('tung/decent.lrt.rds')
out.em <- readRDS('tung/decent.noDE.rds')
data.imp <- t(
  sapply(which(rownames(data.obs) %in% u.gene), function(i) {
    if((i%%500)==1) message(i)
    return(c(sapply(1:ncol(data.obs), function(j)
      rzinb(50, k = exp(-out.lrt$par.DE[i,4]), 
            lambda = exp(out.lrt$par.DE[i,2] + (cell.type-1)[j]*out.lrt$par.DE[i,3]),
            omega = 1/(1+exp(-out.lrt$par.DE[i,1]))))))
  }))
rownames(data.imp) <- rownames(data.obs)[which(rownames(data.obs) %in% u.gene)]

# plot
lib.sf <- colSums(data.obs)
lib.sf <- lib.sf/mean(lib.sf)

vios <- list()
for (g in u.gene) {
  vios[[g]] <- obsVsInf(g, t(apply(data.obs, 1, function(x)x/lib.sf)),
                        data.mimp, cell.type)
}
vios[['ncol']] <- 2; vios[['nrow']] <- 4
do.call(ggarrange, args = vios[!names(vios) %in% c('ENSG00000132406', 'ENSG00000183688')])
# 8.5 9
# 800 850
vios <- vios[names(vios) %in% c('ENSG00000132406', 'ENSG00000183688')]
vios[['ncol']] <- 1; vios[['nrow']] <- 2
do.call(ggarrange, args = vios)
# 5 6
# 500 600


### Type I error control ###
data.split <- data.obs[, cell.type == 1] #
data.split <- data.split[rowSums(data.split) > 10, ] #

smp_size <- floor(0.5 * ncol(data.split))
er.ctrl <- zer.ctrl <- decent.ctrl <- mnc.ctrl <- mst.ctrl <- scde.ctrl <- decent.ctrl <- decent.ctrl.nsp <- list()

sp.obs.ctrl <- data.split[grepl('^ERCC', rownames(data.split)), ]
names(sp.nom.ctrl) <- sp.ref$ERCC.ID
sp.nom.ctrl <- sp.nom.ctrl[match(rownames(sp.obs.ctrl), names(sp.nom.ctrl))]

# spike-in filtering
f <- rowMeans(sp.obs.ctrl) < sp.nom.ctrl & sp.nom.ctrl > 0.05
sp.obs.ctrl <- sp.obs.ctrl[f, ]
sp.nom.ctrl <- sp.nom.ctrl[f]

for (i in 1:20) {
  set.seed(7*i)
  ind <- sample(seq_len(ncol(data.split)), size = smp_size)
  split <- rep(1, ncol(data.split))
  split[ind] <- 2
  
  er.ctrl[[i]] <- runEdgeR(data.split, split)
  zer.ctrl[[i]] <- runZinbwaveEdgeR(data.split, split, J = 1e7)
  mnc.ctrl[[i]] <- runMonocleDE(data.split, split)
  mst.ctrl[[i]] <- runMASTDE(data.split, split)
  scde.ctrl[[i]] <- runSCDE(data.split, split)
  scde.ctrl[[i]]$pval <- 2*pnorm(abs(scde.ctrl[[i]]$Z),lower.tail=F)
  decent.ctrl[[i]] <- decent(data.split, X = ~as.factor(split), tau.global = F, use.spikes = T,
                             spikes = sp.obs, spike.conc = sp.conc)
  decent.ctrl.nsp[[i]] <- decent(data.split, X = ~as.factor(split), tau.global = F)
}
# Again TASC was run separately.
tasc.ctrl <- readRDS('tasc/tung.ctrl.rds')

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
ggarrange(qq(out.ctrl.sp[[11]]$pval[out.ctrl.sp[[11]]$pval > 1e-10], title = 'DECENT sp'),
          qq(out.ctrl.sp[[11]]$pval[out.ctrl.sp[[11]]$pval > 1e-10], title = 'DECENT nsp'),
          qq(mst.ctrl[[11]]$`Pr(>Chisq)`[mst.ctrl[[11]]$`Pr(>Chisq)` > 1e-10], title = 'MAST'),
          qq(scde.ctrl[[11]]$pval[scde.ctrl[[11]]$pval > 1e-10], title = 'SCDE'),
          qq(mnc.ctrl[[11]]$pval[mnc.ctrl[[11]]$pval > 1e-10], title = 'Monocle2'),
          qq(tasc.ctrl[[11]]$group.LRT.PVal[tasc.ctrl[[11]]$group.LRT.PVal > 1e-10], title='TASC'),
          qq(er.ctrl[[11]]$PValue[er.ctrl[[11]]$PValue > 1e-10], title='Z-edgeR'),
          qq(zer.ctrl[[11]]$PValue[zer.ctrl[[11]]$PValue > 1e-10], title='edgeR'),
          ncol = 4, nrow = 2)

p1.fpr.t <-
  annotate_figure(p1, 
                  left = text_grob(expression(Observed~~-log[10](italic(p))), color = "black", size = 10, rot = 90, vjust = 1),
                  bottom = text_grob(expression(Expected~~-log[10](italic(p))), color = "black", size = 10, vjust = -0.3))

# boxplot
p2.fpr.t <- ggboxplot(data = data.frame(FPR=c(sapply(1:20, function(i) mean(out.ctrl.sp[[i]]$pval < 0.05)),
                                        sapply(1:20, function(i) mean(out.ctrl.nsp[[i]]$pval < 0.05)),
                                        sapply(1:20, function(i) mean(mst.ctrl[[i]]$`Pr(>Chisq)` < 0.05)), 
                                        sapply(1:20, function(i) mean(scde.ctrl[[i]]$pval < 0.05)), 
                                        sapply(1:20, function(i) mean(mnc.ctrl[[i]]$pval < 0.05)), 
                                        sapply(1:20, function(i) mean(tasc.ctrl[[i]]$group.LRT.PVal < 0.05)), 
                                        sapply(1:20, function(i) mean(zer.ctrl[[i]]$PValue < 0.05)), 
                                        sapply(1:20, function(i) mean(er.ctrl[[i]]$PValue < 0.05))),
                                  Method=c(rep('DECENT sp', 20), 
                                           rep('DECENT nsp', 20), 
                                           rep('MAST', 20), 
                                           rep('SCDE', 20), 
                                           rep('Monocle2', 20), 
                                           rep('TASC', 20), 
                                           rep('Z-edgeR', 20), 
                                           rep('edgeR', 20))
                                  ), 
                x = 'Method', y = 'FPR', width = 0.7,
                xlab = '', ylab = 'Type I error rate')+ ylim(c(0,0.13))+
  theme_pubclean() + border()+
  theme(axis.text.x = element_text(angle=40, hjust = 1, size = 11)) +
  geom_hline(aes(yintercept = 0.05), color = '#A73030FF', linetype = "dashed", size=0.7)

save(p1.fpr.t, p2.fpr.t, file = './scripts/fpr.plots.tung.rda')

# ggarrange(p1 + theme(plot.margin = margin(l = 5, t = 5)), 
#           p2 + theme(plot.margin = margin(l = 5, t = 10, r= 10, b=30)), 
#           nrow=2, ncol=1, labels = 'auto', font.label = list(size=16))
