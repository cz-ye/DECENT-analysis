require(DECENT)
require(edgeR)
require(ggpubr)
source('func_de_methods.R')

### Read data ###

load('DECENT-sim/simdata_ZINB_BB_Tung6_OD_kb.RData')

data.obs <- data.obs[rowSums(data.obs > 0) >= 3,]
sp.obs <- ercc.obs[rowSums(ercc.obs) != 0, ]
sp.nom <- ercc.true[rowSums(ercc.obs) != 0]

# prepare data for for TASC
sp.len <- read.table('ercc_length.txt', as.is = T, row.names = 1)
sp.len <- sp.len[match(rownames(sp.obs), rownames(sp.len)), ]

sp <- as.data.frame(sp.nom)
sp$len <- sp.len
sp$count <- apply(sp.obs, 1, function(x) paste(x, collapse = ','))
write.table(sp, col.names=F, quote=F, file=paste0('tasc_sim/ercc.txt'), sep='\t')

design <- rbind(paste(rep(1, 500),collapse = ','), paste(ctype-1, collapse = ','))
rownames(design) <- c('intercept', 'group')
write.table(design, col.names=F, quote=F, file=paste0('tasc_sim/x.txt'), sep='\t')

data.obs <- data.obs[rowSums(data.obs)!=0, ]
endo <- apply(data.obs, 1, function(x) paste(x, collapse = ','))
write.table(endo, col.names=F, quote=F, file=paste0('tasc_sim/y.txt'), sep='\t')


### single-cell DE analysis ###
er.table <- runEdgeR(data.obs, ctype)
zer.table <- runZinbwaveEdgeR(data.obs, ctype, J = 1e8)
mnc.table <- runMonocleDE(data.obs, ctype)
mst.table <- runMASTDE(data.obs, ctype)
scde.table <- runSCDE(data.obs, ctype)

# DECENT
decent.table.nsp <- decent(data.obs, ~cell.type, CE.range = range(CE), dir = 'sim/nsp/')
decent.table.nsp.0.5x <- decent(data.obs, ~cell.type, CE.range = 0.5*range(CE), 'sim/nsp_0.5x/')
decent.table.nsp.1.5x <- decent(data.obs, ~cell.type, CE.range = 1.5*range(CE), 'sim/nsp_1.5x/')
decent.table <- decent(data.obs, ~cell.type, use.spikes = T, spikes = sp.obs, spike.conc = sp.nom, tau.est = 'spikes',
                       s.imputed = T, E.imputed = T, dir = 'sim/sp/')
decent.table2 <- decent(data.obs, ~cell.type, use.spikes = T, spikes = sp.obs, spike.conc = sp.nom, dir = 'sim/sp_tau/')

out.em <- readRDS('sim/sp/decent.noDE.rds')
out.lrt <- readRDS('sim/sp/decent.lrt.rds')

# TASC is not a R package. It was run separately and has no tunable parameter. We just provide the result here.
tasc.table <- read.table('tasc/sim2.de', row.names = 1, header = T, sep = '\t')
tasc.table[setdiff(rownames(data.obs), rownames(tasc.table)), ] <- rep(1, 12)

# reference DEG
ref <- rownames(data)[DE.gene]
ref <- intersect(ref, rownames(data.obs))

# calculate and save curves
rank <- list()
rank$er <- rownames(er.table)[order(er.table$PValue)]
rank$zer <- rownames(zer.table)[order(zer.table$PValue)]
rank$mnc <- rownames(mnc.table)[order(mnc.table$pval)]
rank$mst <- mst.table$primerid[order(mst.table$`Pr(>Chisq)`)]
rank$scde <- rownames(scde.table)[order(abs(scde.table$Z), decreasing = T)]
rank$tasc <- rownames(tasc.table)[order(tasc.table$group.LRT.PVal)]
rank$decent <- decent.table$gene[order(decent.table$pvalue)]
rank$decent2 <- decent.table2$gene[order(decent.table2$pvalue)]
rank$decentnsp <- decent.table.nsp$gene[order(decent.table.nsp$pvalue)]
rank$decentnspl <- decent.table.nsp.0.5x$gene[order(decent.table.nsp.0.5x$pvalue)]
rank$decentnsph <- decent.table.nsp.1.5x$gene[order(decent.table.nsp.1.5x$pvalue)]

sim.curves <- list()
sim.curves[['fpr']] <- sim.curves[['tpr']] <- sim.curves[['pr']] <- list()

for (met in c('er', 'zer', 'mnc', 'mst', 'scde', 'tasc', 'decent', 'decent2', 'decentnsp', 'decentnspl', 'decentnsph')){
  sim.curves$fpr[[met]] <-  calcFPR(rank[[met]], ref)
  sim.curves$tpr[[met]] <-  calcTPR(rank[[met]], ref)
  sim.curves$pr[[met]] <-  calcPR(rank[[met]], ref)
}

f <- rownames(data) %in% rownames(data.obs)

p.param <- plotParamEst(out.lrt$par.DE, out.lrt$par.noDE, out.em$est.sf, out.em$CE, pi0[f],
                   mu[f], fc[f, ], disp[f], sf, CE, rowMeans(data.obs))
# 700 500


# Infer pre-dropout distribution by single imputation
data.simp <- readRDS('sim/sp/single.imputed.rds')
data.imp2 <- readRDS('sim/sp/mean.imputed.rds')

plotYEst <- function(data.simp, data.imp, y) {
  p1 <- ggscatter(data = data.frame(x = rowMeans(data.simp==0),
                                    y = rowMeans(y==0)),
                  x = 'x', y = 'y', col = adjustcolor('royalblue4', alpha.f = 0.25), size = 0.4,
                  xlab = 'Estimated proportion',
                  ylab = 'True proportion', title = 'Pre-dropout zeros') + 
    xlim(c(0,1)) + ylim(c(0,1)) +
    geom_abline(slope = 1, intercept = 0, col = 'darkgray', size=0.7) + border() +
    font("title", size = 16) + font("xlab", size = 14) + font("ylab", size = 14)
  p2 <- ggscatter(data = data.frame(x = log10(rowVars(data.simp)),
                                    y = log10(rowVars(y))),
                  x = 'x', y = 'y', col = adjustcolor('royalblue4', alpha.f = 0.25), size = 0.4,
                  xlab = expression(Log[10]~estimated~variance),
                  ylab = expression(Log[10]~true~variance), title = 'Pre-dropout variation') + 
    xlim(log10(c(min(rowVars(y)),max(rowVars(y))))) + ylim(log10(c(min(rowVars(y)),max(rowVars(y))))) +
    geom_abline(slope = 1, intercept = 0, col = 'darkgray', size=0.7) + border() +
    font("title", size = 16) + font("xlab", size = 14) + font("ylab", size = 14)
  set.seed(7)
  ind <- sample(sum(y>0), sum(y>0)*0.05)
  p3 <- ggscatter(data = data.frame(x = log10(c(data.imp[y>0])[ind]),
                                    y = log10(c(y[y>0]))[ind]),
                  x = 'x', y = 'y', col = adjustcolor('royalblue4', alpha.f = 0.1), size = 0.4,
                  xlab = expression(Log[10]~expected~count),
                  ylab = expression(Log[10]~true~count), title = 'Pre-dropout counts') + 
    geom_abline(slope = 1, intercept = 0, col = 'darkgray', size=0.7) + 
    font("title", size = 16) + font("xlab", size = 14) + font("ylab", size = 14)
  ggarrange(p1, p2, p3 + theme(plot.margin = margin(t=5)), 
            ncol = 3, nrow = 1, labels = 'auto', font.label = list(size=20))
}

pdf('fig.2.pdf', width=10, height = 10)
print(plotYEst(data.simp, data.imp2, data[f, ]))
dev.off()

png('fig.2.png', width=640, height = 640)
print(plotYEst(out.lrt$par.noDE, data.imp2, data[f, ]))
dev.off()

# 9 3.5
# 900 350

# main ROC FDR
ngene <- length(sim.curves$fpr$decent)
df <- data.frame(FPR = unlist(sim.curves$fpr[1:7]),
                 TPR = unlist(sim.curves$tpr[1:7]))
df$model <- c(rep('edgeR', ngene), rep('Z-edgeR', ngene), rep('Monocle2', ngene), 
              rep('MAST', ngene), rep('SCDE', ngene), rep('TASC', ngene), rep('DECENT', ngene))
df$model <- factor(df$model, levels = c('DECENT', 'MAST', 'SCDE', 'Monocle2', 'TASC',
                                        'Z-edgeR', 'edgeR'))
proc1 <- ggline(df, "FPR", "TPR", plot_type = 'l', color = "model", palette = 'lancet', size = 0.4,
               numeric.x.axis = T, xlab = 'False positive rate', ylab = 'True positive rate',
               title = 'pROC curve') + scale_x_continuous(breaks = c(0, 0.05, 0.1), limits = c(0, 0.1)) +ylim(0,0.47)
roc1 <- ggline(df, "FPR", "TPR", plot_type = 'l', color = "model", palette = 'lancet', size = 0.4,
               numeric.x.axis = T, xlab = 'False positive rate', ylab = 'True positive rate',
               title = 'ROC curve')

dff <- data.frame(ng = rep(1:ngene, 7),
                  FDR = unlist(sim.curves$pr[1:7]))
dff$FDR <- (1-dff$FDR)
dff$model <- c(rep('edgeR', ngene), rep('Z-edgeR', ngene), rep('Monocle2', ngene), 
               rep('MAST', ngene), rep('SCDE', ngene), rep('TASC', ngene), rep('DECENT', ngene))
dff$model <- factor(df$model, levels = c('DECENT', 'MAST', 'SCDE', 'Monocle2', 'TASC',
                                         'Z-edgeR', 'edgeR'))
fdc1 <- ggline(dff, "ng", "FDR", plot_type = 'l', color = "model", palette = 'lancet', size = 0.3,
               numeric.x.axis = T, xlab = 'Number of discovered DEGs', ylab = 'False discovery rate',
               title = 'FDR curve') + xlim(0, 1500)

ggarrange(proc1, fdc1,
          ncol = 2, nrow = 1, labels = 'auto', font.label = list(size=18), common.legend = T, legend = 'bottom')
# 7 4.3
# 700 430
ggarrange(roc1, ncol=1, nrow = 1, common.legend = T, legend = 'bottom')
# 4.5 5
# 450 500

df <- data.frame(FPR = unlist(sim.curves$fpr[7:11]),
                 TPR = unlist(sim.curves$tpr[7:11]))
df$model <- c(rep('DECENT sp', ngene), rep('DECENT sp tau', ngene), rep('DECENT nsp 1x', ngene), rep('DECENT nsp 0.5x', ngene), 
              rep('DECENT nsp 1.5x', ngene))
df$model <- factor(df$model, levels = c('DECENT sp', 'DECENT sp endo', 'DECENT nsp 1x', 'DECENT nsp 0.5x', 'DECENT nsp 1.5x'))
proc2 <- ggline(df, "FPR", "TPR", plot_type = 'l', color = "model", palette = 'lancet', size = 0.4,
               numeric.x.axis = T, xlab = 'False positive rate', ylab = 'True positive rate',
               title = 'pROC curve') + scale_x_continuous(breaks = c(0, 0.05, 0.1), limits = c(0, 0.1)) +ylim(0,0.47)
proc2

dff <- data.frame(ng = rep(1:ngene, 5),
                  FDR = unlist(sim.curves$pr[7:11]))
dff$FDR <- (1-dff$FDR)
dff$model <- c(rep('DECENT sp', ngene), rep('DECENT sp endo', ngene), rep('DECENT nsp 1x', ngene), rep('DECENT nsp 0.5x', ngene), 
               rep('DECENT nsp 1.5x', ngene))
dff$model <- factor(df$model, levels = c('DECENT sp', 'DECENT sp tau', 'DECENT nsp 1x', 'DECENT nsp 0.5x', 'DECENT nsp 1.5x'))
fdc2 <- ggline(dff, "ng", "FDR", plot_type = 'l', color = "model", palette = 'lancet', size = 0.3,
               numeric.x.axis = T, xlab = 'Number of discovered DEGs', ylab = 'False discovery rate',
               title = 'FDR curve') + xlim(0, 1500)
fdc2
ggarrange(proc2, fdc2,
          ncol = 2, nrow = 1, labels = 'auto', font.label = list(size=18), common.legend = T, legend = 'bottom')
