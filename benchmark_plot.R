require(ggpubr)
require(ZIM)


load('scripts/Tung_benckmark.rda')
load('scripts/Soumillon_benckmark.rda')
load('scripts/Savas_benckmark.rda')
load('scripts/Chen_benckmark.rda')

load('scripts/fpr.plots.soumillon.rda')
load('scripts/fpr.plots.tung.rda')


# Main plot ROC
# ROC
plotPROC <- function(fpr, tpr, ymax, name) {
  ngene <- length(fpr)/7
  df <- data.frame(FPR = fpr, TPR = tpr)
  df$model <- c(rep('edgeR', ngene), rep('Z-edgeR', ngene), rep('Monocle2', ngene), 
                rep('MAST', ngene), rep('SCDE', ngene), rep('TASC', ngene), rep('DECENT', ngene))
  df$model <- factor(df$model, levels = c('DECENT', 'MAST', 'SCDE', 'Monocle2', 'TASC', 
                                          'Z-edgeR', 'edgeR'))
  return(ggline(df, "FPR", "TPR", plot_type = 'l', color = "model", palette = 'lancet', size = 0.4,
                numeric.x.axis = T, xlab = 'False positive rate', ylab = 'True positive rate',
                title = name) + scale_x_continuous(breaks = c(0, 0.05, 0.1), limits = c(0, 0.1)) +ylim(0,ymax))
}
proc.t <- plotPROC(unlist(t.curves$fpr[1:7]), unlist(t.curves$tpr[1:7]), 0.7, 'Tung')
proc.sm <- plotPROC(unlist(sm.curves$fpr[1:7]), unlist(sm.curves$tpr[1:7]), 0.67, 'Soumillon')
proc.sv <- plotPROC(unlist(sv.curves$fpr), unlist(sv.curves$tpr), 0.52, 'Savas')
proc.c <- plotPROC(unlist(c.curves$fpr), unlist(c.curves$tpr), 0.56, 'Chen')

ggarrange(proc.t, proc.sm, proc.sv, proc.c,
          ncol = 2, nrow = 2, labels = 'auto', font.label = list(size=18), common.legend = T, legend = 'bottom')
# 6 7
# 600 700

plotROC <- function(fpr, tpr, name) {
  ngene <- length(fpr)/7
  df <- data.frame(FPR = fpr, TPR = tpr)
  df$model <- c(rep('edgeR', ngene), rep('Z-edgeR', ngene), rep('Monocle2', ngene),
                rep('MAST', ngene), rep('SCDE', ngene), rep('TASC', ngene), rep('DECENT', ngene))
  df$model <- factor(df$model, levels = c('DECENT', 'MAST', 'SCDE', 'Monocle2', 'TASC',
                                          'Z-edgeR', 'edgeR'))
  return(ggline(df, "FPR", "TPR", plot_type = 'l', color = "model", palette = 'lancet', size = 0.4,
                numeric.x.axis = T, xlab = 'False positive rate', ylab = 'True positive rate',
                title = name))
}
roc.t <- plotROC(unlist(t.curves$fpr[1:7]), unlist(t.curves$tpr[1:7]), 'Tung')
roc.sm <- plotROC(unlist(sm.curves$fpr[1:7]), unlist(sm.curves$tpr[1:7]), 'Soumillon')
roc.sv <- plotROC(unlist(sv.curves$fpr), unlist(sv.curves$tpr), 'Savas')
roc.c <- plotROC(unlist(c.curves$fpr), unlist(c.curves$tpr), 'Chen')

ggarrange(roc.t, roc.sm, roc.sv, roc.c,
          ncol = 2, nrow = 2, labels = 'auto', font.label = list(size=18), common.legend = T, legend = 'bottom')

# Main plot FDR
plotFDR <- function(pr, name) {
  ngene <- length(pr)/7
  if (name == 'Tung') {
    dff <- data.frame(ng = rep(1:ngene, 7), FDR = pr)
  } else {
    dff <- data.frame(ng = c(rep(1:ngene, 5), rep(0, ngene), rep(1:ngene, 1)), FDR = pr)
  }
  dff$FDR <- (1-dff$FDR)
  dff$model <- c(rep('edgeR', ngene), rep('Z-edgeR', ngene), rep('Monocle2', ngene), 
                rep('MAST', ngene), rep('SCDE', ngene), rep('TASC', ngene), rep('DECENT', ngene))
  dff$model <- factor(dff$model, levels = c('DECENT', 'MAST', 'SCDE', 'Monocle2', 'TASC', 'Z-edgeR', 'edgeR'))
  return(ggline(dff, "ng", "FDR", plot_type = 'l', color = "model", palette = 'lancet', size = 0.4,
                numeric.x.axis = T, xlab = 'Number of discovered DEGs', ylab = 'False discovery rate',
                title = name) + xlim(0, 1500))
}
fdc.t <- plotFDR(unlist(t.curves$pr[1:7]), 'Tung')
fdc.sm <- plotFDR(unlist(sm.curves$pr[1:7]), 'Soumillon')
fdc.sv <- plotFDR(unlist(sv.curves$pr), 'Savas')
fdc.c <- plotFDR(unlist(c.curves$pr), 'Chen')

ggarrange(fdc.t, fdc.sm, fdc.sv, fdc.c, 
          ncol = 2, nrow = 2, labels = 'auto', font.label = list(size=18), common.legend = T, legend = 'bottom')
# 6 7
# 600 700

# Tung sp vs nsp
ngene <- length(t.curves$fpr$decent)
df.sp <- data.frame(FPR = unlist(t.curves$fpr[7:10]),
                    TPR = unlist(t.curves$tpr[7:10]))
df.sp$model <- c(rep('DECENT', ngene), rep('DECENT nsp 1x', ngene), rep('DECENT nsp 0.5x', ngene),
                 rep('DECENT nsp 1.5x', ngene))
df.sp$model <- factor(df.sp$model, levels = c('DECENT', 'DECENT nsp 1x', 'DECENT nsp 0.5x', 'DECENT nsp 1.5x'))
roc.sp <- ggline(df.sp, "FPR", "TPR", plot_type = 'l', color = "model", palette = 'lancet', size = 0.4,
               numeric.x.axis = T, xlab = 'False positive rate', ylab = 'True positive rate',
               title = 'pROC curve') + scale_x_continuous(breaks = c(0, 0.05, 0.1), limits = c(0, 0.1)) +ylim(0,0.7)

dff.sp <- data.frame(ng = rep(1:ngene, 4),
                     FDR = unlist(t.curves$pr[7:10]))
dff.sp$FDR <- (1-dff.sp$FDR)
dff.sp$model <- c(rep('DECENT', ngene), rep('DECENT nsp 1x', ngene), rep('DECENT nsp 0.5x', ngene),
                  rep('DECENT nsp 1.5x', ngene))
dff.sp$model <- factor(df.sp$model, levels = c('DECENT', 'DECENT nsp 1x', 'DECENT nsp 0.5x', 'DECENT nsp 1.5x'))
fdc.sp <- ggline(dff.sp, "ng", "FDR", plot_type = 'l', color = "model", palette = 'lancet', size = 0.4,
               numeric.x.axis = T, xlab = 'Number of discovered DEGs', ylab = 'False discovery rate',
               title = 'FDR curve') + xlim(0, 1500)

ggarrange(roc.sp, fdc.sp,
          ncol = 2, nrow = 1, labels = 'auto', font.label = list(size=18), common.legend = T, legend = 'bottom')

# 7 4
# 700 400


# Soumillon TMM vs ML
ngene <- length(sm.curves$fpr$decent)
df.n <- data.frame(FPR = unlist(sm.curves$fpr[7:8]), TPR = unlist(sm.curves$tpr[7:8]))
df.n$model <- c(rep('DECENT TMM', ngene), rep('DECENT ML', ngene))
df.n$model <- factor(df.n$model, levels = c('DECENT TMM', 'DECENT ML'))
roc.n <- ggline(df.n, "FPR", "TPR", plot_type = 'l', color = "model", palette = 'lancet', size = 0.4,
               numeric.x.axis = T, xlab = 'False positive rate', ylab = 'True positive rate',
               title = 'pROC curve') + scale_x_continuous(breaks = c(0, 0.05, 0.1), limits = c(0, 0.1)) + ylim(0,0.67)

box.n <- 
  ggboxplot(data = data.frame(logFC = c(sm.logfc.tmm, sm.logfc.ml),
                            normalization = factor(c(rep('TMM', ngene), rep('ML', ngene)), levels=(c('TMM', 'ML')))),
          x='normalization', y='logFC', xlab = '', ylab = 'Log fold-change', width = 0.4, notch = T, outlier.shape= '',
          color = c('#00468BFF', '#ED0000FF')) + geom_hline(yintercept = 0, linetype='dashed', size=0.7) + ylim(c(-2, 2.2))

ggarrange(roc.n, box.n,
          ncol = 2, nrow = 1, labels = 'auto', font.label = list(size=18), common.legend = T, legend = 'bottom')

# 6 3.6
# 600 360

# Soumillon three group
ngene <- length(sm.curves.3g$fpr$decent)
df.3g <- data.frame(FPR = unlist(sm.curves.3g$fpr), TPR = unlist(sm.curves.3g$tpr))
df.3g$model <- c(rep('edgeR', ngene), rep('Z-edgeR', ngene), rep('Monocle2', ngene), 
                 rep('MAST', ngene), rep('DECENT', ngene))
df.3g$model <- factor(df.3g$model, levels = c('DECENT', 'MAST', 'Monocle2', 'Z-edgeR', 'edgeR'))
roc.3g <- ggline(df.3g, "FPR", "TPR", plot_type = 'l', color = "model", palette = 'lancet', size = 0.4,
               numeric.x.axis = T, xlab = 'False positive rate', ylab = 'True positive rate',
               title = 'pROC curve') + scale_x_continuous(breaks = c(0, 0.05, 0.1), limits = c(0, 0.1)) + ylim(0,0.75)

dff.3g <- data.frame(ng = rep(1:ngene, 5), FDR = unlist(sm.curves.3g$pr))
dff.3g$FDR <- (1-dff.3g$FDR)
dff.3g$model <- c(rep('edgeR', ngene), rep('Z-edgeR', ngene), rep('Monocle2', ngene), 
                  rep('MAST', ngene), rep('DECENT', ngene))
dff.3g$model <- factor(dff.3g$model, levels = c('DECENT', 'MAST', 'Monocle2', 'Z-edgeR', 'edgeR'))
fdc.3g <- ggline(dff.3g, "ng", "FDR", plot_type = 'l', color = "model", palette = 'lancet', size = 0.4,
               numeric.x.axis = T, xlab = 'Number of discovered DEGs', ylab = 'False discovery rate',
               title = 'FDR curve') + xlim(0, 1500)

ggarrange(roc.3g, fdc.3g,
          ncol = 2, nrow = 1, labels = 'auto', font.label = list(size=18), common.legend = T, legend = 'bottom')
# 7 4
# 700 400

# # Type I error box plot
# ggarrange(p2.fpr.t, p2.fpr.s, ncol=2, nrow=1,
#           widths = c(8, 6.5), labels = 'auto', font.label = list(size=16))
# # 7 4
# # 700 400
# 
# # Type I error qq plot
# ggarrange(p1.fpr.t + theme(plot.margin = margin(t=10)),
#           p1.fpr.s + theme(plot.margin = margin(t=10, r=75)), 
#           ncol=1, nrow=2, labels = 'auto', font.label = list(size=16))

# Type I error plots
ac <- ggarrange(
  p1.fpr.t + theme(plot.margin = margin(l = 5, t = 5, r = 10)),
  p2.fpr.t + theme(plot.margin = margin(l = 5, t = 10, r = 30, b = 30)), 
  ncol=1, nrow=2, labels=c('a', 'c'), font.label = list(size=16))
bd <- ggarrange(
  p1.fpr.s + theme(plot.margin = margin(l = 5, t = 5)),
  p2.fpr.s + theme(plot.margin = margin(l = 5, t = 10, r = 20, b = 35)),
  ncol=1, nrow=2, labels=c('b', 'd'), font.label = list(size=16))
ggarrange(ac, bd, 
          ncol=2, nrow=1, widths = c(8, 6))

# 11 7.6