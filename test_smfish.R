library(MASS)
library(pscl)
library(AER)
library(ggpubr)

# downloaded from http://linnarssonlab.org/osmFISH/availability/
### Python code for parsing loom file ###
# import loompy
# ds = loompy.connect("/osmFISH_SScortex_mouse_all_cells.loom")
# ds.export(out_file='osmFISH.txt')

### Read data ###
full.table <- read.table('osmFISH.txt', as.is = T, sep = '\t')
data.obs <- as.matrix(full.table[14:46, 5:6476])
mode(data.obs) <- 'numeric'
rownames(data.obs) <- full.table$V2[14:46]
cluster <- as.factor(full.table[3, 5:6476])
names(cluster) <- colnames(data.obs) <- full.table[1, 5:6476]

rm(full.table)

### test over-dispersion ###

# 1. calculate size factor
sf.p4 <- colSums(data.obs[, cluster=='pyramidal L4'])
sf.p4 <- sf.p4/mean(sf.p4)
sf.om <- colSums(data.obs[, cluster=='Oligodendrocyte Mature'])
sf.om <- sf.om/mean(sf.om)
sf.iv <- colSums(data.obs[, cluster=='Inhibitory Vip'])
sf.iv <- sf.iv/mean(sf.iv)

# 2. fit Poisson GLM with size factors as offsets.
pois.p4 <- pois.om <- pois.iv <- list()
for (i in 1:33) {
  pois.p4[[i]] <- glm(formula = data.obs[i, cluster=='pyramidal L4'] ~ 1 + offset(log(sf.p4)), family = poisson)
  pois.om[[i]] <- glm(formula = data.obs[i, cluster=='Oligodendrocyte Mature'] ~ 1 + offset(log(sf.om)), family = poisson)
  pois.iv[[i]] <- glm(formula = data.obs[i, cluster=='Inhibitory Vip'] ~ 1 + offset(log(sf.iv)), family = poisson)
}

# 3. perform score test for overdispersion (Cameron 1990)
ct.pval.p4 <- sapply(1:33, function(i) AER::dispersiontest(pois.p4[[i]], trafo = 2, alternative = 'greater')$p.value)
ct.pval.om <- sapply(1:33, function(i) AER::dispersiontest(pois.om[[i]], trafo = 2, alternative = 'greater')$p.value)
ct.pval.iv <- sapply(1:33, function(i) AER::dispersiontest(pois.iv[[i]], trafo = 2, alternative = 'greater')$p.value)

### test zero-inflation ###

# 1. select three major cell types and calculate size factors
Pyramidal <- grepl('^[P|p]yramidal', cluster)
Oligodendrocyte <- grepl('^Oligodendrocyte', cluster)
Inhibitory <- grepl('^Inhibitory', cluster)
sf.p <- colSums(data.obs[, Pyramidal])/mean(colSums(data.obs[, Pyramidal]))
sf.o <- colSums(data.obs[, Oligodendrocyte])/mean(colSums(data.obs[, Oligodendrocyte]))
sf.i <- colSums(data.obs[, Inhibitory])/mean(colSums(data.obs[, Inhibitory]))

# 2. fit NB and ZINB GLM with size factors as offsets.
nb.p <- nb.o <- nb.i <- zinb.p <- zinb.o <- zinb.i <- list()
for (i in 1:33) {
  nb.p[[i]] <- glm.nb(formula = data.obs[i, Pyramidal] ~ 1 + offset(log(sf.p)))
  nb.o[[i]] <- glm.nb(formula = data.obs[i, Oligodendrocyte] ~ 1 + offset(log(sf.o)))
  nb.i[[i]] <- glm.nb(formula = data.obs[i, Inhibitory] ~ 1 + offset(log(sf.i)))
  zinb.p[[i]] <- pscl::zeroinfl(formula = data.obs[i, Pyramidal] ~ 1 + offset(log(sf.p))| 1, dist = "negbin")
  zinb.o[[i]] <- pscl::zeroinfl(formula = data.obs[i, Oligodendrocyte] ~ 1 + offset(log(sf.o)) | 1, dist = "negbin")
  zinb.i[[i]] <- pscl::zeroinfl(formula = data.obs[i, Inhibitory] ~ 1 + offset(log(sf.i)) | 1, dist = "negbin")
}

# 2. perform likelihood ratio test
lrts.zi.pval.p <- lrts.zi.pval.o <- lrts.zi.pval.i <- rep(NA,33)

for (i in 1:33) {
  lrts.zi.pval.p[i] <- exp(pchisq(zinb.p[[i]]$loglik*2 - nb.p[[i]]$twologlik, df = 1, lower.tail = F, log.p = T))/2
  lrts.zi.pval.o[i] <- exp(pchisq(zinb.o[[i]]$loglik*2 - nb.o[[i]]$twologlik, df = 1, lower.tail = F, log.p = T))/2
  lrts.zi.pval.i[i] <- exp(pchisq(zinb.i[[i]]$loglik*2 - nb.i[[i]]$twologlik, df = 1, lower.tail = F, log.p = T))/2
}

### Plotting ###
qq <- function(pvector, group, title="Quantile-quantile plot of p-values") {
  o = -log10(sort(pvector,decreasing=F))
  e = -log10( 1:length(o)/length(o) )
  data <- data.frame(o = o, e = e, significant = o > -log10(0.05/33))
  data$significant <- factor(data$significant, levels = c(T, F))
  ggscatter(data = data, x = "e", y = "o", color = "significant", palette = c("black", "gray60"),
            ggtheme = theme_bw(), title = title,
            xlab = expression(Expected~~-log[10](italic(p))),
            ylab = expression(Observed~~-log[10](italic(p)))) +
    geom_abline(intercept=0,slope=1, col="#A73030FF", linetype='dashed', size=0.6) + 
    geom_hline(yintercept = -log10(0.05/33), col='orange', linetype='longdash', size=0.6) +
    xlim(c(min(c(o, e)), max(c(o, e)))) + ylim(c(min(c(o, e)), max(c(o, e))))
}

qq0 <- function(pvector, group, title="Quantile-quantile plot of p-values") {
  o = -log10(sort(pvector,decreasing=F))
  e = -log10(c(1:(length(o)/2)/(length(o)/2)/2, rep(0.5, ceiling(length(o)/2))))
  data <- data.frame(o = o, e = e, significant = o > -log10(0.05/33), gene = rownames(data.obs)[order(pvector)])
  data$gene[!data$significant] <- ""
  data$significant <- factor(data$significant, levels = c(T, F))
  ggscatter(data = data, x = "e", y = "o", color = "significant", palette = c("black", "gray60"), 
            label = "gene", repel = T, font.label = list(size=12),
            ggtheme = theme_bw(), title = title,
            xlab = expression(Expected~~-log[10](italic(p))),
            ylab = expression(Observed~~-log[10](italic(p)))) +
    geom_abline(intercept=0,slope=1, col="#7AA6DCFF", linetype='dashed', size=0.6) + 
    geom_hline(yintercept = -log10(0.05/33), col='orange', linetype='longdash', size=0.6) +
    xlim(c(-log10(0.5), max(c(o, e)))) + ylim(c(-log10(0.5), max(c(o, e))))
}

ggarrange(qq(ct.pval.p4, title = 'Pyramidal L4'),
          qq(ct.pval.om, title = 'Oligodendrocyte Mature'),
          qq(ct.pval.iv, title = 'Inhibitory Vip'),
          qq0(lrts.zi.pval.p, title = 'Pyramidal'),
          qq0(lrts.zi.pval.o, title = 'Oligodendrocyte'), 
          qq0(lrts.zi.pval.i, title = 'Inhibitory'), 
          ncol = 3, nrow = 2, labels = c('a', '', '', 'b', '', ''), font.label = list(size=18),
          legend = 'right', common.legend = T)
# 8.5 5.7
# 750 502