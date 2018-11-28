require(DECENT)
require(DECENTNB)
library(ggpubr)
library(matrixStats)
library(VGAM)
source('func_de_methods.R')

### Read data ###
# Zeisel
# available from http://linnarssonlab.org/cortex/
endo.table <- read.table("./expression_mRNA_17-Aug-2014.txt", sep='\t', sep='\t', fill = T, as.is = T)
sp.table <- read.table("./expression_spikes_17-Aug-2014.txt", sep='\t', sep='\t', fill = T, as.is = T)
info <- endo.table[1:10, -c(1, 2)]
rownames(info) <- endo.table[1:10, 2]
endo.table <- endo.table[-1:-11, ]
rownames(endo.table) <- endo.table[, 1]
endo.table <- endo.table[, -c(1,2)]
sp.table <- sp.table[-1:-11, ]
rownames(sp.table) <- sp.table[, 1]
sp.table <- sp.table[, -c(1,2)]

# Lots of cells with extremely large proportion of spike-in counts. Remove those to alleviate.
sp.endo.ratio <- colSums(sp.table)/colSums(endo.table)
endo.table <- endo.table[, sp.endo.ratio < 1]
sp.table <- sp.table[, sp.endo.ratio < 1]
info <- info[sp.endo.ratio < 1, ]

endo.table <- endo.table[rowSums(endo.table) !=0, ]
cut <- median(log(rowSums(endo.table))) - mad(log(rowSums(endo.table)))
endo.table <- endo.table[log(rowSums(endo.table)) > cut, ]
data.z <- endo.table[, info$level1class=='pyramidal CA1']
sp.obs.z <- sp.table[, info$level1class=='pyramidal CA1']

sp.ref <- read.table("cms_095046.txt", sep='\t', row.names=1, header=T)

sp.nom.z <- getSpikeInMolecule(sp.ref$concentration.in.Mix.1..attomoles.ul., 9, 20000)
names(sp.nom.z) <- sp.ref$ERCC.ID
sp.nom.z <- sp.nom.z[match(rownames(sp.obs.z), names(sp.nom.z))]

set.seed(7)
ind <- sample(seq_len(ncol(data.z)), size = floor(0.5 * ncol(data.z)))
cell.type <- rep(1, ncol(data.z))
cell.type[ind] <- 2

out <- decent(data.obs = as.matrix(data.z), X = ~as.factor(cell.type), tau.global = F, tau.init = c(-5, -1),
              use.spikes = T, spikes = as.matrix(sp.obs), spike.conc = sp.ref,
              dir = './zeisel/')
out <- decentnb(data.obs = as.matrix(data.z), X = ~as.factor(cell.type), tau.global = F, tau.init = c(-5, -1),
              use.spikes = T, spikes = as.matrix(sp.obs), spike.conc = sp.ref,
              dir = './zeisel_nb/')

out.nb.z <- readRDS('zeisel_nb/decent.lrt.rds')
out.em.nb.z <- readRDS('zeisel_nb/decent.noDE.rds')
out.zinb.z <- readRDS('zeisel/decent.lrt.rds')
out.em.zinb.z <- readRDS('zeisel/decent.noDE.rds')

# Tung load the data and DECENT output as obatained in 'Tung.R'
data.t <- readRDS('tung/data.obs.rds')
ct <- readRDS('tung/ct.rds')
data.t <- data.t[, ct == 2]
data.t <- data.t[rowSums(data.t) > 10, ]

out.nb.t <- readRDS('tung_nb/decent.lrt.rds')
out.em.nb.t<- readRDS('tung_nb/decent.noDE.rds')
out.zinb.t <- readRDS('tung/ctrl_sp/1/decent.lrt.rds')
out.em.zinb.t<- readRDS('tung/ctrl_sp/1/decent.noDE.rds')

### Tung ###
# 1. Show evidence of overdispersion
var.gene <- rowVars(data.t)
mean.gene <- rowMeans(data.t)
CE <- out.em.nb.t$CE
sf <- out.em.nb.t$est.sf
mu <- exp(out.nb.t$par.noDE[,2])
psi <- exp(out.nb.t$par.noDE[,3])
tau0 <- out.em.nb.t$tau0
tau1 <- out.em.nb.t$tau1
rho <- rep(1, nrow(data.t)) %o% tau0 + log(mu) %o% (tau1*sf)
rho <- 1/(1+exp(-rho))

# mean-variance relationship
the.conc <- 10^(seq(log10(min(mu))-0.2, log10(max(mu))+0.2, len=10000))
rho0 <- rep(1, length(the.conc)) %o% tau0 + log(the.conc) %o% (tau1*sf)
rho0 <- 1/(1+exp(-rho0))

p1.t <- 
  ggscatter(data = data.frame(var=log10(var.gene), mean=log10(mean.gene)), x='mean', y='var', 
          ylab ='Log10 observed gene-wise variance', xlab='Log10 observed gene mean count', size=0.1) + 
  geom_line(data = data.frame(y=log10(mean(CE*sf) * the.conc + var(CE*sf) * the.conc^2 + 
                                        rowMeans((the.conc^2 %o% (CE*(1-CE) * sf^2)) * rho0)),
                              x=log10(mean(CE*sf) * the.conc)), aes(x, y), col='red', size=1)

'limegreen'

# observed vs expected
p2.t <- 
ggscatter(data = data.frame(x=log10(mean(CE*sf)*mu + var(CE*sf)*mu^2 + rowMeans((mu^2 %o% (CE*(1-CE) * sf^2)) * rho)), 
                            y=log10(var.gene)), x='x', y='y', 
          xlab ='Log10 expected variance', ylab='Log10 observed variance', 
          size=0.1, title = expression('Pre-dropout dispersion'~psi[i]~'='~0)) +
  geom_abline(slope = 1, intercept = 0, color = 'gray', size=0.8) + border() +
  xlim(quantile(c(log10(var.pois), log10(var.gene[var.gene!=0])), c(0,1))) + 
  ylim(quantile(c(log10(var.pois), log10(var.gene[var.gene!=0])), c(0,1)))

p3.t <- 
ggscatter(data = data.frame(x=log10(mean(CE*sf)*mu + var(CE*sf)*mu^2 + rowMeans(((mu^2*(1+psi)) %o% (CE*(1-CE) * sf^2)) * rho) +
                                    psi*mean(CE^2*sf^2)*mu^2), 
                            y=log10(var.gene)), x='x', y='y', 
          xlab ='Log10 expected variance', ylab='Log10 observed variance', 
          size=0.1, title = expression('Pre-dropout dispersion'~psi[i]~'='~hat(psi[i]))) +
  geom_abline(slope = 1, intercept = 0, color = 'gray', size=0.8) + border() +
  xlim(quantile(c(log10(var.pois), log10(var.gene[var.gene!=0])), c(0,1))) + 
  ylim(quantile(c(log10(var.pois), log10(var.gene[var.gene!=0])), c(0,1)))


ggarrange(p1.t, p2.t, p3.t, ncol = 3, nrow = 1)

# 2. Examples of zero-inflated gene
prop0.gene <- rowMeans(data.t==0)
plot(log10(mean.gene), prop0.gene, cex=0.3, pch=16)
pi0 <- 1/(1+exp(-out.zinb.t$par.noDE[,1]))

tau0 <- out.em.zinb.t$tau0
tau1 <- out.em.zinb.t$tau1
rho.zinb <- rep(1, nrow(data.t)) %o% tau0 + (1-pi0)*out.zinb.t$par.noDE[,2] %o% (tau1*out.em.zinb.t$est.sf)
rho.zinb <- 1/(1+exp(-rho.zinb))

dbetabinom2 <- function(x,prob,size,rho) {
  a <- prob*(1-rho)/rho ; b <- (1-prob)*(1-rho)/rho
  corr.size <- ifelse(size<x,x,size)
  logl <- lgamma(size+1) - lgamma(x+1) - lgamma(corr.size - x + 1) + lbeta(x+a,corr.size-x+b)-lbeta(a,b)
  return(exp(logl)*(x <= size))
}

decentGOF <- function(g, data.obs, par1, par2, sf1, sf2, CE1, CE2, rho1, rho2, plot=T) {
  zmax <- max(data.obs[g, ])+5
  ncell <- ncol(data.obs)
  set.seed(7)
  y1 <- sapply(1:ncell, function(j) ZIM::rzinb(1, exp(-par1[g, 3]),
                                               exp(par1[g, 2]) * sf1[j],
                                               1/(1+exp(-par1[g, 1]))))
  set.seed(7)
  y2 <- sapply(1:ncell, function(j) ZIM::rzinb(1, exp(-par2[g, ncol(par2)]),
                                                 exp(par2[g, 2]) * sf2[j],
                                                 0))
  # pm.poisb <- sapply(0:zmax, function(k) sum(sapply(1:ncell, function(j) dbinom2(k, out$CE[j], y.direct[j]))))
  pm1 <- sapply(0:zmax, function(k) sum(sapply(1:ncell, function(j) dbetabinom2(k, CE1[j], y1[j], rho1[j]))))
  pm2 <- sapply(0:zmax, function(k) sum(sapply(1:ncell, function(j) dbetabinom2(k, CE2[j], y2[j], rho2[j]))))
  if(plot) {
    plot(0:zmax, table(factor(data.obs[g,], levels=0:zmax)), col=1, lwd=1.5, type = 'l')
    lines(0:zmax, pm1, col=2, lwd=1.5)
    lines(0:zmax, pm2, col=3, lwd=1.5)
  }
  return(list(table(factor(data.obs[g,], levels=0:zmax)), pm1, pm2))
}

t.zf1 <- decentGOF(4644, data.t, out.zinb.t$par.noDE, out.nb.t$par.noDE, out.em.zinb.t$est.sf, out.em.nb.t$est.sf,
                   out.em.zinb.t$CE, out.em.nb.t$CE, rho.zinb, rho, plot = F)
t.zf2 <- decentGOF(4647, data.t, out.zinb.t$par.noDE, out.nb.t$par.noDE, out.em.zinb.t$est.sf, out.em.nb.t$est.sf,
                   out.em.zinb.t$CE, out.em.nb.t$CE, rho.zinb, rho, plot = F)

# Chi-square test
bin.chisq.test <- function(e.tab, o.tab, n.param) {
  e.bin <- combineTable(0:(length(e.tab)-1), e.tab, 5)
  o.bins <- cut(as.numeric(names(o.tab)), breaks=c(-Inf, rev(e.bin[2,])))
  o.bin <- aggregate(o.tab ~ o.bins, FUN = function(x) sum(x))
  E <- rev(e.bin[3, ]); O <- o.bin[, 2]
  stat <- sum((E-O)^2/E)
  pchisq(stat, df = max(1, length(E)-3), lower.tail = F)
}

t1.pval.nb <- bin.chisq.test(t.zf1[[3]], t.zf1[[1]], 2)
t1.pval.zinb <- bin.chisq.test(t.zf1[[2]], t.zf1[[1]], 3)

t2.pval.nb <- bin.chisq.test(t.zf2[[3]], t.zf2[[1]], 2)
t2.pval.zinb <- bin.chisq.test(t.zf2[[2]], t.zf2[[1]], 3)

# Plot
distPlot <- function(obj, title) {
  ggline(data = data.frame(x = rep(0:(length(obj[[1]])-1), 3),
                           y = c(obj[[1]], obj[[3]], obj[[2]]),
                           density = c(rep("Observed", length(obj[[1]])), 
                                       rep("With NB", length(obj[[1]])), 
                                       rep("With ZINB", length(obj[[1]])))),
         x = 'x', y = 'y', ylab = 'Frequency', xlab = 'Value', point.size = 0.01,
         shape = "density", color = "density", palette = c("black", "#FC4E07", "#00AFBB"), title = title)
}
p4.t <- distPlot(t.zf1, 'ENSG00000169715')
p5.t <- distPlot(t.zf2, 'ENSG00000125144')


### Zeisel ###
# 1. Show evidence of overdispersion
var.gene <- rowVars(data.z)
mean.gene <- rowMeans(data.z)
CE <- out.em.nb.z$CE
sf <- out.em.nb.z$est.sf
mu <- exp(out.nb.z$par.noDE[,2])
psi <- exp(out.nb.z$par.noDE[,3])
tau0 <- out.em.nb.z$tau0
tau1 <- out.em.nb.z$tau1
rho <- rep(1, nrow(data.z)) %o% tau0 + log(mu) %o% (tau1*sf)
rho <- 1/(1+exp(-rho))

# mean-varianCE relationship
the.conc <- 10^(seq(log10(min(mu))-0.2, log10(max(mu))+0.2, len=10000))
rho0 <- rep(1, length(the.conc)) %o% tau0 + log(the.conc) %o% (tau1*sf)
rho0 <- 1/(1+exp(-rho0))
p1.z <- 
  ggscatter(data = data.frame(var=log10(var.gene), mean=log10(mean.gene)), x='mean', y='var', 
            ylab ='Log10 observed gene-wise variance', xlab='Log10 observed gene mean count', size=0.1) + 
  geom_line(data = data.frame(y=log10(mean(CE*sf) * the.conc + var(CE*sf) * the.conc^2 + 
                                        rowMeans((the.conc^2 %o% (CE*(1-CE) * sf^2)) * rho0)),
                              x=log10(mean(CE*sf) * the.conc)), aes(x, y), col='red', size=1)
# observed vs expected
var.nb <- mean(CE*sf)*mu + var(CE*sf)*mu^2 + rowMeans(((mu^2*(1+psi)) %o% (CE*(1-CE) * sf^2)) * rho) + psi*mean(CE^2*sf^2)*mu^2
var.pois <- mean(CE*sf)*mu + var(CE*sf)*mu^2 + rowMeans((mu^2 %o% (CE*(1-CE) * sf^2)) * rho)
p2.z <- 
  ggscatter(data = data.frame(x=log10(var.pois), 
                              y=log10(var.gene)), x='x', y='y', 
            xlab ='Log10 expected variance', ylab='Log10 observed variance', 
            size=0.1, title = expression('Pre-dropout dispersion'~psi[i]~'='~0)) +
  geom_abline(slope = 1, intercept = 0, color = 'gray', size=0.8) + border() +
  xlim(quantile(c(log10(var.pois), log10(var.gene[var.gene!=0])), c(0,1))) + 
  ylim(quantile(c(log10(var.pois), log10(var.gene[var.gene!=0])), c(0,1)))

p3.z <- 
  ggscatter(data = data.frame(x=log10(var.nb), 
                              y=log10(var.gene)), x='x', y='y', 
            xlab ='Log10 expected variance', ylab='Log10 observed variance', 
            size=0.1, title = expression('Pre-dropout dispersion'~psi[i]~'='~hat(psi[i]))) +
  geom_abline(slope = 1, intercept = 0, color = 'gray', size=0.8) + border() +
  xlim(quantile(c(log10(var.pois), log10(var.gene[var.gene!=0])), c(0,1))) + 
  ylim(quantile(c(log10(var.pois), log10(var.gene[var.gene!=0])), c(0,1)))


ggarrange(p1.z, p2.z, p3.z, ncol = 3, nrow = 1)

# 2. Examples of zero-inflated gene
prop0.gene <- rowMeans(data.z==0)
plot(log10(mean.gene), prop0.gene, cex=0.3, pch=16)
pi0 <- 1/(1+exp(-out.zinb.z$par.noDE[,1]))

tau0 <- out.em.zinb.z$tau0
tau1 <- out.em.zinb.z$tau1
rho.zinb <- rep(1, nrow(data.z)) %o% tau0 + (1-pi0)*out.zinb.z$par.noDE[,2] %o% (tau1*out.em.zinb.z$est.sf)
rho.zinb <- 1/(1+exp(-rho.zinb))

z.zf1 <- decentGOF(138, data.z, out.zinb.z$par.noDE, out.nb.z$par.noDE, out.em.zinb.z$est.sf, out.em.nb.z$est.sf,
                   out.em.zinb.z$CE, out.em.nb.z$CE, rho.zinb, rho, plot = F) #
z.zf2 <- decentGOF(3043, data.z, out.zinb.z$par.noDE, out.nb.z$par.noDE, out.em.zinb.z$est.sf, out.em.nb.z$est.sf,
                   out.em.zinb.z$CE, out.em.nb.z$CE, rho.zinb, rho, plot = F)


# Chi-square test
z1.pval.nb <- bin.chisq.test(z.zf1[[3]], z.zf1[[1]], 2)
z1.pval.zinb <- bin.chisq.test(z.zf1[[2]], z.zf1[[1]], 3)

z2.pval.nb <- bin.chisq.test(z.zf2[[3]], z.zf2[[1]], 2)
z2.pval.zinb <- bin.chisq.test(z.zf2[[2]], z.zf2[[1]], 3)

# plot
p4.z <- distPlot(z.zf1, 'Xist')
p5.z <- distPlot(z.zf2, 'Lpl')

### combine plots ###
ggarrange(p1.z + theme(plot.margin = margin(10, 0, 0, 15)), p2.z, p3.z, 
          p1.t + theme(plot.margin = margin(10, 0, 0, 15)), p2.t, p3.t, ncol = 3, nrow = 2,
          labels = c('a', '', '', 'b', '', ''), font.label = list(size=18))

# 9.5 7
# 950 700

ggarrange(p4.t , p5.t, 
          p4.z , p5.z, ncol = 2, nrow = 2,
          labels = c('a', '', 'b', ''), font.label = list(size=18), common.legend = T, legend = 'right')

# 8 6
# 750 550

