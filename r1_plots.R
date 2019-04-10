# Visualize filtering cut-off
# read Soumillon and Chen data as in corresponding scripts

s1.cutoff <- gghistogram(data = data.frame(x = log10(rowSums(data.obs))), x = "x", bins = 50, xlab = 'Log10 sum expression') + 
  geom_vline(xintercept = median(log10(rowSums(data.obs))) - mad(log10(rowSums(data.obs))), col = 'red', linetype = 'dashed')

s2.cutoff <- gghistogram(data = data.frame(x = log10(rowSums(data.obs.3g))), x = "x", bins = 50, xlab = 'Log10 sum expression') + 
  geom_vline(xintercept = median(log10(rowSums(data.obs.3g))) - mad(log10(rowSums(data.obs.3g))), col = 'red', linetype = 'dashed')

c.cutoff <- gghistogram(data = data.frame(x = log10(rowSums(sc))), x = "x", bins = 50, xlab = 'Log10 sum expression') + 
  geom_vline(xintercept = log10(100), col = 'red', linetype = 'dashed')

ggarrange(s1.cutoff, s2.cutoff, c.cutoff, ncol = 3, labels = 'auto', font.label = list(size=18))

# 8 3
# 800 300

# distribution of cell-specific tau in Tung data
em.t <- readRDS('tung/sp/decent.noDE.rds')
tau0 <- gghistogram(data = data.frame(x = em.t$tau0), x = "x", bins = 50, xlab = 'tau0 estimates', fill = "#00AFBB", 
                    add = 'median', add.params = list(size = 0.8, linetype = "dashed"))
tau1 <- gghistogram(data = data.frame(x = em.t$tau1), x = "x", bins = 50, xlab = 'tau1 estimates', fill = "#E7B800", 
                    add = 'median', add.params = list(size = 0.8, linetype = "dashed"))
ggarrange(tau0, tau1, ncol = 2, labels = 'auto', font.label = list(size=18))
# 6 3
# 600 300

# Distribution of rho in the 4 datasets
em.sm <- readRDS('soumillon/tmm/decent.noDE.rds')
em.sv <- readRDS('savas/decent.noDE.rds')
em.c.2x <- readRDS('chen/decent.noDE.rds')
  
t.rho <- gghistogram(data = data.frame(x = c(em.t$tau0 + em.t$tau1 %o% log((1-logit_inv(lrt.t$par.noDE[, 1]))*exp(lrt.t$par.noDE[,2])))), 
                     x = "x", bins = 30, xlab = 'Logit rho', fill = "lightgray", add = 'median', title = 'Tung', y = '..density..')
sm.rho <- gghistogram(data = data.frame(x = em.sm$tau0 + log((1-logit_inv(lrt.sm$par.noDE[, 1]))*exp(lrt.sm$par.noDE[,2])) * em.sm$tau1), 
                     x = "x", bins = 30, xlab = 'Logit rho', fill = "lightgray", add = 'median', title = 'Soumillon', y = '..density..')
sv.rho <- gghistogram(data = data.frame(x = em.sv$tau0 + log((1-logit_inv(lrt.sv$par.noDE[, 1]))*exp(lrt.sv$par.noDE[,2])) * em.sv$tau1), 
                     x = "x", bins = 30, xlab = 'Logit rho', fill = "lightgray", add = 'median', title = 'Savas', y = '..density..')
c.rho <- gghistogram(data = data.frame(x = em.c.2x$tau0 + log((1-logit_inv(lrt.c.2x$par.noDE[, 1]))*exp(lrt.c.2x$par.noDE[,2])) * em.c.2x$tau1), 
                     x = "x", bins = 30, xlab = 'Logit rho', fill = "lightgray", add = 'median', title = 'Chen', y = '..density..')
ggarrange(t.rho, sm.rho, sv.rho, c.rho, ncol = 2, nrow = 2, labels = 'auto', font.label = list(size=18))
# 6 5
# 600 500

# How psi changes between models with BB or B capture models
lrt.sm <- readRDS('soumillon/tmm/decent.lrt.rds')
lrt.sm.b <- readRDS('soumillon/tmm_b/decent.lrt.rds')

ggscatter(data = data.frame(x = exp(lrt.sm$par.DE[, 4]), y = exp(lrt.sm.b$par.DE[, 4])), x = 'x', y = 'y', 
          xlab = 'psi estimates with BB capture model', ylab = 'psi estimates with B capture model',
          xlim = range(exp(lrt.sm.b$par.DE[, 4])), size = 0.1) +
  geom_abline(slope = 1, intercept = 0, col = 'red')
# 3.5 3.5
# 350 350

# tau with sp vs endo
em.t.tau <- readRDS('tung/sp_tau/decent.noDE.rds')

s.e.tau1 <- ggscatter(data = data.frame(x = em.t$tau1[em.t.tau$tau1 > -5], y = em.t.tau$tau1[em.t.tau$tau1 > -5]), x = 'x', y = 'y', 
          xlab = 'tau1 estimated using endogenous genes', ylab = 'tau1 estimated using spike-ins',
          size = 0.3,
          add = "reg.line",  # Add regressin line
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
          cor.coeff.args = list(method = "pearson", label.sep = "\n", label.x.npc = "left", label.y.npc = "top"))
s.e.tau0 <- ggscatter(data = data.frame(x = em.t$tau0, y = em.t.tau$tau0), x = 'x', y = 'y', 
          xlab = 'tau0 estimated using endogenous genes', ylab = 'tau0 estimated using spike-ins',
          size = 0.3,
          add = "reg.line",  # Add regressin line
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
          cor.coeff.args = list(method = "pearson", label.sep = "\n", label.x.npc = "left", label.y.npc = "top"))

ggarrange(s.e.tau0, s.e.tau1, ncol = 2, labels = 'auto', font.label = list(size=18))

