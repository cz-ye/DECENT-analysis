require(pROC)
source('func_de_methods.R')

mnc <- mst <- er <- dct <- list()
mnc.auc <- mst.auc <- er.auc <- dct.auc <- list()


for (rep in 1:20) {
  load(paste0('simdata_ZINB_BB_Tung6_OD_kb_010119/simdata_ZINB_BB_Tung6_OD_kb_010119_set', rep, '.RData'))
  f <- rowSums(data.obs > 0) >= 3
  data.obs <- data.obs[f,]
  label <- as.numeric(DE.gene[f])
  
  mnc[[rep]] <- runMonocleDE(data.obs, ctype)
  mnc.auc[[rep]] <- auc(label, mnc[[rep]]$pval, partial.auc=c(0.9,1), partial.auc.correct=T)
  
  mst[[rep]] <- runMASTDE(data.obs, ctype)
  mst.auc[[rep]] <- auc(label, mst[[rep]]$`Pr(>Chisq)`[match(rownames(data.obs), mst[[rep]]$primerid)], partial.auc=c(0.9,1), partial.auc.correct=T)
  
  er[[rep]] <- runEdgeR(data.obs, ctype)
  er.auc[[rep]] <- auc(label, er[[rep]]$PValue, partial.auc=c(0.9,1), partial.auc.correct=T)
  
  sp.obs <- ercc.obs[rowSums(ercc.obs) != 0, ]
  sp.nom <- ercc.true[rowSums(ercc.obs) != 0]
  dct[[rep]] <- decent.table <- decent(data.obs, ~as.factor(ctype), use.spikes = T, spikes = sp.obs, spike.conc = sp.nom, 
                                       tau.est = 'spikes', dir='temp/')
  dct.auc[[rep]] <- auc(label, dct[[rep]]$pval, partial.auc=c(0.9,1), partial.auc.correct=T)
  
  print(paste('rep', rep))
}
