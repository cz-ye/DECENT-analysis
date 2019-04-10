# Wrapper functions for running other DE methods

# version 3.22.2
runEdgeR <- function (data.obs, cell.type) {
  require(edgeR)
  group <- as.factor(cell.type)
  y <- DGEList(counts = data.obs, group = group)
  y <- calcNormFactors(y)
  y <- estimateDisp(y, model.matrix(~group))
  fit <- glmFit(y,design = model.matrix(~group))
  lrt <- glmLRT(fit, coef=2:length(levels(group)))
  return(lrt$table)
}

# version 1.2.0
runZinbwaveEdgeR <- function (data.obs, cell.type, J = 1e12, n.cores = 8) {
  require(edgeR)
  require(zinbwave)
  #  BiocParallel::register(BiocParallel::SerialParam())
  BiocParallel::register(BiocParallel::SnowParam(n.cores))
  input.zw <- SummarizedExperiment(assays = list(counts = as.matrix(data.obs)))
  group <- as.factor(cell.type)
  zw.fit <- zinbwave(input.zw, X = model.matrix(~group), epsilon=J, verbose = T)
  
  y <- DGEList(counts = data.obs, group = group)
  y <- calcNormFactors(y)
  y$weights <- assay(zw.fit, "weights")
  y <- estimateDisp(y, model.matrix(~group))
  
  fit <- glmFit(y,design = model.matrix(~group))
  lrt <- glmWeightedF(fit, coef=2:length(levels(group)))
  return(lrt$table)
}

# version 2.8.0
runMonocleDE <- function (data.obs, cell.type, UMI = T, n.cores = 2) {
  require(monocle)
  pd <- new("AnnotatedDataFrame", 
            data = data.frame(CellType = as.factor(cell.type), row.names = colnames(data.obs)))
  fd <- new("AnnotatedDataFrame", 
            data = data.frame(gene_short_name = rownames(data.obs), row.names = rownames(data.obs)))
  # IF UMI data
  if (UMI){
    MNC <- newCellDataSet(data.obs,
                          phenoData = pd,
                          featureData = fd,
                          expressionFamily=negbinomial.size())
  } else {
    # IF non-UMI data, first convert to TPM outside this function
    MNC <- newCellDataSet(data.obs,
                          phenoData = pd,
                          featureData = fd,
                          lowerDetectionLimit=0.1,
                          expressionFamily=tobit(Lower=0.1))
    rpc_matrix <- relative2abs(MNC) # Census
    MNC <- newCellDataSet(rpc_matrix,
                          phenoData = pd,
                          featureData = fd,
                          lowerDetectionLimit=0.5,
                          expressionFamily=negbinomial.size())
  }
  MNC <- estimateSizeFactors(MNC)
  mnc.de <- differentialGeneTest(MNC, fullModelFormulaStr="~CellType", cores = n.cores, verbose = T)
  return(mnc.de)
}

# version 1.6.1
runMASTDE <- function (data.obs, cell.type, UMI = T) {
  require(MAST)
  group <- as.factor(cell.type)
  if(UMI){
    expr <- log2(edgeR::cpm(data.obs)+1)
  }
  else { # else put tpm in
    expr <- log(data.obs+1)
  }
  sca <- FromMatrix(exprsArray = expr)
  cdr2 <-colSums(assay(sca)>0)
  colData(sca)$cngeneson <- scale(cdr2)
  colData(sca)$condition <- group
  zlmCond <- zlm(~condition + cngeneson, sca)
  summaryCond <- summary(zlmCond, doLRT='condition2')
  summaryDt <- summaryCond$datatable
  fcHurdle <- merge(summaryDt[contrast=='condition2' & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                    summaryDt[contrast=='condition2' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients
  return(fcHurdle)
}

runMASTDE3 <- function (data.obs, cell.type) {
  require(MAST)
  group <- as.factor(cell.type)
  if(UMI){
    require(edgeR)
    expr <- log2(edgeR::cpm(data.obs)+1)
  }
  else { # else put tpm in
    expr <- log(data.obs+1)
  }
  sca <- FromMatrix(exprsArray = expr)
  cdr2 <-colSums(assay(sca)>0)
  colData(sca)$cngeneson <- scale(cdr2)
  colData(sca)$condition <- group
  zlmCond <- zlm(~condition + cngeneson, sca)
  lr_anova <- lrTest(zlmCond, CoefficientHypothesis(c("condition2", "condition3")))
  return(as.data.frame(lr_anova[, 3, ]))
}

# Stable release 1.99.2 is used rather than newer versions due to this issue:
# https://github.com/hms-dbmi/scde/issues/40
runSCDE <- function (data.obs, cell.type, n.cores = 10) {
  cd <- apply(data.obs, 2, function(x) {storage.mode(x) <- 'integer'; x})
  groups <- as.factor(cell.type)
  names(groups) <- colnames(cd)
  library(scde)
  o.ifm <- scde.error.models(counts = cd, groups = groups, n.cores = n.cores, threshold.segmentation = TRUE, save.crossfit.plots = FALSE, save.model.plots = FALSE, 
                             min.count.threshold = 1, min.nonfailed = 10, verbose = 1)
  valid.cells <- o.ifm$corr.a > 0
  o.ifm <- o.ifm[valid.cells, ]
  o.prior <- scde.expression.prior(models = o.ifm, counts = cd, length.out = 400, show.plot = FALSE)
  names(groups) <- row.names(o.ifm)
  ediff <- scde.expression.difference(o.ifm, cd, o.prior, groups, n.randomizations  =  100, n.cores = n.cores, verbose  =  1)
  return(ediff)
}

calcFPR <- function(genes, ref) {
  neg <- setdiff(genes, ref); return(sapply(1:length(genes), function(x) mean(neg %in% genes[1:x])))
}

calcTPR <- function(genes, ref) {
  return(sapply(1:length(genes), function(x) mean(ref %in% genes[1:x])))
}
calcPR <- function(genes, ref) {
  return(sapply(1:length(genes), function(x) mean(genes[1:x] %in% ref)))
}

getSpikeInMolecule <- function(m, n, d){
  #' @param m vector of # attomole per microliter for each spike-in species in the mix
  #' @param n volume of added spike-in mix (nL)
  #' @param d dilution factor
  m*n*(1/d)*1e-3*1e-18*6.02214129*1e23
}

# from NBID 
# https://bitbucket.org/Wenan/nbid/src/8767d8f11cb4142001b6bb303a0d908b31d05acb/R/combineTable.R?at=default&fileviewer=file-view-default
combineTable <- function(value, count, threshold) {
  # combine the table to make sure the count in each cell is above a threshold,
  # e.g., 5
  
  # input
  # value: the value 
  # count: the count of the observed value
  # threshold: a threshold value
  
  # output
  # combined: a 3 * n matrix with count >= threshold
  
  combined = NULL
  pointer = length(value)
  currentBucket = c(NA, Inf, 0) # (leftBound, rightBound, 0)
  while (pointer >= 1) {
    # add to the current bucket
    currentBucket[3] = currentBucket[3] + count[pointer]
    # check whether the bucket count >= threshold
    if (currentBucket[3] >= threshold || pointer == 1) {
      # set the left bound
      currentBucket[1] = value[pointer]
      # add to the combined
      combined = cbind(combined, currentBucket)
      # reset the currentBucket
      currentBucket = c(NA, value[pointer] - 1, 0)
    }
    pointer = pointer - 1
  }
  
  nBin = dim(combined)[2]
  
  # change the last bin to start from 0
  combined[1, nBin] = 0
  
  # make sure the last bin is >= threshold
  if (combined[3, nBin] < threshold) {
    # combine bin 1 and 2
    if (dim(combined)[2] > 1) {
      combined[1, nBin - 1] = combined[1, nBin]
      combined[3, nBin - 1] = combined[3, nBin] + combined[3, nBin - 1]
      combined = combined[, -nBin] # remove the last
    } else {
      print("not enough count with all combined")
    }
  }
  
  return(combined)
}
