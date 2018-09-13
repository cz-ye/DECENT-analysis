## R scripts for the analyses in the DECENT manuscript

### Scripts:
`exam_dropout_model.Rmd`
* Examine dropout model assumptions using ERCC spike-in data. 
* Fig.1, sFig.1, sFig.2

`overdisp_zeroinfl_pre_dropout.R`
* Looking for overdispersion and zero-inflation in the pre-dropout count using scRNA-seq data based on the DECENT model. 
* sFig.3, sFig.4

`test_smfish.R`
* Testing for overdispersion and zero-inflation in the smFISH data.
* sFig.5

`misspecify_eta.R`
* Plotting probability mass function of the observed data given various levels of misspecified eta. 
* sFig.6

`sim.R`
* Simulation studies. 
* Fig.2, Fig.3, sFig.7 

`tung.R`
* Analyses of Tung et al. dataset.
* sFig.8

`soumillon.R`
* Analyses of Soumillon et al. dataset. 

`savas.R`
* Analyses of Savas et al. dataset. 

`chen.R`
* Analyses of Chen et al. dataset. 

`benchmark_plot.R`
* Making benchmarking plots.
* Fig.4, Fig.5, Fig.6, sFig.9, sFig.10, sFig.11

`func_de_methods.R`
* Utility functions.

### Input data:
`data.savas.rds`
* The UMI count matrix of the Savas et al. data.

`ct.savas.rds`
* A vector denoting the two cell types in the Savas et al. data.

`simdata_ZINB_BB_Tung6_OD_kb.RData`
* Simulated dataset.

`ercc_length.txt`
* Lengths of ERCC spike-ins.

`cms_095046.txt`
* Information of ERCC spike-in mix.

* Other public data can be downloaded according to the comments in the scripts.

### Saved data:
`Tung_benckmark.rda`
`Soumillon_benckmark.rda`
`Savas_benckmark.rda`
`Chen_benckmark.rda`
`fpr.plots.tung.rda`
`fpr.plots.soumillon.rda`

### DECENT output:
`sim/`
`tung/`
`tung_nb/`
`zeisel/`
`zeisel_nb/`

### TASC output:
`tasc/`