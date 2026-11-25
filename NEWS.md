# scdRs 0.99.0

## NEW FEATURES

* First version, submitted to Bioconductor: a native R port of scDRS
  (single-cell disease-relevance score; Zhang, Hou et al., Nature
  Genetics 2022, Python scdrs 1.0.4) for `SingleCellExperiment` objects.
  It does not call Python.
* Disease gene sets: `scdrsReadGeneSets()`, `scdrsWriteGeneSets()`,
  `scdrsMungeGeneSets()` (from gene-level GWAS p-values or z-scores) and
  `scdrsGetHomologs()` (mouse-human gene mapping).
* Data preparation: `scdrsNormalize()` (scDRS's filtering and
  normalization of raw counts), `scdrsPreprocess()` and
  `scdrsComputeGeneStats()` (covariate correction and gene statistics).
* Cell scoring: `scdrsScoreCells()`, `scdrsSelectControlGenes()` and
  `scdrsComputeScore()` (several gene sets); score files compatible with
  the Python tool via `scdrsWriteScores()` and `scdrsReadScores()`.
* Downstream analyses: `scdrsTestGroups()` (group association and
  heterogeneity), `scdrsTestCorrelation()`, `scdrsRankGenes()`,
  `scdrsPerformDownstream()`, `scdrsTestGearysC()`,
  `scdrsComputeGearysC()` and `scdrsBuildGraph()` (cell-cell graph as
  built by scanpy).
* Utilities: `scdrsZToP()`, `scdrsPToZ()`, `scdrsMetaAnalyze()` and
  `scdrsTestOverlap()`.
* Every function is tested against the output of the Python scDRS.
  Results agree to float32 precision when the same control gene sets are
  used.
