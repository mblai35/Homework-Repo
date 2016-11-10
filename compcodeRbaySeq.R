# compcodeR.R
# R version 3.2.2 (2015-08-14)
# November 6, 2016. Mallory B. Lai.
# Reviewed by: TODO (Mallory B. Lai) : Find reviewer to proofread
# Practice using compcodeR package for simulating differential gene expression
# data.
#------------------------------------------------------------------------------
source("https://bioconductor.org/biocLite.R")
biocLite("compcodeR")
library(compcodeR)
#------------------------------------------------------------------------------

# Generate synthetic data.
sim <- generateSyntheticData(dataset = "D1", n.vars = 12500,
                                   samples.per.cond = 5, n.diffexp = 1250,
                                   repl.id = 1, seqdepth = 1e7,
                                   fraction.upregulated = 0.5,
                                   between.group.diffdisp = FALSE,
                                   filter.threshold.total = 1,
                                   filter.threshold.mediancpm = 0,
                                   fraction.non.overdispersed = 0,
                                   output.file = "simD1.rds")

# Summarize synthetic data.
summarizeSyntheticDataSet(sim, output.filename = "simD1.html")

# Differential expression analysis
runDiffExp(data.file = "simD1.rds",
           result.extent = "baySeq", Rmdfunction = "baySeq.createRmd",
           output.directory = ".",
           norm.method = "edgeR",
           equaldisp = TRUE, distr.choice = "NB")

