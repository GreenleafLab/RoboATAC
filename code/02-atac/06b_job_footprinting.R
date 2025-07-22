# args ---------------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

dsa_path <- args[1]
motifName <- args[2]
motifDb <- args[3]
outdir <- args[4]
plotdir <- args[5]

# main ---------------------------------------------------------------------------------
library(ChrAccR)
library(tidyverse)
source("../utils/tf_funcs.R")

dir.create(outdir, showWarnings=F, recursive=T)
dir.create(plotdir, showWarnings=F, recursive=T)

dsa <- loadDsAcc(dsa_path)
motifs <- readRDS(motifDb)
samples <- getSampleAnnot(dsa)$sampleName

fps <- getMotifFootprints(dsa, motifName, samples, motifDb=motifs) # takes a long time

motif_out <- gsub("[[:punct:]]", "_", motifName) # clean up motif name for file name
saveRDS(fps, paste0(outdir, "/fps_", motif_out, ".rds"))

plot_motifs(fps, motifName, dsa, "GFP", plotdir=plotdir)
