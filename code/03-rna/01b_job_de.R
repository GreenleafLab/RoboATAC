# Differential analysis of RNAseq counts estimated by kallisto

# args ---------------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

dds_path <- args[1]
id2name_path <- args[2]
cond1 <- args[3] # test level
cond2 <- args[4] # reference level
contrast.col <- args[5]
log2fc.cutoff <- as.numeric(args[6])
padj.cutoff <- as.numeric(args[7])
outdir <- args[8]
plotdir <- args[9]
labels <- eval(parse(text=args[10]))
lfcshrink <- as.logical(args[11])
top.label <- as.integer(args[12])

# main ---------------------------------------------------------------------------------
source("../utils/tf_funcs.R")
print(labels)
suppressWarnings({
  library(DESeq2)
  library(RColorBrewer)
  library(ggpubr)
  library(EnhancedVolcano)
  library(ggrastr)
  library(tidyverse)
  library(data.table)
})


dir.create(outdir, recursive=T, showWarnings=F)
dir.create(plotdir, recursive=T, showWarnings=F)

id2name <- read.table(id2name_path, sep="\t")
dds <- readRDS(dds_path)

p <- plotDDS(dds, cond1, cond2, id2name, contrast.col, log2fc.cutoff=log2fc.cutoff, padj.cutoff=padj.cutoff, 
              labels=labels, outdir=outdir, top.label=top.label, lfcshrink=lfcshrink)

pdf(paste0(plotdir, "/MA_Volcano_plot_",cond1, "_vs_", cond2,"_log2fc",log2fc.cutoff, "_p", padj.cutoff, ".pdf"))
print(p[[1]])
print(p[[2]])
dev.off()
