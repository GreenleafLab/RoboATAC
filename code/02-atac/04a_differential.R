# script modified from Sidney V.
# Import Libraries
suppressWarnings({
  library(ChrAccR)
  library(muLogR)
  library(muRtools)
  library(ggplot2)
  library(chromVAR)
  library(motifmatchr)
  library(JASPAR2018)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(dplyr)
  library(DESeq2)
  library(GenomicFeatures)
  library(ChrAccRAnnotationHg38)
  library(tidyverse)
  library(ggpubr)
  library(ggrastr)
  library(ggrepel)
  library(viridis)
  library(Gviz)
  library(TFBSTools)
  library(reshape2)
})

source("../utils/tf_funcs.R")

outdir <- "../../output/02-atac/04/"
plotdir <- "../../plots/02-atac/04/"
dir.create(plotdir, recursive=T, showWarnings=F)
dir.create(outdir, recursive=T, showWarnings=F)

#########################
######## HEK only #######
#########################

# # redo differential analysis by {plate,tf,dose} ---------------------
dsa_path <- '../../output/02-atac/01/dsa_filt/'
dsa <- loadDsAcc(dsa_path)
dsa@sampleAnnot$Plate_TF_Dose <- paste0("P", dsa@sampleAnnot$PlateNum, "_", dsa@sampleAnnot$TF_Dose)

setConfigElement("regionTypes",c("filteredConsensus"))
setConfigElement("differentialColumns", c("Plate_TF_Dose"))
if (!dir.exists(paste0(outdir, "/ChrAccR_reports_Plate_TF_Dose/"))){
  createReport_differential(dsa, paste0(outdir, "/ChrAccR_reports_Plate_TF_Dose/"))
}

## read motif matching to JASPAR ------------------------------------------------
region <- "filteredConsensus"
motifDb <- "../../output/02-atac/03/JASPAR2020_Motifs.rds"
motif_score_file <- "../../output/02-atac/03/motif.scores.JASPAR2020.rds"

# TF Motif Enrichment on differential peaks by {tf} -------------------------------
dds_file <- "../../output/02-atac/01/ChrAccR_reports/differential_data/diffObj_filteredconsensus.rds"
group1 <- "GFP" # reference level
contrast.col <- "TF"
p.cutoff <- 0.05
log2fc.cutoff <- 0.58
  
all_groups <- unique(dsa@sampleAnnot[contrast.col]) %>% unlist
all_groups <- all_groups[all_groups != group1]
for (group2 in all_groups){
  message(paste0(group2, " vs ", group1))
  command <- sprintf("sbatch -p wjg,sfgf,biochem --time=24:00:00 --mem-per-cpu=150g --out=slurm-logs/differential/slurm-%%j-%s --job-name=diff-%s --wrap \"Rscript 04b_job_differential.R %s %s %s %s %s %s %s %s %s %s %s %s %s\"",
                      paste0(group2, "_vs_", group1), paste0(group2, "_vs_", group1),
                      dsa_path, region, motifDb, motif_score_file, dds_file,
                      group1, group2, contrast.col, p.cutoff, log2fc.cutoff,
                      paste0(outdir, "/deseq_data_", contrast.col), paste0(plotdir, "/de_", contrast.col), paste0(outdir, "/motif_data_", contrast.col))
  print(command)
  system(command)
}

# TF Motif Enrichment on differential peaks by {plate,tf,dose} ------------------------------
dds_file <- paste0(outdir, "/ChrAccR_reports_Plate_TF_Dose/differential_data/diffObj_filteredconsensus.rds")

contrast.col <- "Plate_TF_Dose"
p.cutoff <- 0.05
log2fc.cutoff <- 0.58

# plate 1
group1 <- "P1_GFP_1" # reference level
all_groups <- unique(dsa@sampleAnnot$Plate_TF_Dose) 
all_groups <- all_groups[grepl("^P1", all_groups)] 
all_groups <- all_groups[all_groups != group1]

for (group2 in all_groups){
  message(paste0(group2, " vs ", group1))
  command <- sprintf("sbatch -p wjg,sfgf,biochem --time=4:00:00 --mem-per-cpu=64g --out=slurm-logs/differential/slurm-%%j-%s --job-name=diff-%s --wrap \"Rscript 04b_job_differential.R %s %s %s %s %s %s %s %s %s %s %s %s %s\"",
                      paste0(group2, "_vs_", group1), paste0(group2, "_vs_", group1),
                      dsa_path, region, motifDb, motif_score_file, dds_file,
                      group1, group2, contrast.col, p.cutoff, log2fc.cutoff,
                      paste0(outdir, "/deseq_data_", contrast.col), paste0(plotdir, "/de_", contrast.col), paste0(outdir, "/motif_data_", contrast.col))
  print(command)
  system(command)
}


# plate 2
group1 <- "P2_GFP_1" # reference level
all_groups <- unique(dsa@sampleAnnot$Plate_TF_Dose) 
all_groups <- all_groups[grepl("^P2", all_groups)] 
all_groups <- all_groups[all_groups != group1]

for (group2 in all_groups){
  message(paste0(group2, " vs ", group1))
  command <- sprintf("sbatch -p wjg,sfgf,biochem --time=4:00:00 --mem-per-cpu=64g --out=slurm-logs/differential/slurm-%%j-%s --job-name=diff-%s --wrap \"Rscript 04b_job_differential.R %s %s %s %s %s %s %s %s %s %s %s %s %s\"",
                      paste0(group2, "_vs_", group1), paste0(group2, "_vs_", group1),
                      dsa_path, region, motifDb, motif_score_file, dds_file,
                      group1, group2, contrast.col, p.cutoff, log2fc.cutoff,
                      paste0(outdir, "/deseq_data_", contrast.col), paste0(plotdir, "/de_", contrast.col), paste0(outdir, "/motif_data_", contrast.col))
  print(command)
  system(command)
}

# plate 3
# group1_ls <- c("P3_GFP_1", "P3_neg-ctl_1", "P3_tet-CTCF_1") # reference level
group1_ls <- c("P3_GFP_1") # reference level
all_groups <- unique(dsa@sampleAnnot$Plate_TF_Dose) 
all_groups <- all_groups[grepl("^P3", all_groups)] 

for (group1 in group1_ls) {
  all_groups_sub <- all_groups[all_groups != group1]
  for (group2 in all_groups_sub){
    message(paste0(group2, " vs ", group1))
    command <- sprintf("sbatch -p wjg,sfgf,biochem --time=4:00:00 --mem-per-cpu=64g --out=slurm-logs/differential/slurm-%%j-%s --job-name=diff-%s --wrap \"Rscript 04b_job_differential.R %s %s %s %s %s %s %s %s %s %s %s %s %s\"",
                      paste0(group2, "_vs_", group1), paste0(group2, "_vs_", group1),
                      dsa_path, region, motifDb, motif_score_file, dds_file,
                      group1, group2, contrast.col, p.cutoff, log2fc.cutoff,
                      paste0(outdir, "/deseq_data_", contrast.col), paste0(plotdir, "/de_", contrast.col), paste0(outdir, "/motif_data_", contrast.col))
  print(command)
  system(command)
  }
}
