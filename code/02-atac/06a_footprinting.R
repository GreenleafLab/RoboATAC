library(ChrAccR)
library(tidyverse)
source("../utils/tf_funcs.R")


dsa_path <- "../../output/02-atac/01/dsa_norm_rpkmlog2quantile/"
dsa <- loadDsAcc(dsa_path)

outdir <- "../../output/02-atac/06/"
plotdir <- "../../plots/02-atac/06/"
dir.create(plotdir, recursive=T, showWarnings=F)
dir.create(outdir, recursive=T, showWarnings=F)


# JASPAR 2020 motifs
motifFile <- "../../output/02-atac/03/JASPAR2020_Motifs.rds"
if (!file.exists(motifFile)){
   jaspar.motifs <- getMotifs(database=JASPAR2020::JASPAR2020) 
   saveRDS(jaspar.motifs, motifFile)
}

jaspar.motifs <- readRDS(motifFile)
jaspar.motifs


meta <- read_tsv("atac_sample_sheet.txt")
TF_ls <- meta$TF %>% unique 
TF_ls <- c(TF_ls[!grepl("\\-", TF_ls)], "SPI1", "KLF1", "CTCF", "POU5F1", "PHOX2B")
#TF_ls <- "PHOX2B"

motifNames <- grep(paste0(TF_ls, collapse="|"), names(jaspar.motifs), value=TRUE) # searching for patterns

dir.create("slurm-logs/footprinting/", recursive=T, showWarnings=F)
for (motifName in motifNames){
    motif_out <- gsub("[[:punct:]]", "_", motifName)  # clean up motif name for file name
    out <- paste0(outdir, "/JASPAR2020")
    plot <- paste0(plotdir, "/JASPAR2020")
    fps_file <- paste0(out, "/fps_", motif_out, ".rds")
    if (!file.exists(fps_file)){
        command <- sprintf("sbatch -p wjg,sfgf,biochem --time=20:00:00 --mem-per-cpu=150g --job-name=fp_%s --out=slurm-logs/footprinting/slurm-%%j-%s.out --wrap \"Rscript 06b_job_footprinting.R %s '%s' %s %s %s\"", 
                    motif_out, motif_out, dsa_path, motifName, motifFile, out, plot)
        print(command)
        system(command)
    } else{
        fps <- readRDS(fps_file)
        plot_motifs(fps, motifName, dsa, "GFP", plotdir=plot)
    }
    
}


# Vierstra motifs
url <- "https://jeffgranja.s3.amazonaws.com/ArchR/Annotations/Vierstra_Archetype_Motifs_v2.1.rds"
annoPath <- outdir
motifFile <- file.path(annoPath, basename(url))

if (!file.exists(motifFile)){
   download.file(
        url = url, 
        destfile = file.path(annoPath, basename(url)),
        quiet = FALSE
      )
}
vierstra.motifs <- readRDS(motifFile)
vierstra.motifs

motifNames <- grep("(GATA|SPI|PU1|KLF|SP|TCF|ELF|FOXP|FOXO|IRF|BACH|IKZF|LEF|NR4A1|PRDM|XBP|ALX|MXD|POU5F|OCT|SOX|MYC|CTCF)", names(vierstra.motifs), value=TRUE) # searching for patterns
dir.create("slurm-logs/footprinting/", recursive=T, showWarnings=F)
for (motifName in motifNames){
    motif_out <- gsub("[[:punct:]]", "_", motifName)  # clean up motif name for file name
    out <- paste0(outdir, "/Vierstra")
    plot <- paste0(plotdir, "/Vierstra")
    fps_file <- paste0(out, "/fps_", motif_out, ".rds")
    if (!file.exists(fps_file)){
        command <- sprintf("sbatch -p wjg,sfgf,biochem --time=20:00:00 --mem-per-cpu=150g --job-name=fp_%s --out=slurm-logs/footprinting/slurm-%%j-%s.out --wrap \"Rscript 06b_job_footprinting.R %s '%s' %s %s %s\"", 
                    motif_out, motif_out, dsa_path, motifName, motifFile, out, plot)
        print(command)
        system(command)
    } else{
        fps <- readRDS(fps_file)
        plot_motifs(fps, motifName, dsa, "GFP", plotdir=plot)
    }

}