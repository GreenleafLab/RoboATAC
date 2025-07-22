# Import Libraries
suppressWarnings({
  library(ChrAccR)
  library(muLogR)
  library(muRtools)
  library(GenomicFeatures)
  library(tidyverse)
  library(viridis)
  library(Gviz)
  library(TFBSTools)
  library(reshape2)
})

source("../utils/tf_funcs.R")

outdir <- "../../output/02-atac/02/"
plotdir <- "../../plots/02-atac/02/"
dir.create(plotdir, recursive=T, showWarnings=F)
dir.create(outdir, recursive=T, showWarnings=F)

# color palette ---------------------------------------------------------
# 25 colors
colpal <- read_tsv("../utils/colpal_TF.tsv") %>% 
          column_to_rownames(var="TF")

# functions ---------------------------------------------------------------
plot_pca_all <- function(sample_annot, cm, suffix=""){
  annot <- sample_annot
  annot$PlateNum <- as.character(annot$PlateNum)
  annot$TF <- annot$TF %>% as.factor
  BL_stallion <- colpal[levels(annot$TF), "color"]
  
  f <- paste0(outdir, "/pcaCoord", suffix, ".rds")
  if (file.exists(f)){
    pcaCoord <- readRDS(f)
  } else{
    pcaCoord <- muRtools::getDimRedCoords.pca(t(cm))
    saveRDS(pcaCoord, f)
  }
  # 
  
  pdf(paste0(plotdir, "/PCA", suffix, ".pdf"), width=9, height=7)
  print(muRtools::getDimRedPlot(pcaCoord, annot=annot, colorCol="TF", shapeCol="Plasmid_Dose", colScheme=BL_stallion, ptSize=5) + 
                  theme_classic() + theme(text=element_text(size=20),panel.border = element_rect(color ="black",fill = NA, size = 1)))
  print(muRtools::getDimRedPlot(pcaCoord, annot=annot, colorCol="CellType", ptSize=5) + 
                  theme_classic() + theme(text=element_text(size=20),panel.border = element_rect(color ="black",fill = NA, size = 1)))
  print(muRtools::getDimRedPlot(pcaCoord, annot=annot, colorCol="PlateNum", ptSize=5) + 
                  theme_classic() + theme(text=element_text(size=20),panel.border = element_rect(color ="black",fill = NA, size = 1)))
  print(getDimRedPlot_BL(pcaCoord, annot=annot, colorCol="TF", colScheme=BL_stallion, ptSize=5) + 
          theme_classic() + theme(text=element_text(size=20),panel.border = element_rect(color ="black",fill = NA, size = 1)))
  print(getDimRedPlot_BL(pcaCoord, annot=annot, colorCol="TF", alphaCol="Plasmid_Dose", colScheme=BL_stallion, ptSize=5) + 
          theme_classic() + theme(text=element_text(size=20),panel.border = element_rect(color ="black",fill = NA, size = 1)))
  print(getDimRedPlot_BL(pcaCoord, annot=annot, colorCol="TF", sizeCol="Plasmid_Dose", colScheme=BL_stallion) + 
          theme_classic() + theme(text=element_text(size=20),panel.border = element_rect(color ="black",fill = NA, size = 1)))
  dev.off()
}

plot_pca_celltype <- function(sample_annot, cm, celltype, suffix=""){
  cols <- sample_annot %>% dplyr::filter(CellType==celltype) %>% dplyr::select(sampleName) %>% unlist %>% unname
  
  annot <- sample_annot[cols,]
  annot$PlateNum <- as.character(annot$PlateNum)
  annot$TF <- annot$TF %>% as.factor
  BL_stallion <- colpal[levels(annot$TF), "color"]
  
  f <- paste0(outdir, "/pcaCoord_", celltype, suffix, ".rds")
  if (file.exists(f)){
    pcaCoord <- readRDS(f)
  } else{
    subcm <- cm[,cols]
    pcaCoord <- muRtools::getDimRedCoords.pca(t(subcm))
    saveRDS(pcaCoord, f)
  }
  
  pdf(paste0(plotdir, "/PCA_", celltype, suffix, ".pdf"), width=9, height=7)
  print(muRtools::getDimRedPlot(pcaCoord, annot=annot, colorCol="TF", shapeCol="Plasmid_Dose", colScheme=BL_stallion, ptSize=5) + 
                  theme_classic() + theme(text=element_text(size=20),panel.border = element_rect(color ="black",fill = NA, size = 1)))
  print(muRtools::getDimRedPlot(pcaCoord, annot=annot, colorCol="PlateNum", ptSize=5) + 
                  theme_classic() + theme(text=element_text(size=20), panel.border = element_rect(color ="black",fill = NA, size = 1)))
  dev.off()
}


plot_pca_plate <- function(sample_annot, cm, plate, suffix=""){
  cols <- sample_annot %>% dplyr::filter(PlateNum==plate) %>% dplyr::select(sampleName) %>% unlist %>% unname
  
  annot <- sample_annot[cols,]
  annot$PlateNum <- as.character(annot$PlateNum)
  annot$TF <- annot$TF %>% as.factor
  BL_stallion <- colpal[levels(annot$TF), "color"]

  f <- paste0(outdir, "/pcaCoord_plate", plate, suffix, ".rds")
  if (file.exists(f)){
    pcaCoord <- readRDS(f)
  } else{
    subcm <- cm[,cols]
    pcaCoord <- muRtools::getDimRedCoords.pca(t(subcm))
    saveRDS(pcaCoord, f)
  }

  pdf(paste0(plotdir, "/PCA_plate", plate, suffix, ".pdf"), width=9, height=7)
  print(muRtools::getDimRedPlot(pcaCoord, annot=annot, colorCol="TF", shapeCol="Plasmid_Dose", colScheme=BL_stallion, ptSize=5) + 
                  theme_classic() + theme(text=element_text(size=20),panel.border = element_rect(color ="black",fill = NA, size = 1)))
  print(muRtools::getDimRedPlot(pcaCoord, annot=annot, colorCol="CellType", ptSize=5) + 
                  theme_classic() + theme(text=element_text(size=20),panel.border = element_rect(color ="black",fill = NA, size = 1)))
  dev.off()
}

# PCA (rpkm log2 quantile norm) ---------------------------------------------------------------------
dsa_norm <- loadDsAcc('../../output/02-atac/01/dsa_norm_rpkmlog2quantile/')
cm <- as.data.frame(dsa_norm@counts$filteredConsensus)

# plot all
plot_pca_all(getSampleAnnot(dsa_norm), cm, suffix="_all_rpkmlog2quantile")

# plot by plate
plot_pca_plate(getSampleAnnot(dsa_norm), cm, 1, suffix="_rpkmlog2quantile")
plot_pca_plate(getSampleAnnot(dsa_norm), cm, 2, suffix="_rpkmlog2quantile")
plot_pca_plate(getSampleAnnot(dsa_norm), cm, 3, suffix="_rpkmlog2quantile")


# PCA (quantile norm) ---------------------------------------------------------------------
dsa_norm <- loadDsAcc('../../output/02-atac/01/dsa_norm_quantile/')
cm <- as.data.frame(dsa_norm@counts$filteredConsensus)

# plot all
plot_pca_all(getSampleAnnot(dsa_norm), cm, suffix="_all")

# plot by plate
plot_pca_plate(getSampleAnnot(dsa_norm), cm, 1)
plot_pca_plate(getSampleAnnot(dsa_norm), cm, 2)
plot_pca_plate(getSampleAnnot(dsa_norm), cm, 3)



# PCA (vst norm) ---------------------------------------------------------------------
dsa_norm <- loadDsAcc('../../output/02-atac/01/dsa_norm_vst/')
cm <- as.data.frame(dsa_norm@counts$filteredConsensus)

# plot all
plot_pca_all(getSampleAnnot(dsa_norm), cm, suffix="_all_vst")

# plot by plate
plot_pca_plate(getSampleAnnot(dsa_norm), cm, 1, suffix="_vst")
plot_pca_plate(getSampleAnnot(dsa_norm), cm, 2, suffix="_vst")
plot_pca_plate(getSampleAnnot(dsa_norm), cm, 3, suffix="_vst")

