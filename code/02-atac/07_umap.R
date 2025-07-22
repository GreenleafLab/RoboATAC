library(ChrAccR)
library(umap)
library(tidyverse)
library(muRtools)

outdir <- "../../output/02-atac/07/"
plotdir <- "../../plots/02-atac/07/"
dir.create(plotdir, recursive=T, showWarnings=F)
dir.create(outdir, recursive=T, showWarnings=F)

# color palette ---------------------------------------------------------
# 25 colors
colpal <- read_tsv("../utils/colpal_TF.tsv") %>% 
  column_to_rownames(var="TF")


# functions ---------------------------------------------------------------
plot_umap_coord <- function (coords, annot = NULL, colorCol = NULL, shapeCol = NULL, 
                             colScheme = "[auto]", ptSize = 3, addLabels = FALSE, addDensity = FALSE, 
                             addVoronoi = FALSE, annot.text = NULL, orderCol = NULL, 
                             facetCols = NULL) {
  # adapted from muRtools::getDimRedPlot
  #' @param coords a umap embedding object from umap::umap
  
  if (!is.null(annot)) {
    if (nrow(annot) != nrow(coords$layout)) {
      stop("Non-matching number of rows for dimension reduction coordinates and annotation")
    }
    if (ncol(annot) > 0) {
      if (is.null(colorCol)) {
        colorCol <- 1L
      }
    }
    else {
      stop("annot must have at least one column")
    }
  }
  xLab <- "UMAP1"
  yLab <- "UMAP2"
  df2p <- data.frame(coords$layout)
  colnames(df2p) <- c(xLab, yLab)
  colorNumeric <- FALSE
  if (!is.null(colorCol)) {
    if (length(colorCol) != 1) {
      stop("colorCol must be of length 1")
    }
    if (is.numeric(colorCol)) {
      if (!is.null(colnames(annot))) {
        colorCol <- colnames(annot)[colorCol]
      }
    }
    if (is.character(colorCol)) {
      df2p[, colorCol] <- annot[, colorCol]
    }
    else if (is.numeric(colorCol)) {
      df2p[, "color"] <- annot[, colorCol]
      colorCol <- "color"
    }
    else {
      stop("invalid value for colorCol")
    }
    if (is.numeric(df2p[, colorCol])) {
      colorNumeric <- TRUE
    }
  }
  if (!is.null(shapeCol)) {
    if (length(shapeCol) != 1) {
      stop("shapeCol must be of length 1")
    }
    if (is.numeric(shapeCol)) {
      if (!is.null(colnames(annot))) {
        shapeCol <- colnames(annot)[shapeCol]
      }
    }
    if (is.character(shapeCol)) {
      df2p[, shapeCol] <- annot[, shapeCol]
    }
    else if (is.numeric(shapeCol)) {
      df2p[, "shape"] <- annot[, shapeCol]
      shapeCol <- "shape"
    }
    else if (is.logical(shapeCol) && !shapeCol) {
      shapeCol <- NULL
    }
    else {
      stop("invalid value for shapeCol")
    }
    if (is.numeric(df2p[, shapeCol])) {
      warning("Currently only non-numeric columns are supported for dimRed shapes. --> converting to factor")
      df2p[, shapeCol] <- factor(df2p[, shapeCol])
    }
  }
  if (!is.null(orderCol)) {
    if (length(orderCol) != 1) {
      stop("colorCol must be of length 1")
    }
    if (is.numeric(orderCol)) {
      if (!is.null(colnames(annot))) {
        orderCol <- colnames(annot)[orderCol]
      }
    }
    if (is.character(orderCol) || is.numeric(orderCol)) {
      oo <- order(df2p[, orderCol], decreasing = FALSE, 
                  na.last = FALSE)
    }
    else {
      stop("invalid value for orderCol")
    }
    df2p <- df2p[oo, ]
  }
  if (!is.null(facetCols)) {
    if (length(facetCols) > 2) {
      stop("Too many facet columns")
    }
    for (cn in facetCols) {
      df2p[, cn] <- annot[, cn]
    }
  }
  if (!is.null(rownames(coords))) 
    df2p$observation <- rownames(coords)
  pp <- ggplot(df2p, aes_string(x = xLab, y = yLab, color = colorCol))
  if (!is.null(colScheme)) {
    if (is.character(colScheme) && length(colScheme) == 
        1 && colScheme == "[auto]") {
      pp <- pp + ggAutoColorScale(df2p[, colorCol])
    }
    else {
      if (colorNumeric) {
        pp <- pp + scale_color_gradientn(colours = colScheme, 
                                         na.value = "#C0C0C0")
      }
      else {
        pp <- pp + scale_color_manual(na.value = "#C0C0C0", 
                                      values = colScheme)
      }
    }
  }
  if (addDensity) {
    if (colorNumeric) {
      warning("Currently only non-numeric columns are supported for dimRed coloring in combination with density contours. --> skipping density contours")
    }
    else {
      grpCounts <- table(df2p[, colorCol])
      grps.1sample <- names(grpCounts)[grpCounts < 2]
      df2p2 <- df2p[!(df2p[, colorCol] %in% grps.1sample), 
      ]
      pp <- pp + stat_density2d(data = df2p2, alpha = 0.3)
    }
  }
  if (addVoronoi) {
    if (colorNumeric) {
      warning("Currently only non-numeric columns are supported for dimRed coloring in combination with density contours. --> skipping density contours")
    }
    else {
      require(deldir)
      center.x <- tapply(df2p[, xLab], df2p[, colorCol], 
                         FUN = function(x) {
                           mean(x, na.rm = TRUE)
                         })
      center.y <- tapply(df2p[, yLab], df2p[, colorCol], 
                         FUN = function(x) {
                           mean(x, na.rm = TRUE)
                         })
      voronoi <- deldir(center.x, center.y, rw = c(min(df2p[, 
                                                            xLab], na.rm = TRUE), max(df2p[, xLab], na.rm = TRUE), 
                                                   min(df2p[, yLab], na.rm = TRUE), max(df2p[, 
                                                                                             yLab], na.rm = TRUE)))
      pp <- pp + geom_segment(aes(x = x1, y = y1, xend = x2, 
                                  yend = y2), size = 1, data = voronoi$dirsgs, 
                              linetype = 1, color = "#C0C0C0")
    }
  }
  if (!is.null(facetCols)) {
    pp <- pp + geom_point(data = df2p[, !(colnames(df2p) %in% 
                                            facetCols)], color = "#C0C0C0", size = ptSize, shape = 16)
  }
  if (!is.null(shapeCol)) {
    pp <- pp + geom_point(aes_string(shape = shapeCol), 
                          size = ptSize) + scale_shape_manual(values = c(19, 
                                                                         15, 17, 4, 3, 18, 8, 1, 0, 2, 6))
  }
  else {
    pp <- pp + geom_point(size = ptSize, shape = 16)
  }
  if (addLabels && is.element("observation", colnames(df2p))) {
    pp <- pp + geom_text(aes(label = observation), size = 2)
  }
  if (!is.null(annot.text)) {
    pp <- pp + annotate("text", x = max(df2p[, xLab], na.rm = TRUE), 
                        y = min(df2p[, yLab], na.rm = TRUE), label = annot.text, 
                        parse = FALSE, hjust = 1, vjust = 1, size = 4)
  }
  if (!is.null(facetCols)) {
    if (length(facetCols) == 1) {
      nr <- ceiling(sqrt(length(unique(df2p[, facetCols]))))
      pp <- pp + facet_wrap(as.formula(paste0("~", facetCols)), 
                            nrow = nr)
    }
    else {
      pp <- pp + facet_grid(as.formula(paste(facetCols, 
                                             collapse = "~")))
    }
  }
  
  return(pp)
}

plot_umap_all <- function(sample_annot, cm, suffix="", ndim=30){
  annot <- sample_annot
  annot$PlateNum <- as.character(annot$PlateNum)
  annot$TF <- annot$TF %>% as.factor
  BL_stallion <- colpal[levels(annot$TF), "color"]
  
  f <- paste0(outdir, "/umapCoord", suffix, ".rds")
  if (file.exists(f)){
    umapCoord <- readRDS(f)
  } else{
    pcaCoord <- muRtools::getDimRedCoords.pca(t(cm), components = 1:ndim)
    umapCoord <- umap(pcaCoord)
    saveRDS(umapCoord, f)
  }
  
  pdf(paste0(plotdir, "/UMAP", suffix, ".pdf"), width=9, height=7)
  
  print(plot_umap_coord(umapCoord, annot=annot, colorCol="TF", shapeCol="Plasmid_Dose", colScheme=BL_stallion, ptSize=5) + 
          theme_classic() + theme(text=element_text(size=20),panel.border = element_rect(color ="black",fill = NA, size = 1)))
  print(plot_umap_coord(umapCoord, annot=annot, colorCol="PlateNum", ptSize=5) + 
          theme_classic() + theme(text=element_text(size=20),panel.border = element_rect(color ="black",fill = NA, size = 1)))
  dev.off()
}

plot_umap_plate <- function(sample_annot, cm, plate, suffix="", ndim=30){
  cols <- sample_annot %>% dplyr::filter(PlateNum==plate) %>% dplyr::select(sampleName) %>% unlist %>% unname
  
  annot <- sample_annot[cols,]
  annot$PlateNum <- as.character(annot$PlateNum)
  annot$TF <- annot$TF %>% as.factor
  BL_stallion <- colpal[levels(annot$TF), "color"]
  
  f <- paste0(outdir, "/umapCoord_plate", plate, suffix, ".rds")
  if (file.exists(f)){
    umapCoord <- readRDS(f)
  } else{
    subcm <- cm[,cols]
    pcaCoord <- muRtools::getDimRedCoords.pca(t(subcm), components = 1:ndim)
    umapCoord <- umap(pcaCoord)
    saveRDS(umapCoord, f)
  }
  # 
  
  pdf(paste0(plotdir, "/UMAP_plate", plate, suffix, ".pdf"), width=9, height=7)
  print(plot_umap_coord(umapCoord, annot=annot, colorCol="TF", shapeCol="Plasmid_Dose", colScheme=BL_stallion, ptSize=5) + 
          theme_classic() + theme(text=element_text(size=20),panel.border = element_rect(color ="black",fill = NA, size = 1)))
  dev.off()
}

# UMAP (quantile norm) ---------------------------------------------------------------------
dsa_norm <- loadDsAcc('../../output/02-atac/01/dsa_norm_quantile/')
cm <- as.data.frame(dsa_norm@counts$filteredConsensus)

# plot all
plot_umap_all(getSampleAnnot(dsa_norm), cm, suffix="_all_qnorm_ndim10", ndim=10)
plot_umap_all(getSampleAnnot(dsa_norm), cm, suffix="_all_qnorm_ndim30", ndim=30)
plot_umap_all(getSampleAnnot(dsa_norm), cm, suffix="_all_qnorm_ndim50", ndim=50)


# plot by plate
plot_umap_plate(getSampleAnnot(dsa_norm), cm, 1, suffix = "_qnorm_ndim10", ndim=10)
plot_umap_plate(getSampleAnnot(dsa_norm), cm, 2, suffix = "_qnorm_ndim10", ndim=10)
plot_umap_plate(getSampleAnnot(dsa_norm), cm, 3, suffix = "_qnorm_ndim10", ndim=10)

plot_umap_plate(getSampleAnnot(dsa_norm), cm, 1, suffix = "_qnorm_ndim30", ndim=30)
plot_umap_plate(getSampleAnnot(dsa_norm), cm, 2, suffix = "_qnorm_ndim30", ndim=30)
plot_umap_plate(getSampleAnnot(dsa_norm), cm, 3, suffix = "_qnorm_ndim30", ndim=30)

plot_umap_plate(getSampleAnnot(dsa_norm), cm, 1, suffix = "_qnorm_ndim50", ndim=50)
plot_umap_plate(getSampleAnnot(dsa_norm), cm, 2, suffix = "_qnorm_ndim50", ndim=50)
plot_umap_plate(getSampleAnnot(dsa_norm), cm, 3, suffix = "_qnorm_ndim50", ndim=50)

# UMAP (rpkm log2 quantile norm) ---------------------------------------------------------------------
dsa_norm <- loadDsAcc('../../output/02-atac/01/dsa_norm_rpkmlog2quantile/')
cm <- as.data.frame(dsa_norm@counts$filteredConsensus)

# plot all
plot_umap_all(getSampleAnnot(dsa_norm), cm, suffix="_all_rpkmlog2quantile_ndim10", ndim=10)
plot_umap_all(getSampleAnnot(dsa_norm), cm, suffix="_all_rpkmlog2quantile_ndim30", ndim=30)
plot_umap_all(getSampleAnnot(dsa_norm), cm, suffix="_all_rpkmlog2quantile_ndim50", ndim=50)


# plot by plate
plot_umap_plate(getSampleAnnot(dsa_norm), cm, 1, suffix = "_rpkmlog2quantile_ndim10", ndim=10)
plot_umap_plate(getSampleAnnot(dsa_norm), cm, 2, suffix = "_rpkmlog2quantile_ndim10", ndim=10)
plot_umap_plate(getSampleAnnot(dsa_norm), cm, 3, suffix = "_rpkmlog2quantile_ndim10", ndim=10)

plot_umap_plate(getSampleAnnot(dsa_norm), cm, 1, suffix = "_rpkmlog2quantile_ndim30", ndim=30)
plot_umap_plate(getSampleAnnot(dsa_norm), cm, 2, suffix = "_rpkmlog2quantile_ndim30", ndim=30)
plot_umap_plate(getSampleAnnot(dsa_norm), cm, 3, suffix = "_rpkmlog2quantile_ndim30", ndim=30)

plot_umap_plate(getSampleAnnot(dsa_norm), cm, 1, suffix = "_rpkmlog2quantile_ndim50", ndim=50)
plot_umap_plate(getSampleAnnot(dsa_norm), cm, 2, suffix = "_rpkmlog2quantile_ndim50", ndim=50)
plot_umap_plate(getSampleAnnot(dsa_norm), cm, 3, suffix = "_rpkmlog2quantile_ndim50", ndim=50)

# UMAP (vst norm) ---------------------------------------------------------------------
dsa_norm <- loadDsAcc('../../output/02-atac/01/dsa_norm_vst/')
cm <- as.data.frame(dsa_norm@counts$filteredConsensus)

# plot all
plot_umap_all(getSampleAnnot(dsa_norm), cm, suffix="_all_vst_ndim10", ndim=10)
plot_umap_all(getSampleAnnot(dsa_norm), cm, suffix="_all_vst_ndim30", ndim=30)
plot_umap_all(getSampleAnnot(dsa_norm), cm, suffix="_all_vst_ndim50", ndim=50)


# plot by plate
plot_umap_plate(getSampleAnnot(dsa_norm), cm, 1, suffix = "_vst_ndim10", ndim=10)
plot_umap_plate(getSampleAnnot(dsa_norm), cm, 2, suffix = "_vst_ndim10", ndim=10)
plot_umap_plate(getSampleAnnot(dsa_norm), cm, 3, suffix = "_vst_ndim10", ndim=10)

plot_umap_plate(getSampleAnnot(dsa_norm), cm, 1, suffix = "_vst_ndim30", ndim=30)
plot_umap_plate(getSampleAnnot(dsa_norm), cm, 2, suffix = "_vst_ndim30", ndim=30)
plot_umap_plate(getSampleAnnot(dsa_norm), cm, 3, suffix = "_vst_ndim30", ndim=30)

plot_umap_plate(getSampleAnnot(dsa_norm), cm, 1, suffix = "_vst_ndim50", ndim=50)
plot_umap_plate(getSampleAnnot(dsa_norm), cm, 2, suffix = "_vst_ndim50", ndim=50)
plot_umap_plate(getSampleAnnot(dsa_norm), cm, 3, suffix = "_vst_ndim50", ndim=50)
