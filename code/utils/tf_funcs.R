############################################################################
###################### functions from Z. Shipony ###########################
############################################################################

chromVar_cluster_TF = function(locs,data,genome,motifs,out_dir,prefix,tf_num,height=12,width=8) {
	#browser()
	rowRanges = makeGRangesFromDataFrame(locs[,1:3])
	counts = as.matrix(data)
	colData = DataFrame(Treatment = names(data),row.names=colnames(counts))
	rse <- SummarizedExperiment(assays=SimpleList(counts=counts), rowRanges=rowRanges, colData=colData)	
	counts_filtered	<- addGCBias(rse, genome = genome)
	motif_ix <- matchMotifs(motifs, counts_filtered, genome = genome)
	dev_TF <- computeDeviations(object = counts_filtered, annotations = motif_ix)
	bg <- getBackgroundPeaks(object = counts_filtered)
	expected <- computeExpectations(counts_filtered)
	variability_TF <- computeVariability(dev_TF)
	variabilty_plot_TF <- plotVariability(variability_TF, use_plotly = FALSE,n=8)
	pdf(paste(out_dir,"/",prefix,"_","variability_plot_TF.pdf",sep=""))
	print(variabilty_plot_TF)
	dev.off()
	sample_cor_TF <- getSampleCorrelation(dev_TF)	
	dist_TF <- as.matrix(as.dist(sample_cor_TF))
	TF_colData <- colData(dev_TF)
	TF_colNames <- TF_colData$Treatment
	colnames(dist_TF) <- TF_colNames
	rownames(dist_TF) <- TF_colNames
	#browser()
	
	pheatmap(dist_TF, clustering_distance_rows = as.dist(1-sample_cor_TF),
         clustering_distance_cols = as.dist(1-sample_cor_TF), filename = paste(out_dir,"/",prefix,"_","Heatmap_TF_colnames.pdf",sep=""),
         height=10,width=10)


	Highest_var_TF <- cbind(variability_TF, rownames(variability_TF))
	Highest_var_TF <- Highest_var_TF %>% group_by(name) %>% dplyr::slice(which.max(variability))
	Highest_var_TF <- Highest_var_TF[order(Highest_var_TF$variability, decreasing = TRUE),]
    to_out = Highest_var_TF
	Highest_var_TF <- Highest_var_TF[c(1:tf_num),]
	Dev_Highest_var_TF <- assays(dev_TF)[[1]]
	Dev_Highest_var_TF <- as.data.frame(Dev_Highest_var_TF)
	Dev_Highest_var_TF <- cbind(rownames(variability_TF), variability_TF$name, Dev_Highest_var_TF)
	Dev_Highest_var_TF <- (subset(Dev_Highest_var_TF, Dev_Highest_var_TF$`rownames(variability_TF)` %in% Highest_var_TF$`rownames(variability_TF)`))
	rownames(Dev_Highest_var_TF) <- Dev_Highest_var_TF$`variability_TF$name`
	Dev_Highest_var_TF <- Dev_Highest_var_TF[,c(3:(dim(data)[2]+2))]
	pheatmap(Dev_Highest_var_TF, cluster_rows = T, cluster_cols = FALSE, scale = "row", filename = paste(out_dir,"/",prefix,"_","Top50_TF_Heatmap.pdf",sep=""),height=height,width=width)
	return(list(to_out,Dev_Highest_var_TF))
}

tf_clusters <- function(tf_matrix) {
    clusters = read_tsv("~/scripts/tf_clusters.tsv",col_names=c("Cluster","TFS"))
    cnts = dim(tf_matrix)[1]
    clst_out = data.frame(row=1:cnts,cluster=c(rep(0,cnts)),name=rep("a",cnts), stringsAsFactors=FALSE)
    for(i in 1:cnts) {
        rownames(tf_matrix)[i] = gsub("::","_",rownames(tf_matrix)[i])
        f = length(grep(paste0("\\b",rownames(tf_matrix)[i],"\\b"),clusters$TFS,ignore.case=TRUE)) == 0
        clst_out$cluster[i] = ifelse(f,0,grep(paste0("\\b",rownames(tf_matrix)[i],"\\b"),clusters$TFS,ignore.case=TRUE))
        clst_out$name[i] =    ifelse(f,"Unknown",clusters$Cluster[grep(paste0("\\b",rownames(tf_matrix)[i],"\\b"),clusters$TFS,ignore.case=TRUE)])
    }
    return(clst_out[order(clst_out$cluster),])
}

getMotifs <- function (species = "Homo sapiens", collection = "CORE", database = JASPAR2018::JASPAR2018, ...) 
{
  opts <- list()
  opts["species"] <- species
  opts["collection"] <- collection
  opts <- c(opts, list(...))
  out <- TFBSTools::getMatrixSet(database, opts)
  if (!isTRUE(all.equal(TFBSTools::name(out), names(out)))) 
    names(out) <- paste(names(out), TFBSTools::name(out), 
                        sep = "_")
  return(out)
}

replaceExtension <- function(oldname, new.extension, get.basename = T) {
  if (get.basename == T) { n <- basename(oldname) }
  else { n <- oldname }
  new <- paste0(tools::file_path_sans_ext(n), new.extension)
}

hypergeom_test <- function(clusters.table, motif_matrix, cluster.col = 4) {
  
  compute_hypergeom <- function(i) {
    cluster <- clusters.table[ ,cluster.col] == i
    print(paste0(i, " ", sum(cluster)))
    bg <- clusters.table[ ,cluster.col] != i
    k <- colSums(motif_matrix[cluster, ])
    m <- colSums(motif_matrix[bg, ])
    n <- length(bg) - m
    N <- sum(cluster)
    k[k>m] =  m[k>m]-1
    df <- data.frame(broom::tidy(phyper(k-1, m, n, N, log.p = T, lower.tail = F)))
    df <- cbind(df, matches = k, bg.matches = m, bg.non.matches = n, cluster.size = N, cluster.num = i)
  }
  
  num.clusters <- max(clusters.table[,cluster.col])
  
  do.call(rbind, lapply(1:num.clusters, compute_hypergeom))
}

lt_hypergeom_test <- function(clusters.table, motif_matrix, cluster.col = 4) {
  
  lt_compute_hypergeom <- function(i) {
    cluster <- clusters.table[ ,cluster.col] == i
    print(paste0(i, " ", sum(cluster)))
    bg <- clusters.table[ ,cluster.col] != i
    k <- colSums(motif_matrix[cluster, ])
    m <- colSums(motif_matrix[bg, ])
    n <- length(bg) - m
    N <- sum(cluster)
    k[k>m] =  m[k>m]-1
    # print(paste(k-1, m, n, N))
    df <- data.frame(broom::tidy(phyper(k-1, m, n, N, log.p = T, lower.tail = T)))
    df <- cbind(df, matches = k, bg.matches = m, bg.non.matches = n, cluster.size = N, cluster.num = i)
  }
  
  num.clusters <- max(clusters.table[,cluster.col])
  
  do.call(rbind, lapply(1:num.clusters, lt_compute_hypergeom))
}

plot_motif_volcano <- function(test.out, dat.subset = NULL, dep.n=12, enr.n=12, title.string="") {
  
  max.match <- max(test.out$matches)
  min.match <- min(test.out$matches)
  
  if (!is.null(dat.subset)) {
    test.out <- test.out %>% filter(cluster.num == dat.subset)
  }
  
  dff1 <- test.out %>%
    filter(log2.enrichment <= 0) %>%
    top_n(dep.n, -log10.p)
  dff2 <- test.out %>% 
    filter(log2.enrichment > 0 ) %>%
    top_n(enr.n, -log10.p)
  dff <- rbind(dff1, dff2)
  
  gg <- ggplot(test.out, aes(x = log2.enrichment, y = -log10.p)) +
    geom_point(size = 4, shape = 21, color = "black", stroke = 0.2, aes(fill = matches)) +
    geom_label_repel(data = dff, aes(label = short.name), size = 3, max.overlaps=30) +
    theme_classic() +
    scale_fill_viridis(option="inferno", limits = c(min.match, max.match)) +
    labs(title = title.string) +
    geom_vline(xintercept = 0) +
    geom_hline(yintercept = 0)
  return(gg)
}

add_split_name <- function(df, colname, split.char, idx, new.col.name = "short.name") {
  split.names <- strsplit(as.character(df[[colname]]), split.char)
  df[[new.col.name]] <- as.character(sapply(split.names, "[", idx))
  return(df)
}

cleanGenomicsDataFrame <- function(coord.data.frame) {
  df <- coord.data.frame %>% filter(!grepl("Un", chrom), !grepl("random", chrom), !grepl("EBV", chrom), !grepl("_", chrom))
  return(df)
}


plot_tf_families <- function(tf.string, test.out, color.scheme = "YlOrRd") {

  filt.df <- data.frame(do.call(cbind, 
                                lapply(1:length(test.out), 
                                       function(x) test.out[[x]]))) %>% 
    rownames_to_column(var = "motif.name") %>% 
    reshape2::melt(value.name = "log10.probability", variable.name = "cluster.num") %>% 
    group_by(cluster.num) %>% 
    mutate(q90 = quantile(log10.probability, probs = c(0,0.1, 0.25, 1))[2]) %>% 
    filter(grepl(tf.string, motif.name))
  
  filt.df$log10.probability[filt.df$log10.probability == -Inf] <- -7000
  
  num.motifs <- length(unique(filt.df$motif.name))
  getPalette <- colorRampPalette(brewer.pal(9, color.scheme))
  gg <- ggplot(filt.df, aes(x = motif.name, y = -log10.probability)) +
    geom_bar(aes(fill = motif.name), color = "black", stat = "identity", position = "dodge") +
    scale_fill_manual(values = getPalette(num.motifs)) +
    guides(fill = guide_legend(title = "motif.name")) +
    geom_hline(aes(yintercept = -q90)) +
    theme(axis.text.x = element_blank()) + 
    facet_grid(. ~ cluster.num, scales = "free_y") +
    labs(title = sprintf("%s family transcription factors by cluster", tf.string))
  return(gg)
}

plot_tfs_in_context <- function(test.out, cluster.num, tf.string) {
  clust <- test.out[[cluster.num]]
  tf.fam <- names(clust[grep(tf.string, names(clust))])
  
  df <- data.frame(clust) %>%
    rownames_to_column()
  colnames(df) <- c("motif.name", "log10.probability")
  df <- df %>%
    mutate(family = ifelse(motif.name %in% tf.fam, T, F))
  df$motif.name <- factor(df$motif.name, levels = df$motif.name[order(df$log10.probability)])
  
  gg <- ggplot(df, aes(x = motif.name, y = log10.probability, color = family, alpha = family)) + 
    geom_point() +
    scale_color_manual(values = c("grey", "red")) +
    scale_alpha_manual(values = c(0.4, 1)) +
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank()) +
    labs(title = sprintf("%s family motif enrichments, cluster %s", tf.string, as.character(cluster.num))) +
    theme_classic()
  return(gg)
}

############################################################################
######################## B. Liu custom functions ###########################
############################################################################

tf_motif_de_analysis <- function(dds, dsa, motifs, motif.scores, peak.meta, 
                                 group1, group2, contrast.col, region, p.cutoff=0.05, log2fc.cutoff=0.58,
                                 save_full=T, deseq_dir="deseq_data", plot_dir="plots", motif_dir="motif_data"){
  #-------------------------------------------------------------------------------
  #' tf_motif_de_analysis
  #' 
  #' Plot MA plot of differential ATAC peaks and Volcano plots showing motif enrichment of up- or down-regulated peaks
  #' @param dds            a DESeq object
  #' @param dsa            a ChrAccR object
  #' @param motifs         a PFMatrixList object 
  #' @param motif.scores   output of motifmatchr::matchMotifs with the option \code{out='scores'}
  #' @param peak.meta      TODO remove this argument
  #' @param group1         denominator/reference group in the differential
  #' @param group2         numerator group in the differential 
  #' @param contrast.col   column name in colData(dds) used for differential
  #' @param region         peak set name in dsa to use
  #' @param p.cutoff       cutoff for p value (passing filter if padj < p.cutoff), default to 0.05
  #' @param log2fc.cutoff  cutoff for log2fc value (passing filter if abs(log2fc) > log2fc.cutoff), default to 0.58
  #' @param save_full      if True save the full DESeq results table from the comparison, else only save the rows passing filter, default to True
  #' @param deseq_dir      directory to save the deseq results
  #' @param plot_dir       directory to save the MA and Volcano plots
  #' @param motif_dir      directory to save the motif enrichment data from differential peaks
  #' @return None
  #' @author Betty Liu
  #' @export
  #-------------------------------------------------------------------------------
  ## get deseq results
  message("reading deseq object")
  de <- results(dds, contrast = c(contrast.col,group2,group1))
  near_gene <- findNearestGeneForGr(dsa@coord[[region]])
  de <- cbind(de,near_gene)
  meta <- dsa@coord[[region]]
  de <- as.data.frame(cbind(as.data.frame(meta, row.names=NULL),de))
  de_filt <- de %>% filter((abs(log2FoldChange)>log2fc.cutoff) & (padj<p.cutoff))
  if(nrow(de_filt)==0){
    message(paste0("No differential peaks found for ", group2, "_vs_", group1, ". Skip plotting."))
    return
  }
  if (!dir.exists(deseq_dir)){
    dir.create(deseq_dir)
  }

  message("writing deseq results table")
  if (save_full){
    write_tsv(de, paste0(deseq_dir, "/deseq_",group2, "_vs_", group1, "_", region,"_full.tsv"))  
  } else{
    write_tsv(de_filt, paste0(deseq_dir, "/deseq_",group2, "_vs_", group1, "_", region,"_p", p.cutoff, "_log2fc",log2fc.cutoff,".tsv"))
  }
  # plot MA plot 
  if (!dir.exists(plot_dir)){
    dir.create(plot_dir)
  }
  
  p1 <- rasterise(ggmaplot(de, fdr=p.cutoff, fc=2^log2fc.cutoff, size=2, top=0, xlab = "Log2 mean coverage", 
                           ggtheme = ggplot2::theme_classic()), dpi=300) +
    ggtitle(paste0("MAplot_",group2,"_vs_",group1,"_",region)) 
  ggsave(plot=p1, paste0(plot_dir, "/MAplot_",group2,"_vs_",group1,"_",region,
                         "_p", p.cutoff, "_log2fc",log2fc.cutoff, ".pdf"))
  
  
  # cluster up and down peaks
  de$cluster = 0
  de$cluster[de$log2FoldChange > log2fc.cutoff & de$padj < p.cutoff] = 1 # pos log2fc 
  de$cluster[de$log2FoldChange < -log2fc.cutoff & de$padj < p.cutoff] = 2 # neg log2fc
  
  # get peak meta data
  counts_df <- peak.meta %>% mutate(k.cluster=de$cluster)
  
  # hypergeometric tests of motif enrichment
  motif.mat <- motifMatches(motif.scores)
  enrichments <- hypergeom_test(counts_df, 
                                as.matrix(motif.mat), 
                                cluster.col = 4)
  lt.depletions <- lt_hypergeom_test(counts_df,
                                     as.matrix(motif.mat),
                                     cluster.col = 4)
  enrichments <- add_split_name(enrichments, "names", "_", 2)
  lt.depletions <- add_split_name(lt.depletions, "names", "_", 2)
  
  tt.enrich <- enrichments %>% 
    inner_join(lt.depletions %>% 
                 dplyr::select(depletion.x = x, short.name, depletion.cluster.num = cluster.num), 
               by = c("short.name", c("cluster.num" = "depletion.cluster.num"))) %>%
    mutate(log10.p = pmin(x, depletion.x)/2) %>% # divide by two because need two tailed test
    mutate(log2.enrichment = log2((matches/cluster.size)/(bg.matches/(bg.matches + bg.non.matches))))
  
  if (!dir.exists(motif_dir)){
    dir.create(motif_dir)
  }
  write_tsv(tt.enrich, paste0(motif_dir, "/motif_enrichment_",group2, "_vs_", group1, "_", region, ".tsv"))
  
    # hypergeometric tests of motif family enrichment
  tt.enrich$family = "ss"
  cnt = 1
  for (i in tt.enrich$names) {
    if(is.null(slot(motifs@listData[[grep(i,names(motifs@listData),fixed=T)]],"tags")$family)) {
      tt.enrich$family[cnt] = "NA"
      cnt = cnt + 1
      next
    }
    tt.enrich$family[cnt] = slot(motifs@listData[[grep(i,names(motifs@listData),fixed=T)]],"tags")$family
    cnt = cnt + 1
  }
  
  tt.enrich[1:(nrow(tt.enrich)%/%2),] %>% group_by(family) %>% summarise(x=mean(x),matches=sum(matches),
                                                       bg.matches=sum(bg.matches),bg.non.matches=sum(bg.non.matches),
                                                       cluster.size = mean(cluster.size), cluster.num= mean(cluster.num),
                                                       depletion.x = mean(depletion.x),log10.p = mean(log10.p),
                                                       log2.enrichment = mean(log2.enrichment)) -> tmp
  tmp$short.name = tmp$family
  tt.enrich[(nrow(tt.enrich)%/%2)+1:nrow(tt.enrich),] %>% group_by(family) %>% summarise(x=mean(x),matches=sum(matches),
                                                         bg.matches=sum(bg.matches),bg.non.matches=sum(bg.non.matches),
                                                         cluster.size = mean(cluster.size), cluster.num= mean(cluster.num),
                                                         depletion.x = mean(depletion.x),log10.p = mean(log10.p),
                                                         log2.enrichment = mean(log2.enrichment)) -> tmp1
  tmp1$short.name = tmp1$family
  tt.family = rbind(tmp,tmp1)
  write_tsv(tt.family, paste0(motif_dir, "/motif_enrichment_family_",group2, "_vs_", group1, "_", region,".tsv"))
  
  # plot volcano plots of tf enrichment
  Name = paste0("tf_motif_enrich_",group2,"_vs_",group1,"_",region,
                "_p", p.cutoff, "_log2fc",log2fc.cutoff)
  pdf(sprintf("%s/Volcanos_%s.pdf", plot_dir, Name), useDingbats = F, width = 5, height = 5)
  lapply(1:max(counts_df[,4]), function(Y) {
    plot_motif_volcano(tt.enrich, dat.subset = Y, dep.n = 0, enr.n = 5, title.string = paste0("Cluster",Y,": ",group2,"/",group1))
  }) %>% print
  lapply(1:max(counts_df[,4]), function(Y) {
    plot_motif_volcano(tt.family, dat.subset = Y, dep.n = 0, enr.n = 5, title.string = paste0("Cluster",Y,": ",group2,"/",group1))
  }) %>% print
  dev.off()
}




getPeakSet.snakeATAC.v2 <- function(sampleAnnot, filePrefixCol, genome, dataDir, sampleIdCol=filePrefixCol, 
                                 type="summits_no_fw", unifWidth=500L, replicateCol=NA, replicatePercReq=1.0, 
                                 replicateConsSelect=FALSE, keepOvInfo=FALSE, capMaxPeaks=NA){
  #-------------------------------------------------------------------------------
  #' getPeakSet.snakeATAC.v2
  #' 
  #' Retrieve a consensus set of ATAC peaks from the snakeATAC pipline run
  #' @param sampleAnnot  data.frame specifying the sample annotation table
  #' @param filePrefixCol column name specifying the file prefix for each sample in the sample annotation table
  #' @param genome       genome assembly
  #' @param dataDir      directory where the files are located
  #' @param sampleIdCol  column name or index in the sample annotation table containing unique sample identifiers
  #' @param type         input data type. Currently only "summits_no_fw" (non-overlapping, fixed-width peaks deduced from summits)
  #' @param unifWidth    width of the peaks if the results have uniform peak lengths
  #' @param replicateCol column name specifying the replicate group for cross-checking coverage across replicates
  #' @param replicatePercReq percentile of replicates in a group required to contain a peak in order to keep it.
  #'                     E.g. a value of 1 (default) means that all replicates in a group are required to contain that peak in order
  #'                     to keep it.
  #' @param replicateConsSelect if set, the peak set will also be checked for consistency, i.e. in order to retain a peak
  #'                     it has to be consistently be present or absent in each replicate group (as specified in \code{replicatePercReq} percent of samples)
  #' @param keepOvInfo   keep annotation columns in the elementMetadata of the results specifying whether a consensus peak overlaps with a
  #'                     peak in each sample
  #' @param capMaxPeaks  the max number of peaks used from each sample, sorted by score
  #' @return \code{GRanges} object containing consensus peak set
  #' @author Fabian Mueller, edited by Betty Liu
  #' @export
  #-------------------------------------------------------------------------------
	library(GenomicRanges)
  if (!is.element(type, c("summits_no_fw", "summits_filt_no_fw"))){
		logger.error(c("Unsupported import type:", type))
	}
	if (replicatePercReq > 1 || replicatePercReq < 0){
		logger.error(c("Invalid value for replicatePercReq. Must be in [0,1]"))
	}

	sampleIds <- as.character(sampleAnnot[,sampleIdCol])
	rownames(sampleAnnot) <- sampleIds

	inputFns <- c()
	if (is.element(type, c("summits_no_fw", "summits_filt_no_fw"))){
		if (nchar(dataDir) > 0){
			fns <- NULL
			if (type=="summits_no_fw"){
				fns <- paste0(sampleAnnot[,filePrefixCol], "_summits.bed")
			} else if (type=="summits_filt_no_fw"){
				fns <- paste0(sampleAnnot[,filePrefixCol], "_peaks.narrowPeak")
			}
			inputFns <- file.path(dataDir, fns)
		} else {
			inputFns <- sampleAnnot[,filePrefixCol]
		}
		names(inputFns) <- sampleIds
	}
	if (!is.character(inputFns)) inputFns <- as.character(inputFns)
	if (!all(file.exists(inputFns))){
		missingSamples <- sampleIds[!file.exists(inputFns)]
		logger.error(c("Missing input files for samples:", paste(missingSamples, collapse=", ")))
	}

	peakFun <- NULL
	res <- NULL
	if (type=="summits_no_fw"){
		peakFun <- function(fn, sid){
			rr <- rtracklayer::import(fn, format="BED")
			rr <- setGenomeProps(rr, genome, onlyMainChrs=TRUE, silent=TRUE)
			rr <- rr[isCanonicalChrom(as.character(seqnames(rr)))]
			# scale scores to their percentiles
			scs <- elementMetadata(rr)[,"score"]
			elementMetadata(rr)[,"score_norm"] <- ecdf(scs)(scs)
			elementMetadata(rr)[,"sampleId"] <- sid

			#extend peaks around the summit
			rr <- trim(promoters(rr, upstream=ceiling(unifWidth/2), downstream=ceiling(unifWidth/2)+1)) #extend each summit on each side by half the width
			rr <- rr[width(rr)==median(width(rr))] #remove too short regions which might have been trimmed

      # sort by score and take the top capMaxPeaks peaks 
      if (!is.na(capMaxPeaks)){
        cap <- min(length(rr), capMaxPeaks)
        rr <- rr[rev(order(elementMetadata(rr)[,"score"]))[1:cap]]
      }
			return(rr)
		}
	} else if (type=="summits_filt_no_fw"){
		peakFun <- function(fn, sid){
			rr <- readMACS2peakFile(fn)
			rr <- setGenomeProps(rr, genome, onlyMainChrs=TRUE, silent=TRUE)
			rr <- rr[isCanonicalChrom(as.character(seqnames(rr)))]
			elementMetadata(rr)[,"calledPeakStart"] <- start(rr)
			elementMetadata(rr)[,"calledPeakEnd"] <- end(rr)
			start(rr) <- end(rr) <- elementMetadata(rr)[,"summit"]

			scs <- elementMetadata(rr)[,"negLog10qval"]
			elementMetadata(rr)[,"score_norm"] <- ecdf(scs)(scs)
			elementMetadata(rr)[,"sampleId"] <- sid
			# only retain peaks with q-value < 0.01
			rr <- rr[elementMetadata(rr)[,"negLog10qval"] > -log10(0.01)]

			rr <- trim(promoters(rr, upstream=ceiling(unifWidth/2), downstream=ceiling(unifWidth/2)+1)) #extend each summit on each side by half the width
			rr <- rr[width(rr)==median(width(rr))] #remove too short regions which might have been trimmed

      # sort by score and take the top capMaxPeaks peaks 
      if (!is.na(capMaxPeaks)){
        cap <- min(length(rr), capMaxPeaks)
        rr <- rr[rev(order(elementMetadata(rr)[,"negLog10qval"]))[1:cap]]
      }
      
			return(rr)
		}
	}

	logger.start("Reading peak sets")
		nSamples <- length(sampleIds)
		peakGrl <- GRangesList(lapply(1:nSamples, FUN=function(i){
			sid <- sampleIds[i]
			logger.status(c("sample:", sid, paste0("(", i, " of ", nSamples, ")")))
			return(peakFun(inputFns[sid], sid))
		}))
		names(peakGrl) <- sampleIds
	logger.completed()

	grps <- NULL
	if (is.element(replicateCol, colnames(sampleAnnot)) && replicatePercReq>0){
		grps <- sampleAnnot[,replicateCol]
	}

	logger.start("Retrieving consensus peak set")
		res <- getConsensusPeakSet(peakGrl, mode="no_by_score", grouping=grps, groupAgreePerc=replicatePercReq, groupConsSelect=replicateConsSelect, scoreCol="score_norm", keepOvInfo=keepOvInfo)
	logger.completed()
	return(res)
}


plot_motifs <- function(fps, motif, dsa, ctl_group, plotdir="plots/footprinting"){
  #-------------------------------------------------------------------------------
  #' plot_motifs
  #' 
  #' Plot motif footprints 
  #' @param fps     output from ChrAccR::getMotifFootprints
  #' @param motif   name of the TF motif
  #' @param dsa     ChrAccR object
  #' @param ctl_group the control group to normalize footprint to 
  #' @param plotdir  directory to save the output plots
  #' @return None
  #' @author Betty Liu
  #' @export
  #-------------------------------------------------------------------------------
    fp <- fps[[motif]]$footprintDf
    meta <- dsa@sampleAnnot
    fp <- fp %>% mutate(TF=meta[fp$sampleId,"TF"], 
                        Dose=meta[fp$sampleId,"Plasmid_Dose"] %>% as.factor)

    tmp <- fp[fp$TF==ctl_group,]  %>% group_by(pos) %>% summarise(avg_ctl_countNormBiasCor = mean(countNormBiasCor), 
                                                               avg_ctl_countNorm = mean(countNorm),
                                                               avg_ctl_Tn5biasNorm = mean(Tn5biasNorm)) %>% as.data.frame
    fp <- fp %>% dplyr::mutate(avg_ctl_countNormBiasCor = lapply(fp$pos, function(p){tmp[tmp$pos==p, "avg_ctl_countNormBiasCor"]}) %>% unlist,
                               avg_ctl_countNorm = lapply(fp$pos, function(p){tmp[tmp$pos==p, "avg_ctl_countNorm"]}) %>% unlist)

    fp <- fp %>% group_by(TF, Dose, pos) %>% summarise(countNorm = mean(countNorm),
                                                       countNormBiasCor = mean(countNormBiasCor),
                                                       countNormBiasCor_ctlNorm = mean(countNormBiasCor / avg_ctl_countNormBiasCor),
                                                       countNorm_ctlNorm = mean(countNorm / avg_ctl_countNorm),
                                                       Tn5biasNorm = mean(Tn5biasNorm))

    # clean up motif name for file name
    motif_out <- gsub("[[:punct:]]", "_", motif)
    dir.create(plotdir, recursive=T, showWarnings=F)
        
    ggplot(fp) + aes(x=pos, y=countNorm, col=Dose) + geom_line() + theme_classic() + facet_grid(TF ~ .)
    ggsave(paste0(plotdir, "/", motif_out,"_countNorm.pdf"), width=3, height=1*length(unique(fp$TF)))

    ggplot(fp) + aes(x=pos, y=countNormBiasCor, col=Dose) + geom_line() + theme_classic() + facet_grid(TF ~ .)
    ggsave(paste0(plotdir, "/", motif_out, "_countNormBiasCor.pdf"), width=3, height=1*length(unique(fp$TF)))

    ggplot(fp) + aes(x=pos, y=countNormBiasCor_ctlNorm, col=Dose) + geom_line(linewidth=0.5) + theme_classic() + facet_grid(TF ~ .)
    ggsave(paste0(plotdir, "/", motif_out, "_countNormBiasCor_ctlNorm.pdf"), width=3, height=1*length(unique(fp$TF)))

    ggplot(fp) + aes(x=pos, y=countNorm_ctlNorm, col=Dose) + geom_line() + theme_classic() + facet_grid(TF ~ .)
    ggsave(paste0(plotdir, "/", motif_out, "_countNorm_ctlNorm.pdf"), width=3, height=1*length(unique(fp$TF)))
    
    ggplot(fp) + aes(x=pos, y=Tn5biasNorm, col=Dose) + geom_line() + theme_classic() + facet_grid(TF ~ .)
    ggsave(paste0(plotdir, "/", motif_out, "_Tn5biasNorm.pdf"), width=3, height=1*length(unique(fp$TF)))
    
    pdf(paste0(plotdir, "/", motif_out, "_chraccrplot.pdf"))
    print(fps[[motif]]$plot)
    invisible(dev.off())
}


plotDDS <- function(dds, cond1, cond2, id2name, contrast.col, 
                    log2fc.cutoff=0.58, padj.cutoff=0.05, labels=c(""), 
                    lfcshrink=T, top.label=50, outdir="output/"){
  # generates MA plot and volcano plot based on deseq object
  # dds: DESeq object
  # cond1: name of condition 1 to compare (numerator in log2FC)
  # cond2: name of condition 2 to compare (denominator in log2FC)
  # id2name: a dataframe that has gene_id as rownames, and columns "external_gene_name" and "description" for conversion
  # contrast.col: the name of the variable to contrast for deseq results
  # log2fc.cutoff: default to 0.58 (=FC of 1.5)
  # padj.cutoff: default to 0.05
  # labels: select labels to display in addition to the top X labels, default to None
  # lfcshrink: whether to use lfcshrink, default to TRUE
  # top.label: top features to label, MA plot sorts by log2FC and volcano sorts by padj, default to 50
  
  savepath <- paste0(outdir, "/deseq_dds_",cond1,"_vs_", cond2, ".rds")
  if (!file.exists(savepath)){
    message(paste0(cond1, " vs ", cond2))
    if (lfcshrink){
        res <- lfcShrink(dds, coef=paste0(contrast.col,"_", cond1, "_vs_", cond2))
      } else{
        res <- results(dds, contrast=c(contrast.col, cond1, cond2))
      }
      
      
      # map gene id to gene name
      res$gene_name <- id2name[rownames(res),"external_gene_name"]
      res$description <- id2name[rownames(res),"description"]
      
      res_filt <- res %>% as.data.frame %>% dplyr::filter(padj<padj.cutoff & abs(log2FoldChange)>log2fc.cutoff) %>% arrange(padj)
      dir.create(outdir, showWarnings=F, recursive=T)
      write.table(res_filt, file=paste0(outdir, "/deseq_",cond1,"_vs_", cond2,"_log2fc",log2fc.cutoff, "_p", padj.cutoff,".txt"), sep='\t',quote=F)
      
      saveRDS(res, savepath)
  } else{
    message(paste0("reading deseq results from ", savepath))
    res <- readRDS(savepath)
  }

  #plot MA plot
  p <- c()
  p[[1]] <- ggmaplot(res, fdr=padj.cutoff, fc=2^log2fc.cutoff, size=2, top=top.label, label.select=labels,
                ggtheme = ggplot2::theme_classic(), genenames = res$gene_name, 
                select.top.method = c("fc")) %>% rasterise(dpi=300)
  vol.labels <- c(labels,arrange(res %>% as.data.frame, padj) %>% 
                        head(top.label) %>% dplyr::select(gene_name) %>% unlist)
  p[[2]] <- EnhancedVolcano(res,
                            lab = res$gene_name,
                            x = 'log2FoldChange',
                            y = 'padj',
                            pCutoff = padj.cutoff,
                            FCcutoff = log2fc.cutoff,
                            selectLab = vol.labels) %>% rasterise(dpi=300)
  return(p)
}

getDimRedPlot_BL <- function (coords, annot = NULL, colorCol = NULL, shapeCol = NULL, alphaCol = NULL, sizeCol = NULL, 
          colScheme = "[auto]", ptSize = 3, addLabels = FALSE, addDensity = FALSE, 
          addVoronoi = FALSE, annot.text = NULL, orderCol = NULL, 
          facetCols = NULL){
  #-------------------------------------------------------------------------------
  #' getDimRedPlot_BL (adapted from muRtools::getDimRedPlot)
  #' 
  #' Plot reduced dimension coordinates
  #' @param coords  dimension reduction coordinates
  #' @param annot   annotation matrix with the same number of rows as coord
  #' @param colorCol     name or index in the annotation matrix (annot) that should be used for coloring the points. if colorCol not supplied but annot is supplied, it defaults to the first annotation column
  #' @param shapeCol    name or index in the annotation matrix (annot) that should be used for point shapes. Note that currently only supports 1 of shapeCol/alphaCol/sizeCol.
  #' @param alphaCol   name or index in the annotation matrix (annot) that should be used for point alpha (transparency)
  #' @param sizeCol   name or index in the annotation matrix (annot) that should be used for point size, this overrides ptSize.
  #' @param colScheme	  color sheme to be used in coloring the points. can be a character vector with the supplied colors. Alternatively, if it is a one-element character vector "[auto]" the color scheme will be selected automatically using muRtools::ggAutoColorScale. If NULL, ggplots default color scheme will be used.
  #' @param ptSize    size of the points in the scatterplot
  #' @param addLabels	 should observation labels be added to each point
  #' @param addDensity	 should Gaussian Kernel density estimation be performed and the contour lines plotted for each color group
  #' @param addVoronoi	 should a Voronoi tessalation grid (based on colorCol) be added to the plot
  #' @param annot.text  optional text to be added in the lower right corner of the plot
  #' @param orderCol	 name or index in the annotation matrix (annot) that should be used for ordering the points. If not NULL Points will be ordered increasingly by their value, i.e. higher-valued points are plottet over lower-valued points
  #' @param facetCols	 name (string) of columns to be used for faceting the resulting plot. Each facet will contain all the points not in the facet as grey points.
  
  #' @return a ggplot2 object containing the dimension reduction plot
  #' @author Betty Liu
  #' @export
  #-------------------------------------------------------------------------------
  if (!is.null(annot)) {
    if (nrow(annot) != nrow(coords)) {
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
    # if (ncol(annot) > 1) {
    #   if (is.null(shapeCol)) {
    #     shapeCol <- 2L
    #   }
    # }
  }
  xLab <- colnames(coords)[1]
  yLab <- colnames(coords)[2]
  df2p <- data.frame(coords)
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
  if (!is.null(alphaCol)) {
    if (length(alphaCol) != 1) {
      stop("alphaCol must be of length 1")
    }
    if (is.numeric(alphaCol)) {
      if (!is.null(colnames(annot))) {
        alphaCol <- colnames(annot)[alphaCol]
      }
    }
    if (is.character(alphaCol)) {
      df2p[, alphaCol] <- annot[, alphaCol]
    }
    else if (is.numeric(alphaCol)) {
      df2p[, "alpha"] <- annot[, alphaCol]
      alphaCol <- "alpha"
    }
    else if (is.logical(alphaCol) && !alphaCol) {
      alphaCol <- NULL
    }
    else {
      stop("invalid value for alphaCol")
    }
    if (!is.numeric(df2p[, alphaCol])) {
      stop("Error: Currently only numeric columns are supported for dimRed shapes")
    }
  }
  
  if (!is.null(sizeCol)) {
    if (length(sizeCol) != 1) {
      stop("sizeCol must be of length 1")
    }
    if (is.numeric(sizeCol)) {
      if (!is.null(colnames(annot))) {
        sizeCol <- colnames(annot)[sizeCol]
      }
    }
    if (is.character(sizeCol)) {
      df2p[, sizeCol] <- annot[, sizeCol]
    }
    else if (is.numeric(sizeCol)) {
      df2p[, "alpha"] <- annot[, sizeCol]
      sizeCol <- "alpha"
    }
    else if (is.logical(sizeCol) && !sizeCol) {
      sizeCol <- NULL
    }
    else {
      stop("invalid value for sizeCol")
    }
    if (!is.numeric(df2p[, sizeCol])) {
      stop("Error: Currently only numeric columns are supported for dimRed shapes")
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
  if (!is.null(rownames(coords))){ 
    df2p$observation <- rownames(coords)}
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
    } else {
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
  } else if (!is.null(alphaCol)) {
    pp <- pp + geom_point(aes_string(alpha = alphaCol, stroke=0), 
                          size = ptSize) + scale_alpha_continuous(range = c(0.2, 1))
  } else if (!is.null(sizeCol)) {
    pp <- pp + geom_point(aes_string(size = sizeCol, stroke=0)) + scale_size_continuous(range=c(1,6))
  } else {
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
  percVar <- attr(coords, "percVar")
  if (!is.null(percVar)) {
    pp <- pp + xlab(paste0(pp$labels$x, " (", round(percVar[pp$labels$x], 
                                                    2), "%)"))
    pp <- pp + ylab(paste0(pp$labels$y, " (", round(percVar[pp$labels$y], 
                                                    2), "%)"))
  }
  return(pp)
}

