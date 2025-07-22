
# RoboATAC
![](img/roboatac-logo.png) 

This repository accompanies the preprint Liu et al, bioRxiv 2025b.

## Code

* `code`
    * [`utils`](https://github.com/GreenleafLab/RoboATAC/tree/main/code/utils): helper functions for motif analysis, ChromBPNet motif clustering, and color palette
    * [`00-techdev`](https://github.com/GreenleafLab/RoboATAC/tree/main/code/00-techdev)
        * `01`: comparing OmniATAC and RoboATAC QC metrics from benchmarking datasets
        * `02`: script to check fraction of reads in peaks
    * [`01-preprocessing`](https://github.com/GreenleafLab/RoboATAC/tree/main/code/01-preprocessing)
        * Using the **SE** branch of the [snakeATAC](https://github.com/GreenleafLab/snakeATAC_singularity/) preprocessing pipeline
        * `Snakefile.py`: snakemake file to process single-ended Ultima reads of RoboATAC libraries
    * [`02-atac`](https://github.com/GreenleafLab/RoboATAC/tree/main/code/02-atac)
        * `01`: create ChrAccR object, consensus peak calling, read normalization, chromVAR
        * `02`: dimensionality reduction with PCA
        * `03`: matching sequences in consensus peaks to JASPAR2020 motifs
        * `04`: differential peak analysis
        * `05`: correlating ATAC and RNA differentials
        * `06`: TF footprinting
        * `07`: dimensionality reduction with UMAP 
        * `08`: analyze peak type compositions of differential peak sets
        * `09`: motif scores and motif counts within differential peak sets
        * `10`: linear and Hill fits of peak dose response to determine peak sensitivity group
        * `11`: correlating motif scores and motif counts with Hill-fitted parameters
        * `12`: calling nucleosome position with NucleoATAC
        * `13`: calculate motif distance to nucleosomes
        * `14`: in silico marginalization with ChromBPNet models
        * `15`: ChromBPNet model performance evaluation, multinomial logistic regression models
        * `16`: hit dose analysis, PWM and pileups of different hit dose sets
        * `17`: overlap with ENCODE ChIP-seq data 
        * `18`: ChromHMM annotations
        * `19`: motif pattern distribution in the genome
    * [`03-rna`](https://github.com/GreenleafLab/RoboATAC/tree/main/code/03-rna)
        * `00`: snakemake pipeline to preprocess RNA data with kallisto
        * `01`: PCA, differential analysis
        * `02`: plot average TPMs for overexpressed TFs
        * `03`: GO term enrichment for differential gene sets
    * [`04-chrombpnet`](https://github.com/GreenleafLab/RoboATAC/tree/main/code/04-chrombpnet)
        * snakemake pipeline to prepare input regions, train ChromBPNet models, interpret models, discovery motifs, and identify hit instances
    * [`05-bravo`](https://github.com/GreenleafLab/RoboATAC/tree/main/code/05-bravo)
        * scripts and device configuration file for running RoboATAC on an Agilent Bravo liquid handling robot (NGS Option B layout)


## Citation
If you use this data or code, please cite: Liu et al, bioRxiv 2025.


