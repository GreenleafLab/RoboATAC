# Configuration file 

import os
import sys
import glob
from collections import defaultdict

################## User Inputs ###################
SPECIES_GENOME = "hg38"
FASTQ_DIR = "../../data/fastq"
METADATA_FILE = "./meta.txt"
FASTQ_SCREEN_CONF = "./fastq_screen.conf"
OUTPUT_DIR = "../../output/01-preprocessing/"

# unfortunately not everything in the shared lab resources folders have consistent naming
# 	so always double check if these files exist if you are working outside of hg19/hg38/mm9
BEDS_DICT = {"mm9": "../../data/snakeatac_resources/mm9/mm9_tss.ensembl.bed",
	     "hg19": "../../data/snakeatac_resources/hg19/hg19_tss.ensembl.bed",
	     "hg38": "../../data/snakeatac_resources/hg38/hg38.tss.bed"}
BEDS = {"TSS": BEDS_DICT[SPECIES_GENOME]}
TSS = '../../data/snakeatac_resources/hg38/hg38.tss.bed' if SPECIES_GENOME=="hg38" \
	else '../../data/snakeatac_resources/%s/%s_tss.ensembl.bed' % (SPECIES_GENOME, SPECIES_GENOME)

REFERENCE_FILE = '../../data/snakeatac_resources/genomes/hg38-no-haps/hg38' if SPECIES_GENOME=="hg38" \
	else '../../data/snakeatac_resources/%s/%s' % (SPECIES_GENOME, SPECIES_GENOME)
CHROM_SIZES = '../../data/snakeatac_resources/genomes/gSizes/%s.all.genomsize' % (SPECIES_GENOME)
GENOME_SIZE_DICT = {"mm9": 1.87e9, "sacCer3": 1.2e7, "hg19": 2.7e9, "hg38": 3.0e9}
EFFECTIVE_GENOME_SIZE = GENOME_SIZE_DICT[SPECIES_GENOME]
BLACKLIST = "../../data/snakeatac_resources/%s/%s.blacklist.bed" % (SPECIES_GENOME, SPECIES_GENOME) 
################## End User Inputs ###############

# Adding paths to useful tools for running snakeATAC 
# 	DON'T change this section, the directories are fixed inside the singularity container file system
PICARD_JAR = '/usr/local/bin/picard.jar'
SNAKE_DIR = '/usr/local/snakeATAC/'
ATAC_TOOLS = os.path.join(SNAKE_DIR, 'atac_tools')

# metadata file
def make_meta(filename):
    r1_files = sorted([os.path.abspath(f) for f in glob.glob(FASTQ_DIR+"/*.fastq.gz")])
    if (len(r1_files) < 1):
        sys.exit("No fastqs found.")
    sample_labels = [os.path.basename(r1_file).split("-Z")[0] for r1_file in r1_files]
    with open(filename, 'w') as outfile:
        outfile.write("\t".join(["Name", "Read1"]) + "\n")
        for sample_label, r1_file in zip(sample_labels, r1_files):
            outfile.write("\t".join([sample_label, r1_file]) + "\n")

# make the meta data file
if __name__ == "__main__":
    make_meta(METADATA_FILE)
