#!/bin/bash
test=$1

ml biology samtools bedtools


if [ "$test" = "dry" ]; then
    snakemake -n --profile profile --configfile inputs.yaml
else
    snakemake --profile profile --configfile inputs.yaml
fi

## getting the K562 pretrained bias model (bulk K562 ATAC)
# wget https://storage.googleapis.com/chrombpnet_data/input_files/bias_models/ATAC/ENCSR868FGK_bias_fold_0.h5 -O output/model_bias/ENCSR868FGK_bias_fold_0.h5
