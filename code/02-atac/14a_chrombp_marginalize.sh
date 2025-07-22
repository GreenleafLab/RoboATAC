#!/bin/bash
set -euo pipefail

eval "$(conda shell.bash hook)"
base_dir=../../

# TF=LEF1
for TF in LEF1 KLF4 IRF4 TCF3 ALX4 SP4 PRDM1 NR4A1 OCT4 SOX2; do
    ########################
    # step 1: convert the modisco CWMs to PFMs for gimme cluster
    ########################
    conda activate chrombpnet
    config_all=14b_modisco_all.config
    config=tmp_modisco.config
    cat $config_all | grep $TF > $config
    out_file1=$base_dir/output/02-atac/14/HEK293T_$TF/gimme/pos.pfm
    out_file2=$base_dir/output/02-atac/14/HEK293T_$TF/gimme/neg.pfm
    mkdir -p $(dirname $out_file1)
    mkdir -p $(dirname $out_file2)

    python ../utils/modisco_to_pfm.py -c $config -o $out_file1 -p pos_patterns
    python ../utils/modisco_to_pfm.py -c $config -o $out_file2 -p neg_patterns

    ########################
    # step 2: use gimme cluster to cluster patterns
    ########################
    conda activate gimme
    #t=1
    #t=0.9
    t=0.8 
    #t=0.3 # cluster threshold
    gimme cluster $out_file1 $(dirname $out_file1)/pos_$t -t $t -N 16
    gimme cluster $out_file2 $(dirname $out_file2)/neg_$t -t $t -N 16
done

########################
# step 2b: look at the gimme cluster outputs and do some manual reclustering (e.g. head-head, head-tail, tail-tail separation)
# also rename the consensus motifs based on best matches from looking at modisco reports
########################

########################
# step 3: merge patterns based on gimme clusters
########################
TF=TCF3
conda activate chrombp
outdir=$base_dir/output/02-atac/14/HEK293T_$TF/merge_modisco/
key=$base_dir/output/02-atac/14/HEK293T_$TF/gimme/pos_0.8/cluster_key_manual_BL.txt
modisco_dir=$base_dir/output/04-chrombpnet/output/models/fold_0/
contribs_dir=$base_dir/output/04-chrombpnet/output/models/fold_0/

mkdir -p $outdir
python ../utils/merge_modisco.py --out-dir $outdir --model-head counts --cluster-key $key --modisco-dir $modisco_dir --contribs-dir $contribs_dir --batch 1
python ../utils/merge_modisco.py --out-dir $outdir --model-head counts --cluster-key $key --modisco-dir $modisco_dir --contribs-dir $contribs_dir --batch 2
python ../utils/merge_modisco.py --out-dir $outdir --model-head counts --cluster-key $key --modisco-dir $modisco_dir --contribs-dir $contribs_dir --batch 3
python ../utils/merge_modisco.py --out-dir $outdir --model-head counts --cluster-key $key --modisco-dir $modisco_dir --contribs-dir $contribs_dir --batch 4

########################
# step 4: merge into a single modisco object
########################
conda activate chrombp
merged_outdir=$base_dir/output/02-atac/14/HEK293T_$TF/merge_modisco/
compiled_outdir=$base_dir/output/02-atac/14/HEK293T_$TF/modisco_compiled/
mkdir -p $compiled_outdir
python ../utils/compile_modisco_obj.py --out-dir $compiled_outdir --merged-modisco-dir $merged_outdir

########################
# step 5: hit calling on this merged modisco object using global peak set
########################
# use finemo-gpu, run through chrombp_smk 

########################
# step 6: marginalize using hit calling results
########################
# run 14c_marginalize.ipynb interactively