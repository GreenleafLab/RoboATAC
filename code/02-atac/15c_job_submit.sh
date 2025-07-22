conda activate chrombpnet
module load cuda/11.2.0 
module load cudnn/8.1.1.33

TF_LS=("SPI1" "KLF1" "KLF4" "SOX2" "OCT4" "TCF3" "IRF4" "ALX4" "NR4A1" "PRDM1" "ELF1" "LEF1" "SP4")
BASEDIR=../../

for TF in "${TF_LS[@]}"; do
    # TF=SPI1
    CHR=("chr1" "chr3" "chr6")
    CMD="python 15b_job_model_eval_gpu.py --TF $TF --chr ${CHR[@]}"
    LOGDIR="slurm-logs/chrombp_eval_gpu"
    OUTDIR="$BASEDIR/output/02-atac/15/"
    OUTPUT_FILE="${OUTDIR}/TF_obs_peaks_test_preds_$TF.pkl"
    # Check if the file exists
    if [ ! -f "$OUTPUT_FILE" ]; then
        echo "File $OUTPUT_FILE does not exist. Submitting job."
        sbatch -p wjg,sfgf,gpu,biochem,owners --gpus 1 --mem=8g --time=1:00:00 --out=$LOGDIR/slurm-%j-$TF-test --job-name=eval_$TF --wrap="$CMD"
    else 
        echo "File $OUTPUT_FILE exists. Skipping."
    fi
done

