
TF_LS=("SPI1" "KLF1" "KLF4" "SOX2" "OCT4" "TCF3" "IRF4" "ALX4" "PRDM1" "ELF1" "LEF1" "SP4")
BASEDIR=../../

for TF in "${TF_LS[@]}"; do
    sbatch -p wjg,sfgf,biochem,owners --mem-per-cpu=8g --time=2:00:00 --out=slurm-logs/chrombp_lr_seq_hm/slurm-%j-$TF.out --job-name=lr_seq_$TF --wrap "python 15d_job_lr_seq_hm.py $TF"
done