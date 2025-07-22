#!/usr/bin/bash
snakemake --unlock -ns Snakefile.py
snakemake -p -j 500 -s Snakefile.py --cluster "sbatch -p sfgf,wjg,biochem -n 1 -t 48:00:00 --mem-per-cpu 32g"