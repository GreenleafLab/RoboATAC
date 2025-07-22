# %% [markdown]
# ## Eval model performance across TF dosage
# Adapted from Surag Nair

# conda activate chrombpnet
# module load cuda/11.2.0 
# module load cudnn/8.1.1.33

# %%
import pandas as pd
import numpy as np
import tensorflow as tf
import keras
import pyBigWig
import pyfaidx
import scipy.stats
import matplotlib.pyplot as plt
import sys
import os
import glob
import json
import seaborn as sns
import pickle
import argparse
import chrombpnet.training.utils.one_hot as one_hot
import chrombpnet.training.utils.losses as losses
from tensorflow.keras.models import load_model
from tensorflow.keras.utils import get_custom_objects
import chrombpnet.evaluation.make_bigwigs.bigwig_helper as bigwig_helper
from collections import Counter

# %%
def load_model_wrapper(model_hdf5):
    # read .h5 model
    custom_objects={"multinomial_nll":losses.multinomial_nll, "tf": tf}    
    get_custom_objects().update(custom_objects)    
    model=load_model(model_hdf5, compile=False)
    print("got the model")
    model.summary()
    return model

def get_seq(peaks_df, genome, width):
    """
    fetches sequence from a given genome.
    adapted from chrombpnet.evaluation.make_bigwigs.bigwig_helper.get_seq()
    """
    vals = []
    peaks_used = []
    for i, r in peaks_df.iterrows():
        sequence = str(genome[r['chr']][(r['start']+r['summit'] - width//2):(r['start'] + r['summit'] + width//2)])
        if len(sequence) == width:
            vals.append(sequence)
            peaks_used.append(True)
        else:
            peaks_used.append(False)
    return np.array(vals), one_hot.dna_to_one_hot(vals), np.array(peaks_used)

# insert replacement into the sequence between motif_start and motif_end
def mod_sequence(s, rep, motif_start, motif_end):
    return s[:motif_start] + rep + s[motif_end:]

def softmax(x):
    norm_x = x - np.mean(x)
    return np.exp(norm_x)/np.sum(np.exp(norm_x))

def parse_inputs():
    parser = argparse.ArgumentParser(description="Parse TF and chromosomes.")
    
    parser.add_argument("--TF", type=str, default="SPI1", help="A transcription factor, e.g., SPI1")
    
    parser.add_argument(
        "--chr",
        type=str,
        nargs='+',  # Accepts one or more arguments
        default=["chr1", "chr3", "chr6"],
        help="A list of chromosomes, e.g., chr1 chr3 chr6"
    )

    parser.add_argument("--group", type=str, default="test", help="test or train, used to name outputs")
    
    args = parser.parse_args()
    return args


# %%
# user inputs
args = parse_inputs()
print(args)
TF = args.TF
TEST_CHRS = args.chr
group=args.group

base_dir = "../../"
models_dir = base_dir + "/output/04-chrombpnet/output/models/fold_0/"
celltype = "HEK293T"
control = "GFP_d100"

genome = f"{base_dir}/data/chrombp_resources/hg38.fa"
chromsizes = f"{base_dir}/data/chrombp_resources/hg38.chrom.sizes"
control_peaks_path = f"{base_dir}/output/02-atac/01/consensus_peaks_HEK293T_10col.bed"

outdir = f"{base_dir}/output/02-atac/15/"
plotdir = f"{base_dir}/output/02-atac/15/"
os.makedirs(outdir, exist_ok=True)
os.makedirs(plotdir, exist_ok=True)

NARROWPEAK_SCHEMA = ["chr", "start", "end", "1", "2", "3", "4", "5", "6", "summit"]
INP_LEN, OUT_LEN = 2114, 1000

n_bgd = 1000 # number of backgrounds for prediction
batch_size = 1000 # number of bgds per batch, used to speed up processing

SAMPLES = ["GFP_d100", f"{TF}_d005", f"{TF}_d025", f"{TF}_d050", f"{TF}_d075", f"{TF}_d100"]
hg38 = pyfaidx.Fasta(genome)


# %% [markdown]
# ### Load models

# %%
models = {}
for m in SAMPLES:
    models[m] = load_model_wrapper(glob.glob(f"{models_dir}/{celltype}_*{m}/models/chrombpnet_nobias.h5")[0])


# %% [markdown]
# ## Pre-existing peak sets
# 
# Take peak sets that we called (open nonsensitive, closed nonsensitive, saturating sensitive, nonsaturating sensitive) and see avg prediction for those. 

# %%
# or read the full set of motif containing peaks
obs_peaks = pd.read_csv(f"{base_dir}/output/02-atac/10/{TF}_motif_containing_peak_meta.tsv", sep="\t")

# %%
allpeaks = pd.read_csv(f"{base_dir}/output/02-atac/01/consensus_peaks_HEK293T.bed", sep='\t',
                            names=['chr','start','end','name','x1','x2'])
allpeaks['start'] = allpeaks['start'].astype(int)
allpeaks['end'] = allpeaks['end'].astype(int)
allpeaks['summit'] = (allpeaks['end'] - allpeaks['start'])//2
allpeaks.index = allpeaks['name']
# %%
obs_peaks = obs_peaks.join(allpeaks, how='left')

# %%
obs_peaks_test = obs_peaks[obs_peaks['chr'].isin(TEST_CHRS)]

# %%
seqs_full, obs_peaks_test_seqs, regions_used = get_seq(obs_peaks_test, hg38, INP_LEN)

# %%
# takes 6hr on cpu, 20min on gpu
for s in SAMPLES:
    preds = models[s].predict(obs_peaks_test_seqs, verbose=True, batch_size=128)
    with open(f"{outdir}/TF_obs_peaks_{group}_preds_{TF}_{s}.pkl", 'wb') as f:
        pickle.dump(preds, f)
    del preds

obs_peaks_test_preds = {}
# pickle dump
for s in SAMPLES:
    with open(f"{outdir}/TF_obs_peaks_{group}_preds_{TF}_{s}.pkl", 'rb') as f:
        obs_peaks_test_preds[s] = pickle.load(f)

with open(f"{outdir}/TF_obs_peaks_{group}_preds_{TF}.pkl", 'wb') as f:
    pickle.dump(obs_peaks_test_preds, f)
