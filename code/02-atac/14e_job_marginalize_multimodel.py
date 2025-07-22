# # Marginalize motifs of interest to determine its effect on predicted accessibility
# adapted from Nair et al biorxiv 2023: https://github.com/kundajelab/scATAC-reprog/blob/master/src/analysis/20210608_SMC_BPNet/ZEB_effect.ipynb 
# and https://github.com/kundajelab/chrombpnet/blob/master/chrombpnet/evaluation/make_bigwigs/predict_to_bigwig.py

import numpy as np
import pandas as pd
import keras
import tensorflow as tf
import pyfaidx
import math
import tqdm
import sys
import glob
import os
import argparse
import seaborn as sns
import matplotlib.pyplot as plt
from deeplift.dinuc_shuffle import dinuc_shuffle
from tensorflow.keras.models import load_model
from tensorflow.keras.utils import get_custom_objects
import chrombpnet.training.utils.losses as losses
import chrombpnet.training.utils.one_hot as one_hot

# will run into some cuda warning during import but it's fine, we don't need to run this on a gpu

def load_model_wrapper(model_hdf5):
    # read .h5 model
    custom_objects={"multinomial_nll":losses.multinomial_nll, "tf": tf}    
    get_custom_objects().update(custom_objects)    
    model=load_model(model_hdf5, compile=False)

# will run into some cuda warning during import but it's fine, we don't need to run this on a gpu

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

def revcomp(s):
    # get reverse complement of a DNA sequence
    m = {"A": "T", "C": "G", "G": "C", "T": "A"}
    comp = [m[i] for i in s]
    return("".join(comp[::-1]))

# insert replacement into the sequence between motif_start and motif_end
def mod_sequence(s, rep, motif_start):
    # s: sequence string
    # rep: replacement string
    # motif_start: integer, position in sequence string to insert the left side of motif into 
    motif_end = motif_start + len(rep)
    return s[:motif_start] + rep + s[motif_end:]

def softmax(x):
    norm_x = x - np.mean(x)
    return np.exp(norm_x)/np.sum(np.exp(norm_x))

def parse_inputs():
    parser = argparse.ArgumentParser(description="Parse in silico marginalization arguments.")

    # Define the positional arguments
    parser.add_argument("--motifs", type=str, nargs="*", default=[], help="The motifs to insert separated by space (e.g., GGAA TTCC)")
    parser.add_argument("--motif_starts", type=int, nargs="*", default=[], help="The positions in input sequence to insert each motif (align to motif start) separated by space (e.g., 11 20)")
    parser.add_argument("--name", type=str, default="", help="The name identifier (e.g., ETS_HH_9bp)")
    parser.add_argument("--exp_ls", default=[], nargs="*", type=str, help="List of models to run")
    parser.add_argument("--fold_ls", default=[0], nargs="*", type=int, help="List of model folds to run")
    parser.add_argument("--outdir", type=str, default="", help="Path to output folder")
    parser.add_argument("--plotdir", type=str, default="", help="Path to plot folder")

    # Parse the arguments
    args = parser.parse_args()

    assert len(args.motifs) == len(args.motif_starts)

    # Return the parsed arguments as a list
    return [args.motifs, args.motif_starts, args.name, args.exp_ls, args.fold_ls, args.outdir, args.plotdir]

# ## User Inputs

# user inputs
base_dir = "../../"
models_dir = base_dir + "/output/04-chrombpnet/output/models/"

celltype = "HEK293T"
control = "GFP_d100"

genome = f"{base_dir}/data/chrombp_resources/hg38.fa"
chromsizes = f"{base_dir}/data/chrombp_resources/hg38.chrom.sizes"
control_peaks_path = f"{base_dir}/output/02-atac/01/consensus_peaks_HEK293T_10col.bed"


NARROWPEAK_SCHEMA = ["chr", "start", "end", "1", "2", "3", "4", "5", "6", "summit"]

n_bgd = 1000 # number of backgrounds for prediction

# ## parse inputs
motifs, motif_starts, name, exp_ls, fold_ls, outdir, plotdir = parse_inputs()
os.makedirs(outdir, exist_ok=True)
os.makedirs(plotdir, exist_ok=True)
fold_ls = [str(f) for f in fold_ls]

if os.path.exists(f"{plotdir}/tracks/track_insert_{name}_{n_bgd}bgd.pdf"):
    print(f"{plotdir}/tracks/track_insert_{name}_{n_bgd}bgd.pdf already exists. Skipping.")
    sys.exit(0)


# read genome
hg38 = pyfaidx.Fasta(genome)

# load models
model_control = load_model_wrapper(glob.glob(f"{models_dir}/fold_{fold_ls[0]}/{celltype}_*{control}/models/chrombpnet_nobias.h5")[0])

inputlen = model_control.input_shape[1]
outputlen = model_control.output_shape[0][1]
print(inputlen)
print(outputlen)

# ## insert motifs

np.random.seed(16)

# generate random background sequences from shuffling some control peaks
regions_df =  pd.read_csv(control_peaks_path,
                        sep='\t',
                        names=NARROWPEAK_SCHEMA)

# sample n_bgd * 1.5 number of peaks, buffer for peaks that get clipped and thrown out in the next step
idx = np.random.choice(range(regions_df.shape[0]), size=int(n_bgd*1.5), replace=False)
regions_df = regions_df.loc[idx].reset_index()

# extract sequences with a width of inputlen centered around the given peak summits
seqs_full, seqs, regions_used = get_seq(regions_df, hg38, inputlen)

# sample n_bgd of given peaks
idx = np.random.choice(range(len(seqs)), size=n_bgd, replace=False)
bgd = seqs[idx]
bgd_seqs_full = seqs_full[idx]
bgd_hits = regions_df.loc[idx].reset_index()

# dinucleotide shuffling of the selected bgd peaks
bgd_seqs_full = np.array([dinuc_shuffle(str(i), rng=np.random.RandomState(16)) for i in bgd_seqs_full])

# one hot encoding 
bgd = one_hot.dna_to_one_hot(bgd_seqs_full)

# run predictions for experiments
output = []
output_logits = []
os.makedirs(f"{plotdir}/boxplots", exist_ok=True)
os.makedirs(f"{plotdir}/boxplots_no_outlier", exist_ok=True)
os.makedirs(f"{plotdir}/tracks", exist_ok=True)

for experiment in exp_ls:
    exp_output = []
    exp_output_logits = []
    for f in fold_ls:
        print("---------------------------")
        print(f"{name} {experiment} fold{f}")
        # load model
        model_kd = load_model_wrapper(glob.glob(f"{models_dir}/fold_{f}/{celltype}_*{experiment}/models/chrombpnet_nobias.h5")[0])
        inputlen = model_kd.input_shape[1]
        outputlen = model_kd.output_shape[0][1]
        
        # predict over random backgrounds
        
        base_pred_counts_kd = []
        base_pred_logits_kd = []
        
        cur_raw_seqs = bgd
        cur_pred = model_kd.predict([cur_raw_seqs], verbose=False)

        base_pred_counts_kd = np.vstack(cur_pred[1]).ravel()
        base_pred_logits_kd = np.vstack(cur_pred[0]).squeeze().mean(axis=0)

        # replace with canonical motif and predict
        synth_pred_counts_kd = []
        synth_pred_logits_kd = []
        cur_raw_seqs = []
        for j in tqdm.tqdm(range(n_bgd)):
            tmp = bgd_seqs_full[j]

            # add each motif iteratively
            for k in range(len(motifs)):
                motif = motifs[k]
                motif_start = motif_starts[k]    
                tmp = mod_sequence(tmp, motif, motif_start)

            cur_raw_seqs.append(tmp)

        cur_raw_seqs = one_hot.dna_to_one_hot(cur_raw_seqs)
        cur_pred = model_kd.predict([cur_raw_seqs], verbose=False)

        synth_pred_counts_kd = np.vstack(cur_pred[1]).ravel()
        synth_pred_logits_kd = np.vstack(cur_pred[0]).squeeze().mean(axis=0)

        df = pd.DataFrame([base_pred_counts_kd, synth_pred_counts_kd], 
                        index=[f"base_pred_counts_{experiment}_fold{f}", f"synth_pred_counts_{experiment}_fold{f}"]).transpose()
        df2 = pd.DataFrame([base_pred_logits_kd, synth_pred_logits_kd],
                        index=[f"base_pred_avglogits_{experiment}_fold{f}", f"synth_pred_avglogits_{experiment}_fold{f}"]).transpose()

        exp_output.append(df)
        exp_output_logits.append(df2)
    
    # summarise all folds
    allfolds = pd.concat(exp_output, axis=1)
    allfolds_logits = pd.concat(exp_output_logits, axis=1)
    allfolds[f"base_pred_counts_{experiment}"] = allfolds.loc[:,allfolds.columns.str.contains(f"base_pred_counts_{experiment}")].mean(axis=1)
    allfolds[f"synth_pred_counts_{experiment}"] = allfolds.loc[:,allfolds.columns.str.contains(f"synth_pred_counts_{experiment}")].mean(axis=1)
    allfolds_logits[f"base_pred_avglogits_{experiment}"] = allfolds_logits.loc[:,allfolds_logits.columns.str.contains(f"base_pred_avglogits_{experiment}")].mean(axis=1)
    allfolds_logits[f"synth_pred_avglogits_{experiment}"] = allfolds_logits.loc[:,allfolds_logits.columns.str.contains(f"synth_pred_avglogits_{experiment}")].mean(axis=1)
    
    # append to output
    output.append(allfolds)
    output_logits.append(allfolds_logits)

# save concatenated output
df = pd.concat(output, axis=1)
df.to_csv(f"{outdir}/insert_{name}_{n_bgd}bgd.csv", sep=",")

df2 = pd.concat(output_logits, axis=1)
df2.to_csv(f"{outdir}/insert_{name}_{n_bgd}bgd_avglogits.csv", sep=",")

# plot per TF_dose avg across folds--------------
# plot fold change
l = df.shape[0]
x = np.array([[i]*l for i in exp_ls]).flatten().tolist()
y = np.array([(df[f"synth_pred_counts_{i}"] - df[f"base_pred_counts_{i}"])/np.log(2) for i in exp_ls]).flatten().tolist()

plt.figure(figsize=(10, 10))
p = sns.boxplot(x=x, y=y)
plt.axhline(ls='--')
p.set_xticklabels(exp_ls,size = 16)
plt.ylabel("Log2 Fold change", fontsize=16)
plt.title(f"Motif {name} Insertion", fontsize=20)
plt.savefig(f"{plotdir}/boxplots/insert_{name}_{n_bgd}bgd.pdf", dpi=150)

# plot fold change no outliers
plt.figure(figsize=(10, 10))
p = sns.boxplot(x=x, y=y, showfliers=False)
plt.axhline(ls='--')
p.set_xticklabels(exp_ls,size = 16)
plt.ylabel("Log2 Fold change", fontsize=16)
plt.title(f"Motif {name} Insertion", fontsize=20)
plt.savefig(f"{plotdir}/boxplots_no_outlier/insert_{name}_{n_bgd}bgd_nooutlier.pdf", dpi=150)


# plot tracks
fig, ax = plt.subplots(len(exp_ls)+1, 1, figsize=(20,18))
cols = ["black", "blue", "green", "orange", "red", "magenta"]
for i in range(len(exp_ls)):
    n = exp_ls[i]
    avg_probs = softmax(df2[f"base_pred_avglogits_{n}"])
    avg_counts = np.exp(df[f"base_pred_counts_{n}"]).mean()
    avg_profile_base = avg_probs * avg_counts

    avg_probs = softmax(df2[f"synth_pred_avglogits_{n}"])
    avg_counts = np.exp(df[f"synth_pred_counts_{n}"]).mean()
    avg_profile_synth = avg_probs * avg_counts

    sns.lineplot(avg_profile_synth - avg_profile_base, label=n, color=cols[i], ax=ax[0])

    sns.lineplot(avg_profile_base, label=f"bgd", color="gray", ax=ax[i+1])
    sns.lineplot(avg_profile_synth, label=n, color=cols[i], ax=ax[i+1])
    ax[i+1].set_ylabel("predicted counts")

ax[0].set_title(f"inserted {name}")
ax[0].set_ylabel("predicted counts (synth - base)")
plt.xlabel("bp position")

# highlight the motif locations
for a in ax:
    for j in range(len(motifs)):
        pos_to_insert = motif_starts[j] - (inputlen - outputlen)//2 
        a.axvspan(pos_to_insert, pos_to_insert + len(motifs[j]), ymin=0, ymax=1, alpha=0.2, color='red')


plt.tight_layout()
plt.savefig(f"{plotdir}/tracks/track_insert_{name}_{n_bgd}bgd.pdf", dpi=150)

# plot per fold -------------
# plot fold change
l = df.shape[0]
toplot = []
for f in fold_ls: 
    x = np.array([[i]*l for i in exp_ls]).flatten().tolist()
    y = np.array([(df[f"synth_pred_counts_{i}_fold{f}"] - df[f"base_pred_counts_{i}_fold{f}"])/np.log(2) for i in exp_ls]).flatten().tolist()
    toplot.append(pd.DataFrame([x,y,[f]*len(y)], index=["exp", "log2fc", "fold"]).transpose())

toplot = pd.concat(toplot, axis=0)

plt.figure(figsize=(10, 10))
p = sns.boxplot(data=toplot, x="exp", y="log2fc", hue="fold")
plt.axhline(ls='--')
p.set_xticklabels(exp_ls,size = 16)
plt.ylabel("Log2 Fold change", fontsize=16)
plt.title(f"Motif {name} Insertion", fontsize=20)
plt.savefig(f"{plotdir}/boxplots/insert_{name}_{n_bgd}bgd_perfold.pdf", dpi=150)
