# %% [markdown]
# ## Eval model performance across TF dosage
# Adapted from Surag Nair

# %%
import pandas as pd
import numpy as np
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
import chrombpnet.training.utils.one_hot as one_hot
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

# %%
# user inputs
base_dir = "../../"
models_dir = base_dir + "/output/04-chrombpnet/output/models/fold_0/"
celltype = "HEK293T"
control = "GFP_d100"
# TF = "KLF1" # example
TF = sys.argv[1]

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

TEST_CHRS = ['chr1','chr3','chr6']
SAMPLES = ["GFP_d100", f"{TF}_d005", f"{TF}_d025", f"{TF}_d050", f"{TF}_d075", f"{TF}_d100"]
hg38 = pyfaidx.Fasta(genome)

# %% [markdown]
# ### Normalizing factors

# %%
import subprocess

DEPTH_VECTOR = []

counts_file = f"{base_dir}/output/01-preprocessing/output/bams/qc/compiled_counts.txt"
for s in SAMPLES:
    cmd = f"cat {counts_file}| grep {s} | awk '{{sum += $5}} END {{print sum}}'"
    DEPTH_VECTOR.append(int(subprocess.check_output(cmd, text=True, shell=True).strip()))

DEPTH_VECTOR = np.array(DEPTH_VECTOR)/1e6
print(DEPTH_VECTOR)

# %% [markdown]
# ## Pre-existing peak sets
# 
# Take peak sets that I called (open nonsensitive, closed nonsensitive, saturating sensitive, nonsaturating sensitive and see avg prediction for those. 

# %%
obs_peaks = pd.read_csv(f"{base_dir}/output/02-atac/10/{TF}_motif_containing_peak_meta.tsv", sep="\t")
obs_peaks.head()


# %%
allpeaks = pd.read_csv(f"{base_dir}/output/02-atac/01/consensus_peaks_HEK293T.bed", sep='\t',
                            names=['chr','start','end','name','x1','x2'])
allpeaks['start'] = allpeaks['start'].astype(int)
allpeaks['end'] = allpeaks['end'].astype(int)
allpeaks['summit'] = (allpeaks['end'] - allpeaks['start'])//2
allpeaks.index = allpeaks['name']
allpeaks.head()

# %%
obs_peaks = obs_peaks.join(allpeaks, how='left')

# %%
obs_peaks_test = obs_peaks[obs_peaks['chr'].isin(TEST_CHRS)]
obs_peaks_test.shape

# %%
# pickle load
with open(f"{outdir}/TF_obs_peaks_test_preds_{TF}.pkl", 'rb') as f:
    obs_peaks_test_preds = pickle.load(f)

# %%
obs_peaks_test_preds_500bp_cpm = []
for i,s in enumerate(SAMPLES):
    center_cts = (scipy.special.softmax(obs_peaks_test_preds[s][0], axis=-1) * (np.exp(obs_peaks_test_preds[s][1])-1))[:, 250:750].sum(-1)
    obs_peaks_test_preds_500bp_cpm.append(center_cts/DEPTH_VECTOR[i])

obs_peaks_test_preds_500bp_cpm =  np.array(obs_peaks_test_preds_500bp_cpm).T

# %% [markdown]
# Same as above but with actual data.

# %%
obs_peaks_test_obs_500bp_cpm = []

for j in range(len(SAMPLES)):
    s = SAMPLES[j]
    print(s)
    with pyBigWig.open(glob.glob(f"{models_dir}/{celltype}_*{s}/auxiliary/data_unstranded.bw")[0]) as bw:
        counts = []
        reg = obs_peaks_test.reset_index()
        for i in range(reg.shape[0]):
            start = reg.loc[i, "start"] + reg.loc[i, "summit"] - 500//2
            end = start + 500
            val = np.nansum(bw.values(reg.loc[i, "chr"], start, end))
            counts.append(val)
        obs_peaks_test_obs_500bp_cpm.append(counts/DEPTH_VECTOR[j])

obs_peaks_test_obs_500bp_cpm =  np.array(obs_peaks_test_obs_500bp_cpm).T


# %% [markdown]
# ### heatmap of obs and pred, plus ChromHMM annotation

# %%
chromhmm = pd.read_csv(f"{base_dir}/output/02-atac/18/{TF}_chromHMM_annot_1to1_peak2state.csv", sep=",")
chromhmm.index = chromhmm.peakname
print(chromhmm.shape)
chromhmm.head()

# %%
obs_peaks_test = pd.merge(obs_peaks_test, chromhmm, left_index=True, right_index=True, how="left")
obs_peaks_test.chromhmm_state[np.where((obs_peaks_test.chromhmm_state).isnull())[0]] = "none"
obs_peaks_test.head()

# %%
obs_peaks_test["chromhmm_state_id"] = obs_peaks_test.chromhmm_state.str.split("_").str[0]
obs_peaks_test.chromhmm_state_id[obs_peaks_test.chromhmm_state_id=="none"] = "-1"
obs_peaks_test["chromhmm_state_broad"] = obs_peaks_test.chromhmm_state.str.replace(r'[0-9\W_]+', '', regex=True)

from sklearn.preprocessing import LabelEncoder
label_encoder = LabelEncoder()
numeric_labels = label_encoder.fit_transform(obs_peaks_test.chromhmm_state_broad)
obs_peaks_test['chromhmm_state_id_broad'] = numeric_labels
obs_peaks_test.head()

# %%
from matplotlib.colors import ListedColormap
from scipy.cluster.hierarchy import linkage, dendrogram

idx = np.argsort(obs_peaks_test['group']).values
## save ordered peaks
# obs_peaks_test.iloc[idx].to_csv(f"{outdir}/{TF}_obs_peaks_test_sorted.csv")

toplot = obs_peaks_test_obs_500bp_cpm[idx]
# Map groups to colors
sorted_groups = obs_peaks_test['group'][idx]
unique_groups = np.unique(sorted_groups)

sorted_groups2 = obs_peaks_test['chromhmm_state_broad'][idx]
unique_groups2 = np.unique(list(sorted_groups2))

# Perform hierarchical clustering on the subset of rows for the specific group
out = []
og_idx = []
for gr in unique_groups:
    group_data = toplot[sorted_groups == gr]
    linkage_matrix = linkage(group_data, method='ward')
    dendro = dendrogram(linkage_matrix, no_plot=True)
    row_order = dendro['leaves']
    og_idx.extend(idx[sorted_groups == gr][row_order])
    out.extend(group_data[row_order])

cmap = ListedColormap(sns.color_palette("Set1", n_colors=len(unique_groups)).as_hex())
group_colors = pd.Series(sorted_groups).map(dict(zip(unique_groups, cmap.colors)))
cmap2 = ListedColormap(sns.color_palette("Set2", n_colors=len(unique_groups2)).as_hex())
group_colors2 = pd.Series(sorted_groups2).map(dict(zip(unique_groups2, cmap2.colors)))

row_colors = pd.DataFrame({
    "sensitivity_group": group_colors,
    "chromhmm_state": group_colors2
})
np.random.seed(16)
g = sns.clustermap(out, row_cluster=False, col_cluster=False, 
                row_colors=[group_colors, group_colors2], cmap="viridis", figsize=(8, 8), 
                xticklabels=[0, 0.05, 0.25, 0.5, 0.75, 1], yticklabels=False,
                vmin = np.percentile(out, 5),
                vmax = np.percentile(out, 95),
                )
legend_labels = [plt.Line2D([0], [0], marker='o', color=c, linestyle='') for c in cmap.colors]
plt.legend(legend_labels, unique_groups, title="Sensitivity Group", loc="center left", bbox_to_anchor=(4, 0.5))
legend_labels2 = [plt.Line2D([0], [0], marker='o', color=c, linestyle='') for c in cmap2.colors]
plt.legend(legend_labels, unique_groups2, title="ChromHMM State", loc="center left", bbox_to_anchor=(0.5, 4))
# Add a label to the color bar
colorbar = g.ax_heatmap.collections[0].colorbar
colorbar.set_label("Observed CPM", rotation=270, labelpad=20)  # Rotation and padding for better alignment
g.ax_heatmap.set_xlabel(f"{TF} dose", fontsize=12)


g.savefig(f"{plotdir}/{TF}_obs_cpm_heatmap.pdf", dpi=300)
g.savefig(f"{plotdir}/{TF}_obs_cpm_heatmap.png", dpi=600)

# %%
toplot = obs_peaks_test_preds_500bp_cpm[og_idx]

# using the same color scale as obs
g = sns.clustermap(toplot, row_cluster=False, col_cluster=False, 
                row_colors=list(group_colors), cmap="viridis", figsize=(8, 8), 
                xticklabels=[0, 0.05, 0.25, 0.5, 0.75, 1], yticklabels=False,
                vmin = np.percentile(out, 5), # keep the same vmin vmax as the preds color scale
                vmax = np.percentile(out, 95),
                )
legend_labels = [plt.Line2D([0], [0], marker='o', color=c, linestyle='') for c in cmap.colors]
plt.legend(legend_labels, unique_groups, title="Group", loc="center left", bbox_to_anchor=(4, 0.5))
# Add a label to the color bar
colorbar = g.ax_heatmap.collections[0].colorbar
colorbar.set_label("Predicted CPM", rotation=270, labelpad=20)  # Rotation and padding for better alignment
g.ax_heatmap.set_xlabel(f"{TF} dose", fontsize=12)

g.savefig(f"{plotdir}/{TF}_preds_cpm_heatmap.pdf", dpi=300)
g.savefig(f"{plotdir}/{TF}_preds_cpm_heatmap.png", dpi=600)

# %%
# or using a diff color scale than obs
g = sns.clustermap(toplot, row_cluster=False, col_cluster=False, 
                row_colors=list(group_colors), cmap="viridis", figsize=(8, 8), 
                xticklabels=[0, 0.05, 0.25, 0.5, 0.75, 1], yticklabels=False,
                vmin = np.percentile(toplot, 5),
                vmax = np.percentile(toplot, 95),
                )
legend_labels = [plt.Line2D([0], [0], marker='o', color=c, linestyle='') for c in cmap.colors]
plt.legend(legend_labels, unique_groups, title="Group", loc="center left", bbox_to_anchor=(4, 0.5))
# Add a label to the color bar
colorbar = g.ax_heatmap.collections[0].colorbar
colorbar.set_label("Predicted CPM", rotation=270, labelpad=20)  # Rotation and padding for better alignment
g.ax_heatmap.set_xlabel(f"{TF} dose", fontsize=12)

g.savefig(f"{plotdir}/{TF}_preds_cpm_heatmap_diffscale.pdf", dpi=300)
g.savefig(f"{plotdir}/{TF}_preds_cpm_heatmap_diffscale.png", dpi=600)

# %%
# generate one heatmap per sensitivity group
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np

# Create a PDF file to save the plots
with PdfPages(f'{plotdir}/{TF}_obs_cpm_heatmap_pergroup.pdf') as pdf:        
    for gr in unique_groups:
        mask = (sorted_groups == gr)
        subout = np.array(out)[mask]
        np.random.seed(16)
        g = sns.clustermap(subout, row_cluster=False, col_cluster=False, 
                        row_colors=[group_colors[mask]], 
                        cmap="viridis", figsize=(8, 8), 
                        xticklabels=[0, 0.05, 0.25, 0.5, 0.75, 1], yticklabels=False,
                        vmin = np.percentile(out, 5),
                        vmax = np.percentile(out, 95),
                        )
        legend_labels = [plt.Line2D([0], [0], marker='o', color=c, linestyle='') for c in cmap.colors]
        plt.legend(legend_labels, unique_groups, title="Sensitivity Group", loc="center left", bbox_to_anchor=(0.2, -2))
        # Add a label to the color bar
        colorbar = g.ax_heatmap.collections[0].colorbar
        colorbar.set_label("Observed CPM", rotation=270, labelpad=20)  # Rotation and padding for better alignment
        g.ax_heatmap.set_xlabel(f"{TF} dose", fontsize=12)

        pdf.savefig()
        plt.close()


# %%
# pred 
# Create a PDF file to save the plots
toplot = obs_peaks_test_preds_500bp_cpm[og_idx]

with PdfPages(f'{plotdir}/{TF}_preds_cpm_heatmap_pergroup.pdf') as pdf:        
    for gr in unique_groups:
        mask = (sorted_groups == gr)
        subout = np.array(toplot)[mask]
        np.random.seed(16)
        g = sns.clustermap(subout, row_cluster=False, col_cluster=False, 
                        row_colors=[group_colors[mask]], 
                        cmap="viridis", figsize=(8, 8), 
                        xticklabels=[0, 0.05, 0.25, 0.5, 0.75, 1], yticklabels=False,
                        vmin = np.percentile(toplot, 5),
                        vmax = np.percentile(toplot, 95),
                        )
        legend_labels = [plt.Line2D([0], [0], marker='o', color=c, linestyle='') for c in cmap.colors]
        plt.legend(legend_labels, unique_groups, title="Sensitivity Group", loc="center left", bbox_to_anchor=(0.2, -2))
        # Add a label to the color bar
        colorbar = g.ax_heatmap.collections[0].colorbar
        colorbar.set_label("Predicted CPM", rotation=270, labelpad=20)  # Rotation and padding for better alignment
        g.ax_heatmap.set_xlabel(f"{TF} dose", fontsize=12)

        pdf.savefig()
        plt.close()


# %% [markdown]
# ### same heatmap but sorted by chromhmm

# %%
from matplotlib.colors import ListedColormap
from scipy.cluster.hierarchy import linkage, dendrogram

idx = np.argsort(obs_peaks_test['chromhmm_state_broad']).values
## save ordered peaks
# obs_peaks_test.iloc[idx].to_csv(f"{outdir}/{TF}_obs_peaks_test_sorted.csv")

toplot = obs_peaks_test_obs_500bp_cpm[idx]

# Map groups to colors
# sorted_groups2 = obs_peaks_test['group'][idx]
# unique_groups2 = np.unique(sorted_groups2)

sorted_groups = obs_peaks_test['chromhmm_state_broad'][idx]
unique_groups = np.unique(list(sorted_groups))

# Perform hierarchical clustering on the subset of rows for the specific group
out = []
og_idx = []
for gr in unique_groups:
    group_data = toplot[sorted_groups == gr]
    linkage_matrix = linkage(group_data, method='ward')
    dendro = dendrogram(linkage_matrix, no_plot=True)
    row_order = dendro['leaves']
    og_idx.extend(idx[sorted_groups == gr][row_order])
    out.extend(group_data[row_order])

cmap = ListedColormap(sns.color_palette("tab20", n_colors=len(unique_groups)).as_hex())
group_colors = pd.Series(sorted_groups).map(dict(zip(unique_groups, cmap.colors)))
# cmap2 = ListedColormap(sns.color_palette(sns.color_palette(), n_colors=len(unique_groups2)).as_hex())
# group_colors2 = pd.Series(sorted_groups2).map(dict(zip(unique_groups2, cmap2.colors)))

np.random.seed(16)
g = sns.clustermap(out, row_cluster=False, col_cluster=False, 
                #row_colors=[group_colors, group_colors2], 
                row_colors=[group_colors], 
                cmap="viridis", figsize=(8, 8), 
                xticklabels=[0, 0.05, 0.25, 0.5, 0.75, 1], yticklabels=False,
                vmin = np.percentile(out, 5),
                vmax = np.percentile(out, 95),
                )
# legend_labels2 = [plt.Line2D([0], [0], marker='o', color=c, linestyle='') for c in cmap2.colors]
# plt.legend(legend_labels2, unique_groups2, title="Sensitivity Group", loc="center left", bbox_to_anchor=(4, 0.5))
legend_labels = [plt.Line2D([0], [0], marker='o', color=c, linestyle='') for c in cmap.colors]
plt.legend(legend_labels, unique_groups, title="ChromHMM State", loc="center left", bbox_to_anchor=(0.2, -2))
# Add a label to the color bar
colorbar = g.ax_heatmap.collections[0].colorbar
colorbar.set_label("Observed CPM", rotation=270, labelpad=20)  # Rotation and padding for better alignment
g.ax_heatmap.set_xlabel(f"{TF} dose", fontsize=12)

plt.tight_layout()

g.savefig(f"{plotdir}/{TF}_obs_cpm_heatmap_sortbychromhmm.pdf", dpi=300)
g.savefig(f"{plotdir}/{TF}_obs_cpm_heatmap_sortbychromhmm.png", dpi=600)

# %%
toplot = obs_peaks_test_preds_500bp_cpm[og_idx]

# using a diff color scale than obs
g = sns.clustermap(toplot, row_cluster=False, col_cluster=False, 
                row_colors=list(group_colors), cmap="viridis", figsize=(8, 8), 
                xticklabels=[0, 0.05, 0.25, 0.5, 0.75, 1], yticklabels=False,
                vmin = np.percentile(toplot, 5),
                vmax = np.percentile(toplot, 95),
                )
legend_labels = [plt.Line2D([0], [0], marker='o', color=c, linestyle='') for c in cmap.colors]
plt.legend(legend_labels, unique_groups, title="ChromHMM State", loc="center left", bbox_to_anchor=(0.2, -2))
# Add a label to the color bar
colorbar = g.ax_heatmap.collections[0].colorbar
colorbar.set_label("Predicted CPM", rotation=270, labelpad=20)  # Rotation and padding for better alignment
g.ax_heatmap.set_xlabel(f"{TF} dose", fontsize=12)

g.savefig(f"{plotdir}/{TF}_preds_cpm_heatmap_sortbychromhmm_diffscale.pdf", dpi=300)
g.savefig(f"{plotdir}/{TF}_preds_cpm_heatmap_sortbychromhmm_diffscale.png", dpi=600)

# %%
# generate one heatmap per chromhmm state
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np

# Create a PDF file to save the plots
with PdfPages(f'{plotdir}/{TF}_obs_cpm_heatmap_sortbychromhmm_perhmm.pdf') as pdf:        
    for gr in unique_groups:
        mask = (sorted_groups == gr)
        subout = np.array(out)[mask]
        np.random.seed(16)
        g = sns.clustermap(subout, row_cluster=False, col_cluster=False, 
                        row_colors=[group_colors[mask]], 
                        cmap="viridis", figsize=(8, 8), 
                        xticklabels=[0, 0.05, 0.25, 0.5, 0.75, 1], yticklabels=False,
                        vmin = np.percentile(out, 5),
                        vmax = np.percentile(out, 95),
                        )
        legend_labels = [plt.Line2D([0], [0], marker='o', color=c, linestyle='') for c in cmap.colors]
        plt.legend(legend_labels, unique_groups, title="ChromHMM State", loc="center left", bbox_to_anchor=(0.2, -2))
        # Add a label to the color bar
        colorbar = g.ax_heatmap.collections[0].colorbar
        colorbar.set_label("Observed CPM", rotation=270, labelpad=20)  # Rotation and padding for better alignment
        g.ax_heatmap.set_xlabel(f"{TF} dose", fontsize=12)

        pdf.savefig()


# %%
# pred 
# Create a PDF file to save the plots
toplot = obs_peaks_test_preds_500bp_cpm[og_idx]

with PdfPages(f'{plotdir}/{TF}_preds_cpm_heatmap_sortbychromhmm_perhmm.pdf') as pdf:        
    for gr in unique_groups:
        mask = (sorted_groups == gr)
        subout = np.array(toplot)[mask]
        np.random.seed(16)
        g = sns.clustermap(subout, row_cluster=False, col_cluster=False, 
                        row_colors=[group_colors[mask]], 
                        cmap="viridis", figsize=(8, 8), 
                        xticklabels=[0, 0.05, 0.25, 0.5, 0.75, 1], yticklabels=False,
                        vmin = np.percentile(toplot, 5),
                        vmax = np.percentile(toplot, 95),
                        )
        legend_labels = [plt.Line2D([0], [0], marker='o', color=c, linestyle='') for c in cmap.colors]
        plt.legend(legend_labels, unique_groups, title="ChromHMM State", loc="center left", bbox_to_anchor=(0.2, -2))
        # Add a label to the color bar
        colorbar = g.ax_heatmap.collections[0].colorbar
        colorbar.set_label("Predicted CPM", rotation=270, labelpad=20)  # Rotation and padding for better alignment
        g.ax_heatmap.set_xlabel(f"{TF} dose", fontsize=12)

        pdf.savefig()
        

# %% [markdown]
# ## ChromBPNet class prediction

# %%
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import label_binarize
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import classification_report, confusion_matrix
from sklearn.metrics import precision_recall_curve, auc

# %%
allpr = pd.DataFrame()
allauc = pd.DataFrame()

# %% [markdown]
# ### chrombpnet features only

# %%
suffix = "chrombpfeatonly"

X = obs_peaks_test_preds_500bp_cpm
y = obs_peaks_test.group 

# Split the dataset into training and testing sets
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=16)

# fit elastic net LR model

model = LogisticRegression(multi_class='multinomial', penalty='elasticnet', solver='saga', l1_ratio=0.5, max_iter=1000)
model.fit(X_train, y_train)

y_pred = model.predict(X_test)

# Evaluate the model
print("Confusion Matrix:")
print(confusion_matrix(y_test, y_pred))
print("\nClassification Report:")
print(classification_report(y_test, y_pred))

# Get predicted probabilities
y_scores = model.predict_proba(X_test)  # Probability estimates for the positive class

# Binarize the output for calculating precision-recall
classes = ['saturating sensitive', 'nonsaturating sensitive',
       'open nonsensitive', 'closed nonsensitive']
y_test_bin = label_binarize(y_test, classes=np.unique(y))
n_classes = y_test_bin.shape[1]

# Plot Precision-Recall curve for each class
plt.figure(figsize=(7.5, 4))

for i in range(n_classes):
    precision, recall, _ = precision_recall_curve(y_test_bin[:, i], y_scores[:, i])
    pr_auc = auc(recall, precision)  # Calculate AUC for PR curve
    allpr = pd.concat([allpr, pd.DataFrame({"feat": suffix, "group": classes[i], "precision": precision, "recall": recall})], ignore_index=True)
    allauc = pd.concat([allauc, pd.DataFrame({"feat": suffix, "group": [classes[i]], "aucpr": [pr_auc]})], ignore_index=True)

    plt.plot(recall, precision, label=f'Class {classes[i]} (AUC={pr_auc:.2f})')

plt.xlabel('Recall')
plt.ylabel('Precision')
plt.title(f'Precision-Recall Curve, {suffix}')
plt.axhline(0.5, linestyle='--', color='grey', label='Random model')
plt.legend(loc="lower left", bbox_to_anchor=[1, 0.4])
plt.grid()
plt.tight_layout()
plt.savefig(f"{plotdir}/{TF}_elasticnetLR_AUPR_{suffix}.pdf")


# %% [markdown]
# ### chromHMM features only

# %%
suffix="chromhmmfeatonly"

X = np.column_stack([obs_peaks_test.chromhmm_state_id, obs_peaks_test.chromhmm_state_id_broad])
y = obs_peaks_test.group 

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=16)

# fit elastic net LR model
model = LogisticRegression(multi_class='multinomial', penalty='elasticnet', solver='saga', l1_ratio=0.5, max_iter=1000)
model.fit(X_train, y_train)

y_pred = model.predict(X_test)

# Evaluate the model
print("Confusion Matrix:")
print(confusion_matrix(y_test, y_pred))
print("\nClassification Report:")
print(classification_report(y_test, y_pred))

# Get predicted probabilities
y_scores = model.predict_proba(X_test)  # Probability estimates for the positive class

# Binarize the output for calculating precision-recall
classes = ['saturating sensitive', 'nonsaturating sensitive',
       'open nonsensitive', 'closed nonsensitive']
y_test_bin = label_binarize(y_test, classes=np.unique(y))
n_classes = y_test_bin.shape[1]

# Plot Precision-Recall curve for each class
plt.figure(figsize=(7.5, 4))

for i in range(n_classes):
    precision, recall, _ = precision_recall_curve(y_test_bin[:, i], y_scores[:, i])
    pr_auc = auc(recall, precision)  # Calculate AUC for PR curve
    allpr = pd.concat([allpr, pd.DataFrame({"feat": suffix, "group": classes[i], "precision": precision, "recall": recall})], ignore_index=True)
    allauc = pd.concat([allauc, pd.DataFrame({"feat": suffix, "group": [classes[i]], "aucpr": [pr_auc]})], ignore_index=True)

    plt.plot(recall, precision, label=f'Class {classes[i]} (AUC={pr_auc:.2f})')

plt.xlabel('Recall')
plt.ylabel('Precision')
plt.title('Precision-Recall Curve for Multiclass Classification')
plt.axhline(0.5, linestyle='--', color='grey', label='Random model')
plt.legend(loc="lower left", bbox_to_anchor=[1, 0.4])
plt.grid()
plt.tight_layout()
plt.savefig(f"{plotdir}/{TF}_elasticnetLR_AUPR_{suffix}.pdf")


# %% [markdown]
# ### chrombpnet + chromHMM features

# %%
suffix = "chrombp+chromhmmfeat"

X = np.column_stack([obs_peaks_test_preds_500bp_cpm, obs_peaks_test.chromhmm_state_id, obs_peaks_test.chromhmm_state_id_broad])
y = obs_peaks_test.group 

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=16)

# fit elastic net LR model
model = LogisticRegression(multi_class='multinomial', penalty='elasticnet', solver='saga', l1_ratio=0.5, max_iter=1000)
model.fit(X_train, y_train)

y_pred = model.predict(X_test)

# Evaluate the model
print("Confusion Matrix:")
print(confusion_matrix(y_test, y_pred))
print("\nClassification Report:")
print(classification_report(y_test, y_pred))

# Get predicted probabilities
y_scores = model.predict_proba(X_test)  # Probability estimates for the positive class

# Binarize the output for calculating precision-recall
classes = ['saturating sensitive', 'nonsaturating sensitive',
       'open nonsensitive', 'closed nonsensitive']
y_test_bin = label_binarize(y_test, classes=np.unique(y))
n_classes = y_test_bin.shape[1]

# Plot Precision-Recall curve for each class
plt.figure(figsize=(7.5, 4))

for i in range(n_classes):
    precision, recall, _ = precision_recall_curve(y_test_bin[:, i], y_scores[:, i])
    pr_auc = auc(recall, precision)  # Calculate AUC for PR curve
    allpr = pd.concat([allpr, pd.DataFrame({"feat": suffix, "group": classes[i], "precision": precision, "recall": recall})], ignore_index=True)
    allauc = pd.concat([allauc, pd.DataFrame({"feat": suffix, "group": [classes[i]], "aucpr": [pr_auc]})], ignore_index=True)

    plt.plot(recall, precision, label=f'Class {classes[i]} (AUC={pr_auc:.2f})')

plt.xlabel('Recall')
plt.ylabel('Precision')
plt.title('Precision-Recall Curve for Multiclass Classification')
plt.axhline(0.5, linestyle='--', color='grey', label='Random model')
plt.legend(loc="lower left", bbox_to_anchor=[1, 0.4])
plt.grid()
plt.tight_layout()
plt.savefig(f"{plotdir}/{TF}_elasticnetLR_AUPR_{suffix}.pdf")


# %% [markdown]
# #### check the performance of train data

# %%
from sklearn.metrics import classification_report, confusion_matrix
y_pred = model.predict(X_train)

# Evaluate the model
print("Confusion Matrix:")
print(confusion_matrix(y_train, y_pred))
print("\nClassification Report:")
print(classification_report(y_train, y_pred))


from sklearn.metrics import precision_recall_curve, auc
# Get predicted probabilities
y_scores = model.predict_proba(X_train)  # Probability estimates for the positive class

# Binarize the output for calculating precision-recall
classes = ['saturating sensitive', 'nonsaturating sensitive',
       'open nonsensitive', 'closed nonsensitive']
y_train_bin = label_binarize(y_train, classes=np.unique(y))
n_classes = y_train_bin.shape[1]

# Plot Precision-Recall curve for each class
plt.figure(figsize=(7.5, 4))

for i in range(n_classes):
    precision, recall, _ = precision_recall_curve(y_train_bin[:, i], y_scores[:, i])
    pr_auc = auc(recall, precision)  # Calculate AUC for PR curve
    
    plt.plot(recall, precision, label=f'Class {classes[i]} (AUC={pr_auc:.2f})')

plt.xlabel('Recall')
plt.ylabel('Precision')
plt.title('Precision-Recall Curve for Multiclass Classification')
plt.axhline(0.5, linestyle='--', color='grey', label='Random model')
plt.legend(loc="lower left", bbox_to_anchor=[1, 0.4])
plt.grid()
plt.tight_layout()
plt.savefig(f"{plotdir}/{TF}_elasticnetLR_AUPR_chrombp+chromhmmfeat_train.pdf")

# %% [markdown]
# ### plot PR curves together

# %%
# Set the color palette for the groups
palette = sns.color_palette(n_colors=allpr['group'].nunique())

# Create a color dictionary mapping group to colors
color_dict = {group: palette[i] for i, group in enumerate(classes)}
alpha_dict = {"chrombpfeatonly": 1, "chrombp+chromhmmfeat": 0.6, "chromhmmfeatonly": 0.3, }

# Create the PR curve plot
plt.figure(figsize=(10, 5))

# Loop through each feature to plot the PR curves
for feat, data in allpr.groupby('feat'):
    alpha_value = alpha_dict[feat]
    for group, group_data in data.groupby('group'):
          # Example: Adjust alpha based on feat
        auc_value = allauc.loc[(allauc['feat'] == feat) & (allauc['group'] == group), 'aucpr'].values[0]
        plt.plot(group_data['recall'], group_data['precision'],
                 label=f'{feat} - {group} (AUC={auc_value:.2f})', 
                 color=color_dict[group], 
                 alpha=alpha_value)

# Formatting the plot
plt.title('Precision-Recall Curves')
plt.xlabel('Recall')
plt.ylabel('Precision')
plt.legend(loc="lower left", bbox_to_anchor=[1, 0.4])
plt.grid()
plt.tight_layout()
plt.savefig(f"{plotdir}/{TF}_elasticnetLR_AUPR_combined.pdf")

# %%
plt.figure(figsize=(5,7.5))

for idx, group in enumerate(allauc['group'].unique()):
    # Filter data for current group
    group_data_raw = allauc[allauc['group'] == group]
    
    for feat in alpha_dict.keys():
        group_data = group_data_raw[group_data_raw["feat"] == feat]
        # Plot bars for the current group
        bars = plt.bar(group_data['feat'] + f' ({group})', group_data['aucpr'], 
                    color=[color_dict[group]]*len(group_data),
                    alpha=alpha_dict[feat])
        # Add AUC text on top of bars
        for bar in bars:
            yval = bar.get_height()
            plt.text(bar.get_x() + bar.get_width()/2, yval + 0.01, round(yval, 2), 
                    ha='center', va='bottom')

# Adding labels and title
plt.title('AUC-PR')
plt.xlabel('Group')
plt.ylabel('AUC-PR')
plt.xticks(rotation=90)
plt.tight_layout()
plt.savefig(f"{plotdir}/{TF}_elasticnetLR_AUPR_barplot_combined.pdf")

# %%



