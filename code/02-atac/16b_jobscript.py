# %% [markdown]
# # Hit dose analysis

# %%
# imports
import glob
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import seaborn as sns
import sys
sys.path.insert(0, '../../software/finemo_gpu/src/finemo')
from finemo.data_io import *
from finemo.evaluation import *
from scipy.stats import entropy
from matplotlib import rcParams
import logomaker
rcParams['pdf.fonttype'] = 42  # Embed fonts (Type 3 or Type 42) otherwise can't edit text in illustrator

# %%
refdict = {"SPI1": "ETS_single__merged_pattern_0",
            # "ELF1": "ETS_single_1__merged_pattern_0",
            "ELF1": "ETS_single_2__merged_pattern_0",
            "KLF1": "SP.KLF_single_1__merged_pattern_0",
            "KLF4": "SP.KLF_single_1__merged_pattern_0",
            "ALX4": "ALX.HD__merged_pattern_0",
            "IRF4": "IRF_single__merged_pattern_0",
            "LEF1": "LEF_single__merged_pattern_0",
            "OCT4": "POU_single__merged_pattern_0",
            "SOX2": "SOX_single__merged_pattern_0",
            "SP4": "SP.KLF_single__merged_pattern_0",
            "TCF3": "bHLH_single__merged_pattern_0"}

def load_hits(path, filter=False, filter_keyword=None):
    # path: path to hits.tsv output from finemo-gpu output
    hits = pd.read_csv(path, sep="\t")
    hits["hit_name"] = hits["chr"] + "_" + hits["start"].astype(str) + "_" + hits["end"].astype(str)
    if (filter == True) and (filter_keyword != None):
        hits = hits[hits.motif_name.str.contains(filter_keyword)]
    return(hits)
# %%
# user inputs
# TF = "IRF4"
base_dir = "../../"
dose_ls = ["d005", "d025", "d050", "d075", "d100"]
outdir = f"{base_dir}/output/02-atac/16/"
plotdir = f"{base_dir}/plots/02-atac/16/"
os.makedirs(outdir, exist_ok=True)
os.makedirs(plotdir, exist_ok=True)

for TF in refdict.keys():
    print(TF)
    keyword = refdict[TF]
    hits_paths = {d: glob.glob(f"{base_dir}/output/04-chrombpnet/output/models/fold_0/*{TF}_{d}/finemo_out/{TF}/hits/hits.tsv")[0] for d in dose_ls}
    hits_paths

    # %%
    # check to make sure all hits paths exist
    assert all([os.path.exists(f) for f in list(hits_paths.values())])

    hits_file = f"{outdir}/{TF}_{keyword}_hit_dose.csv"
    
    if not os.path.exists(hits_file):
        # %% [markdown]
        # ## read motif hits

        # %%
        hits_all = load_hits(hits_paths["d100"], filter=True, filter_keyword=keyword)
        hits_all["hit_d100"] = 1
        print(hits_all.shape)
        hits_all.head()

        # %% [markdown]
        # ## group motifs by when they get called by finemo as hits
        # start with all hits called at highest dose, then annotate if these are also hits in each of the lower doses.

        # %%
        for dose in ["d075", "d050", "d025", "d005"]:
            print(dose)
            hits = load_hits(hits_paths[dose], filter=True, filter_keyword=keyword)
            hits[f"hit_{dose}"] = 1
            hits_all = pd.merge(hits_all, hits, on=["chr", "start", "end", "start_untrimmed", "end_untrimmed", "motif_name", "strand", "peak_name",  "peak_id", "hit_name"], 
                                how="left", suffixes=["", f"_{dose}"])


        # %%
        hits_all["hit_dose"] = "d100"
        hits_all.loc[hits_all.hit_d075==1, "hit_dose"] = "d075"
        hits_all.loc[(hits_all.hit_dose=="d075") & (hits_all.hit_d050==1), "hit_dose"] = "d050"
        hits_all.loc[(hits_all.hit_dose=="d050") & (hits_all.hit_d025==1), "hit_dose"] = "d025"
        hits_all.loc[(hits_all.hit_dose=="d025") & (hits_all.hit_d005==1), "hit_dose"] = "d005"

        # %%
        hits_all.hit_dose.value_counts()

        # %%
        hits_all.to_csv(hits_file)

    else: 
        # %%
        hits_all = pd.read_csv(hits_file)
        hits_all.head()

    # %% [markdown]
    # ### generate an upset plot

    # %%
    # outer merge instead of left merge
    hits_all = load_hits(hits_paths["d100"], filter=True, filter_keyword=keyword)
    hits_all["hit_d100"] = True
    for dose in ["d075", "d050", "d025", "d005"]:
        print(dose)
        hits = load_hits(hits_paths[dose], filter=True, filter_keyword=keyword)
        hits[f"hit_{dose}"] = True
        hits_all = pd.merge(hits_all, hits, on=["chr", "start", "end", "hit_name"], 
                            how="outer", suffixes=["", f"_{dose}"])
    hits_all.fillna(False, inplace=True)

    # %%
    hits_all.head()

    # %%
    hits_all.set_index(["hit_d005", "hit_d025", "hit_d050", "hit_d075", "hit_d100"], inplace=True)
    hits_all.head()

    # %%
    from upsetplot import UpSet, plot
    UpSet(hits_all, subset_size="count",sort_by="cardinality").plot()
    plt.savefig(f"{plotdir}/upset_{TF}_{keyword}_hit_dose.pdf")


    # %% [markdown]
    # ## plot pileups for each hit dose group

    # %%
    hits = pd.read_csv(hits_file)
    hits["is_revcomp"] = (hits.strand=="-")
    hits.head()

    # %%
    genome = f"{base_dir}/data/chrombp_resources/hg38.fa"
    hg38 = pyfaidx.Fasta(genome)

    def get_seq_polar(hits, genome, pad=0):
        """
        fetches sequence from a given genome.
        """
        vals = []
        for i in range(hits.shape[0]):
            r = hits[i]
            sequence_raw = genome[r['chr'][0]][(r['start'][0]-pad):(r['end'][0]+pad)]
            if r['is_revcomp'][0]:
                sequence = str(sequence_raw.reverse.complement.seq)
            else:
                sequence = str(sequence_raw)
            vals.append(sequence)
        return vals

    def get_seq(hits, genome, pad=0):
        """
        fetches sequence from a given genome.
        """
        vals = []
        for i in range(hits.shape[0]):
            r = hits.iloc[i]
            sequence_raw = genome[r['chr']][(r['start']-pad):(r['end']+pad)]
            if r['is_revcomp']:
                sequence = str(sequence_raw.reverse.complement.seq)
            else:
                sequence = str(sequence_raw)
            vals.append(sequence)
        return vals

    # %%
    from matplotlib.colors import ListedColormap

    # pad=0 # number of bases around each side of motif to include in pileup
    for pad in [0,5]:
        for m in set(hits["hit_dose"]):
            plot = f"{plotdir}/{TF}_{keyword}_cwms/{m}"
            os.makedirs(plot, exist_ok=True)
            subhits = hits[hits.hit_dose==m]
            # sample 1000 motifs for speed
            np.random.seed(16)
            idx = np.random.randint(0, len(subhits)-1, size=1000)
            subhits = subhits.iloc[idx]
            sequences = get_seq(subhits, hg38, pad=pad)

            # Define colors for nucleotides
            nucleotide_colors = {
                "A": "green",
                "T": "red",
                "G": "orange",
                "C": "blue",
                "N": "gray",
            }

            # Create a dynamic colormap for only the nucleotides present in the sequences
            unique_nucleotides = sorted(set("".join(sequences)))
            present_mapping = {nuc: idx for idx, nuc in enumerate(unique_nucleotides)}
            present_colors = [nucleotide_colors[nuc] for nuc in unique_nucleotides]
            colormap = ListedColormap(present_colors)

            # Convert sequences to a numeric matrix
            numeric_matrix = np.array([[present_mapping[nuc] for nuc in seq] for seq in sequences])

            # cluster the sequences
            from scipy.cluster.hierarchy import linkage, dendrogram
            linked = linkage(numeric_matrix, method="average")
            dend = dendrogram(linked, orientation="left", no_plot=True)
            # plt.gca().invert_yaxis()
            # plt.xticks([])
            # plt.yticks([])
            seq_order = dend["leaves"]

            # Plot the heatmap
            plt.figure(figsize=(8, 4))
            sns.heatmap(
                numeric_matrix[seq_order],
                annot=False,  # Set to True to annotate with letters
                cmap=colormap,  # Use custom colors
                cbar=False,  # Hide the color bar
            )

            # plt.yticks(ticks=range(len(seq_order)),labels=seq_order)
            plt.yticks([])
            plt.xticks([])
            plt.xlabel("Position")
            plt.ylabel("Sequence")
            plt.title(f"DNA Sequence Pileup {TF} Hit Dose {m}")
            plt.tight_layout()
            plt.savefig(f"{plot}/pileup_pad{pad}.pdf")
            plt.show()


    # %% [markdown]
    # ## plot PWMs for each hit dose group  

    # %%
    from Bio import motifs
    from Bio.Seq import Seq
    def dna_sequences_to_pwm(sequences):
        """
        Converts a list of DNA sequences to a PWM.

        Args:
            sequences (list): A list of DNA sequences (strings).

        Returns:
            Bio.motifs.matrix.PositionWeightMatrix: The PWM.
        """
        # Create a Motif object from the sequences
        motif = motifs.create([Seq(s) for s in sequences])
        # Calculate the PWM
        pwm = motif.counts.normalize(pseudocounts=0.5) # Adding pseudocounts to avoid zero probabilities

        return pwm

    # %%
    pad = 5
    for m in set(hits["hit_dose"]):
        subhits = hits[hits.hit_dose==m]
        # sample 5000 motifs for speed
        np.random.seed(16)
        idx = np.random.randint(0, len(subhits)-1, size=5000)
        subhits = subhits.iloc[idx]
        sequences = get_seq(subhits, hg38, pad=pad)

        pwm = dna_sequences_to_pwm(sequences)
        print(pwm)

        # Create a Pandas DataFrame from the PWM data
        pwm_df = pd.DataFrame(pwm)

        logo = logomaker.Logo(pwm_df)
        logo.style_glyphs()
        logo.style_spines(visible=False)
        logo.style_spines(spines=['left', 'bottom'], visible=True)
        logo.ax.set_xticks(range(len(pwm_df)))
        logo.ax.set_xticklabels(range(1, len(pwm_df) + 1))
        logo.ax.set_xlabel('Position')
        logo.ax.set_ylabel('Information Content (bits)')
        logo.fig.show()
            
        plt.savefig(f"{plotdir}/{TF}_{keyword}_cwms/{m}/pwm_pad{pad}_raw.pdf", dpi=150)

        # normalize for information content at each position
        pwm_df_norm = (pwm_df.transpose() * (2-entropy(pwm_df.transpose(), base=2))).transpose()
        pwm_df_norm.to_csv(f"{outdir}/{TF}_{keyword}_{m}_pwm_pad{pad}_norm.csv")

        logo = logomaker.Logo(pwm_df_norm)
        logo.style_glyphs()
        logo.style_spines(visible=False)
        logo.style_spines(spines=['left', 'bottom'], visible=True)
        logo.ax.set_xticks(range(len(pwm_df)))
        logo.ax.set_xticklabels(range(1, len(pwm_df) + 1))
        logo.ax.set_ylim(0, 2)
        logo.ax.set_xlabel('Position')
        logo.ax.set_ylabel('Information Content (bits)')
        logo.fig.show()
        plt.savefig(f"{plotdir}/{TF}_{keyword}_cwms/{m}/pwm_pad{pad}_norm.pdf", dpi=150)

    # %%


    # %% [markdown]
    # ## plot CWMs for each hit dose group

    # %% [markdown]
    # reading in hits again in the polar format to use finemo cwm functions

    # %%
    # regions_path = glob.glob(f"../../output/04-chrombpnet/output/models/fold_0/HEK293T_*{TF}_d100/finemo_out/intermediate_inputs.npz")[0]
    # peaks_path = glob.glob(f"../../output/04-chrombpnet/output/models/fold_0/HEK293T_*{TF}_d100/finemo_out/{TF}/hits/peaks_qc.tsv")[0]

    # # takes 5m
    # motif_width = 30
    # sequences, contribs = load_regions_npz(regions_path)
    # peaks = pl.scan_csv(peaks_path, has_header=True, separator="\t").collect()

    # if len(contribs.shape) == 3:
    #     regions = contribs * sequences
    # elif len(contribs.shape) == 2:
    #     regions = contribs[:,None,:] * sequences

    # %%
    # hits_path = hits_file
    # HITS_DTYPES = {
    #     "chr": pl.Utf8,
    #     "start": pl.UInt32,
    #     "end": pl.UInt32,
    #     "start_untrimmed": pl.UInt32,
    #     "end_untrimmed": pl.UInt32,
    #     "motif_name": pl.Utf8,
    #     "hit_coefficient": pl.Float32,
    #     "hit_coefficient_global": pl.Float32,
    #     "hit_correlation": pl.Float32,
    #     "hit_importance": pl.Float32,
    #     "strand": pl.Utf8,
    #     "peak_name": pl.Utf8,
    #     "peak_id": pl.UInt32,
    # }
    # hits_df = (
    #         pl.scan_csv(hits_path, separator=',', quote_char=None, dtypes=HITS_DTYPES)
    #         .with_columns(pl.lit(1).alias("count"))
    #     )

    # %%
    # hits = hits_df.collect()
    # hits = hits.join(peaks, on="peak_name").with_columns(is_revcomp=(pl.col("strand") == '-'))
    # hits.head()

    # %%
    # cwms = {}
    # for m in set(hits["hit_dose"]):
    #     cwms[m] = {"hits_fc": get_cwms(regions, hits.filter(pl.col("hit_dose")==m), motif_width)}
    # cwms["all_d100_hits"] = {"hits_fc": get_cwms(regions, hits, motif_width)}
    # cwms

    # %%
    # def plot_cwms_bl(cwms, out_dir, alphabet=LOGO_ALPHABET, colors=LOGO_COLORS, font=LOGO_FONT):
    #     for m, v in cwms.items():
    #         motif_dir = os.path.join(out_dir, m)
    #         os.makedirs(motif_dir, exist_ok=True)
    #         for cwm_type, cwm in v.items():
    #             output_path = os.path.join(motif_dir, f"{cwm_type}.png")

    #             fig, ax = plt.subplots(figsize=(10,2))

    #             plot_logo(ax, cwm, alphabet, colors=colors, font_props=font)

    #             for name, spine in ax.spines.items():
    #                 spine.set_visible(False)
                
    #             plt.savefig(output_path, dpi=100)
    #             plt.close(fig)

    # plot_cwms_bl(cwms, out_dir=f"{plotdir}/{TF}_{keyword}_cwms")


