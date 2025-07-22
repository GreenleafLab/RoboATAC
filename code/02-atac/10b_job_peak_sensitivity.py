## Peak sensitivity analysis (ATAC)

# imports 
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import seaborn as sns
import pandas as pd
import itertools as it
import argparse
from sklearn.decomposition import PCA
from tqdm import tqdm
from scipy.odr import odrpack
import os
import glob

def init_grid(grid_dimensions: tuple, fig_dimensions: tuple, caption_height: int = 0, **kwargs) -> tuple:
    """ Initializes a grid of figures of correct dimension.

    To use this, please call plt.subplot(*grid_dimensions, next(counter)) before each figure to
    include in the grid

    Args:
        grid_dimensions: The dimension of grid (num rows, num columns)
        fig_dimensions: The dimension of each figure in the grid (width, height)
        caption_height: The height of the caption
        **kwargs: Additional arguments to pass to plt.figure

    Returns:
        tuple of (grid_dimensions, fig_dimensions, counter)

    """
    counter = it.count(1)
    plt.figure(figsize=(grid_dimensions[1] * fig_dimensions[0],
                        grid_dimensions[0] * fig_dimensions[1] + caption_height), **kwargs)
    return grid_dimensions, fig_dimensions, counter

def plot_peaks(peakset,pdf,t,fit,func=None, showparam=False, titletexts=None, cpm=False):
    # func is only needed if fit=="hill"
    g, d, c = init_grid((len(peakset)//5+1,5), (7,7))
    for peak in peakset:
        plt.subplot(*g, next(c))
        df = pdf[pdf["index"] == peak]
        x = np.linspace(0, 1)
        offset = df.loc[df.dose==0, "normalized peak count"].mean()
        scale = df.loc[df.dose==1, "normalized peak count"].mean() - offset

        titletext = peak
        if titletexts is not None:
            titletext = titletexts[peak]
        matplotlib.rc("font",size=20)
        sns.scatterplot(data = df, x = "dose", y = "normalized peak count", s=150).set(title=titletext)

        # then choose one of the three below to plot a line
        if fit == "nofit":
            sns.lineplot(data = df, x = "dose", y = "normalized peak count")
        elif fit == "linear":
            plt.plot(x, [np.polyval([t.loc[peak,"m"], t.loc[peak,"b"]], i) * scale + offset for i in x], linewidth=5)
            if showparam:
                plt.title(titletext + "\n" + f"m={round(t.loc[peak, 'm'],3)}, b={round(t.loc[peak, 'b'],3)}")
        elif fit == "hill":
            plt.plot(x, [func(i, t.loc[peak, "h"], t.loc[peak, "ka"]) * scale + offset for i in x], linewidth=5)
            if showparam:
                plt.title(titletext + "\n" + f"h={round(t.loc[peak, 'h'],3)} (p={round(t.loc[peak,'h_pval'],3)}), " +  
                          f"ka={round(t.loc[peak, 'ka'],3)} (p={round(t.loc[peak,'ka_pval'],3)})")
        elif fit == "hillamp":
             plt.plot(x, [func(i, t.loc[peak, "h"], t.loc[peak, "ka"], t.loc[peak, "amp"]) + offset for i in x])
             if showparam:
                plt.title(titletext + "\n" + f"h={round(t.loc[peak, 'h'],3)} (p={round(t.loc[peak,'h_pval'],3)}), " +  
                          f"ka={round(t.loc[peak, 'ka'],3)} (p={round(t.loc[peak,'ka_pval'],3)})" + "\n" + "amp={round(t.loc[peak, 'ka'], 3)}" )
        
        # ylab
        if cpm:
            label="ATAC CPM"
        else:
            label="Scaled ATAC CPM (0-1)"

        plt.tight_layout()
        plt.ylabel(label)
        plt.xlabel("TF Plasmid Dose ng/uL")
        
def parse_inputs():
    parser = argparse.ArgumentParser(description="Parse arguments.")
    # Define the positional arguments
    parser.add_argument("--TF", default=[], nargs="*", type=str, help="List of TFs to run")

    # Parse the arguments
    args = parser.parse_args()

    # Return the parsed arguments
    return args.TF


##############################
# inputs
##############################
qnorm_rpkm_df = pd.read_csv("../../output/02-atac/01/cpm.tsv", sep=" ") # just use the cpm normalized peak count matrix

group1 = "GFP_1" # reference level
region = "filteredConsensus"

TF_ls = parse_inputs()

plotdir = "../../plots/02-atac/10/"
outdir = "../../output/02-atac/10/"
os.makedirs(plotdir, exist_ok=True)
os.makedirs(outdir, exist_ok=True)

##############################
# read meta data
##############################
metadata = pd.read_csv("atac_sample_sheet.txt", sep="\t")
metadata.loc[metadata.TF=="GFP","Plasmid_Dose"] = 0 # set control plasmid dose to zero
dosage_dic = {metadata.loc[i, "sampleName"]: metadata.loc[i, "Plasmid_Dose"] 
            if metadata.loc[i, "TF"]!="GFP" else 0 for i in metadata.index}
metadata.index = metadata.sampleName
    
# loop through TFs
#for TF in ["SPI1", "KLF1", "KLF4", "SOX2", "OCT4", "ALX4", "IRF4", "SP4", "TCF3"]:
#for TF in ["ALX4", "IRF4", "SP4", "TCF3", "OCT4"]:
# for TF in ["IRF4", "SP4", "TCF3", "OCT4", "ALX4"]:
for TF in TF_ls:
    print(TF)
    group2 = TF+"_1"

    ##############################
    # read differential peaks
    ##############################
    de = pd.read_csv(glob.glob(f"../../output/02-atac/04/deseq_data_Plate_TF_Dose/deseq_*{group2}_vs_*{group1}_{region}_full.tsv")[0], sep = "\t")
    peaks = set(de.loc[(de["padj"] < 0.05) & (de["log2FoldChange"] > 0.58), "name"]) # only look at peaks that increase in accessibility
    de.index = de.name
    de["peakloc"] = de.seqnames + ": " + de.start.astype(str) + "-" + de.end.astype(str)

    ##############################
    # read cpm data
    ##############################
    subset = qnorm_rpkm_df.loc[peaks, qnorm_rpkm_df.columns.str.contains(TF) | qnorm_rpkm_df.columns.str.contains("GFP")]
    subset.columns = subset.columns.str.lstrip('X')
    subset[pd.isna(subset)]=0


    ##############################
    # plot a heatmap of cpm values of DE peaks
    ##############################
    col_order = metadata.loc[subset.columns,].sort_values(by="Plasmid_Dose").index.tolist()
    # g = sns.clustermap(subset[col_order], col_cluster = False, z_score=0, figsize=(14,14), rasterized=True)
    # _ = plt.setp(g.ax_heatmap.xaxis.get_majorticklabels(), fontsize=8)
    # _ = plt.setp(g.ax_heatmap.yaxis.get_majorticklabels(), fontsize=8)
    # os.system("mkdir -p plots/sensitivity")
    # for ax in g.ax_row_dendrogram, g.ax_col_dendrogram:
    #     ax.set_rasterized(True)
    # g.savefig("plots/sensitivity/heatmap_cpm_de_peaks_" + TF + ".pdf", dpi=100)

    ##############################
    # flatten cpm table
    ##############################
    r_save = subset[col_order].rename(columns = dosage_dic)
    pdf_raw = r_save.reset_index().melt(id_vars = "index").rename(columns={"variable": "dose", "value": "normalized peak count"})
    r_save = r_save.subtract(r_save[0].mean(axis = 1), axis =0) # control counts at 0
    r_save = r_save.divide(r_save[1].mean(axis = 1), axis = 0) # highest dosage counts at 1
    pdf = r_save.reset_index().melt(id_vars = "index").rename(columns={"variable": "dose", "value": "normalized peak count"})

    pdf.to_csv(f"{outdir}/pdf_{TF}_log2fc0.58.csv")
    pdf_raw.to_csv(f"{outdir}/pdf_raw_{TF}_log2fc0.58.csv")

    ##############################
    # linear fits
    ##############################
    print("starting linear fit")
    linear_fits = {}
    for peak, df in pdf.groupby("index"):
        out, residuals, _, _, _ = np.polyfit(df["dose"].astype(float).values, 
                                                df["normalized peak count"].values, 1, full = True)
        linear_fits[peak] = {"m": out[0], "b": out[1], "res_linear": residuals[0]}
    linear_fit_df = pd.DataFrame(linear_fits).T
    print("DONE linear fit")

    ##############################
    # hill fits
    ##############################
    from scipy.optimize import curve_fit
    import scipy

    def func(x, a, b):
        return (x**a/(b**a + x**a))

    def f_wrapper_for_odr(beta, x): # parameter order for odr
        return func(x, *beta)

    def calc_square_error(x, y, a, b):
        error = 0
        for i in range(len(x)): 
            xi = x[i]
            yi = y[i]
            error += (func(xi, a, b) - yi)**2
        return error

    print("starting hill fit")
    hill_fits = {}

    for peak, df in tqdm(pdf.groupby("index")):
        df["shifted value"] = df["normalized peak count"] - df.loc[df["dose"] == 0, "normalized peak count"].mean()
        xdata = df["dose"].astype(float).values
        ydata = df["shifted value"].values

        # Find guess 
        ka_guess = .025 if df.loc[df["dose"] == .05, "shifted value"].mean() > 0.5 else .3
        
        try:
            popt, pcov = curve_fit(func, xdata, ydata, p0=[ 1, ka_guess], bounds = ([0.01, 0], 
                                                                                        [10, 1]))
            a, b = popt[0], popt[1]
            # determine the significance of the fit
            model = odrpack.Model(f_wrapper_for_odr)
            data = odrpack.Data(xdata,ydata)
            myodr = odrpack.ODR(data, model, beta0=popt,  maxit=0)
            myodr.set_job(fit_type=2)

            parameterStatistics = myodr.run()
            # df_e = len(xdata) - len(popt) # degrees of freedom, error
            #cov_beta = parameterStatistics.cov_beta # parameter covariance matrix from ODR
            #sd_beta = parameterStatistics.sd_beta * parameterStatistics.sd_beta
            # ci = []
            # t_df = scipy.stats.t.ppf(0.975, df_e)
            # ci = []
            # for i in range(len(popt)):
            #     ci.append([popt[i] - t_df * parameterStatistics.sd_beta[i], popt[i] + t_df * parameterStatistics.sd_beta[i]])

            tstat_beta = popt / [sd if sd>0 else 1e-3 for sd in parameterStatistics.sd_beta] # coeff t-statistics
            pstat_beta = (1.0 - scipy.stats.t.cdf(np.abs(tstat_beta), df_e)) * 2.0   # coef. p-values

            hill_fits[peak] = {"h": a, "ka": b, "res_hill": calc_square_error(xdata, ydata, a, b), "h_pval": pstat_beta[0], "ka_pval": pstat_beta[1]}
        except:
            print(peak)
            continue

    hill_fits_df = pd.DataFrame(hill_fits).T
    print("DONE hill fit")

    ##############################
    # concatenate fit results
    ##############################
    t = pd.concat([linear_fit_df, hill_fits_df], axis = 1)
    t = t.dropna(axis=0)
    t.to_csv(f"{outdir}/sensitivity_{TF}_log2fc0.58_fit.csv")