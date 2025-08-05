import numpy as np
import pandas as pd
from scipy.stats import f_oneway
from ENT3C.core import utils
from itertools import combinations
import matplotlib
import matplotlib.colors as mcolors
from matplotlib.backends.backend_pdf import PdfPages

matplotlib.use("Qt5Agg")
import matplotlib.pyplot as plt


def locate_largest_euclidean_diff(config_file, group1, group2):
    (
        config_df,
        SUB_M_SIZE_FIX,
        PHI_MAX,
        CHRSPLIT,
        phi,
        NormM,
        CHROMOSOMES,
        RESOLUTIONS,
        BR,
        FNs,
        OUT_DIR,
        OUT_PREFIX,
        entropy_out_FN,
        similarity_out_FN,
        LOG_FN,
    ) = utils.check_config(config_file)

    ENT3C_OUT = pd.read_csv(f"{entropy_out_FN}", sep="\t", dtype={"ChrNr": str})
    ENT3C_OUT["ChrNr"] = ENT3C_OUT["ChrNr"].astype(str)
    ENT3C_OUT["GROUP"] = ENT3C_OUT["Name"].str.extract(r"(.*?)_BR")[0]

    n_group1 = ENT3C_OUT[ENT3C_OUT["GROUP"] == group1]["Name"].unique().shape[0]
    n_group2 = ENT3C_OUT[ENT3C_OUT["GROUP"] == group2]["Name"].unique().shape[0]
    cmap1 = plt.cm.get_cmap("Greys", n_group1).reversed()
    cmap1 = [mcolors.rgb2hex(cmap1(j)) for j in np.linspace(0, 0.5, n_group1)]
    cmap2 = plt.cm.get_cmap("Oranges", n_group2).reversed()
    cmap2 = [mcolors.rgb2hex(cmap2(j)) for j in np.linspace(0, 0.5, n_group2)]
    # print(list(zip(ENT3C_OUT[ENT3C_OUT["GROUP"] == group1]["Name"].unique(), cmap1)))

    colormap = {}
    for name, color in zip(
        ENT3C_OUT[ENT3C_OUT["GROUP"] == group1]["Name"].unique(), cmap1
    ):
        colormap[name] = color
    for name, color in zip(
        ENT3C_OUT[ENT3C_OUT["GROUP"] == group2]["Name"].unique(), cmap2
    ):
        colormap[name] = color

    ENT3C_OUT["Color"] = ENT3C_OUT["Name"].map(colormap)

    for Resolution in RESOLUTIONS:
        with PdfPages(
            f"{OUT_DIR}/{OUT_PREFIX}_distances_{Resolution / 1e3}kb.pdf",
        ) as pdf:
            for i, ChrNr in enumerate(CHROMOSOMES):
                fig, axs = plt.subplots(
                    nrows=2, ncols=1, figsize=(8, 4), constrained_layout=False
                )

                fig.suptitle(
                    f"{Resolution / 1e3}kb Largest Differences across {group1} and {group2} Chr{ChrNr}"
                )

                upper_ax, lower_ax = axs

                DF = ENT3C_OUT[
                    (ENT3C_OUT["ChrNr"] == ChrNr)
                    & (ENT3C_OUT["Resolution"] == Resolution)
                ]
                # print(DF)
                # print(DF["Name"].unique())

                # ax = axs[0]
                # for name, group_df in DF.groupby("Name"):
                #    x = range(0, len(group_df["S"]))
                #    ax.plot(x, group_df["S"], label=name)
                #### zscores
                DF.loc[:, "S"] = DF.groupby("Name")["S"].transform(
                    lambda S: (S - S.mean()) / S.std()
                )
                for name, group_df in DF.groupby("Name"):
                    x = range(0, len(group_df["S"]))
                    COLOR = group_df["Color"].unique()[0]
                    upper_ax.plot(
                        x, group_df["S"], label=name, color=COLOR, linewidth=0.3
                    )

                G1 = DF[DF["GROUP"] == group1][["S", "Name", "START"]].reset_index(
                    drop=True
                )
                G1 = G1.pivot(index="START", columns="Name", values="S")

                G2 = DF[DF["GROUP"] == group2][["S", "Name", "START"]].reset_index(
                    drop=True
                )
                G2 = G2.pivot(index="START", columns="Name", values="S")

                G1 = G1.to_numpy().T  # (n1, n_bins)
                G2 = G2.to_numpy().T  # (n2, n_bins)

                G1 = G1[:, None, :]  # (n1,1,n_bins)
                G2 = G2[None, :, :]  # (1,n2,n_bins)
                DIFF = np.abs(G1 - G2)
                # print(DIFF.shape)  # (n1, n2, n_bins)

                DIFFERENCES = DF.drop_duplicates("binNrStart")[["START", "END"]].copy()
                DIFFERENCES["meanS_Euclidean"] = np.mean(DIFF, axis=(0, 1))

                G1 = np.squeeze(G1)
                G1 = np.mean(G1, axis=0)
                G2 = np.squeeze(G2)
                G2 = np.mean(G2, axis=0)

                DIFF2 = abs(G1 - G2)
                DIFFERENCES["mean_EuclideanD"] = DIFF2

                x_positions = DIFFERENCES.index
                x_labels = DIFFERENCES["START"].astype(str)

                lower_ax.plot(
                    x_positions,
                    DIFFERENCES["mean_EuclideanD"],
                    # DIFFERENCES["meanS_Euclidean"],
                    label=name,
                    color="black",
                    linewidth=1,
                )

                tick_positions = x_positions[::16]
                tick_labels = x_labels[::16]

                lower_ax.set_xticks(tick_positions)
                lower_ax.set_xticklabels(tick_labels, rotation=45, fontsize=6)

                upper_ax.set_xticks(tick_positions)
                upper_ax.set_xticklabels(tick_labels, rotation=45, fontsize=6)

                upper_ax.autoscale(enable=True, axis="both", tight=True)
                lower_ax.autoscale(enable=True, axis="both", tight=True)

                lower_ax.xaxis.label.set_fontsize(5)
                lower_ax.tick_params(axis="y", which="major", labelsize=8)
                lower_ax.set_title(
                    # f"Euclidean distance between average S of {group1} and {group2}",
                    f"Average of pairwise Euclidean distance between S of {group1} and {group2}",
                    fontsize=9,
                )
                lower_ax.set_ylim(0, 1)

                upper_ax.xaxis.label.set_fontsize(5)
                upper_ax.tick_params(axis="y", which="major", labelsize=8)
                upper_ax.set_title("Zscores ENT3C", fontsize=9)

                lower_ax.grid(True)
                upper_ax.grid(True)

                handles, labels = upper_ax.get_legend_handles_labels()
                labels, handles = zip(
                    *sorted(zip(labels, handles), key=lambda t: t[0].lower())
                )
                # legend does not fit
                fig.subplots_adjust(
                    left=0.1, right=0.85, top=0.8, bottom=0.2, hspace=0.8
                )

                leg = upper_ax.legend(
                    handles,
                    labels,
                    loc="upper right",
                    bbox_to_anchor=(1, 0.9),
                    borderaxespad=0.0,
                    handlelength=1,  # legend line thickness
                    fontsize="x-small",  # smaller legend text
                    frameon=False,
                    bbox_transform=fig.transFigure,
                )

                for line in leg.get_lines():
                    line.set_linewidth(1)

                pdf.savefig(fig)
                plt.close(fig)

                DIFFERENCES = DIFFERENCES.sort_values(
                    "mean_EuclideanD", ascending=False
                ).reset_index(drop=True)

    return DIFFERENCES
