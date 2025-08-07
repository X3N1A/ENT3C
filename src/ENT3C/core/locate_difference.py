import numpy as np
import pandas as pd
from ENT3C.core import utils
import matplotlib.colors as mcolors
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import sys


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
    print(BR)
    if BR:
        ENT3C_OUT["GROUP"] = ENT3C_OUT["Name"].str.extract(r"(.*)_BR")[0]
    else:
        ENT3C_OUT["GROUP"] = ENT3C_OUT["Name"]

    print(ENT3C_OUT)

    if ENT3C_OUT[(ENT3C_OUT["GROUP"] == group1) | (ENT3C_OUT["GROUP"] == group2)].empty:
        sys.exit(
            f"Error: {entropy_out_FN} does not contain {group1} and {group2} samples."
        )

    n_group1 = ENT3C_OUT[ENT3C_OUT["GROUP"] == group1]["Name"].unique().shape[0]
    n_group2 = ENT3C_OUT[ENT3C_OUT["GROUP"] == group2]["Name"].unique().shape[0]

    cmap = mcolors.LinearSegmentedColormap.from_list(
        "black_to_grey", ["black", "grey"], n_group1
    )
    cmap1 = [mcolors.to_hex(cmap(i)) for i in np.linspace(0, 1, n_group1)]

    cmap = mcolors.LinearSegmentedColormap.from_list(
        "red_to_magenta", ["darkred", "red", "magenta"], n_group2
    )
    # cmap = plt.cm.get_cmap("rainbow_r", n_group2).reversed()
    cmap2 = [mcolors.to_hex(cmap(i)) for i in np.linspace(0, 1, n_group2)]

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
    EUCLIDEAN = pd.DataFrame()

    eucl_out_FN = f"{OUT_DIR}/{OUT_PREFIX}_Eucl_{group1}vs{group2}.csv"

    for Resolution in RESOLUTIONS:
        with PdfPages(
            f"{OUT_DIR}/{OUT_PREFIX}_Eucl_{int(Resolution / 1e3)}kb_{group1}vs{group2}.pdf",
        ) as pdf:
            for i, ChrNr in enumerate(CHROMOSOMES):
                fig, axs = plt.subplots(
                    nrows=3, ncols=1, figsize=(8, 4), constrained_layout=False
                )

                fig.suptitle(
                    f"{Resolution / 1e3}kb Largest Differences across {group1} and {group2} Chr{ChrNr}",
                    fontsize=9,
                )

                upper_ax, mid_ax, lower_ax = axs
                DF = ENT3C_OUT[
                    (ENT3C_OUT["ChrNr"] == ChrNr)
                    & (ENT3C_OUT["Resolution"] == Resolution)
                    & ((ENT3C_OUT["GROUP"] == group1) | (ENT3C_OUT["GROUP"] == group2))
                ]
                DF = DF.copy()
                #### centering: mean = 0
                # DF.loc[:, "S_centered"] = DF.groupby("START")["S"].transform(
                #    lambda S: (S - S.mean())
                # )
                #### zscores: mean = 0 and std=1
                DF.loc[:, "S_zscore"] = DF.groupby("Name")["S"].transform(
                    lambda S: (S - S.mean()) / S.std()
                )

                for name, group_df in DF.groupby("Name"):
                    COLOR = group_df["Color"].unique()[0]
                    upper_ax.plot(
                        # group_df["START"],
                        group_df["START"],
                        group_df["S"],
                        label=name,
                        color=COLOR,
                        linewidth=0.3,
                        alpha=0.9,
                    )
                    mid_ax.plot(
                        group_df["START"],
                        group_df["S_zscore"],
                        label=name,
                        color=COLOR,
                        linewidth=0.3,
                        alpha=0.4,
                    )

                G1 = DF[DF["GROUP"] == group1][
                    ["S_zscore", "Name", "START"]
                ].reset_index(drop=True)
                G1 = G1.pivot(index="START", columns="Name", values="S_zscore")

                G2 = DF[DF["GROUP"] == group2][
                    ["S_zscore", "Name", "START"]
                ].reset_index(drop=True)
                G2 = G2.pivot(index="START", columns="Name", values="S_zscore")

                G1 = G1.to_numpy().T  # (n1, n_bins)
                G2 = G2.to_numpy().T  # (n2, n_bins)

                G1 = G1[:, None, :]  # (n1,1,n_bins)
                G2 = G2[None, :, :]  # (1,n2,n_bins)
                DIFF = np.abs(G1 - G2)
                # print(DIFF.shape)  # (n1, n2, n_bins)

                DIFFERENCES = DF.drop_duplicates("binNrStart")[
                    ["Resolution", "ChrNr", "START", "END"]
                ].copy()
                DIFFERENCES["meanS_Euclidean"] = np.mean(DIFF, axis=(0, 1))

                G1 = np.squeeze(G1)
                G1 = np.mean(G1, axis=0)
                G2 = np.squeeze(G2)
                G2 = np.mean(G2, axis=0)

                DIFF2 = abs(G1 - G2)

                x_positions = DIFFERENCES["START"]
                x_labels = (
                    DIFFERENCES["START"].astype(str)
                    + "-"
                    + DIFFERENCES["END"].astype(str)
                )
                lower_ax.plot(
                    DIFFERENCES["START"],
                    DIFFERENCES["meanS_Euclidean"],
                    color="black",
                    linewidth=1,
                )
                tick_positions = x_positions[::16]
                tick_labels = x_labels[::16]

                for ax in [upper_ax, mid_ax, lower_ax]:
                    ax.autoscale(enable=True, axis="both", tight=True)
                    ax.set_xticks(tick_positions)
                    ax.xaxis.label.set_fontsize(3)
                    ax.tick_params(axis="y", which="major", labelsize=8)
                    ax.grid(True)

                for ax in [upper_ax, mid_ax]:
                    ax.tick_params(axis="x", labelbottom=True)
                    ax.set_xticklabels([])

                lower_ax.set_xticklabels(tick_labels, rotation=90, fontsize=6.5)

                lower_ax.set_title(
                    f"Euclidean distance between average z-scores of S over {group1} and {group2}",
                    fontsize=7,
                )
                # lower_ax.set_ylim(0, 2)

                # upper_ax.set_ylabel(r"$z-score$ of centered $S$", fontsize=5)
                upper_ax.set_title(r"raw $S$", fontsize=8)
                mid_ax.set_title(r"$z$-scores $S$", fontsize=8)

                handles, labels = upper_ax.get_legend_handles_labels()
                labels, handles = zip(
                    *sorted(zip(labels, handles), key=lambda t: utils.natural_key(t[0]))
                )

                leg = upper_ax.legend(
                    handles,
                    labels,
                    loc="upper right",
                    bbox_to_anchor=(1, 0.9),
                    borderaxespad=0.0,
                    handlelength=1,
                    handleheight=1,
                    fontsize="x-small",
                    frameon=False,
                    bbox_transform=fig.transFigure,
                )

                for line in leg.get_lines():
                    line.set_linewidth(1)

                # legend does not fit
                fig.subplots_adjust(
                    left=0.1, right=0.85, top=0.9, bottom=0.3, hspace=0.35
                )

                pdf.savefig(fig)
                plt.close(fig)

                EUCLIDEAN = pd.concat([EUCLIDEAN, DIFFERENCES])

        EUCLIDEAN = EUCLIDEAN.sort_values(
            "meanS_Euclidean", ascending=False
        ).reset_index(drop=True)

        print(EUCLIDEAN)

    EUCLIDEAN.to_csv(
        eucl_out_FN,
        index=False,
        sep="\t",
    )

    return EUCLIDEAN
