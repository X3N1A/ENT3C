import numpy as np
import matplotlib.pyplot as plt
import itertools
import pandas as pd
from ENT3c.core import utils


def get_similarity(config):
    (
        SUB_M_SIZE_FIX,
        PHI_MAX,
        CHRSPLIT,
        phi,
        NormM,
        CHROMOSOMES,
        RESOLUTIONS,
        BR,
        entropy_out_FN,
        similarity_out_FN,
    ) = utils.check_config(config)

    ENT3C_OUT = pd.read_csv(f"{entropy_out_FN}", sep="\t")
    ENT3C_OUT["ChrNr"] = ENT3C_OUT["ChrNr"].astype(str)

    SAMPLES = set(ENT3C_OUT["Name"])
    if any("BR" in element for element in SAMPLES):
        Biological_replicates = True
    else:
        Biological_replicates = False

    Similarity = pd.DataFrame(
        {
            "Resolution": pd.Series(dtype="int"),
            "ChrNr": pd.Series(dtype="str"),
            "Sample1": pd.Series(dtype="str"),
            "Sample2": pd.Series(dtype="str"),
            "Q": pd.Series(dtype="float"),
        }
    )

    cols = int(np.ceil(np.sqrt(len(CHROMOSOMES))))
    rows = int(np.ceil(len(CHROMOSOMES) / cols))
    if len(SAMPLES) > 1:
        comparisons = list(itertools.combinations(SAMPLES, 2))
        for Resolution in RESOLUTIONS:
            fig, axs = plt.subplots(
                nrows=rows, ncols=cols, figsize=(cols * 8, rows * 4)
            )  # Width, height in inches
            fig.suptitle(f"{Resolution}kb")

            if isinstance(axs, np.ndarray):
                axs = axs.flatten()
            else:
                axs = np.array([axs])

            for i, ChrNr in enumerate(CHROMOSOMES):
                ax = axs[i]
                ax.set_title(f"Chr{ChrNr}")

                plotted = set()

                for comp in comparisons:
                    S1 = ENT3C_OUT[
                        (ENT3C_OUT["Name"] == comp[0])
                        & (ENT3C_OUT["ChrNr"] == ChrNr)
                        & (ENT3C_OUT["Resolution"] == Resolution)
                    ]["S"]
                    S1 = np.array(S1, dtype=np.float64)
                    S2 = ENT3C_OUT[
                        (ENT3C_OUT["Name"] == comp[1])
                        & (ENT3C_OUT["ChrNr"] == ChrNr)
                        & (ENT3C_OUT["Resolution"] == Resolution)
                    ]["S"]
                    S2 = np.array(S2, dtype=np.float64)

                    non_nan_idx = ~np.isnan(S1) & ~np.isnan(S2)
                    Q = np.corrcoef(S1[non_nan_idx], S2[non_nan_idx])[0, 1]

                    new_row = pd.DataFrame(
                        {
                            "Resolution": [Resolution],
                            "ChrNr": [str(ChrNr)],
                            "Sample1": [comp[0]],
                            "Sample2": [comp[1]],
                            "Q": [Q],
                        }
                    )

                    Similarity = pd.concat([Similarity, new_row], ignore_index=True)

                    if comp[0] not in plotted:
                        ax.plot(S1, label=f"{comp[0]}")
                        plotted.add(comp[0])
                        i += 1
                    if comp[1] not in plotted:
                        ax.plot(S2, label=f"{comp[1]}")
                        plotted.add(comp[1])
                        i += 1

                if Biological_replicates:
                    Q_BR = Similarity[
                        (Similarity["ChrNr"] == ChrNr)
                        & (Similarity["Resolution"] == Resolution)
                        & (
                            Similarity["Sample1"].str.split("_").str[0]
                            == Similarity["Sample2"].str.split("_").str[0]
                        )
                    ]["Q"].mean()

                    Q_NR = Similarity[
                        (Similarity["ChrNr"] == ChrNr)
                        & (Similarity["Resolution"] == Resolution)
                        & (
                            Similarity["Sample1"].str.split("_").str[0]
                            != Similarity["Sample2"].str.split("_").str[0]
                        )
                    ]["Q"].mean()

                    title_str = (
                        rf"Chr{ChrNr} {Resolution / 1e3}kb"
                        + "\n"
                        + rf"$\overline{{Q}}_{{BR}} = {Q_BR:.2f}$ "
                        + rf"$\overline{{Q}}_{{NR}} = {Q_NR:.2f}$"
                    )
                    ax.set_title(title_str, fontsize=25)
                else:
                    title_str = rf"Chr{ChrNr} {Resolution / 1e3}kb"
                    ax.set_title(title_str, fontsize=25)

            if ChrNr == CHROMOSOMES[-1]:
                plotted = [p.replace("_", " ") for p in plotted]
                ax.legend(plotted, loc="upper right")
                ax.tick_params(labelsize=15)
                ax.title.set_fontsize(15)
                ax.xaxis.label.set_fontsize(15)
                ax.yaxis.label.set_fontsize(15)
                ax.autoscale(enable=True, axis="both", tight=True)

            for i in range(len(CHROMOSOMES), len(axs)):
                axs[i].axis("off")

            plt.tight_layout()
            plt.savefig(
                f"{config['OUT_DIR'].iloc[0]}/{config['OUT_PREFIX'].iloc[0]}_{Resolution}_ENT3C_signals.svg",
                bbox_inches="tight",
            )  # bbox_inches trims whitespace
    print("Output similarity table:")
    print(Similarity)

    Similarity.to_csv(
        similarity_out_FN,
        index=False,
        sep="\t",
    )

    return Similarity
