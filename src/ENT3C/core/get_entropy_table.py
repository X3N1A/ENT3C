import numpy as np
import pandas as pd
import logging

from ENT3C.core import utils
from ENT3C.core.vN_entropy import vN_entropy


def get_entropy(config_file):
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

    print(f"Inputs: {', '.join(config_df['NAME'])}")
    print(
        "Apply cooler weights? no."
        if config_df["NormM"][0] == 0
        else f"Apply cooler weights? yes. Name in cooler:{config_df['WEIGHTS_NAME'][0]}"
    )
    print(f"CHRSPLIT: {CHRSPLIT}")
    print(f"Sub matrix size PHI: {SUB_M_SIZE_FIX}")
    print(f"Window shift phi: {phi}")
    print(f"Chromosomes: {CHROMOSOMES}")
    print(f"Resolutions: {RESOLUTIONS}")

    print(f"Output Directory: {config_df['OUT_DIR'][0]}")
    print(f"Output Prefix: {config_df['OUT_PREFIX'][0]}")

    if utils.overwrite(entropy_out_FN):
        header_written = False
        with open(entropy_out_FN, "w") as f:
            pass
        ENT3C_OUT = pd.DataFrame()
        print(f"Generating new file: {entropy_out_FN}")

        logging.basicConfig(
            filename=LOG_FN,
            level=logging.INFO,
            format="%(asctime)s - %(levelname)s - %(message)s",
        )
        for Resolution in RESOLUTIONS:
            Resolution = int(Resolution)
            for ChrNr in CHROMOSOMES:
                ######## locate common empty bins ########
                EXCLUDE = set()  # use set() to keep unique values

                for FN in FNs:
                    BIN_TABLE, M = utils.load_cooler(
                        FN,
                        ChrNr,
                        Resolution,
                        config_df["NormM"][0],
                        config_df["WEIGHTS_NAME"][0],
                    )

                    if NormM == 0:
                        EXCLUDE.update(BIN_TABLE[BIN_TABLE["CONTACT"].isna()]["binNr"])
                    elif NormM == 1:
                        EXCLUDE.update(
                            BIN_TABLE[
                                BIN_TABLE["CONTACT"].isna() | BIN_TABLE["weight"].isna()
                            ]["binNr"]
                        )

                #########################################
                # N = len(FNs)
                # cols = int(np.ceil(np.sqrt(N)))
                # rows = int(np.ceil(N / cols))
                # fig, axs = plt.subplots(rows, cols, figsize=(cols * 4, rows * 4))
                # axs = axs.flatten()
                # cmap = mpl.colormaps.get_cmap("magma")
                # cmap.set_bad(color="black")
                # CELL = config[(config["DATA_PATH"] + "/" + config["FILE"] == FN)]["NAME"].iloc[0]
                # axs[i].set_title(FN)
                # i = i + 1
                # plt.tight_layout()
                # plt.show()
                ######### compute entropies ##############
                EXCLUDE = set(sorted(EXCLUDE))

                for f in range(0, len(FNs)):
                    BIN_TABLE, M = utils.load_cooler(
                        FNs[f],
                        ChrNr,
                        Resolution,
                        config_df["NormM"][0],
                        config_df["WEIGHTS_NAME"][0],
                    )

                    INCLUDE = set(range(0, M.shape[0]))
                    INCLUDE = INCLUDE - EXCLUDE

                    BIN_TABLE = BIN_TABLE.iloc[list(INCLUDE), :].reset_index(drop=True)
                    M = M[np.meshgrid(sorted(INCLUDE), sorted(INCLUDE), indexing="ij")]

                    S, SUB_M_SIZE, WN, phi, BIN_TABLE_NEW = vN_entropy(
                        M, SUB_M_SIZE_FIX, CHRSPLIT, PHI_MAX, phi, BIN_TABLE
                    )

                    N = len(S)
                    DICT = {
                        "Name": [config_df["NAME"].iloc[f]] * N,
                        "ChrNr": [str(ChrNr)] * N,
                        "Resolution": [Resolution] * N,
                        "n": [SUB_M_SIZE] * N,
                        "PHI": [WN] * N,
                        "phi": [phi] * N,
                        "binNrStart": BIN_TABLE_NEW.iloc[:, 0],
                        "binNrEnd": BIN_TABLE_NEW.iloc[:, 1],
                        "START": BIN_TABLE_NEW["start"],
                        "END": BIN_TABLE_NEW["end"],
                        "S": S,
                    }

                    OUT_CURRENT = pd.DataFrame(DICT)  # stores DFs in list
                    ENT3C_OUT = pd.concat([ENT3C_OUT, OUT_CURRENT], ignore_index=True)

                    OUT_CURRENT.to_csv(
                        entropy_out_FN,
                        index=False,
                        sep="\t",
                        na_rep="NaN",
                        mode="a",
                        header=not header_written,
                    )
                    header_written = True

                    print(
                        f"Chr{ChrNr} {Resolution / 1e3}kb for {FNs[f]} written to {entropy_out_FN}"
                    )

                    logging.info(
                        f"Chr{ChrNr} {Resolution / 1e3}kb for {FNs[f]} written to {entropy_out_FN}"
                    )

    else:
        ENT3C_OUT = pd.read_csv(
            entropy_out_FN, sep="\t", na_values="NaN", dtype={"ChrNr": str}
        )

    print("ENT3C Output table:")
    print(ENT3C_OUT)

    return ENT3C_OUT
