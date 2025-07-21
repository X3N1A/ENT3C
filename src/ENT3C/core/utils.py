import numpy as np
import cooler as cl
import pandas as pd
import os
import logging


def check_config(config_path):
    with open(config_path, "r") as file:
        config_df = pd.read_json(file)

    FILES = config_df["FILES"]
    # pandas aligns the rows by index.
    # FILES[::2] has indices 0,2,4,.. and FILES[1::2] has 1,2,3
    config_files = pd.DataFrame(
        {
            "FILE": FILES[::2].reset_index(drop=True),
            "NAME": FILES[1::2].reset_index(drop=True),
        }
    )

    del config_df["FILES"]

    config_df = config_df[::2].reset_index(drop=True)
    config_df = pd.concat([config_df, config_files], axis=1)

    if "BR" in config_df["NAME"].astype(str):
        BR = True
    else:
        BR = False

    config_df["OUT_DIR"] = os.path.join(config_df["OUT_DIR"][0], "PYTHON")

    if pd.isna(config_df["OUT_PREFIX"]).all():
        config_df["OUT_PREFIX"] = "%dkb" % (config_df["Resolution"][0] / 1e3)
    # print("%dkb" % (config_df["Resolution"][0]/1e3))

    if not os.path.exists(config_df["OUT_DIR"][0]):
        os.makedirs(config_df["OUT_DIR"][0])

    config_df["ChrNr"] = config_df["ChrNr"].astype(str)

    if "," in config_df["ChrNr"][0]:
        CHROMOSOMES = config_df["ChrNr"][0].split(",")
    else:
        CHROMOSOMES = [config_df["ChrNr"][0]]

    config_df["Resolution"] = config_df["Resolution"].astype(str)

    if "," in str(config_df["Resolution"][0]):
        RESOLUTIONS = [int(float(r)) for r in config_df["Resolution"][0].split(",")]
    else:
        RESOLUTIONS = [int(config_df["Resolution"][0])]

    FNs = config_df["DATA_PATH"] + "/" + config_df["FILE"]
    FNs = FNs.tolist()

    SUB_M_SIZE_FIX = config_df["SUB_M_SIZE_FIX"].unique()
    PHI_MAX = config_df["PHI_MAX"].unique()
    CHRSPLIT = config_df["CHRSPLIT"].unique()
    OUT_DIR = config_df["OUT_DIR"].unique()[0]
    OUT_PREFIX = config_df["OUT_PREFIX"].unique()[0]

    os.makedirs(OUT_DIR, exist_ok=True)
    phi = config_df["phi"].unique()
    NormM = config_df["NormM"].unique()
    if SUB_M_SIZE_FIX.size > 1:
        raise ValueError("SUB_M_SIZE_FIX defined multiple times in config file")
    elif PHI_MAX.size > 1:
        raise ValueError("PHI_MAX defined multiple times in config file")
    elif CHRSPLIT.size > 1:
        raise ValueError("CHRSPLIT defined multiple times in config file")
    elif phi.size > 1:
        raise ValueError("phi defined multiple times in config file")
    else:
        SUB_M_SIZE_FIX = SUB_M_SIZE_FIX[0]
        PHI_MAX = PHI_MAX[0]
        CHRSPLIT = CHRSPLIT[0]
        phi = phi[0]
        NormM = NormM[0]

    entropy_out_FN = f"{OUT_DIR}/{OUT_PREFIX}_ENT3C_OUT.csv"
    similarity_out_FN = f"{OUT_DIR}/{OUT_PREFIX}_ENT3C_similarity.csv"

    LOG_FN = f"{OUT_DIR}/{OUT_PREFIX}_logfile.log"

    logging.basicConfig(
        filename=LOG_FN,
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
    )

    logging.info(f"Inputs: {', '.join(config_df['NAME'])}")
    logging.info(
        "Apply cooler weights? no."
        if config_df["NormM"][0] == 0
        else f"Apply cooler weights? yes. Name in cooler:{', '.join(config_df['WEIGHTS_NAME'][0])}"
    )
    logging.info(f"CHRSPLIT: {CHRSPLIT}")
    logging.info(f"Sub matrix size PHI: {SUB_M_SIZE_FIX}")
    logging.info(f"Window shift phi: {phi}")
    logging.info(f"Chromosomes: {CHROMOSOMES}")
    logging.info(f"Resolutions: {RESOLUTIONS}")

    return (
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
    )


def load_cooler(FN, ChrNr, Resolution, NormM, weights_name):
    if "mcool" in FN:
        clr = cl.Cooler("%s::resolutions/%d" % (FN, Resolution))
    else:
        clr = cl.Cooler(FN)

    R = clr.binsize

    assert R == Resolution, f"Expected {Resolution} resolution in {FN}, found {R}"

    # Table selectors (chroms, bins, pixels)
    BIN_TABLE = clr.bins().fetch("chr" + str(ChrNr))

    BIN_TABLE = BIN_TABLE.reset_index()
    BIN_TABLE.rename(
        columns={
            "index": "BINS_ALL",
            "chrom": "chrs",
            "start": "START",
            "end": "END",
            "weight": "weights",
        }
    )
    BIN_TABLE["binNr"] = range(0, len(BIN_TABLE))
    if NormM:
        M = clr.matrix(balance=True).fetch("chr" + str(ChrNr))
    else:
        M = clr.matrix(balance=False).fetch("chr" + str(ChrNr))

    INCLUDE = np.where(~np.all(np.isnan(M), axis=1) & ~np.all(M == 0, axis=1))[0]

    BIN_TABLE["CONTACT"] = np.full((BIN_TABLE.shape[0],), np.nan)

    BIN_TABLE.loc[INCLUDE, "CONTACT"] = 1

    return BIN_TABLE, M


def get_cell_line(sample):
    cell = sample.split("_BR")[0]
    return cell
