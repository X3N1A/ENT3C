import numpy as np
import cooler as cl
import pandas as pd
import os


def check_config(config):
    FILES = config["FILES"]
    # pandas aligns the rows by index.
    # FILES[::2] has indices 0,2,4,.. and FILES[1::2] has 1,2,3
    config_files = pd.DataFrame(
        {
            "FILE": FILES[::2].reset_index(drop=True),
            "NAME": FILES[1::2].reset_index(drop=True),
        }
    )

    del config["FILES"]

    config = config[::2].reset_index(drop=True)
    config = pd.concat([config, config_files], axis=1)

    if "BR" in config["NAME"].astype(str):
        BR = True
    else:
        BR = False

    config["OUT_DIR"] = os.path.join(config["OUT_DIR"][0], "PYTHON")

    if pd.isna(config["OUT_PREFIX"]).all():
        config["OUT_PREFIX"] = "%dkb" % (config["Resolution"][0] / 1e3)
    # print("%dkb" % (config["Resolution"][0]/1e3))

    if not os.path.exists(config["OUT_DIR"][0]):
        os.makedirs(config["OUT_DIR"][0])

    config["ChrNr"] = config["ChrNr"].astype(str)

    if "," in config["ChrNr"][0]:
        CHROMOSOMES = config["ChrNr"][0].split(",")
    else:
        CHROMOSOMES = [config["ChrNr"][0]]

    config["Resolution"] = config["Resolution"].astype(str)

    if "," in str(config["Resolution"][0]):
        RESOLUTIONS = [int(r) for r in config["Resolution"][0].split(",")]
    else:
        RESOLUTIONS = [int(config["Resolution"][0])]

    FNs = config["DATA_PATH"] + "/" + config["FILE"]
    FNs = FNs.tolist()

    SUB_M_SIZE_FIX = config["SUB_M_SIZE_FIX"].unique()
    PHI_MAX = config["PHI_MAX"].unique()
    CHRSPLIT = config["CHRSPLIT"].unique()
    phi = config["phi"].unique()
    NormM = config["NormM"].unique()
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

    entropy_out_FN = (
        f"{config['OUT_DIR'].iloc[0]}/{config['OUT_PREFIX'].iloc[0]}_ENT3C_OUT.csv"
    )

    similarity_out_FN = f"{config['OUT_DIR'].iloc[1]}/{config['OUT_PREFIX'].iloc[0]}_ENT3C_similarity.csv"

    return (
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
