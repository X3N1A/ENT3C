import os
import pandas as pd
import argparse
import logging
from ENT3c.core import get_similarity, get_entropy
from ENT3c.core.utils import check_config


def run_get_entropy(config):
    get_entropy(config)


def run_get_similarity(config):
    get_similarity(config)


def run_all(config):
    get_entropy(config)
    get_similarity(config)


def main():
    # top-level parser: main ArgumentParser instance

    parser = argparse.ArgumentParser(
        prog="ENT3C",
        description="Compute the similarity between genomic contact matrices.",
        epilog="For more info, visit https://github.com/X3N1A/ENT3C",
    )

    subparsers = parser.add_argument(
        dest="command", required=True
    )  # subcommand stored in args.command

    parser_get_entropy = subparsers.add_parser("get_entropy")
    parser_get_entropy.add_argument(
        "--config",
        required=True,
        help="Path to json configuration file. Generates an entropy output file.",
    )

    parser_get_similarity = subparsers.add_parser("get_similarity")
    parser_get_similarity.add_argument(
        "--config",
        required=True,
        help="Path to json configuration file. Generates a similarity output file from an entropy input file <entropy_out_FN>.",
    )
    # parser_get_similarity.add_argument(
    #    "--entropy_out_FN", required=True, help="Path to entropy table output file."
    # )

    parser_run_all = subparsers.add_parser("run_all")
    parser_run_all.add_argument(
        "--config",
        required=True,
        help="Path to json configuration file. Generates both entropy and similarity output files.",
    )

    args = parser.parse_args()

    if not args.command:
        parser.print_help()

    config_path = os.path.abspath(args.config)

    with open(config_path, "r") as file:
        config = pd.read_json(file)

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
    ) = check_config(config)

    LOG_FN = f"{config['OUT_DIR'].iloc[0]}/{config['OUT_PREFIX'].iloc[0]}_logfile.log"

    logging.basicConfig(
        filename=LOG_FN,
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
    )

    logging.info(f"{args}")
    logging.info(f"Inputs: {', '.join(config['NAME'])}")
    logging.info(
        "Apply cooler weights? no."
        if config["NormM"][0] == 0
        else f"Apply cooler weights? yes. Name in cooler:{', '.join(config['WEIGHTS_NAME'][0])}"
    )
    logging.info(f"CHRSPLIT: {CHRSPLIT}")
    logging.info(f"Sub matrix size PHI: {SUB_M_SIZE_FIX}")
    logging.info(f"Window shift phi: {phi}")
    logging.info(f"Chromosomes: {CHROMOSOMES}")
    logging.info(f"Resolutions: {RESOLUTIONS}")

    if args.command == "get_entropy":
        print("Using config file:", args.config)
        get_entropy(config)
        print(f"entropy table in: {entropy_out_FN}")

    if args.command == "get_similarity":
        print("Using config file:", args.config)
        print("Similarity accordinf to entropy file:", {entropy_out_FN})
        get_similarity(config)
        print(f"similarity table in: {similarity_out_FN}")

    if args.command == "run_all":
        print("Using config file:", args.config)
        get_entropy(config)
        get_similarity(config)
        print(f"entropy table in: {entropy_out_FN}")
        print(f"similarity table in: {similarity_out_FN}")

    print(f"log file in: {LOG_FN}")


if __name__ == "__main__":
    main()
