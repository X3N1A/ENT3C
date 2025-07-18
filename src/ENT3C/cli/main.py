import os
import argparse
import logging
from ENT3C.core import get_similarity, get_entropy
from ENT3C.core.utils import check_config


def run_get_entropy(config_file):
    get_entropy(config_file)


def run_get_similarity(config_file):
    get_similarity(config_file)


def run_all(config_file):
    get_entropy(config_file)
    get_similarity(config_file)


def main():
    # top-level parser: main ArgumentParser instance

    parser = argparse.ArgumentParser(
        prog="ENT3C",
        description="Compute the similarity between genomic contact matrices.",
        epilog="For more info, visit https://github.com/X3N1A/ENT3C",
    )

    subparsers = parser.add_subparsers(
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

    config_file = os.path.abspath(args.config)

    print(config_file)
    print(args)

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
    ) = check_config(config_file)

    LOG_FN = f"{OUT_DIR}/{OUT_PREFIX}_logfile.log"

    logging.basicConfig(
        filename=LOG_FN,
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
    )

    logging.info(f"{args}")
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

    if args.command == "get_entropy":
        config_file = os.path.abspath(args.config)
        print("Using config file:", config_file)
        get_entropy(config_file)
        print(f"entropy table in: {entropy_out_FN}")

    if args.command == "get_similarity":
        config_file = os.path.abspath(args.config)
        print("Using config file:", config_file)
        print("Similarity accordinf to entropy file:", {entropy_out_FN})
        get_similarity(config_file)
        print(f"similarity table in: {similarity_out_FN}")

    if args.command == "run_all":
        config_file = os.path.abspath(args.config)

        print("Using config file:", config_file)
        get_entropy(config_file)
        get_similarity(config_file)
        print(f"entropy table in: {entropy_out_FN}")
        print(f"similarity table in: {similarity_out_FN}")

    LOG_FN = f"{OUT_DIR}/{OUT_PREFIX}_logfile.log"
    logging.info(args)
    print(f"log file in: {LOG_FN}")


if __name__ == "__main__":
    main()
