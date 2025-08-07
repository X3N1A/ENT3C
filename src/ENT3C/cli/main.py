import os
import argparse
import logging
from ENT3C.core import get_similarity, get_entropy, locate_largest_euclidean_diff
from ENT3C.core.utils import check_config
from ENT3C.__version__ import __version__


def run_get_entropy(config_file):
    ENT3C_OUT = get_entropy(config_file)
    return ENT3C_OUT


def run_get_similarity(config_file):
    Similarity = get_similarity(config_file)
    return Similarity


def run_all(config_file):
    ENT3C_OUT = get_entropy(config_file)
    Similarity = get_similarity(config_file)
    return ENT3C_OUT, Similarity


def run_compare_groups(config_file, group1, group2):
    EUCLIDEAN = locate_largest_euclidean_diff(config_file, group1, group2)
    return EUCLIDEAN


def main():
    # top-level parser: main ArgumentParser instance

    parser = argparse.ArgumentParser(
        prog="ENT3C",
        description="Compute the similarity between genomic contact matrices.",
        epilog="""
        Usage:
        ENT3C <command> --config=<path/to/config.json> [options]

        Commands:
            get_entropy        Generates entropy output file <entropy_out_FN> .
            get_similarity            Generates similarity output file <similarity_out_FN> from <entropy_out_FN>.
            run_all            Generates <entropy_out_FN> and <similarity_out_FN>.
            compare_groups     Compare signal groups (requires --group1 and --group2 options)

        Global Options:
            --config=<path>    Path to config JSON file (required for all commands)

        <compare_groups> Options:
        --group1=<GROUP>        First group name, must correspond to what comes before _BR* in config file.
        --group2=<GROUP>        Second group name, must correspond to what comes before _BR* in config file.

        Examples:
            ENT3C run_all --config=configs/myconfig.json
            ENT3C get_entropy --config=configs/myconfig.json
            ENT3C get_similarity --config=configs/myconfig.json
            ENT3C compare_groups --config=configs/myconfig.json --group1=H1-hESC --group2=K562
            
        For more info, please visit https://github.com/X3N1A/ENT3C.
        """,
    )

    parser.add_argument(
        "--version", action="version", version=f"%(prog)s {__version__}"
    )

    subparsers = parser.add_subparsers(
        dest="command",
        required=True,
    )  # subcommand stored in args.command

    parser_get_entropy = subparsers.add_parser(
        "get_entropy",
        description="Generates entropy output file <entropy_out_FN> according to the specifications in configuration file.",
        epilog=("Example usage:\nENT3C get_entropy --config config.json"),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser_get_entropy.add_argument(
        "--config",
        required=True,
        help="Path to json configuration file.",
    )

    parser_get_similarity = subparsers.add_parser(
        "get_similarity",
        description="Generates a similarity output file <similarity_out_FN> from an entropy input file <entropy_out_FN>.",
        epilog=("Example usage:\nENT3C get_similarity --config config.json"),
    )
    parser_get_similarity.add_argument(
        "--config",
        required=True,
        help="Path to json configuration file.",
    )
    # parser_get_similarity.add_argument(
    #    "--entropy_out_FN", required=True, help="Path to entropy table output file."
    # )

    parser_run_all = subparsers.add_parser(
        "run_all",
        description="Generates both entropy and similarity output files <entropy_out_FN> and <similarity_out_FN>.",
        epilog=("Example usage:\nENT3C run_all --config config.json"),
    )
    parser_run_all.add_argument(
        "--config",
        required=True,
        help="Path to json configuration file.",
    )

    parser_run_compare_groups = subparsers.add_parser(
        "compare_groups",
        description="Identifies the regions of lowest similarity between two groups as the largest Euclidean distance between the z-score-transformed signals.",
        epilog=(
            "Example usage:\n"
            "ENT3C compare_groups --config config.json --group1 Hffc6 --group2 A549\n"
            "\n"
            "This command compares group 'Hffc6' and group 'A549' using the specified config file. The config file contains for example:"
            "\n"
            "HFFc6_BR1, HFFc6_BR2, Hffc6_BR3 and A549_BR1, A549_BR2"
            "\n",
            "Outputs saved as: <OUT_DIR>/<OUTPUT_PREFIX>_Eucl_HFFc6vsA549.csv and"
            "\n"
            "<OUT_DIR>/{OUT_PREFIX}_Eucl_<Resolution / 1e3>kb_{group1}vs{group2}.pdf",
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser_run_compare_groups.add_argument(
        "--config",
        required=True,
        help="Path to json configuration file. Generates both entropy and similarity output files.",
    )
    parser_run_compare_groups.add_argument(
        "--group1",
        required=True,
        help="Indication of group1: must correspond to what comes before _BR* in config file.",
    )
    parser_run_compare_groups.add_argument(
        "--group2",
        required=True,
        help="Indication of group 2: must correspond to what comes before _BR* in config file.",
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
        LOG_FN,
    ) = check_config(config_file)

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

    if args.command == "compare_groups":
        config_file = os.path.abspath(args.config)
        group1 = args.group1
        group2 = args.group2

        eucl_out_FN1 = f"{OUT_DIR}/{OUT_PREFIX}_Eucl_{group1}vs{group2}.csv"
        eucl_out_FN2 = (
            f"{OUT_DIR}/{OUT_PREFIX}_Eucl_<Resolution>kb{group1}vs{group2}.csv"
        )

        print("Using config file:", config_file)
        print("Comparing", group1, " with ", group2)

        locate_largest_euclidean_diff(config_file, group1, group2)
        print(f"comparison table in: {eucl_out_FN1}")
        print(f"figures in: {eucl_out_FN2}")

    LOG_FN = f"{OUT_DIR}/{OUT_PREFIX}_logfile.log"
    logging.info(args)
    print(f"log file in: {LOG_FN}")


if __name__ == "__main__":
    main()
