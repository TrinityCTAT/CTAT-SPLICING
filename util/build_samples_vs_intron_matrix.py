#!/usr/bin/env python3

import sys, os, re
import csv
from collections import defaultdict
import argparse
import logging


logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s : %(levelname)s : %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger(__name__)


def main():

    parser = argparse.ArgumentParser(
        description="build sample-vs-intron count matrix from *.introns files",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "--output_matrix",
        type=str,
        required=True,
        help="output prefix for .count.matrix and .TPM.matrix",
    )

    parser.add_argument(
        "--intron_files",
        type=str,
        nargs="+",
        required=True,
        help="space-delimited list of ${sample_name}.*.introns files to build into matrices (need at least 2)",
    )

    args = parser.parse_args()

    output_matrix = args.output_matrix
    intron_files = args.intron_files

    if len(intron_files) < 2:
        print(
            "Error, need at least 2 introns files specified to build the matrices",
            file=sys.stderr,
        )

    intron_counts_matrix_data = defaultdict(lambda: defaultdict(int))
    intron_ids = set()

    for intron_file in intron_files:
        sample_name = os.path.basename(intron_file)
        sample_name = sample_name.split(".")[0]

        logger.info("Parsing {}".format(intron_file))
        with open(intron_file, "rt") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:

                intron = row["intron"]
                splice_pair = row["splice_pair"]
                splice_flag = row["splice_flag"]

                count = int(row["count"])

                intron_token = "^".join([intron, splice_pair, splice_flag])

                intron_ids.add(intron_token)
                intron_counts_matrix_data[sample_name][intron_token] += count

    logger.info("-writing output intron count matrix {}".format(output_matrix))

    with open(output_matrix, "wt") as counts_ofh:

        sample_names = intron_counts_matrix_data.keys()

        print("\t" + "\t".join(sample_names), file=counts_ofh)

        for intron_id in intron_ids:
            count_vals = [intron_id]

            for sample_name in sample_names:
                intron_id_sample_count = intron_counts_matrix_data[sample_name][
                    intron_id
                ]
                count_vals.append(str(intron_id_sample_count))

            print("\t".join(count_vals), file=counts_ofh)

    return


if __name__ == "__main__":
    main()
