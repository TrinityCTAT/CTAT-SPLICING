#!/usr/bin/env python

import sys, os, re
import pysam
from collections import defaultdict


MIN_INTRON_LEN = 20


def main():

    usage = "usage: {} input.bam sample_name\n\n".format(sys.argv[0])
    if len(sys.argv) < 3:
        sys.stderr.write(usage)
        sys.exit(1)

    input_bam_filename = sys.argv[1]
    sample_name = sys.argv[2]

    intron_counter = defaultdict(int)

    samfile = pysam.AlignmentFile(input_bam_filename, "rb")
    for read in samfile.fetch():
        blocks = read.get_blocks()

        num_blocks = len(blocks)
        if num_blocks < 2:
            continue

        chr_name = samfile.get_reference_name(read.reference_id)

        for i in range(num_blocks - 1):
            intron_lend = blocks[i][1] + 1
            intron_rend = blocks[i + 1][0]

            if intron_rend - intron_lend >= MIN_INTRON_LEN:
                intron_name = "{}:{}-{}".format(chr_name, intron_lend, intron_rend)
                intron_counter[intron_name] += 1

    for intron, count in intron_counter.items():
        print("\t".join([sample_name, intron, str(count)]))

    sys.exit(0)


if __name__ == "__main__":
    main()
