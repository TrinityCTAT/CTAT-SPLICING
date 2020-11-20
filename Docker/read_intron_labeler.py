#!/usr/bin/env python

import sys, os, re
import pysam
from collections import defaultdict


MIN_INTRON_LEN = 20


def main():

    usage = "usage: {} input.bam\n\n".format(sys.argv[0])
    if len(sys.argv) < 2:
        sys.stderr.write(usage)
        sys.exit(1)

    input_bam_filename = sys.argv[1]

    intron_counter = defaultdict(int)

    samfile = pysam.AlignmentFile(input_bam_filename, "rb")
    for read in samfile.fetch():
        blocks = read.get_blocks()

                
        num_blocks = len(blocks)
        if num_blocks < 2:
            continue

        chr_name = samfile.get_reference_name(read.reference_id)
        read_name = read.query_name

        
        for i in range(num_blocks - 1):
            intron_lend = blocks[i][1] + 1
            intron_rend = blocks[i + 1][0]

            if intron_rend - intron_lend >= MIN_INTRON_LEN:
                intron_name = "{}:{}-{}".format(chr_name, intron_lend, intron_rend)
                print("\t".join([read_name, intron_name]))
                

    sys.exit(0)


if __name__ == "__main__":
    main()
