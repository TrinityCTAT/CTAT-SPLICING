#!/usr/bin/env python

import sys, os, re
import pysam
from collections import defaultdict

def main():

    usage = "usage: {} input.bam region\n\tregion format is chr?:lend-rend\n\n".format(sys.argv[0])
    if len(sys.argv) < 3:
        sys.stderr.write(usage)
        sys.exit(1)

    input_bam_filename = sys.argv[1]
    region = sys.argv[2]

    intron_counter = defaultdict(int)
    
    samfile = pysam.AlignmentFile(input_bam_filename, "rb")
    for read in samfile.fetch():
        blocks = read.get_blocks()
        num_blocks = len(blocks)
        if num_blocks == 1:
            continue

        chr_name = samfile.get_reference_name(read.reference_id)
        
        for i in range(num_blocks):
            intron_lend = blocks[i][1] + 1 + 1
            intron_rend = blocks[i+1][0] + 1 - 1

            intron_name = "{}:{}-{}".format(chr_name, intron_lend, intron_rend)

            intron_counter[intron_name] += 1


    for intron, count in intron_counter.items():
        print("\t".join(intron, count))
        

    sys.exit(0)


if __name__=='__main__':
    main()
