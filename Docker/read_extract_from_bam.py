#!/usr/bin/env python

import sys, os, re
import pysam
from collections import defaultdict


MIN_INTRON_LEN = 20


def main():

    usage = "usage: {} input.bam read.accs.file output.bam\n\n".format(sys.argv[0])
    if len(sys.argv) < 4:
        sys.stderr.write(usage)
        sys.exit(1)

    input_bam_filename = sys.argv[1]
    input_read_accs_file = sys.argv[2]
    output_bam_filename = sys.argv[3]
    
    reads_want = set()
    with open(input_read_accs_file, 'rt') as fh:
        for line in fh:
            line = line.rstrip()
            reads_want.add(line)

    
    intron_counter = defaultdict(int)

    samfile = pysam.AlignmentFile(input_bam_filename, "rb")
    samoutfile = pysam.AlignmentFile(output_bam_filename, 'wb', template=samfile)

    for read in samfile.fetch():
        read_name = read.query_name
                
        if read_name in reads_want:
            samoutfile.write(read)
                

    sys.exit(0)


if __name__ == "__main__":
    main()
