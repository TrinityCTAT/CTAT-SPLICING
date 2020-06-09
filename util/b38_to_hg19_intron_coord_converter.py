#!/usr/bin/env python

import sys, os, re
from pyliftover import LiftOver


# get hg38ToHg19.over.chain.gz from http://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz

# for pyliftover docs, see: https://github.com/konstantint/pyliftover

def main():

    usage = "\n\n\tusage: {} cancer_introns.b38.annot_ready.tsv hg38ToHg19.over.chain.gz > cancer_introns.b37.annot_ready.tsv\n\n".format(sys.argv[0])

    if len(sys.argv) < 3:
        print(usage, file=sys.stderr)
        sys.exit(1)

    cancer_introns_file = sys.argv[1]
    hg_chain_file = sys.argv[2]

    lo = LiftOver('hg38ToHg19.over.chain.gz')

    with open(cancer_introns_file, 'rt') as fh:
        header = next(fh)
        header = header.rstrip()
        print(header)
        for line in fh:
            line = line.rstrip()
            vals = line.split("\t")
            intron = vals[0]
            chr, coordset = intron.split(":")
            (lend, rend) = coordset.split("-")
            lend = int(lend)
            rend = int(rend)
            
            new_lend = lo.convert_coordinate(chr, lend-1)
            #print("new_lend: {}".format(str(new_lend)))
            new_rend = lo.convert_coordinate(chr, rend-1)
            #print("new_rend: {}".format(str(new_rend)))
            if new_lend and new_rend:
                                
                new_lend_chr = new_lend[0][0]
                new_lend_coord = new_lend[0][1] + 1

                new_rend_chr = new_rend[0][0]
                new_rend_coord = new_rend[0][1] + 1
                
                if new_lend_chr != new_rend_chr or new_lend_chr != chr:
                    sys.stderr.write("-failed conversion of {}".format(line) + "  --> {} {}, {} {}\n".format(new_lend_chr, new_lend_coord, new_rend_chr, new_rend_coord))
                    continue
                
                if new_lend_coord > new_rend_coord:
                    (new_lend_coord, new_rend_coord) = (new_rend_coord, new_lend_coord)
                    
                new_intron_feature = "{}:{}-{}".format(chr, new_lend_coord, new_rend_coord)
                vals[0] = new_intron_feature
                print("\t".join(vals))

    sys.exit(0)

if __name__=='__main__':
    main()
    
