#!/usr/bin/env python

import sys, os, re
import pysam
from collections import defaultdict
import argparse




def main():

    parser = argparse.ArgumentParser(description="extract read alignemnts for gene and supporting cancer introns", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("--igv_introns_bed", type=str, required=True, help="igv introns bed file containing cancer introns")
    parser.add_argument("--bam", type=str, required=True, help="coordinate-sorted bam file containing all genome alignments")
    parser.add_argument("--output_prefix", type=str, required=True, help="prefix for output bam files")

    args = parser.parse_args()

    igv_introns_bed_file = args.igv_introns_bed
    bam_file = args.bam
    output_prefix = args.output_prefix

    regions_to_cancer_introns = parse_cancer_intron_regions(igv_introns_bed_file)

    samfile = pysam.AlignmentFile(bam_file, "rb")
    
    introns_only_samfile_name = output_prefix + ".cancer_intron_reads.bam"
    introns_only_samfile_obj = pysam.AlignmentFile(introns_only_samfile_name, 'wb', template=samfile)
    gene_reads_samfile_name = output_prefix + ".gene_reads.bam"
    gene_reads_samfile_obj = pysam.AlignmentFile(gene_reads_samfile_name, 'wb', template=samfile)



    for region, cancer_introns in regions_to_cancer_introns.items():

        #print(region)
        #print("\t" + str(cancer_introns))

        # get the region reads
        region_chr, region_coords = region.split(":")
        region_lend, region_rend = region_coords.split("-")
        
        gene_reads = samfile.fetch(region_chr, int(region_lend), int(region_rend))
        for read in gene_reads:
            gene_reads_samfile_obj.write(read)

            if read_has_cancer_intron(read, region_chr, cancer_introns):
                introns_only_samfile_obj.write(read)



    sys.exit(0)


def read_has_cancer_intron(read, chromosome_name, cancer_introns):

    
    blocks = read.get_blocks()

    num_blocks = len(blocks)
    if num_blocks < 2:
        return False


    for i in range(num_blocks - 1):
        intron_lend = blocks[i][1] + 1
        intron_rend = blocks[i + 1][0]
        
        intron_name = "{}:{}-{}".format(chromosome_name, intron_lend, intron_rend)

        if intron_name in cancer_introns:
            return True

    return False




def parse_cancer_intron_regions(igv_introns_bed_filename):

    regions_to_cancer_introns = defaultdict(set)

    with open(igv_introns_bed_filename, 'rt') as fh:
        for line in fh:
            line = line.rstrip()
            chr, lend, rend, annot, score, strand = line.split("\t")
            m = re.search("viewport=([^; ]+)", annot)
            if m:
                viewport = m.group(1)
                intron = "{}:{}-{}".format(chr, int(lend) + 1, rend)
                regions_to_cancer_introns[viewport].add(intron)

    return regions_to_cancer_introns



        

if __name__ == "__main__":
    main()
