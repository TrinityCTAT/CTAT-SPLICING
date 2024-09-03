#!/usr/bin/env python3

import sys, os, re
import pysam
from collections import defaultdict
import logging
import gzip
import intervaltree as itree
import argparse

logging.basicConfig(level=logging.INFO, 
                    format='%(asctime)s : %(levelname)s : %(message)s',
                    datefmt='%H:%M:%S')
logger = logging.getLogger(__name__)


MIN_MAPPING_QUALITY = 60

OK_SPLICES = ("GT--AG", "GC--AG", "AT--AC", # forward strand
              "CT--AC", "CT--GC", "GT--AT" # reverse strand
              )


def main():

    parser = argparse.ArgumentParser(description="annotate introns", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("--genome_lib_dir", type=str, default=os.environ["CTAT_GENOME_LIB"], required=False, help="CTAT genome lib")

    parser.add_argument("--bam", type=str, required = False, help="bam file")

    parser.add_argument("--introns", type=str, required=False, help="introns file")

    args = parser.parse_args()

    genome_lib_dir = args.genome_lib_dir
    bam_file = args.bam
    introns_file = args.introns

    if genome_lib_dir is None or not os.path.exists(genome_lib_dir):
        exit("Error, must set --genome_lib_dir.  Use -h for usage info")
    if not (bam_file is not None or introns_file is not None):
        exit("Error, must set --bam or --introns, run with -h for usage info")

        
    usage = "\n\nusage: {} alignments.bam genome.fasta\n\n".format(sys.argv[0])

    if len(sys.argv) < 3:
        exit(usage)

    
    genome_fasta_filename = os.path.join(genome_lib_dir, "ref_genome.fa")
    
    exon_coords_file = os.path.join(genome_lib_dir, "ref_annot.gtf.mini.sortu")
    intron_db_file = os.path.join(genome_lib_dir, "ctat_splicing_lib/known_intron_annotations.tsv.gz")

    exon_itree = get_exon_itree(exon_coords_file)

    intron_db = get_intron_db(intron_db_file)
    
    fasta_reader = pysam.FastaFile(genome_fasta_filename)

    intron_counter = defaultdict( lambda: defaultdict(int) )
    
    if (bam_file is not None):
        evaluate_introns_from_bam_file(bam_file, intron_counter)

    elif (introns_file):
        capture_introns_from_file(introns_file, intron_counter)
    else:
        raise RuntimeError("shouldn't get here")


    ######
    logger.info("reporting intron annotations")

    # print header
    print("\t".join(["intron","splice_pair", "splice_flag", "count", "gene", "TCGA|Cryptic", "GTEx"]))
    
    for chrom, chrom_icounter in intron_counter.items():

        if "_" in chrom: # only main chromosomes
            continue
        
        chrom_seq = fasta_reader.fetch(chrom)
        
        for intron, count in sorted(chrom_icounter.items()):

            chromval, coords_val = intron.split(":")
            lend, rend = coords_val.split("-")
                        
            lend = int(lend)
            rend = int(rend)

            intron_key = intron

            cryptic_or_known = intron_db[intron_key] if intron_key in intron_db else "Cryptic"
                        
            left_dinuc = chrom_seq[lend-1:lend+1]
            right_dinuc = chrom_seq[rend-1-1:rend]

            splice_tok = f"{left_dinuc}--{right_dinuc}"

            splice_flag = "OK" if splice_tok in OK_SPLICES else "NON"

            gene_annots = get_overlapping_gene(chrom, lend, rend, exon_itree)
                        
            print("\t".join([intron_key,
                             splice_tok,
                             splice_flag,
                             str(count),
                             ";".join(list(gene_annots)),
                             cryptic_or_known ] ) )
            


def evaluate_introns_from_bam_file(bam_filename, intron_counter):

    logger.info("-searching bam file")

    bam_reader = pysam.AlignmentFile(bam_filename, "rb")
    for read in bam_reader.fetch():
        if read.mapping_quality < MIN_MAPPING_QUALITY:
            continue

        if read.is_secondary:
            continue

        chrom = bam_reader.get_reference_name(read.reference_id)

        introns = bam_reader.find_introns([read])

        for coordpair in introns.keys():
            lend, rend = coordpair
            lend += 1

            intron_key = f"{chrom}:{lend}-{rend}"
            
            intron_counter[chrom][intron_key] += 1
                

    return



            
def capture_introns_from_file(introns_file, intron_counter):

    with open(introns_file, "rt") as fh:
        for line in fh:
            m = re.match("(chr[^\\:]+):(\d+\\-\d+)\s", line)
            if m is None:
                raise RuntimeError("Error, couldn't parse intron from {} in chr\d+:\d+-\d+ format".format(line))
            chrom = m.group(1)
            coords = m.group(2)
            token = f"{chrom}:{coords}"
            #print("Token: {}".format(token))
            intron_counter[chrom][token] += 1

    return
        
    

def get_intron_db(known_introns_file):

    intron_db = dict()

    logger.info("building intron db")
    
    with gzip.open(known_introns_file, "rt") as fh:
        header = next(fh)
        assert re.match("intron\tgenes", header), "Error, {} lacks expected header".format(header)
        
        for intron_info in fh:
            intron_info = intron_info.rstrip()
            intron_vals = intron_info.split("\t")
            intron = intron_vals[0]
            intron_info = "\t".join(intron_vals[2:])
            intron_db[intron] = intron_info
    
    return intron_db


def get_exon_itree(exon_coords_file):

    exon_itree = defaultdict(lambda: itree.IntervalTree())
    
    with open(exon_coords_file, "rt") as fh:
        for line in fh:
            line = line.rstrip()
            vals = line.split("\t")
            gene_info = vals[8].split()[1].replace("\"", "")
            chrom = vals[0]
            lend = int(vals[3])
            rend = int(vals[4])

            exon_itree[chrom][lend:rend+1] = gene_info

    return exon_itree
            



def get_overlapping_gene(chrom, lend, rend, exon_itree):

    ## try ranges of coordinate extensions to find nearest overlapping exon.

    # get overlapping genes
    overlapping_genes = set()
    
    ranges = (0, 100, 500, 1000, 5000, 10000)

    for i_range in ranges:
        for overlapping_gene in exon_itree[chrom][ (lend - i_range) : (rend + i_range)]:
            overlapping_genes.add(overlapping_gene.data)        
    
        if len(overlapping_genes) > 0:
            return overlapping_genes

    # no overlaps found
    overlapping_genes.add("NA")
    return overlapping_genes



if __name__=='__main__':
    main()
