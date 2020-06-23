#!/usr/bin/env python3

import sys, os, re
import subprocess
import argparse
import logging


import logging
FORMAT = "%(asctime)-15s: %(levelname)s %(module)s.%(name)s.%(funcName)s %(message)s"
logger = logging.getLogger(__file__)
logging.basicConfig(stream=sys.stderr, format=FORMAT, level=logging.INFO)

UTILDIR = os.path.join(os.path.dirname(__file__), "util")



def main():

    parser = argparse.ArgumentParser()
    
    parser.add_argument("--cancer_introns_tsv", dest="cancer_introns_tsv", type=str, required=True,
                        help="cancer introns database tsv file")
    
    parser.add_argument("--genome_lib_dir", dest="genome_lib_dir", type=str, required=True,
                        help="CTAT genome lib build dir")
    args=parser.parse_args()
    
    genome_lib_dir = args.genome_lib_dir
    cancer_introns_tsv_file = args.cancer_introns_tsv

    
    ensure_sorted_gene_bed(genome_lib_dir)

    index_cancer_db(cancer_introns_tsv_file, genome_lib_dir)

    logger.info("done")

    sys.exit(0)



def index_cancer_db(cancer_introns_tsv_file, genome_lib_dir):

    logger.info("indexing cancer introns database")

    cmd = str(os.path.join(UTILDIR, "index_cancer_introns_ctat_genome_db.pl") +
              " --cancer_introns_tsv {} ".format(cancer_introns_tsv_file) +
              " --ctat_genome_lib {} ".format(genome_lib_dir) )

    subprocess.check_call(cmd, shell=True)

    return




def ensure_sorted_gene_bed(genome_lib_dir):


    refGene_sorted_bed_file = os.path.join(genome_lib_dir, "refGene.sort.bed")

    if not os.path.exists(refGene_sorted_bed_file + ".gz.tbi"):

        logger.info("prepping gene regions bed")

        refGene_bed_file = os.path.join(genome_lib_dir, "refGene.bed")

        cmd = str(os.path.join(UTILDIR, "gencode_gtf_to_bed.pl") +
                  " " + os.path.join(genome_lib_dir, "ref_annot.gtf") +
                  " > {}".format(refGene_bed_file) )
        subprocess.check_call(cmd, shell=True)
        

        cmd = "sort -k 1,1 -k2,2g -k3,3g < {} > {}".format(refGene_bed_file, refGene_sorted_bed_file)
        subprocess.check_call(cmd, shell=True)

        cmd = "bgzip {}".format(refGene_sorted_bed_file)
        subprocess.check_call(cmd, shell=True)

        cmd = "tabix -p bed {}.gz".format(refGene_sorted_bed_file)
        subprocess.check_call(cmd, shell=True)

    return



if __name__=='__main__':
    main()
