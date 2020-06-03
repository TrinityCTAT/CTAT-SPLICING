#!/usr/bin/env python

import sys, os, re, io
import pandas as pd
import pyranges as pr
import gzip as gz
import logging
import argparse
from collections import defaultdict


if sys.version_info[0] != 3:
    print("This script requires Python 3")
    exit(1)


logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

utildir=os.path.dirname(os.path.realpath(__file__))


def main():


    parser = argparse.ArgumentParser(description="capture gene to intron usage stats", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("--ctat_genome_lib", dest="ctat_genome_lib", type=str, required=True, help="ctat genome lib build dir")
    parser.add_argument("--tab_gz_files_list_file", dest="tab_gz_files_list_file", type=str, required=True, help="file containing lists of SJ.tab.gz files")
    parser.add_argument("--output_file_name", dest="output_file_name", type=str, required=True, help="name of output file")
    parser.add_argument("--db_class", dest="db_class", type=str, required=True, help="database class: ie. GTEx or TCGA")
    
    args = parser.parse_args()
    
    ctat_genome_lib = args.ctat_genome_lib
    tab_gz_files_list_file = args.tab_gz_files_list_file
    output_file_name = args.output_file_name
    db_class = args.db_class


    ## get target list info.
    targets_list_file = os.path.join(ctat_genome_lib, "ref_annot.gtf.mini.sortu")
    chr_intron_bounds = populate_intron_bounds(targets_list_file)

    ofh = open(output_file_name, 'wt')
    
    column_header = [
        "class", "sample", "genes",
        "Chromosome", "Start", "End",
        "strandval", "intron_motif", "annot_status",
        "unique_mappings", "multi_mappings", "max_spliced_align_overhang" ]
    
    ofh.write("\t".join(column_header) + "\n")
    
    with open(tab_gz_files_list_file) as fh:
        counter = 0
        for gz_filename in fh:
            gz_filename = gz_filename.rstrip()
            counter += 1
            logger.info("-[{}] processing {}".format(counter, gz_filename))
            
            map_introns(gz_filename, chr_intron_bounds, db_class, ofh)


    logger.info("Done.")
    
    sys.exit(0)
    

def map_introns(tab_gz_filename : str, chr_intron_bounds : dict, db_class : str, ofh : io.IOBase):

    
    """  from the STAR manual
    4.4 Splice junctions.
    SJ.out.tab contains high confidence collapsed splice junctions in tab-delimited format. Note that
    STAR defines the junction start/end as intronic bases, while many other software define them as
    exonic bases. The columns have the following meaning:
    column 1: chromosome
    column 2: first base of the intron (1-based)
    column 3: last base of the intron (1-based)
    column 4: strand (0: undefined, 1: +, 2: -)
    column 5: intron motif: 0: non-canonical; 1: GT/AG, 2: CT/AC, 3: GC/AG, 4: CT/GC, 5:
    AT/AC, 6: GT/AT
    column 6: 0: unannotated, 1: annotated (only if splice junctions database is used)
    column 7: number of uniquely mapping reads crossing the junction
    column 8: number of multi-mapping reads crossing the junction
    column 9: maximum spliced alignment overhang
    """

    """ example row
               0      1     2    3  4  5   6    7   8
             chr1  14830  14969  2  2  1  18  233  38
    """
    

    """
    column_header = [
        "class", "sample", "genes",
        "Chromosome", "Start", "End",
        "strandval", "intron_motif", "annot_status",
        "unique_mappings", "multi_mappings", "max_spliced_align_overhang" ]
    """

    with gz.open(tab_gz_filename, 'rt') as fh:
        for line in fh:
            line = line.rstrip()
            vals = line.split("\t")
            chr = vals[0]
            intron_lend = vals[1]
            intron_rend = vals[2]

            intron_tok_A = "{}:{}".format(chr, intron_lend)
            intron_tok_B = "{}:{}".format(chr, intron_rend)
            
            genesA = None
            genesB = None

            if intron_tok_A in chr_intron_bounds:
                genesA = chr_intron_bounds[intron_tok_A]

            if intron_tok_B in chr_intron_bounds:
                genesB = chr_intron_bounds[intron_tok_B]

            if genesA and genesB:

                if genesA == genesB:
                    genes_entry = ",".join(list(genesA))
                else:
                    genes_entry = ",".join(list(genesA)) + "--" + ",".join(list(genesB))

                
                ## output record:
                sample_name = os.path.basename(tab_gz_filename).replace(".SJ.out.tab.gz", "")
                
                ofh.write("\t".join([db_class, sample_name, genes_entry, *vals]) + "\n")
                
                                    



    return



def populate_intron_bounds(targets_list_file : str) -> dict:

    logger.info("-reading targets list: {}".format(targets_list_file))
    
    chr_intron_bounds = defaultdict(set)

    with open(targets_list_file) as fh:
        for line in fh:
            line = line.rstrip()
            vals = line.split("\t")
            assert(len(vals) >= 9)
            chr = vals[0]
            lend = int(vals[3])
            rend = int(vals[4])
            info = vals[8]

            m = re.search("gene_id \"([^\"]+)\"", info)
            if m:
                gene_id = m.group(1)
            else:
                raise RuntimeError("Error, no gene id extracted from line: {}".format(line))

            lend_token = ":".join([chr, str(lend-1)])
            rend_token = ":".join([chr, str(rend+1)])

            chr_intron_bounds[lend_token].add(gene_id)
            chr_intron_bounds[rend_token].add(gene_id)


    return chr_intron_bounds



if __name__ == '__main__':
    main()
