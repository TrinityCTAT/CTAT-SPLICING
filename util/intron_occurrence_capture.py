#!/usr/bin/env python

import sys, os, re, io
import pandas as pd
import pyranges as pr
import gzip as gz
import logging
import argparse
from collections import defaultdict
import subprocess

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
    parser.add_argument("--debug", "-d", dest="DEBUG", action='store_true', default=False)
        
    args = parser.parse_args()
    
    ctat_genome_lib = args.ctat_genome_lib
    tab_gz_files_list_file = args.tab_gz_files_list_file
    output_file_name = args.output_file_name
    db_class = args.db_class

    if args.DEBUG:
        logger.setLevel(logging.DEBUG)


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
        for filename_pair in fh:
            filename_pair = filename_pair.rstrip()
            splice_tab_gz_file, chimeric_out_introns_file = filename_pair.split("\t")
            
            counter += 1
            logger.info("-[{}] processing {}".format(counter, splice_tab_gz_file))
            
            introns_dict = map_introns_from_splice_tab(splice_tab_gz_file, chr_intron_bounds)

            logger.info("-[{}] processing {}".format(counter, chimeric_out_introns_file))
            introns_dict = supplement_introns_from_chimeric_junctions_file(chimeric_out_introns_file, introns_dict, chr_intron_bounds)
            
            ## output record:
            sample_name = os.path.basename(splice_tab_gz_file).replace(".SJ.out.tab.gz", "")
            for intron in introns_dict.values():
                ofh.write("\t".join([db_class, sample_name, intron.genes,
                                     intron.chromosome, intron.lend, intron.rend,
                                     intron.strand,
                                     intron.intron_motif, intron.annotated_flag,
                                     str(intron.uniq_mapped), str(intron.multi_mapped),
                                     intron.max_splice_overhang]) + "\n")
                
            


    logger.info("Done.")
    
    sys.exit(0)
    


class Intron:

    def __init__(self, chromosome, lend, rend, strand, intron_motif, annotated_flag,
                 uniq_mapped, multi_mapped, max_splice_overhang, genes):

        self.chromosome = chromosome
        self.lend = lend
        self.rend = rend
        self.strand = strand
        self.intron_motif = intron_motif
        self.annotated_flag = annotated_flag
        self.uniq_mapped = uniq_mapped
        self.multi_mapped = multi_mapped
        self.max_splice_overhang = max_splice_overhang
        self.genes = genes

    def __repr__(self):
        return("^".join([self.chromosome, self.lend, self.rend, str(self.uniq_mapped), str(self.multi_mapped)]))

        

def map_introns_from_splice_tab(tab_gz_filename : str,
                                chr_intron_bounds : dict) -> dict:
        
    
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


    introns_dict = dict()
        
    with gz.open(tab_gz_filename, 'rt') as fh:
        for line in fh:
            line = line.rstrip()
            vals = line.split("\t")
            chr = vals[0]
            intron_lend = vals[1]
            intron_rend = vals[2]

            strand = '+' if vals[3] == "1" else '-'

            intron_motif = vals[4]
            annotated_flag = vals[5]
            uniq_mapped = vals[6]
            multi_mapped = vals[7]
            max_splice_overhang = vals[8]
            
            
            intron_tok_A = "{}:{}:{}".format(chr, intron_lend, strand)
            intron_tok_B = "{}:{}:{}".format(chr, intron_rend, strand)
            
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

                    
                intron_tok = "{}:{}-{}".format(chr, intron_lend, intron_rend)
                intron_obj = Intron(chr, intron_lend, intron_rend, strand,
                                    intron_motif, annotated_flag, int(uniq_mapped),
                                    int(multi_mapped), max_splice_overhang, genes_entry)

                introns_dict[intron_tok] = intron_obj


    return introns_dict



def supplement_introns_from_chimeric_junctions_file(chimeric_out_introns_file : str,
                                                    introns_dict : dict,
                                                    chr_intron_bounds : dict) -> dict:

    with open(chimeric_out_introns_file) as fh:
        for line in fh:
            if line[0] == "#":
                continue

            line = line.rstrip()
            if line == "":
                continue
            vals = line.split("\t")
            if len(vals) != 3:
                raise RuntimeError("Error, couldn't parse line in to three fields: {}".format(line))
            intron, uniq_map, multi_map = vals
            uniq_map = int(uniq_map)
            multi_map = int(multi_map)
            if intron in introns_dict:
                intron_obj = introns_dict[intron]
                intron_obj.uniq_mapped += uniq_map
                intron_obj.multi_mapped += multi_map
                logger.info("-supplementing existing intron: " + str(intron_obj) + " with uniq: {}, multi: {}".format(uniq_map, multi_map))
            else:
                # see if intron has known splice sites.
                intron_obj = try_make_intron_obj(intron, chr_intron_bounds, uniq_map, multi_map)
                if intron_obj is not None:
                    introns_dict[intron] = intron_obj
                    logger.info("-supplementing NEW intron: " + str(intron_obj))
            
    return introns_dict



def try_make_intron_obj(intron : str, chr_intron_bounds : dict, uniq_map : int, multi_map : int) -> Intron:
    chr, coords = intron.split(':')
    lend, rend = coords.split('-')

    for orient in ('+', '-'):
        left_splice_token = ":".join([chr, lend, orient])
        right_splice_token = ":".join([chr, rend, orient])

        if (left_splice_token in chr_intron_bounds and
            right_splice_token in chr_intron_bounds):

            genes_left = chr_intron_bounds[left_splice_token]
            genes_right = chr_intron_bounds[right_splice_token]
            genes = genes_left.union(genes_right)
            genes = ",".join(list(genes))
            
            intron_obj = Intron(chr, lend, rend, orient, "-1", "1", uniq_map, multi_map, "-1", genes)

            return intron_obj

    # no intron could be created
    return None



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
            orient = vals[6]
            info = vals[8]

            m = re.search("gene_id \"([^\"]+)\"", info)
            if m:
                gene_id = m.group(1)
            else:
                raise RuntimeError("Error, no gene id extracted from line: {}".format(line))

            lend_token = ":".join([chr, str(lend-1), orient])
            rend_token = ":".join([chr, str(rend+1), orient])
            
            chr_intron_bounds[lend_token].add(gene_id)
            chr_intron_bounds[rend_token].add(gene_id)


    return chr_intron_bounds



if __name__ == '__main__':
    main()
