#!/usr/bin/env python

import sys, os, re
import pandas as pd
import pyranges as pr
import gzip as gz
import logging
import argparse


if sys.version_info[0] != 3:
    print("This script requires Python 3")
    exit(1)


logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


utildir=os.path.dirname(os.path.realpath(__file__))


def main():

    parser = argparse.ArgumentParser(description="capture gene to intron usage stats", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("--targets_list_file", dest="targets_list_file", type=str, required=True, help="targets list file")
    parser.add_argument("--tab_gz_files_list_file", dest="tab_gz_files_list_file", type=str, required=True, help="file containing lists of SJ.tab.gz files")
    parser.add_argument("--output_file_name", dest="output_file_name", type=str, required=True, help="name of output file")

    args = parser.parse_args()
    
    targets_list_file = args.targets_list_file
    tab_gz_files_list_file = args.tab_gz_files_list_file
    output_file_name = args.output_file_name
    
    logger.info("-reading targets list: {}".format(targets_list_file))
    df = pd.read_csv(targets_list_file, sep="\t")
    df.columns = ["gene_name", "Chromosome", "Start", "End", "Strand", "gene_symbol", "gene_type"]

    df_pr = pr.PyRanges(df)


    df_concat = None

    with open(tab_gz_files_list_file) as fh:
        counter = 0
        for gz_filename in fh:
            gz_filename = gz_filename.rstrip()
            counter += 1
            logger.info("-[{}] processing {}".format(counter, gz_filename))
            
            df = map_introns(df_pr, gz_filename)
            df['sample_file'] = os.path.basename(gz_filename)
            
            if df_concat is None:
                df_concat = df
            else:
                df_concat = df_concat.append(df)

    logger.info("-* writing final output to {} *".format(output_file_name))
    df_concat.to_csv(output_file_name, sep="\t")
    
    sys.exit(0)
    

def map_introns(df_pr : pr.PyRanges, tab_gz_filename : str) -> pd.DataFrame:
    
    introns_df = pd.read_csv(tab_gz_filename, sep="\t", header=None)
    

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
          0      1      2  3  4  5   6    7   8
          0  chr1  14830  14969  2  2  1  18  233  38
    """
    

    introns_df.columns = ["Chromosome", "Start", "End",
                          "strandval", "intron_motif", "annot_status",
                          "unique_mappings", "multi_mappings", "max_spliced_align_overhang"]
    

    introns_pr = pr.PyRanges(introns_df)

    overlaps_pr = df_pr.join(introns_pr)

    overlaps_df = overlaps_pr.as_df()

    return overlaps_df



if __name__ == '__main__':
    main()
