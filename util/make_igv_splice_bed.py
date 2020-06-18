#!/usr/bin/env python


import os
import sys
import json
import math
import argparse
import pandas as pd



def main():

     ## Input Arguments
    args_parser = argparse.ArgumentParser(
        description = "Creates the IGV Report for spilce data."
    )

    args_parser.add_arguments("--all_introns", type=str, required=True, help="all introns file")
    args_parser.add_argument("--cancer_introns", type=str, required=True,  help="cancer introns file.")
    args_parser.add_argument("--genome_lib_dir", type=str, requried=True, help="path to ctat genome lib") 
    
    args_parsed = args_parser.parse_args()


    all_introns_file = args_parsed.all_introns         ## example: ../testing/__expected_output/ctat.introns.b38
    cancer_introns_file = args_parsed.cancer_introns   ## example:  ../testing/__expected_output/ctat.cancer.introns.b38
    genome_lib_dir = args_parsed.genome_lib_dir        ## for b38, use: /seq/RNASEQ/__ctat_genome_lib_building/Apr2020/GRCh38_gencode_v22_CTAT_lib_Apr032020.plug-n-play/ctat_genome_lib_build_dir


    ## what we want for output is a bed file that:
    ## -contains a row for each entry of 'all_introns'
    ## -the cancer introns should contain a bed field entry for 'display_in_table=true;'
    ## -the cancer introns also will contain additional annotation info for TCGA and GTEx occurrence info as provided in the cancer.introns file.

    ## an example entry might look like this:
    # (column headers provided below only here for ease of reference below)
    #   chr     start           end             annotations       score    strand
    ##  chr7    55019365        55155829        uniquely_mapped=74;multi_mapped=0;gene=EGFR;viewport=chr7:55013358-55207969;TCGA=GBM:28:16.57,LGG:9:1.73,STAD:1:0.25,HNSC:1:0.18;GTEx=NA;variant_name=EGFRvIII;display_in_table=true     74      +
    #
    # The 'score' can be computed as uniquely_mapped + multi_mapped
    #
    # The viewport annotation attribute can be computed based on finding the range for the gene span of the gene (or genes) which you can access from file: ${CTAT_GENOME_LIB}/ref_annot.gtf.gene_spans
    # Note, it might be easiest to cross-reference the 'all_introns' file with this ref_annot.gtf.gene_spans  file to get the viewport coordinates because of the gene identifier formatting. The ENSG-style gene ids would be best to cross-reference here.
    #
    # Again, all 'all_introns' entries end up in the BED file, but only the 'cancer_introns' have the annotations: TCGA, GTEx, variant_name, and display_in_table.
    
    

    ## earlier code as starting point below:
    
    
    # Read in the cander intron data 
    dt = pd.read_table(cancer_introns_file)

    # split chr:start:stop into seperate columns 
    split1 = dt["intron"].str.split(":",n = 1, expand = True) 
    split2 = split1[1].str.split("-",n = 1, expand = True) 
    split1["START"] = pd.to_numeric(split2[0]) -1 
    split1["END"] = pd.to_numeric(split2[1])
    split1.drop(columns =[1], inplace = True) 

    # Add the name column to the data table 
    bed_file = split1
    # Make the name column for the bed file 
    dt['uniq_mapped_str'] = 'uniquely_mapped=' + dt['uniq_mapped'].astype(str)
    dt['multi_mapped_str'] = 'multi_mapped=' + dt['multi_mapped'].astype(str)
    name = dt['uniq_mapped_str'] + ";" + dt['multi_mapped_str']

    # insert them into the bed file 
    bed_file.insert(loc = 3, column = "NAME", value = name)
    bed_file.insert(loc = 4, column = "uniquely_mapped", value = dt['uniq_mapped'])
    bed_file.insert(loc = 5, column = "strand", value = dt['strand'])

    # Sort the BED File 
    bed_file.sort_values(by=['START','END'], inplace=True, ascending=True)

    self.bed_file = bed_file

    return(self)



if __name__=='__main__':
    main()
