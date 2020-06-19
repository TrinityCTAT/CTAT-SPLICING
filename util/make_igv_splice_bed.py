#!/usr/bin/env python


import os
import sys
import json
import math
import argparse
import pandas as pd
import logging
logging.basicConfig(level=logging.INFO, 
                    format='%(asctime)s : %(levelname)s : %(message)s',
                    datefmt='%H:%M:%S')
logger = logging.getLogger(__name__)

def split_intron(dt):
    # split chr:start:stop into seperate columns 
    split1 = dt["intron"].str.split(":",n = 1, expand = True) 
    split2 = split1[1].str.split("-",n = 1, expand = True) 
    split1["START"] = pd.to_numeric(split2[0]) -1 
    split1["END"] = pd.to_numeric(split2[1])
    split1.drop(columns =[1], inplace = True) 
    split1.rename({0: 'CHR'}, axis=1, inplace=True)
    split1["CHR"] = split1["CHR"].apply(str)
    split1.set_index(['CHR', 'START', 'END'])

    return(split1)


def combineColumns(df, cols):
    df = df[cols].apply(lambda row: ';'.join(row.values.astype(str)), axis=1)
    return(df)



class BEDfile:

    def __init__(self, args):

        import warnings
        self.args = args

        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Seperate Arguments and add to Object 
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        self.all_introns_file = args.all_introns         ## example: ../testing/__expected_output/ctat.introns.b38
        self.cancer_introns_file = args.cancer_introns   ## example:  ../testing/__expected_output/ctat.cancer.introns.b38
        self.genome_lib_dir = args.genome_lib_dir        ## for b38, use: /seq/RNASEQ/__ctat_genome_lib_building/Apr2020/GRCh38_gencode_v22_CTAT_lib_Apr032020.plug-n-play/ctat_genome_lib_build_dir
        self.output_bed = args.output_bed

    def createBedFile(self):
        logger.info(" Creating the BED File.")
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Read in the intron and cander intron data 
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        dt = pd.read_table(self.all_introns_file)
        cancer_dt = pd.read_table(self.cancer_introns_file)

        # set up the bed file 
        bed_file = split_intron(dt)
        cancer_local = split_intron(cancer_dt)
        # uniquely_mapped=74;multi_mapped=0;gene=EGFR;viewport=chr7:55013358-55207969;TCGA=GBM:28:16.57,LGG:9:1.73,STAD:1:0.25,HNSC:1:0.18;GTEx=NA;variant_name=EGFRvIII;display_in_table=true     74      +

        #~~~~~~~~~~~~~~~~~~~~~
        # Create the viewport
        #~~~~~~~~~~~~~~~~~~~~~
        gene_spans = os.path.join(self.genome_lib_dir, "ref_annot.gtf.gene_spans")
        gene_spans_df = pd.read_table(gene_spans, header = None)
        split_gene = dt['genes'].str.split("^",n = 1, expand = True) 

        df_gene = split_gene.merge(gene_spans_df, 
                            how = "left", 
                            left_on = 1, 
                            right_on = 0)

        #~~~~~~~~~~~~~~~~~~~~~
        # Make the name column for the bed file 
        # Put the columns together
        #~~~~~~~~~~~~~~~~~~~~~
        dt['uniq_mapped_str']  = 'uniquely_mapped=' + dt['uniq_mapped'].astype(str)
        dt['multi_mapped_str'] = 'multi_mapped=' + dt['multi_mapped'].astype(str)
        dt['gene']             = 'gene=' + split_gene[0].astype(str)

        dt['viewport']         = df_gene['1_y'] + ":" + df_gene[2].astype(str) + "-" + df_gene[3].astype(str)

        # concatenate to make name column
        name = dt['uniq_mapped_str'] + ";" + dt['multi_mapped_str'] + ";" + dt['gene']

        # insert them into the bed file 
        bed_file.insert(loc = 3, column = "NAME", value = name)
        bed_file.insert(loc = 4, column = "total_mapped", value = dt['uniq_mapped'] + dt['multi_mapped'])
        bed_file.insert(loc = 5, column = "strand", value = dt['strand'])
        bed_file = bed_file.assign(viewport="viewport=" + dt['viewport']) 
        
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Edit the Cancer introns 
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # covert the dataframe into a series 
        ## The Keys to use, (column names)
        keys = list(cancer_local.columns.values)
        bed_file = bed_file.set_index(keys)
        i1 = bed_file.index
        cancer_local = cancer_local.set_index(keys)
        i2 = cancer_local.index
        
        # ~ means not in
        cancer_temp_df = bed_file.loc[i1[i1.isin(i2)]] 
        temp_df = pd.DataFrame()
        
        # convert NaN'ss into NA's
        cancer_dt = cancer_dt.fillna('NA')

        # Add the new columns into the temp dataframe 
        temp_df['NAME'] = cancer_temp_df["NAME"]
        temp_df['TCGA'] = list('TCGA=' + cancer_dt['TCGA_sample_counts'].astype(str))
        temp_df['GTEx'] = list('GTEx=' + cancer_dt['GTEx_sample_counts'].astype(str))
        temp_df['variant_name'] = list('variant_name=' + cancer_dt['variant_name'].astype(str))
        
        # combine the new columns into one, then replace NAME with the newly created column 
        cols = ["NAME", "TCGA", "GTEx", "variant_name"]
        replace_NAME = combineColumns(df = temp_df, cols = cols)
        cancer_temp_df = cancer_temp_df.assign(NAME=replace_NAME)


        # Now update the existing be file to include the new NAME changes for the cancer introns 
        bed_file.update(cancer_temp_df)

        # tack on viewport to just the cancer intron names:
        cancer_temp_df = bed_file.loc[i1[i1.isin(i2)]]
        cols = ["NAME", "viewport"]
        replace_NAME = combineColumns(df = cancer_temp_df, cols = cols)
        cancer_temp_df = cancer_temp_df.assign(NAME=replace_NAME) 
        bed_file.update(cancer_temp_df)
        bed_file = bed_file.drop('viewport', axis=1)
        bed_file = bed_file.reset_index()
        
        # Sort the BED File 
        bed_file.sort_values(by=['START','END'], inplace=True, ascending=True)

        self.bed_file = bed_file
        return(self)
    
    def saveBedFile(self):
        logger.info("Saving Bed File as {}".format(self.output_bed))
        ofh = open(self.output_bed, "wt") 
        
        # Convert the bed file pandas Data Frame to a csv string format 
        text = self.bed_file.to_csv(index=False, header=None, sep="\t")
        ofh.write(text) # Write to the temporary file 
        ofh.close()
                    
                    
def main():

    ## Input Arguments
    args_parser = argparse.ArgumentParser(
        description = "Creates the IGV Report for spilce data."
    )

    args_parser.add_argument("--all_introns",    type=str, required=True,  help="all introns file")
    args_parser.add_argument("--cancer_introns", type=str, required=True,  help="cancer introns file.")
    args_parser.add_argument("--genome_lib_dir", type=str, required=True,  help="path to ctat genome lib") 
    args_parser.add_argument("--output_bed", type=str, required=True, help='output bed filename')

    args = args_parser.parse_args()

    # Create teh object
    bed_file = BEDfile(args)
    # Create the bed file 
    bed_file = bed_file.createBedFile()
    # Save the Bed file 
    bed_file.saveBedFile()


if __name__=='__main__':
    main()
