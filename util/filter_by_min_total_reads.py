#!/usr/bin/env python3

import sys, os, re
import pandas as pd
import argparse


def main():

    args_parser = argparse.ArgumentParser(
        description = "filter results by total number of reads"
        )

    args_parser.add_argument("--cancer_intron_candidates",
                             type=str, required=True,
                             help="cancer intron candidates file")

    args_parser.add_argument("--min_total_reads", type=int, required=True,
                             help="min number of total reads required")

    args = args_parser.parse_args()


    cancer_introns_file = args.cancer_intron_candidates
    min_total_reads = args.min_total_reads



    data = pd.read_table(cancer_introns_file)
    data['total_reads'] = data['uniq_mapped'] + data['multi_mapped'] 

    filtered_data = data[ data['total_reads'] >= min_total_reads ]

    filtered_data = filtered_data.sort_values('total_reads', ascending=False)

    filtered_data = filtered_data.drop('total_reads', axis=1)

    sys.stdout.write(filtered_data.to_csv(sep="\t", index=False, na_rep="NA"))
    
    

    sys.exit(0)



if __name__=='__main__':
    main()
