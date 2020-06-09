#!/usr/bin/env python

import os, re, sys
import sqlite3
import logging
import argparse

if sys.version_info[0] != 3:
    print("This script requires Python 3")
    exit(1)


logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

COMMIT_EVERY=1000


def main():

    parser = argparse.ArgumentParser(description="loads introns into intron sqlite3 db", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("--sqlite3_db", dest="sqlite3_db", type=str, required=True, help="sqlite3_db name")
    parser.add_argument("--input", dest="input", type=str, required=True, help="input data table", nargs='+')

    args = parser.parse_args()

    sqlite3_dbname = args.sqlite3_db
    input_files = args.input

    GTEx_sample_to_tissue_type = parse_GTEx_sample_types()
    
    ## ----------
    ##  samples
    bulk_samples_ofh = open("bulk.{}.samples.tsv".format(sqlite3_dbname), 'wt')
    # fields: sample_name, db_class, sample_type, total_uniq_count, total_multi_count, total_count
    samples = dict()

    
    ## --------------
    ## intron feature
    bulk_intron_feature_ofh = open("bulk.{}.intron_feature.tsv".format(sqlite3_dbname), 'wt')
    # fields: intron, chromosome, start, end, strand, intron_motif, annot_status, genes
    intron_features = set()

    ## -----------------
    ## intron_occurrence
    bulk_intron_occurrence_ofh = open("bulk.{}.intron_occurrence.tsv".format(sqlite3_dbname), 'wt')
    # fields: intron, sample, unique_mappings, multi_mappings, all_mappings, max_spliced_align_overhang, norm_unique_mappings, norm_multi_mappings, norm_all_mappings
    
    
    ## populate data
    for input_file in input_files:
        logger.info("-processing file: " + input_file)
        with open(input_file, 'rt') as fh:
            header = next(fh)
            counter = 0
            for line in fh:

                counter += 1
                if counter % 1000 == 0:
                    sys.stderr.write("\r[{}]  ".format(counter))
                
                line = line.rstrip()
                (classname, sample, genes,
                 Chromosome, Start, End, strandval,
                 intron_motif, annot_status,
                 unique_mappings, multi_mappings, max_spliced_align_overhang) = line.split("\t")
                
                unique_mappings = int(unique_mappings)
                multi_mappings = int(multi_mappings)


                ## sample:     # fields: sample_name, db_class, sample_type, total_uniq_count, total_multi_count, total_count
                if sample not in samples:

                    ## determine if tumor or normal
                    TN = 'N'

                    if classname == "TCGA":
                        TN = 'T'
                        
                        if sample[-3:] == "-NT":
                            TN = 'N'
                        
                        sample_type = sample.split("-")[0]

                    elif classname == "GTEx":

                        assert(sample in GTEx_sample_to_tissue_type)
                        sample_type = GTEx_sample_to_tissue_type[sample]

                    else:
                        raise RuntimeError("Error, not recognizing class type: {}".format(classname))
                                           
                                        
                    sample_struct = samples[sample] = { 'db_class' : classname,
                                                        'sample_type' : sample_type,
                                                        'total_uniq_count' : 0,
                                                        'total_multi_count' : 0,
                                                        'total_count' : 0,
                                                        'TN' : TN}
                else:
                    sample_struct = samples[sample]

                sample_struct['total_uniq_count'] += unique_mappings
                sample_struct['total_multi_count'] += multi_mappings
                sample_struct['total_count'] += unique_mappings + multi_mappings


                ## intron feature: fields: intron, chromosome, start, end, strand, intron_motif, annot_status, genes
                intron_feature_key = "{}:{}-{}".format(Chromosome, Start, End)
                if intron_feature_key not in intron_features:

                    bulk_intron_feature_ofh.write("\t".join([intron_feature_key,
                                                             Chromosome,
                                                             Start,
                                                             End,
                                                             strandval,
                                                             intron_motif,
                                                             annot_status,
                                                             genes ]) + "\n")
                    intron_features.add(intron_feature_key)
                    

                ## intron occurrence: fields: intron, sample, unique_mappings, multi_mappings, all_mappings,
                ##                    max_spliced_align_overhang, norm_unique_mappings, norm_multi_mappings, norm_all_mappings 

                intron_occurrence = { 'intron' : intron_feature_key,
                                      'sample' : sample,
                                      'unique_mappings' : unique_mappings,
                                      'multi_mappings' : multi_mappings,
                                      'all_mappings' : unique_mappings + multi_mappings,
                                      'max_spliced_align_overhang' : max_spliced_align_overhang }

                bulk_intron_occurrence_ofh.write("\t".join([intron_feature_key,
                                                            sample,
                                                            str(unique_mappings),
                                                            str(multi_mappings),
                                                            str(unique_mappings + multi_mappings),
                                                            max_spliced_align_overhang
                                                            ]) + "\n")
        
        sys.stderr.write("  ok\n")


    ## write samples table data
    for sample, sample_struct in samples.items():
        bulk_samples_ofh.write( "\t".join([sample,
                                           sample_struct['db_class'],
                                           sample_struct['sample_type'],
                                           str(sample_struct['total_uniq_count']),
                                           str(sample_struct['total_multi_count']),
                                           str(sample_struct['total_count']),
                                           sample_struct['TN']] ) + "\n")
        

    logger.info("Writing bulk load files")
    
    sys.exit(0)




def parse_GTEx_sample_types():

    gtex_sample_to_tissue = dict()

    gtex_sample_file = os.path.join(os.path.dirname(__file__), "gtex_sample_info.tsv")
    with open(gtex_sample_file, 'rt') as fh:
        header = next(fh)
        for line in fh:
            line = line.rstrip()
            (sample_name, tissue_type) = line.split("\t")
            gtex_sample_to_tissue[sample_name] = tissue_type

    return gtex_sample_to_tissue


if __name__ == '__main__':
    main()

    
