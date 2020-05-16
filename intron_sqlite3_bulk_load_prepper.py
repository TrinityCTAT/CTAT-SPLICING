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
                    sample_struct = samples[sample] = { 'db_class' : classname,
                                                        'total_uniq_count' : 0,
                                                        'total_multi_count' : 0,
                                                        'total_count' : 0 }
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
                                                            max_spliced_align_overhang,
                                                            "-1", # norm unique mappings
                                                            "-1", # norm multi mappings
                                                            "-1" # norm all mappings
                                                            ]) + "\n")
        
        sys.stderr.write("  ok\n")


    ## write samples table data
    for sample, sample_struct in samples.items():
        bulk_samples_ofh.write( "\t".join([sample,
                                           sample_struct['db_class'],
                                           "NA",
                                           str(sample_struct['total_uniq_count']),
                                           str(sample_struct['total_multi_count']),
                                           str(sample_struct['total_count']) ] ) + "\n")
        

    logger.info("Writing bulk load files")
    
    sys.exit(0)



def build_tables():
    
    conn = sqlite3.connect(sqlite3_dbname)

    c = conn.cursor()

    ## samples table:
    c.execute("CREATE TABLE samples \
                 (sample_name TEXT, \
                  db_class TEXT, \
                  sample_type TEXT, \
                  total_uniq_count INT, \
                  total_multi_count INT, \
                  total_count INT)")

    
    ## intron_feature table
    c.execute("CREATE TABLE intron_feature \
                 (intron TEXT, \
                  chromosome TEXT, \
                  start INT, \
                  end INT, \
                  strand INT, \
                  intron_motif INT, \
                  annot_status INT, \
                  genes TEXT)")

    ## intron_occurrence table
    c.execute("CREATE TABLE intron_occurrence \
                  (intron TEXT, \
                   sample TEXT, \
                   unique_mappings INT, \
                   multi_mappings INT, \
                   all_mappings INT, \
                   max_spliced_align_overhang INT, \
                   norm_unique_mappings REAL, \
                   norm_multi_mappings REAL, \
                   norm_all_mappings REAL) ")
    
    conn.commit()

    


if __name__ == '__main__':
    main()

    
