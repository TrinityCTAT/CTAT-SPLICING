#!/usr/bin/env python

import sys, os, re
import typing
import sqlite3
import logging
import collections
from collections import defaultdict
import argparse
import time

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

"""
## -----------------------
## intron_occurrence table

sqlite> select * from intron_occurrence limit 1;
                    intron = chr1:15039-15795
                    sample = GTEX-O5YT-1026-SM-3MJGF
           unique_mappings = 3
            multi_mappings = 64
              all_mappings = 67
max_spliced_align_overhang = 36
      norm_unique_mappings = -1.0
       norm_multi_mappings = -1.0
         norm_all_mappings = -1.0


## --------------------
## intron_feature table

sqlite> select * from intron_feature limit 1;
      intron = chr1:15039-15795
  chromosome = chr1
       start = 15039
         end = 15795
      strand = 2
intron_motif = 2
annot_status = 1
       genes = WASH7P^ENSG00000227232.5



## --------------
### samples table

sqlite> select * from samples limit 1;
      sample_name = GTEX-O5YT-1026-SM-3MJGF
         db_class = GTEx
      sample_type = NA
 total_uniq_count = 21408796
total_multi_count = 4092132
      total_count = 25500928

"""


def main():

    parser = argparse.ArgumentParser(description="examines intron feature for tumor type enrichment", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("--sqlite3_db", dest="sqlite3_db", type=str, required=True, help="sqlite3_db name")
    parser.add_argument("--intron_features_file", dest="intron_features_file", type=str, required=True, help="file containing list of intron features")
    parser.add_argument("--output_file", dest="output_file", type=str, required=True, help="output filename")

    args = parser.parse_args()
    sqlite3_dbname = args.sqlite3_db
    intron_features_file = args.intron_features_file
    output_filename = args.output_file

    conn = sqlite3.connect(sqlite3_dbname)
    c = conn.cursor()


    ofh = open(output_filename, 'wt')

    ## get counts of samples according to tissue type
    query = "select db_class, sample_type, count(*) from samples group by db_class, sample_type"
    c.execute(query)
    rows = c.fetchall()
    sample_type_counts = defaultdict(int)
    for row in rows:
        (db_class, sample_type, count) = row
        sample_type = "^".join([db_class, sample_type])
        sample_type_counts[sample_type] += count
        base_sample_type = "^".join([db_class, "ALL"])
        sample_type_counts[base_sample_type] += count
        

    with open(intron_features_file) as fh:
        for intron_feature in fh:
            intron_feature = intron_feature.rstrip()
            logger.info(intron_feature)
            start_time = time.time()
            examine_intron_feature_for_enrichment(intron_feature, c, sample_type_counts, ofh)
            end_time = time.time()
            seconds = int(end_time - start_time)
            logger.info("-took {} seconds".format(seconds))

    
    ofh.close()

    sys.exit(0)



def examine_intron_feature_for_enrichment(intron_feature : str,
                                          c : sqlite3.Cursor,
                                          sample_type_counts : collections.defaultdict,
                                          ofh : typing.TextIO):
    
    query = str("select s.sample_name, s.db_class, s.sample_type, s.total_uniq_count, s.total_multi_count, s.total_count, "
                + " io.intron, io.unique_mappings, io.multi_mappings, io.all_mappings "
                + " from samples as s, intron_occurrence as io "
                + " where s.sample_name = io.sample "
                + "       and io.intron = ? " #\"{}\" ".format(intron_feature)
                + "       and not (db_class = \"TCGA\" and TN = \"N\") ")  # skip the tumor normals for now... analyze separately.
                

    #logger.info(query)

    #sys.exit(1)

    
    c.execute(query, (intron_feature,))

    rows = c.fetchall()


    sample_type_counter = dict()
    def add_to_counter(sample_grp_type : str, category : str):
        if sample_grp_type not in sample_type_counter:
            sample_type_counter[sample_grp_type] = defaultdict(int)
        sample_type_counter[sample_grp_type][category] += 1
        return
    

    for row in rows:
        (sample_name, db_class, sample_type, sample_total_uniq_count, sample_total_multi_count, sample_total_count,
         intron_token, intron_unique_mappings, intron_multi_mappings, intron_all_mappings) = row


        ## simple counting intron occurrences according to sample type
        specific_sample_type = "^".join([db_class, sample_type])
        base_sample_type = "^".join([db_class, "ALL"])

        if intron_unique_mappings != "0":
            add_to_counter(specific_sample_type, "intron_unique_mappings")
            add_to_counter(base_sample_type, "intron_unique_mappings")

        if intron_multi_mappings != "0":
            add_to_counter(specific_sample_type, "intron_multi_mappings")
            add_to_counter(base_sample_type, "intron_multi_mappings")

        if intron_all_mappings != "0":
            add_to_counter(specific_sample_type, "intron_all_mappings")
            add_to_counter(base_sample_type, "intron_all_mappings")


        ## generate normalized counts
        norm_unique_mappings = int(intron_unique_mappings) / int(sample_total_uniq_count) * 1e6
        norm_multi_mapings = int(intron_multi_mappings) / int(sample_total_multi_count) * 1e6
        norm_all_mappings = int(intron_all_mappings) / int(sample_total_count) * 1e6
        
        print("\t".join(["INTRON_NORM_VALS", intron_token, sample_name, "{:.4f}".format(norm_unique_mappings),
                         "{:.4f}".format(norm_multi_mapings), "{:.4f}".format(norm_all_mappings)]))
        
    
    ## report on findings.
    categories = ("intron_unique_mappings", "intron_multi_mappings", "intron_all_mappings")
    for sample_grp_type, counter_dict in sample_type_counter.items():
        output = ["INTRON_SAMPLE_TYPE_COUNTER", sample_grp_type]
        num_sample_counts = sample_type_counts[sample_grp_type]
        
        for cat in categories:
            num_cat = counter_dict.get(cat, 0)
            frac_cat = num_cat / num_sample_counts
            output += [str(num_cat), "{:.4f}".format(frac_cat)]
            
        print("\t".join(output))
    
    
    

if __name__ == '__main__':
    main()






