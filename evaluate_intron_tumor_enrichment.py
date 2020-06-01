#!/usr/bin/env python

import sys, os, re
import typing
import sqlite3
import logging
import collections
from collections import defaultdict
import argparse
import time
import scipy.stats as stats
import statistics

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def main():

    parser = argparse.ArgumentParser(description="examines intron feature for tumor type enrichment", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("--sqlite3_db", dest="sqlite3_db", type=str, required=True, help="sqlite3_db name")
    parser.add_argument("--intron_features_file", dest="intron_features_file", type=str, required=True, help="file containing list of intron features")
    parser.add_argument("--output_file", dest="output_file", type=str, required=True, help="output filename")
    parser.add_argument("--pseudocount", dest="pseudocount", type=float, required=False, default=0.1, help="pseudocount used in ratio computations: (A + pseudo)/(B + pseudo)") 
    
    args = parser.parse_args()
    sqlite3_dbname = args.sqlite3_db
    intron_features_file = args.intron_features_file
    output_filename = args.output_file
    pseudocount = args.pseudocount

    conn = sqlite3.connect(sqlite3_dbname)
    c = conn.cursor()


    ## get counts of samples according to tissue type
    query = "select db_class, sample_type, count(*) from samples group by db_class, sample_type"
    c.execute(query)
    rows = c.fetchall()
    sample_type_counts = defaultdict(int)
    gtex_counts = list()
    tcga_counts = list()
    for row in rows:
        (db_class, sample_type, count) = row
        sample_type = "^".join([db_class, sample_type])
        sample_type_counts[sample_type] += count
        base_sample_type = "^".join([db_class, "ALL"])
        sample_type_counts[base_sample_type] += count

        if db_class == "TCGA":
            tcga_counts.append(count)

        if db_class == "GTEx":
            gtex_counts.append(count)

    mean_tcga_count = round(statistics.mean(tcga_counts))
    mean_gtex_count = round(statistics.mean(gtex_counts))
    

    ofh = open(output_filename, 'wt')

    with open(intron_features_file) as fh:
        for intron_feature in fh:
            intron_feature = intron_feature.rstrip()
            logger.info(intron_feature)
            start_time = time.time()
            examine_intron_feature_for_enrichment(intron_feature, c, sample_type_counts, ofh, pseudocount, mean_tcga_count, mean_gtex_count)
            end_time = time.time()
            seconds = int(end_time - start_time)
            logger.info("-took {} seconds".format(seconds))

    
    ofh.close()

    sys.exit(0)



def examine_intron_feature_for_enrichment(intron_feature : str,
                                          c : sqlite3.Cursor,
                                          sample_type_counts : collections.defaultdict,
                                          ofh : typing.TextIO,
                                          pseudocount : int,
                                          mean_tcga_count : int,
                                          mean_gtex_count : int) -> None:

    query = str("select db_class, sample_type, uniq_count, uniq_pct, all_count, all_pct " +
                " from intron_sample_type_counts " +
                " where intron = ? ")

    c.execute(query, (intron_feature,))
    
    rows = c.fetchall()

    class db_sample_obj:

        def __init__(self):
            # null settings by default
            self.db_class = "NONE"
            self.sample_type = "NONE"
            self.uniq_count = 0
            self.uniq_pct = 0.0
            self.all_count = 0
            self.all_pct = 0.0
        


    tcga_all_obj = db_sample_obj()
    gtex_all_obj = db_sample_obj()

    tcga_top_obj = db_sample_obj()
    gtex_top_obj = db_sample_obj()
    

    for row in rows:
        (db_class, sample_type, uniq_count, uniq_pct, all_count, all_pct) = row        

        obj = db_sample_obj()
        obj.db_class = db_class
        obj.sample_type = sample_type
        obj.uniq_count = uniq_count
        obj.uniq_pct = uniq_pct
        obj.all_count = all_count
        obj.all_pct = all_pct

        assert(obj.db_class in ("TCGA", "GTEx", "NONE"))

        if obj.sample_type == "ALL":
            if obj.db_class == "TCGA":
                tcga_all_obj = obj
            elif obj.db_class == "GTEx":
                gtex_all_obj = obj

        elif obj.db_class == "TCGA":
            if tcga_top_obj.db_class == "NONE" or tcga_top_obj.all_pct < obj.all_pct:
                tcga_top_obj = obj

        elif obj.db_class == "GTEx":
            if gtex_top_obj.db_class == "NONE" or gtex_top_obj.all_pct < obj.all_pct:
                gtex_top_obj = obj





    ## examine tumor enrichment stats

    # examine enrichment based on ALL (cumulative) samples.

    count_tcga_all = sample_type_counts["TCGA^ALL"]
    count_gtex_all = sample_type_counts["GTEx^ALL"]

    count_tcga_top = sample_type_counts["TCGA^{}".format(tcga_top_obj.sample_type)]
    if count_tcga_top == 0:
        count_tcga_top = mean_tcga_count

    count_gtex_top = sample_type_counts["GTEx^{}".format(gtex_top_obj.sample_type)]
    if count_gtex_top == 0:
        count_gtex_top = mean_gtex_count
    

    #################
    ## All comparison
        
    tcga_all_enrichment = ( (tcga_all_obj.all_count + pseudocount) / (count_tcga_all + pseudocount) ) / ( (gtex_all_obj.all_count + pseudocount) / (count_gtex_all + pseudocount) )


    oddsratio, pvalue = stats.fisher_exact([[tcga_all_obj.all_count, count_tcga_all - tcga_all_obj.all_count],
                                            [gtex_all_obj.all_count, count_gtex_all - gtex_all_obj.all_count]],
                                           alternative='greater')

    # write sql for ALL comparison
    sql = str("insert into tumor_vs_normal (intron, tumor_sample_type, normal_sample_type, tumor_yes, tumor_no, normal_yes, normal_no, odds_ratio) " +
              "values (\"{}\",\"ALL\", \"ALL\",".format(intron_feature) +
              "{},".format(tcga_all_obj.all_count) +
              "{},".format(count_tcga_all - tcga_all_obj.all_count) + 
              "{},".format(gtex_all_obj.all_count) +
              "{},".format(count_gtex_all - gtex_all_obj.all_count) +
              "{:.4}, ".format(tcga_all_enrichment) + 
              "{:.4});".format(pvalue) ) 
    print(sql, file=ofh)


    #################
    ## top comparison

    tcga_top_enrichment = ( (tcga_top_obj.all_count + pseudocount) / (count_tcga_top + pseudocount) ) / ( (gtex_top_obj.all_count + pseudocount) / (count_gtex_top + pseudocount) )


    oddsratio, pvalue = stats.fisher_exact([[tcga_top_obj.all_count, count_tcga_top - tcga_top_obj.all_count],
                                            [gtex_top_obj.all_count, count_gtex_top - gtex_top_obj.all_count] ],
                                           alternative='greater')

    # write sql for top entry comparisons
    sql = str("insert into tumor_vs_normal (intron, tumor_sample_type, normal_sample_type, tumor_yes, tumor_no, normal_yes, normal_no, odds_ratio) " +
              "values (\"{}\", \"{}\", \"{}\", ".format(intron_feature, tcga_top_obj.sample_type, gtex_top_obj.sample_type) +
              "{},".format(tcga_top_obj.all_count) +
              "{},".format(count_tcga_top - tcga_top_obj.all_count) + 
              "{},".format(gtex_top_obj.all_count) +
              "{},".format(count_gtex_top - gtex_top_obj.all_count) +
              "{:.4}, ".format(tcga_top_enrichment) +
              "{:.4});".format(pvalue))
    print(sql, file=ofh)
    
    return
    
if __name__ == '__main__':
    main()






