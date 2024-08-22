#!/usr/bin/env python

import sys, os, re
import sqlite3
import argparse

def main():

    parser = argparse.ArgumentParser(description="generate data file for cancer intron annotation",
                                      formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument("--sqlite3_db", dest="sqlite3_db", type=str, required=True, help="sqlite3_db name")
    
    
    args = parser.parse_args()

    sqlite3_dbname = args.sqlite3_db


    conn = sqlite3.connect(sqlite3_dbname)
    c = conn.cursor()

    # print header:
    print("\t".join(["intron", "genes", "TCGA_sample_counts", "GTEx_sample_counts"]))

    counter = 0
    
    query = "select intron, genes from intron_feature"
    c.execute(query)
    rows = c.fetchall()
    for row in rows:
        intron_feature, genes = row

        write_intron_feature_annotation(c, intron_feature, genes)
        counter += 1
        if counter % 100 == 0:
            sys.stderr.write("\r[{}]  ".format(counter))

    sys.stderr.write("\nDone parsing.\n")

    sys.exit(0)


def write_intron_feature_annotation(c : sqlite3.Cursor, intron_feature : str, genes : str) -> None:
    

    query = str("select sample_type, all_map_sample_count, all_map_sample_pct " +
                " from intron_sample_type_counts " +
                " where intron = ? and db_class = ? and sample_type != 'total' " +
                " order by all_map_sample_pct desc")

    TCGA_vals = list()
    c.execute(query, (intron_feature, "TCGA") )
    rows = c.fetchall()
    for row in rows:
        (sample_type, sample_count, sample_fraction) = row
        sample_info = ":".join([sample_type, str(sample_count), "{:.2f}".format(sample_fraction*100)])
        TCGA_vals.append(sample_info)
    if not TCGA_vals:
        TCGA_vals = ["NA"]
        
    
    GTEx_vals = list()
    c.execute(query, (intron_feature, "GTEx") )
    rows = c.fetchall()
    for row in rows:
        (sample_type, sample_count, sample_fraction) = row
        sample_info = ":".join([sample_type, str(sample_count), "{:.2f}".format(sample_fraction*100)])
        GTEx_vals.append(sample_info)
    if not GTEx_vals:
        GTEx_vals = ["NA"]


    print("\t".join([intron_feature, genes, ",".join(TCGA_vals), ",".join(GTEx_vals)]))
    
    return
        
    


if __name__=='__main__':
    main()
