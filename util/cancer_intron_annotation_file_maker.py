#!/usr/bin/env python

import sys, os, re
import sqlite3
import argparse

def main():

    parser = argparse.ArgumentParser(description="generate data file for cancer intron annotation",
                                      formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument("--sqlite3_db", dest="sqlite3_db", type=str, required=True, help="sqlite3_db name")
    parser.add_argument("--cancer_introns", dest="cancer_introns", type=str, required=True, help="selected cancer introns file")
    parser.add_argument("--intron_feature_names", dest="intron_feature_names", type=str, required=True, help="intron feature names file")
    
    
    args = parser.parse_args()

    sqlite3_dbname = args.sqlite3_db
    cancer_introns_file = args.cancer_introns
    intron_feature_names_file = args.intron_feature_names

    intron_feature_names = parse_intron_feature_names(intron_feature_names_file)


    conn = sqlite3.connect(sqlite3_dbname)
    c = conn.cursor()

    counter = 0
    with open(cancer_introns_file) as fh:
        header = next(fh)
        for line in fh:
            vals = line.split("\t")
            genes = vals[1]
            intron_feature = vals[2]
            
            intron_feature_name = "NA"
            if intron_feature in intron_feature_names:
                intron_feature_name = intron_feature_names[intron_feature]
                del intron_feature_names[intron_feature]

            write_intron_feature_annotation(c, intron_feature, intron_feature_name)
            counter += 1
            if counter % 100 == 0:
                sys.stderr.write("\r[{}]  ".format(counter))

    sys.stderr.write("\nDone parsing.\n")
    
    # capture any known/annotated introns not defined as cancer introns based on empirical data
    for intron in intron_feature_names:
        write_intron_feature_annotation(c, intron, intron_feature_names.get(intron, "NA"))

    
    sys.exit(0)



def parse_intron_feature_names(intron_feature_names_file : str) -> dict:

    intron_names_dict = dict()

    with open(intron_feature_names_file) as fh:
        for line in fh:
            line = line.rstrip()
            (intron_feature, intron_name) = line.split("\t")
            intron_names_dict[intron_feature] = intron_name

    return intron_names_dict



def write_intron_feature_annotation(c : sqlite3.Cursor, intron_feature : str, intron_feature_name : str) -> None:
    

    query = "select genes from intron_feature where intron = ?"
    c.execute(query, (intron_feature,))
    genes, = c.fetchone()


    query = str("select sample_type, all_count, all_pct " +
                " from intron_sample_type_counts " +
                " where intron = ? and db_class = ? and sample_type != 'ALL' " +
                " order by all_pct desc")

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


    print("\t".join([intron_feature, genes, ",".join(TCGA_vals), ",".join(GTEx_vals), intron_feature_name]))
    
    return
        
    


if __name__=='__main__':
    main()
