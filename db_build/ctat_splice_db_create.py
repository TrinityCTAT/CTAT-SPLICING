#!/usr/bin/env python

import sys, os, re
import sqlite3
import argparse
import subprocess
import logging


logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)





def main():

    parser = argparse.ArgumentParser(description="ctat intron database builder", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument("--sqlite3_db", dest="sqlite3_db", type=str, required=True, help="sqlite3_db name")

    parser.add_argument("--create", dest='create', action='store_true', default=False, help="create a new database")

    parser.add_argument("--index", dest='index', type=str, required=False, nargs='+', help="tables to index")

    args = parser.parse_args()

    sqlite3_dbname = args.sqlite3_db

    
    if not (args.create or args.index):
        sys.stderr.write("\n\n\tMust select --create or --index <tablename>  ..... nothing to do here.\n\n")
        sys.exit(1)


    if args.create and os.path.exists(sqlite3_dbname):
        sys.stderr.write("\n\n\tDatabase file {} already exists. Please rename it or provide a different database anme via --sqlite3_db\n\n".format(sqlite3_dbname))
        sys.exit(1)


    conn = sqlite3.connect(sqlite3_dbname)
    c = conn.cursor()


    if args.create:

        logger.info("-creating database: {}".format(sqlite3_dbname))
                
        c.execute("CREATE TABLE samples " +
                  " (sample_name TEXT, " +
                  "  db_class TEXT, " +
                  "  sample_type TEXT, " +
                  "  total_uniq_count INT, " +
                  "  total_multi_count INT, " +
                  "  total_count INT, " +
                  "  TN char)")

        
        c.execute("CREATE TABLE intron_feature " +
                  " (intron TEXT, " +
                  "  chromosome TEXT, " +
                  "  start INT, " +
                  "  end INT, " +
                  "  strand INT, " +
                  "  intron_motif INT, " +
                  "  annot_status INT, " +
                  "  genes TEXT)")
        
        c.execute("CREATE TABLE intron_sample_type_counts " +
                  " (intron TEXT, " +
                  "  db_class TEXT, " +
                  "  sample_type TEXT, " +
                  "  uniq_count INT, " +
                  "  uniq_pct REAL, " +
                  "  multi_count INT, " +
                  "  multi_pct REAL, " +
                  "  all_count INT, " +
                  "  all_pct REAL)")
    
        c.execute("CREATE TABLE intron_occurrence " +
                  " (intron TEXT, " +
                  "  sample TEXT, " +
                  "  unique_mappings INT, " +
                  "  multi_mappings INT, " +
                  "  all_mappings INT, " +
                  "  max_spliced_align_overhang INT)")


        c.execute("CREATE TABLE tumor_vs_normal " +
                  " (intron TEXT, " +
                  "  tumor_sample_type TEXT, " +
                  "  normal_sample_type TEXT, " +
                  "  tumor_yes INT, " +
                  "  tumor_no INT, " +
                  "  normal_yes INT, " +
                  "  normal_no INT, "
                  "  enrichment REAL, " +
                  "  odds_ratio REAL, " +
                  "  pvalue REAL)")


        conn.commit()

        ## done database creation


    ## table indexing
    if args.index:

        for tablename in args.index:

            if tablename == 'samples':

                logger.info("-indexing table: samples")
                
                c.execute("CREATE UNIQUE INDEX samples_table_idx_sample_name ON samples(sample_name)")
                conn.commit()

            if tablename == 'intron_feature':

                logger.info("-indexing table: intron_feature")

                c.execute("CREATE UNIQUE INDEX intron_feature_idx_intron ON intron_feature (intron)")
                conn.commit()

            if tablename == 'intron_occurrence':

                logger.info("-indexing table: intron_occurrence")

                c.execute("CREATE UNIQUE INDEX intron_occurrence_idx_intron_sample ON intron_occurrence(intron, sample)")

                c.execute("CREATE INDEX intron_occurrence_idx_intron ON intron_occurrence(intron)")


            if tablename == 'intron_sample_type_counts':

                logger.info("-indexing table: intron_sample_type_counts")

                c.execute("CREATE UNIQUE INDEX intron_sample_type_counts_idx_key ON intron_sample_type_counts(intron, db_class, sample_type)")
                
                c.execute("CREATE INDEX intron_sample_type_counts_idx_db_class_sample_type ON intron_sample_type_counts(db_class, sample_type)")
                
                c.execute("CREATE INDEX intron_sample_type_counts_idx_intron ON intron_sample_type_counts(intron)")
                    
                conn.commit()

            
            if tablename == "tumor_vs_normal":

                logger.info("-indexing table: tumor_vs_normal")
                
                c.execute("CREATE UNIQUE INDEX tumor_vs_normal_idx_intron_tumor_normal ON tumor_vs_normal (intron, tumor_sample_type, normal_sample_type)")
                
                c.execute("CREATE INDEX tumor_vs_normal_idx_intron ON tumor_vs_normal(intron)")

                conn.commit()


    logger.info("-done")
    
    sys.exit(0)



if __name__=='__main__':
    main()

    
