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

    if os.path.exists(sqlite3_dbname):
        raise RuntimeError("Error, database {} already exists. Remove it or use a different database name for creation step.".format(sqlite3_dbname))



    ## populate data
    for input_file in input_files:
        with open(input_file, 'rt') as fh:
            header = next(fh)
            counter = 0
            for line in fh:
                line = line.rstrip()
                vals = line.split("\t")

                query = "INSERT INTO introns VALUES (?,?,?,?,?,?,?,?,?,?,?,?)"
                c.execute(query, vals)

                counter += 1
                
                if counter % COMMIT_EVERY == 0:
                    conn.commit()
                    logger.info("[{}]".format(counter))
                    
            # commit last ones.
            if counter % COMMIT_EVERY != 0:
                conn.commit()

        logger.info("-done loading in {}".format(input_file))


    sys.exit(0)




def build_tables():
    
    conn = sqlite3.connect(sqlite3_dbname)

    c = conn.cursor()

    ## samples table:
    c.execute("CREATE TABLE samples \
                 (sample_name TEXT,
                  db_class TEXT,
                  sample_type TEXT,
                  total_uniq_count INT,
                  total_multi_count INT,
                  total_count INT)")

    
    ## intron_feature table
    c.execute("CREATE TABLE intron_feature \
                 (intron TEXT, \
                  chromosome TEXT,
                  start INT,
                  end INT,
                  strand INT,
                  intron_motif INT,
                  annot_status INT,
                  genes TEXT)")

    ## intron_occurrence table
    c.execute("CREATE TABLE intron_occurrence \
                  (intron TEXT,
                   sample TEXT,
                   unique_mappings INT, \
                   multi_mappings INT, \
                   all_mappings INT,
                   max_spliced_align_overhang INT, \
                   norm_unique_mappings REAL,
                   norm_multi_mappings REAL,
                   norm_all_mappings REAL) ")
    
    conn.commit()

    


if __name__ == '__main__':
    main()

    
