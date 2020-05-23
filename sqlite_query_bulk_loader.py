#!/usr/bin/env python

import sys, os, re
import sqlite3
import argparse
import subprocess

def main():

    parser = argparse.ArgumentParser(description="bulk database loader", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--sqlite3_db", dest="sqlite3_db", type=str, required=True, help="sqlite3_db name")
    parser.add_argument("--sql_file", dest="sql_file", type=str, required=True, help="file containing sql queries")
    parser.add_argument("--transaction_size", dest="transaction_size", type=int, required=True, help="transaction size")

    args = parser.parse_args()

    sqlite3_dbname = args.sqlite3_db
    sql_file = args.sql_file
    transaction_size = int(args.transaction_size)

    sys.stderr.write("-examining count of sql queries...\n")
    dbsize = subprocess.check_output("wc -l {}".format(sql_file), shell=True).decode().split(" ")[0]

    print("dbsize: {}".format(dbsize))
    dbsize = int(dbsize)
    
    conn = sqlite3.connect(sqlite3_dbname)
    c = conn.cursor()


    with open(sql_file) as fh:

        counter = 0
        for line in fh:
            line = line.rstrip()
            query = re.sub(";$", "", line)
            c.execute(query)
            counter += 1
            if counter % transaction_size == 0:
                conn.commit()
                pct_done = (counter / dbsize * 100)
                sys.stderr.write("\r[{} = {:3f}%]   ".format(counter, pct_done))


        # get last ones
        if counter % transaction_size != 0:
            conn.commit()
            sys.stderr.write("\r[{} = 100%]  ".format(counter))


    conn.close()
    
    sys.stderr.write("\n\nDone.\n")
    
    sys.exit(0)



if __name__ == '__main__':
    main()
