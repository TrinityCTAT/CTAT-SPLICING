#!/usr/bin/env python

import sys, os, re


def main():

    usage = "\n\n\tusage: {} intron.stats\n\n".format(sys.argv[0])
    if len(sys.argv) < 2:
        sys.stderr.write(usage)
        sys.exit(1)

    stats_file = sys.argv[1]

    with open(stats_file) as fh:
        for line in fh:
            line = line.rstrip()
            vals = line.split("\t")
            generate_sql(vals)
    
    sys.exit(0)


def generate_sql(vals : list) -> None:

    if vals[0] == "INTRON_SAMPLE_TYPE_COUNTER":
        (table_type, db_class_tissue_type,
         uniq_count, uniq_pct,
         multi_count, multi_pct,
         all_count, all_pct) = vals

        db_class, tissue_type = db_class_tissue_type.split("^")
        
        print(str("insert intron_sample_type_counts (db_class, sample_type, " +
                  " uniq_count, uniq_pct, all_count, all_pct) " +
                  " values (\"{}\", \"{}\", {}, {}, {}, {}, {} )".format(
                      db_class, tissue_type, uniq_count, uniq_pct,
                      multi_count, multi_pct, all_count, all_pct) ) )

    elif vals[0] == "INTRON_NORM_VALS":
        # intron_occurrence  (intron TEXT,  sample TEXT,  unique_mappings INT, multi_mappings INT,all_mappings INT,max_spliced_align_overhang INT,norm_unique_mappings REAL,norm_multi_mappings REAL, norm_all_mappings REAL);
        (table_type, intron, sample, norm_unique_mappings, norm_multi_mappings, norm_all_mappings) = vals

        print(str("update intron_occurrence " +
                  " set norm_unique_mappings = {}, ".format(norm_unique_mappings) +
                  "     norm_multi_mappings = {}, ".format(norm_multi_mappings) +
                  "     norm_all_mappings = {} ".format(norm_all_mappings) +
                  " where intron = \"{}\" ".format(intron) +
                  "       and sample_type = \"{}\" ".format(sample) ) )


    return



if __name__ == '__main__':
    main()
