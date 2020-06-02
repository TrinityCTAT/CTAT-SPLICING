#!/usr/bin/env python

import sys, os, re

INTRON = None  ## //FIXME: put intron into the file being read for each line.


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

    global INTRON
    
    if vals[0] == "INTRON_SAMPLE_TYPE_COUNTER":
        (table_type, db_class_tissue_type,
         uniq_count, uniq_pct,
         multi_count, multi_pct,
         all_count, all_pct) = vals

        db_class, tissue_type = db_class_tissue_type.split("^")

        if INTRON is None:
            raise RuntimeError("Error, intron not identified!")
        
        print(str("insert into intron_sample_type_counts (intron, db_class, sample_type, " +
                  " uniq_count, uniq_pct, multi_count, multi_pct, all_count, all_pct) " +
                  " values (\"{}\", \"{}\", \"{}\", ".format( INTRON, db_class, tissue_type) +
                  " {}, {}, {}, {}, {}, {} );".format(uniq_count, uniq_pct, multi_count, multi_pct, all_count, all_pct) ) )


    elif vals[0] == "INTRON_NORM_VALS":
        # intron_occurrence  (intron TEXT,  sample TEXT,  unique_mappings INT, multi_mappings INT,all_mappings INT,max_spliced_align_overhang INT,norm_unique_mappings REAL,norm_multi_mappings REAL, norm_all_mappings REAL);
        (table_type, intron, sample, norm_unique_mappings, norm_multi_mappings, norm_all_mappings) = vals

        print(str("update intron_occurrence " +
                  " set norm_unique_mappings = {}, ".format(norm_unique_mappings) +
                  "     norm_multi_mappings = {}, ".format(norm_multi_mappings) +
                  "     norm_all_mappings = {} ".format(norm_all_mappings) +
                  " where intron = \"{}\" ".format(intron) +
                  "       and sample = \"{}\"; ".format(sample) ) )

        
        INTRON = intron
        
    return



if __name__ == '__main__':
    main()
