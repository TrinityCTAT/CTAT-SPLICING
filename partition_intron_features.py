#!/usr/bin/env python

import sys, os, re


def main():

    usage = "\n\n\tusage: {} intron_features.list.file  num_features_per_interval output_directory files_per_bin\n\n\n".format(sys.argv[0]) 
    
    if len(sys.argv) < 5:
        sys.stderr.write(usage)
        sys.exit(1)

    intron_features_list_file = sys.argv[1]
    num_features_per_interval = int(sys.argv[2])
    output_directory = sys.argv[3]
    files_per_bin = int(sys.argv[4])

    if not os.path.exists(output_directory):
        os.makedirs(output_directory)


    files_per_bin = 100

    feature_counter = 0
    file_counter = 0
    
    ofh = None
    with open(intron_features_list_file) as fh:
        for line in fh:
            if feature_counter % num_features_per_interval == 0:
                ## make a new file
                if ofh is not None:
                    ofh.close()
                new_interval_bin = int(file_counter / files_per_bin)
                new_interval_outdir = os.path.join(output_directory, "ibin_{}".format(new_interval_bin))
                if not os.path.exists(new_interval_outdir):
                    os.makedirs(new_interval_outdir)
                new_interval_filename = os.path.join(new_interval_outdir, "introns.{}.txt".format(feature_counter))
                ofh = open(new_interval_filename, 'wt')
                sys.stderr.write("-writing to {}\n".format(new_interval_filename))
                file_counter += 1
                
            ofh.write(line)
            
            feature_counter += 1
            
    sys.stderr.write("Done.\n")
    
    sys.exit(0)


if __name__ == '__main__':
    main()




