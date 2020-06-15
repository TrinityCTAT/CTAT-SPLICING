#!/usr/bin/env python




#~~~~~~~~~~~~~~~~~~~~~~~~~~
# Import packages 
#~~~~~~~~~~~~~~~~~~~~~~~~~~
from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
import os
import sys
import json
import math
import argparse
import logging
from urllib.request import urlopen
from igv_reports import fasta, ideogram, datauri, tracks, feature, bam, vcf, utils
from igv_reports.varianttable import VariantTable
from igv_reports.bedtable import BedTable

logging.basicConfig(level=logging.INFO, 
                    format='%(asctime)s : %(levelname)s : %(message)s',
                    datefmt='%H:%M:%S')
logger = logging.getLogger(__name__)

import pandas as pd
import tempfile




class RNAseqSplice:

    def __init__(self, args_parsed):

        logger.info(" Creating the IGV Report.")


        import warnings

        # arguments added as object attribute 
        self.args = args_parsed


        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Seperate Arguments and add to Object 
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        self.cancer_introns = args_parsed.cancer_introns
        self.flanking       = args_parsed.flanking
        self.fasta          = args_parsed.fasta
        self.ref_bam        = args_parsed.bam
        self.ref_bed        = args_parsed.bed

        if os.path.exists(args_parsed.output_dir):
            self.output_dir = args_parsed.output_dir
        else:
            os.mkdir(args_parsed.output_dir)
            self.output_dir = args_parsed.output_dir

        self.output         = os.path.join(self.output_dir, args_parsed.output)


    def createBedFile(self):
        '''
        Fucntion to take in the cancer inton file and create the BED input file 
        to create the IGV report.
        '''
        cancer_introns = self.cancer_introns

        # Read in the cander intron data 
        dt = pd.read_table(cancer_introns)

        # split chr:start:stop into seperate columns 
        split1 = dt["intron"].str.split(":",n = 1, expand = True) 
        split2 = split1[1].str.split("-",n = 1, expand = True) 
        split1[2] = split2[0]
        split1[3] = split2[1]
        split1.drop(columns =[1], inplace = True) 
        
        # Add the name column to the data table 
        bed_file = split1
        bed_file.insert(loc = 3, column = "NAME", value = dt["variant_name"])
        self.bed_file = bed_file

        return(self)
        







    def CTAT_splicing_create_report(self):




        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Create the output bed file
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
        bed_output =os.path.join(self.output_dir, "Introns.bed")
        file = open(bed_output, "w") 

        # Convert the bed file pandas Data Frame to a csv string format 
        text = self.bed_file.to_csv(index=False, header=None, sep="\t")
        file.write(text) # Write to the temporary file 
        file.close()

        

        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # create the IGV table 
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        table = BedTable(bed_output)
        table_json = table.to_JSON()

        session_dict = {}


        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Set the TRACKS 
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        self.tracks = [self.ref_bam, self.ref_bed, bed_output]

        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Create file readers for tracks.  This is done outside the loop so initialization happens onc
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        trackreaders = []
        if self.tracks is not None:
            for track in self.tracks:
                print(track)
                reader = utils.getreader(track)
                print(reader)
                trackreaders.append({
                    "track": track,
                    "reader": reader
                })



        # loop through variants creating an igv.js session for each one
        for tuple in table.features:
            feature = tuple[0]
            unique_id = tuple[1]

            # Define a genomic region around the variant
            chr = feature.chr
            position = int(math.floor((feature.start + feature.end) / 2)) + 1   # center of region in 1-based coordinates
            start = int (math.floor(feature.start - int(self.flanking) / 2))
            end = int (math.ceil(feature.end + int(self.flanking) / 2))
            region = {
                "chr": chr,
                "start": start,
                "end": end
            }

            # Fasta
            data = fasta.get_data(self.fasta, region)
            fa = '>' + chr + ':' + str(start) + '-' + str(end) + '\n' + data
            fasta_uri = datauri.get_data_uri(fa)
            fastaJson = {
                "fastaURL": fasta_uri,
            }

            # Initial locus, +/- 20 bases
            initial_locus = chr + ":" + str(position - 20) + "-" + str(position + 20)
            session_json = {
                "locus": initial_locus,
                "reference": fastaJson,
                "tracks": []
            }

            for tr in trackreaders:

                track = tr["track"]
                reader = tr["reader"]
                trackObj = tracks.get_track_json_dict(track)
                data = reader.slice(region)
                trackObj["url"] = datauri.get_data_uri(data)
                if(trackObj["type"] == "alignment"):
                    trackObj["height"] = 500

                    # Sort TODO -- do this only for SNV
                    # if (trackObj["type"]) == "alignment":
                    #     trackObj["sort"] = {
                    #         "option": "NUCLEOTIDE",
                    #         "locus": chr + ":" + str(variant.pos - 1)
                    #     }

                session_json["tracks"].append(trackObj)


            # Build the session data URI

            session_string = json.dumps(session_json);

            session_uri = datauri.get_data_uri(session_string)

            session_dict[str(unique_id)] = session_uri







        session_dict = json.dumps(session_dict)

        
        template_file = None
        if None == template_file:
            template_file = os.path.dirname(sys.modules['igv_reports'].__file__) + '/templates/variant_template.html'

        output_file = self.output

        # standalone = args.standalone
        standalone = None
        with open(template_file, "r") as f:
            data = f.readlines()

            with open(output_file, "w") as o:

                for i, line in enumerate(data):

                    if standalone and line.find("<script") and line.find(".js") > 0:
                        print(inline_script(line, o))

                    else:
                        j = line.find('"@TABLE_JSON@"')
                        if j >= 0:
                            line = line.replace('"@TABLE_JSON@"', table_json)

                        j = line.find('"@SESSION_DICTIONARY@"')
                        if j >= 0:
                            line = line.replace('"@SESSION_DICTIONARY@"', session_dict)

                        o.write(line)




def main():
    ## Input Arguments
    args_parser = argparse.ArgumentParser(
        description = "Creates the IGV Report for spilce data."
        )

    args_parser.add_argument("cancer_introns",  help="TSV File that holds the splicing information.")
    args_parser.add_argument("fasta",           help="reference fasta file, required")
    
    args_parser.add_argument("--bam",           help="reference BAM file")
    args_parser.add_argument("--bed",           help="reference BED file")

    args_parser.add_argument("--output",        help="output file name",       default="igvjs_viewer.html")
    args_parser.add_argument("--output_dir",    help="output directory",       required=False,     default="IGV_OUTPUT")
    args_parser.add_argument("--flanking",      help="genomic region to include either side of variant", default=10000)
    args_parser.add_argument("--debug",         help="sets debug mode for logger", action="store_true")
  

    return args_parser

if __name__ == "__main__":
    
    args_parser = main()
    args_parsed = args_parser.parse_args()


    if args_parsed.debug:
        logger.setLevel(logging.DEBUG)



    splice_obj = RNAseqSplice(args_parsed)

    splice_obj.createBedFile()

    splice_obj.CTAT_splicing_create_report()











