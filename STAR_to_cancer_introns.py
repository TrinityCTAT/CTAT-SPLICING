#!/usr/bin/env python

import sys, os, re
import argparse
import subprocess


scriptdir=os.path.dirname(__file__)
utildir=scriptdir + "/util"
sys.path.append(utildir)
import intron_occurrence_capture as ioc



def main():

    parser = argparse.ArgumentParser(description="capture gene to intron usage stats", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    ctat_genome_lib = os.environ.get('CTAT_GENOME_LIB', None)
    parser.add_argument("--ctat_genome_lib", dest="ctat_genome_lib", type=str, required=False, default=ctat_genome_lib, help="ctat genome lib build dir")
    parser.add_argument("--SJ_tab_file", dest="SJ_tab_file", type=str, required=True, help="STAR SJ.out.tab file")
    parser.add_argument("--chimJ_file", dest="chimJ_file", type=str, required=False, default=None, help="STAR Chimeric.out.junction file")
    parser.add_argument("--output_prefix", dest="output_prefix", type=str, required=True, help="prefix for all output files")

    parser.add_argument("--vis", action='store_true', default=False, help="Generate igv html ctat splicing visualization (requires --bam_file to be set)")
    parser.add_argument("--bam_file", dest="bam_file", type=str, required=False, default=None, help="STAR generated BAM file")
        
    args = parser.parse_args()

    ctat_genome_lib = args.ctat_genome_lib

    if not os.path.exists(ctat_genome_lib):
        raise RuntimeError("Error, must set --ctat_genome_lib ")
    
    SJ_tab_file = args.SJ_tab_file
    chimJ_file = args.chimJ_file
    output_prefix = args.output_prefix
    bam_file = args.bam_file
    VIS_flag = args.vis

    if VIS_flag and not bam_file:
        raise RuntimeError("Error, if --vis, must specify --bam_file ")
    
    
    if not os.path.exists(SJ_tab_file):
        raise RuntimeError("Error, cannot locate expected splice junction tab file: {} ".format(SJ_tab_file))


    targets_list_file = os.path.join(ctat_genome_lib, "ref_annot.gtf.mini.sortu")
    chr_intron_bounds = ioc.populate_intron_bounds(targets_list_file)
    introns_dict = ioc.map_introns_from_splice_tab(SJ_tab_file, chr_intron_bounds)

    if chimJ_file is not None:
        if not os.path.exists(chimJ_file):
            raise RuntimeError("Error, cannot locate expected chimeric Junctiom out file: {} ".format(chimJ_file))
        
        # must make splice file:
        chimJ_introns_file = output_prefix + "." + os.path.basename(chimJ_file) + ".introns.tmp"
        cmd = str(os.path.join(utildir, "STAR_chimeric_junctions_to_introns.pl") +
                  " -J {} > {}".format(chimJ_file, chimJ_introns_file))
        subprocess.check_call(cmd, shell=True)
        
        introns_dict = ioc.supplement_introns_from_chimeric_junctions_file(chimJ_introns_file, introns_dict, chr_intron_bounds)

    introns_output_file = output_prefix + ".introns"
    with open(introns_output_file, 'wt') as ofh:
        # write header
        ofh.write("\t".join(["intron", "strand", "genes", "uniq_mapped", "multi_mapped"]) + "\n")
        
        for intron in introns_dict.values():
            ofh.write("\t".join(["{}:{}-{}".format(intron.chromosome, intron.lend, intron.rend),
                                 intron.strand,
                                 intron.genes,
                                 str(intron.uniq_mapped), str(intron.multi_mapped)]) + "\n")

    # annotate for cancer introns.
    cmd = str(os.path.join(utildir, "annotate_cancer_introns.pl") +
              " --introns_file {} ".format(introns_output_file) +
              " --ctat_genome_lib {} ".format(ctat_genome_lib) +
              " --intron_col 0 " +
              " > {} ".format(output_prefix + ".cancer.introns") )

    subprocess.check_call(cmd, shell=True)
    


    if VIS_flag:

        # Create the IGV Reports 
        CTAT_Splicing_IGV = str(os.path.join(utildir,"CTAT_Splicing_IGV.py"))
        ref_genome = str(os.path.join(ctat_genome_lib, "ref_genome.fa"))
        ref_bed = str(os.path.join(ctat_genome_lib, "ref_annot.sorted.bed.gz"))
        introns_cancer_output_file = str(output_prefix + ".cancer.introns")

        IGV_cmd  = " ".join([   CTAT_Splicing_IGV,
                                introns_cancer_output_file,
                                ref_genome,
                                "--bam ", bam_file,
                                "--bed ", ref_bed, 
                                "--flanking 10000",
                                "--output igvjs_viewer.html"
                                ])

        subprocess.check_call(IGV_cmd, shell=True)

    sys.exit(0)



if __name__=='__main__':
    main()
