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
    cancer_introns_file = output_prefix + ".cancer.introns"
    cmd = str(os.path.join(utildir, "annotate_cancer_introns.pl") +
              " --introns_file {} ".format(introns_output_file) +
              " --ctat_genome_lib {} ".format(ctat_genome_lib) +
              " --intron_col 0 " +
              " > {} ".format(cancer_introns_file) )
    
    subprocess.check_call(cmd, shell=True)
    


    if VIS_flag:
        
        # generate the intron/junctions bed needed by igv
        igv_introns_bed_file = introns_output_file + ".for_IGV.bed"
        cmd = str(os.path.join(utildir, "make_igv_splice_bed.py") +
                  " --all_introns {} ".format(introns_output_file) +
                  " --cancer_introns {} ".format(cancer_introns_file) +
                  " --genome_lib_dir {} ".format(ctat_genome_lib) +
                  " --output_bed {} ".format(igv_introns_bed_file) )

        subprocess.check_call(cmd, shell=True)
        

        igv_tracks_config_file = write_igv_config(output_prefix, ctat_genome_lib,
                                                  igv_introns_bed_file, bam_file,
                                                  os.path.join(utildir, "misc/igv.tracks.json"))
        
                  
        # Create the IGV Reports
        cmd = str("create_report {} ".format(igv_introns_bed_file) +
                  " {} ".format(os.path.join(ctat_genome_lib, "ref_genome.fa")) +
                  " --type junction " +
                  " --output {}.ctat-splicing.igv.html ".format(output_prefix) +
                  " --track-config {} ".format(igv_tracks_config_file) +
                  " --info-columns TCGA GTEx variant_name " +
                  " --title 'CTAT_Splicing: my sample name' ")
        
        subprocess.check_call(cmd, shell=True)
        
        
    sys.exit(0)



def write_igv_config(output_prefix, ctat_genome_lib, igv_introns_bed_file, bam_file, template_json_file):

    json_template_text = subprocess.check_output("cat {}".format(template_json_file), shell=True).decode()

    json_template_text = json_template_text.replace("__IGV_SPLICE_BED_FILE__", igv_introns_bed_file)

    ref_annotations_file = os.path.join(ctat_genome_lib, "ref_annot.sorted.bed.gz")
    json_template_text = json_template_text.replace("__REF_GENE_STRUCTURE_ANNOTATIONS__", ref_annotations_file)


    igv_track_config_file = os.path.join(output_prefix + ".igv.tracks")

    with open(igv_track_config_file, 'wt') as ofh:
        ofh.write(json_template_text)

    return igv_track_config_file
                                                    


if __name__=='__main__':
    main()
