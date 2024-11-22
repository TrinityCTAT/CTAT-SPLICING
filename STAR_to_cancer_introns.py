#!/usr/bin/env python

import sys, os, re
import argparse
import subprocess
import logging
import pandas as pd

VERSION = "0.0.3"


logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s : %(levelname)s : %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger(__name__)

scriptdir = os.path.dirname(os.path.abspath(__file__))

sys.path.insert(0, os.path.join(scriptdir, "PyLib"))
from Pipeliner import *

utildir = scriptdir + "/util"
sys.path.append(utildir)
import intron_occurrence_capture as ioc


def main():

    parser = argparse.ArgumentParser(
        description="capture gene to intron usage stats",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    ctat_genome_lib = os.environ.get("CTAT_GENOME_LIB", None)
    parser.add_argument(
        "--ctat_genome_lib",
        dest="ctat_genome_lib",
        type=str,
        required=False,
        default=ctat_genome_lib,
        help="ctat genome lib build dir",
    )
    parser.add_argument(
        "--SJ_tab_file",
        dest="SJ_tab_file",
        type=str,
        required=True,
        help="STAR SJ.out.tab file",
    )
    parser.add_argument(
        "--chimJ_file",
        dest="chimJ_file",
        type=str,
        required=False,
        default=None,
        help="STAR Chimeric.out.junction file",
    )
    parser.add_argument(
        "--output_prefix",
        dest="output_prefix",
        type=str,
        required=True,
        help="prefix for all output files",
    )
    parser.add_argument(
        "--min_total_reads",
        dest="min_total_reads",
        type=int,
        required=False,
        default=5,
        help="minimum reads supporting cancer intron",
    )

    parser.add_argument(
        "--vis",
        action="store_true",
        default=False,
        help="Generate igv html ctat splicing visualization (requires --bam_file to be set)",
    )
    parser.add_argument(
        "--bam_file",
        dest="bam_file",
        type=str,
        required=False,
        default=None,
        help="STAR generated BAM file",
    )
    parser.add_argument(
        "--sample_name",
        dest="vis_sample_name",
        type=str,
        required=False,
        default="",
        help="sample name for vis title",
    )

    args = parser.parse_args()

    ctat_genome_lib = args.ctat_genome_lib

    if not os.path.exists(ctat_genome_lib):
        raise RuntimeError("Error, must set --ctat_genome_lib ")

    SJ_tab_file = args.SJ_tab_file
    chimJ_file = args.chimJ_file
    output_prefix = args.output_prefix
    bam_file = args.bam_file
    VIS_flag = args.vis
    min_total_reads = args.min_total_reads
    vis_sample_name = args.vis_sample_name


    # check for splicing info installation in ctat genome lib:
    splicing_db = os.path.join(ctat_genome_lib, "cancer_splicing_lib/cancer_splicing.idx")
    if not os.path.exists(splicing_db):
        logger.error(f"Splicing database doesn't appear to be installed in the ctat genome lib. Cannot locate: {splicing_db}\nBe sure to follow the ctat splicing setup instructions")
        sys.exit(1)
    

    
    if VIS_flag and not bam_file:
        raise RuntimeError("Error, if --vis, must specify --bam_file ")

    if not os.path.exists(SJ_tab_file):
        raise RuntimeError(
            "Error, cannot locate expected splice junction tab file: {} ".format(
                SJ_tab_file
            )
        )

    chckpts_dir = output_prefix + ".chckpts"
    if not os.path.exists(chckpts_dir):
        os.makedirs(chckpts_dir)

    pipeliner = Pipeliner(chckpts_dir)

    introns_output_file = output_prefix + ".introns"
    introns_output_file_chckpt = os.path.join(chckpts_dir, "introns.ok")

    
    if not os.path.exists(introns_output_file_chckpt):

        targets_list_file = os.path.join(ctat_genome_lib, "ref_annot.gtf.mini.sortu")
        chr_intron_bounds = ioc.populate_intron_bounds(targets_list_file)
        introns_dict = ioc.map_introns_from_splice_tab(SJ_tab_file, chr_intron_bounds)

        if chimJ_file is not None:
            if not os.path.exists(chimJ_file):
                raise RuntimeError(
                    "Error, cannot locate expected chimeric Junctiom out file: {} ".format(
                        chimJ_file
                    )
                )

            # must make splice file:
            chimJ_introns_file = (
                output_prefix + "." + os.path.basename(chimJ_file) + ".introns.tmp"
            )
            cmd = str(
                os.path.join(utildir, "STAR_chimeric_junctions_to_introns.pl")
                + " -J {} > {}".format(chimJ_file, chimJ_introns_file)
            )
            subprocess.check_call(cmd, shell=True)

            introns_dict = ioc.supplement_introns_from_chimeric_junctions_file(
                chimJ_introns_file, introns_dict, chr_intron_bounds
            )

        with open(introns_output_file, "wt") as ofh:
            # write header
            ofh.write(
                "\t".join(["intron", "strand", "genes", "uniq_mapped", "multi_mapped"])
                + "\n"
            )

            for intron in introns_dict.values():
                ofh.write(
                    "\t".join(
                        [
                            "{}:{}-{}".format(
                                intron.chromosome, intron.lend, intron.rend
                            ),
                            intron.strand,
                            intron.genes,
                            str(intron.uniq_mapped),
                            str(intron.multi_mapped),
                        ]
                    )
                    + "\n"
                )

        # done, add checkpoint
        subprocess.check_call("touch {}".format(introns_output_file_chckpt), shell=True)

    # annotate for cancer introns.
    cancer_introns_file_prelim = output_prefix + ".cancer.introns.prelim"
    cmd = str(
        os.path.join(utildir, "annotate_cancer_introns.pl")
        + " --introns_file {} ".format(introns_output_file)
        + " --ctat_genome_lib {} ".format(ctat_genome_lib)
        + " --intron_col 0 "
        + " > {} ".format(cancer_introns_file_prelim)
    )

    pipeliner.add_commands([Command(cmd, "prelim_introns.ok")])

    # filter for min support
    cancer_introns_file = output_prefix + ".cancer.introns"
    cmd = str(
        os.path.join(utildir, "filter_by_min_total_reads.py")
        + " --cancer_intron_candidates {}".format(cancer_introns_file_prelim)
        + " --min_total_reads {} ".format(min_total_reads)
        + " > {} ".format(cancer_introns_file)
    )

    pipeliner.add_commands([Command(cmd, "introns_filtered.ok")])

    pipeliner.run()

    cancer_introns = pd.read_csv(cancer_introns_file, sep="\t")
    num_cancer_introns = len(cancer_introns)
    logger.info(f"-found {num_cancer_introns} cancer introns") 
    if num_cancer_introns == 0:
        # nothing more to do here.
        sys.exit(0)
    
    
    if VIS_flag:

        # generate the intron/junctions bed needed by igv
        igv_introns_bed_file = introns_output_file + ".for_IGV.bed"
        cmd = str(
            os.path.join(utildir, "make_igv_splice_bed.py")
            + " --all_introns {} ".format(introns_output_file)
            + " --cancer_introns {} ".format(cancer_introns_file)
            + " --genome_lib_dir {} ".format(ctat_genome_lib)
            + " --output_bed {} ".format(igv_introns_bed_file)
        )

        pipeliner.add_commands([Command(cmd, "intron_igv_bed.ok")])
        pipeliner.run()

        igv_tracks_config_file = write_igv_config(
            output_prefix,
            ctat_genome_lib,
            igv_introns_bed_file,
            bam_file,
            os.path.join(utildir, "misc/igv.tracks.json"),
            pipeliner,
        )

        # Create the IGV Reports
        cmd = str(
            "create_report {} ".format(igv_introns_bed_file)
            + " {} ".format(os.path.join(ctat_genome_lib, "ref_genome.fa"))
            + " --type junction "
            + " --output {}.ctat-splicing.igv.html ".format(output_prefix)
            + " --track-config {} ".format(igv_tracks_config_file)
            + " --info-columns gene variant_name uniquely_mapped multi_mapped TCGA GTEx "
            + " --title 'CTAT_Splicing: {}' ".format(vis_sample_name)
        )

        pipeliner.add_commands([Command(cmd, "igv_create_html.ok")])
        pipeliner.run()

    logger.info("done.")

    sys.exit(0)


def write_igv_config(
    output_prefix,
    ctat_genome_lib,
    igv_introns_bed_file,
    bam_file,
    template_json_file,
    pipeliner,
):

    json_template_text = subprocess.check_output(
        "cat {}".format(template_json_file), shell=True
    ).decode()

    json_template_text = json_template_text.replace(
        "__IGV_SPLICE_BED_FILE__", igv_introns_bed_file
    )

    ref_annotations_file = os.path.join(ctat_genome_lib, "refGene.sort.bed.gz")
    json_template_text = json_template_text.replace(
        "__REF_GENE_STRUCTURE_ANNOTATIONS__", ref_annotations_file
    )

    (
        gene_reads_bam_file,
        cancer_intron_reads_bam_file,
    ) = get_gene_and_cancer_intron_reads_bam_files(
        output_prefix, igv_introns_bed_file, bam_file, pipeliner
    )

    json_template_text = json_template_text.replace(
        "__RNASEQ_GENE_ALIGNMENTS__", gene_reads_bam_file
    )
    json_template_text = json_template_text.replace(
        "__RNASEQ_CANCER_INTRON_ALIGNMENTS__", cancer_intron_reads_bam_file
    )

    igv_track_config_file = os.path.join(output_prefix + ".igv.tracks")

    with open(igv_track_config_file, "wt") as ofh:
        ofh.write(json_template_text)

    return igv_track_config_file


def get_gene_and_cancer_intron_reads_bam_files(
    output_prefix, igv_introns_bed_file, bam_file, pipeliner
):

    # extract the aligned reads

    cmd = str(
        os.path.join(utildir, "igv_read_alignment_extractor.py")
        + " --igv_introns_bed {} ".format(igv_introns_bed_file)
        + " --bam {} ".format(bam_file)
        + " --output_prefix {} ".format(output_prefix)
    )

    pipeliner.add_commands([Command(cmd, "reads_alignments_extracted.ok")])

    gene_reads_tmp_bam = output_prefix + ".gene_reads.bam"
    cancer_intron_reads_tmp_bam = output_prefix + ".cancer_intron_reads.bam"

    cancer_intron_reads_bam = output_prefix + ".cancer_intron_reads.sorted.bam"

    pipeliner.add_commands(
        [
            Command(
                "samtools sort -o {} {}".format(
                    cancer_intron_reads_bam, cancer_intron_reads_tmp_bam
                ),
                "sort_cancer_intron_reads.ok",
            )
        ]
    )

    gene_reads_bam = output_prefix + ".gene_reads.sorted.bam"
    pipeliner.add_commands(
        [
            Command(
                "samtools sort -o {} {}".format(gene_reads_bam, gene_reads_tmp_bam),
                "sort_gene_reads.ok",
            )
        ]
    )

    ## downsample the reads
    max_coverage = 50
    gene_reads_bam_sifted = sift_bam(gene_reads_bam, max_coverage, pipeliner)

    ## index final bams
    index_bam(gene_reads_bam_sifted, pipeliner)
    index_bam(cancer_intron_reads_bam, pipeliner)

    return (gene_reads_bam_sifted, cancer_intron_reads_bam)


def sift_bam(bam_file, max_coverage, pipeliner):

    sifted_bam_tmp_file, count = re.subn(".bam$", ".sifted.bam.tmp", bam_file)
    if count != 1:
        raise RuntimeError("Error changing extension of .bam")

    bamsifter_prog = os.path.join(scriptdir, "bamsifter/bamsifter")
    cmd = str(
        bamsifter_prog
        + " -c {} ".format(max_coverage)
        + " -i 50 "
        + " -o {} ".format(sifted_bam_tmp_file)
        + " --keep_secondary "
        + " {} ".format(bam_file)
    )

    pipeliner.add_commands(
        [Command(cmd, os.path.basename(sifted_bam_tmp_file) + ".ok")]
    )

    sifted_bam_file, count = re.subn(".bam$", ".sifted.bam", bam_file)
    if count != 1:
        raise RuntimeError("Error changing extension of .bam")

    pipeliner.add_commands(
        [
            Command(
                "samtools sort -o {} {}".format(sifted_bam_file, sifted_bam_tmp_file),
                os.path.basename(sifted_bam_tmp_file) + "sorted.ok",
            )
        ]
    )

    return sifted_bam_file


def index_bam(bam_file, pipeliner):

    pipeliner.add_commands(
        [
            Command(
                "samtools index {} ".format(bam_file),
                os.path.basename(bam_file) + ".indexed.ok",
            )
        ]
    )


if __name__ == "__main__":
    main()
