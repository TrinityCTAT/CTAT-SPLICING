

test_b38:
	../STAR_to_cancer_introns.py  --SJ_tab_file SJ.out.tab.b38 --chimJ_file Chimeric.out.junction.b38 --output_prefix ctat --ctat_genome_lib ${CTAT_GENOME_LIB} 


test_w_vis_b38:
	../STAR_to_cancer_introns.py  --SJ_tab_file SJ.out.tab.b38 --chimJ_file Chimeric.out.junction.b38 --output_prefix ctat --ctat_genome_lib ${CTAT_GENOME_LIB} --vis --bam alignments.b38.sorted.bam


test_w_vis_hg19:
	../STAR_to_cancer_introns.py  --SJ_tab_file SJ.out.tab.hg19 --chimJ_file Chimeric.out.junction.hg19 --output_prefix ctat --ctat_genome_lib ${CTAT_GENOME_LIB} --vis --bam alignments.hg19.sorted.bam


clean:
	rm -rf ./ctat.*





