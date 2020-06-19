

test:
	../STAR_to_cancer_introns.py  --SJ_tab_file SJ.out.tab.b38 --chimJ_file Chimeric.out.junction.b38 --output_prefix ctat --ctat_genome_lib ${CTAT_GENOME_LIB} 


test_w_vis:
	../STAR_to_cancer_introns.py  --SJ_tab_file SJ.out.tab.b38 --chimJ_file Chimeric.out.junction.b38 --output_prefix ctat --ctat_genome_lib ${CTAT_GENOME_LIB} --vis --bam alignments.b38.sorted.bam

clean:
	rm -f ./ctat.introns ./ctat.cancer.introns ./ctat.Chimeric.out.junction.b38.introns.tmp





