workflow targeted {

  String? docker = "trinityctat/ctat_splicing"
  File bam
  String sample_name
  String region

  
  call extract_introns {
    input:
      docker=docker,
      bam=bam,
      sample_name=sample_name,
      region=region
  }


  output {
    File extracted_introns_tsv = extract_introns.extracted_introns_tsv
    File extracted_region_bam = extract_introns.extracted_region_bam
  }

}


task extract_introns {

  String docker
  File bam
  String sample_name
  String region

  String extracted_region_bam_filename = "${sample_name}.${region}.bam"
  String extracted_introns_tsv_filename = "${sample_name}.${region}.tsv"
  
  command <<<

    set -e

    samtools view -b ${bam} ${region} -o ${extracted_region_bam_filename}

    intron_counter.py ${sample_name} ${extracted_region_bam_filename} > ${extracted_introns_tsv_filename}

  >>>

  runtime {
    docker: docker
  }

  
  output {
    File extracted_regions_bam = "${extracted_region_bam_filename}"
    File extracted_introns_tsv = "${extracted_introns_tsv_filename}"
  }

  
}
