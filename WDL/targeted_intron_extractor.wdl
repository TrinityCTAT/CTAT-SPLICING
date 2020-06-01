workflow targeted_intron_extraction {

  String? docker = "trinityctat/ctat_splicing"
  File bam
  File bam_index
  String sample_name
  String region

  
  call extract_introns {
    input:
      docker=docker,
      bam=bam,
      bam_index=bam_index,
      sample_name=sample_name,
      region=region
  }


  output {
    File extracted_introns_tsv = extract_introns.extracted_introns_tsv
    File extracted_region_bam = extract_introns.extracted_region_bam
    File extracted_region_bam_index = extract_introns.extracted_region_bam_index
  }
  
}


task extract_introns {

  String docker
  File bam
  File bam_index
  String sample_name
  String region

  String region_adj = sub(region, ":", "-")
  
  String extracted_region_bam_filename = "${sample_name}.${region_adj}.bam"
  String extracted_introns_tsv_filename = "${sample_name}.${region_adj}.tsv"
  
  command <<<

    set -e

    samtools view -b ${bam} ${region} -o ${extracted_region_bam_filename}
    
    samtools index  ${extracted_region_bam_filename} 
    
    intron_counter.py ${extracted_region_bam_filename} ${sample_name} > ${extracted_introns_tsv_filename}

  >>>

  runtime {
    docker: docker
  }

  
  output {
    File extracted_region_bam = "${extracted_region_bam_filename}"
    File extracted_region_bam_index = "${extracted_region_bam_filename}.bai"
    File extracted_introns_tsv = "${extracted_introns_tsv_filename}"
  }

  
}
