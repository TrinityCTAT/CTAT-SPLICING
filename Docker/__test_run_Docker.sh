#!/bin/bash

set -ve

if [ -z ${CTAT_GENOME_LIB} ]; then
    echo "Error, must have CTAT_GENOME_LIB env var set"
    exit 1
fi


VERSION=`cat VERSION.txt`

docker run -v `pwd`/../:/data -v ${CTAT_GENOME_LIB}:/ctat_genome_lib --rm trinityctat/ctat_splicing:${VERSION} \
       /usr/local/src/CTAT-SPLICING/STAR_to_cancer_introns.py \
       --SJ_tab_file /data/testing/SJ.out.tab.b38 \
       --chimJ_file /data/testing/Chimeric.out.junction.b38 \
       --vis --bam /data/testing/alignments.b38.sorted.bam \
       --output_prefix /data/testing/docker_test \
       --ctat_genome_lib /ctat_genome_lib 

