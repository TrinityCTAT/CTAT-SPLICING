#!/bin/bash

set -ex

VERSION=`cat VERSION.txt`

singularity build ctat_splicing.v${VERSION}.simg docker://trinityctat/ctat_splicing:$VERSION

singularity exec -e ctat_splicing.v${VERSION}.simg env

cp ctat_splicing.v${VERSION}.simg ctat_splicing.vLATEST.simg

