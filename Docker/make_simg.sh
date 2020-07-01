#!/bin/bash

VERSION=`cat VERSION.txt`

singularity build ctat_splicing.v${VERSION}.simg docker://trinityctat/ctat_splicing:$VERSION




