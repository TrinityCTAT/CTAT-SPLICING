#!/bin/bash

set -e

VERSION=`cat VERSION.txt`

docker build -t trinityctat/ctat_splicing:$VERSION .
docker build -t trinityctat/ctat_splicing:latest .


