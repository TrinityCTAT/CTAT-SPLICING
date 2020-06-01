#!/bin/bash

set -e

VERSION=`cat VERSION.txt`

docker push trinityctat/ctat_splicing:$VERSION 
docker push trinityctat/ctat_splicing:latest 


