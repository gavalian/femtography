#!/bin/bash

SCRIPT_FILE=$0
CURRENT_DIR=`pwd`
echo 'script file name : ' $SCRIPT_FILE
echo 'script directory : ' $CURRENT_DIR
SOURCE_DIR=`dirname $0`
echo 'starting docker from directory : ' $CURRENT_DIR
MOUNT_DIR=$SOURCE_DIR/codes
#docker run  -v /Users/gavalian/Work/Software/project-6a.0.0/Femtography/codes:/root/workspace/codes -it pawelsznajder/partons-test
docker run  -v $CURRENT_DIR/codes:/root/workspace/codes -it pawelsznajder/partons-test
