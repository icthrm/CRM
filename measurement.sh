#!/bin/bash
path_mrc1=$1
path_mrc_ref=$2
radius_mask=$3
length_image=$4
peakval=$5

BIN_PATH="/data/projects/constraint-reconstruction/data";
#echo "measurement(${path_mrc1},${path_mrc_ref},${radius_mask},${length_image},1,${peakval})"
# matlab -nosplash -nodesktop -r "measurement('${path_mrc1}','${path_mrc_ref}',${radius_mask},${length_image},1,${peakval});exit"
matlab -nosplash -nodesktop -r "addpath('$BIN_PATH');measurement('${path_mrc1}','${path_mrc_ref}',${radius_mask},${length_image},1);exit"
