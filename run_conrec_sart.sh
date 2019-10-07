#!/bin/bash

groups=$1
number=$2
proj_abs_range=$3
interval=$4
snr=$5

echo "total groups: $groups"
echo "number of particles each group: $number"
echo "projection range: -$proj_abs_range +$proj_abs_range"
echo "projection angle interval: $interval"
echo "SNR: $snr"

root="./${groups}_${number}_${proj_abs_range}_${interval}_${snr}"

mkdir $root

for k in `seq 1 $groups`
do
    mpirun -n 30 ./bin/cnstrec -d "$root/group$k" -o "$root/group$k""_conrec_sart.mrc" -r 65 -g 0,0,0,200 -m SART,20,0.2
done
 
