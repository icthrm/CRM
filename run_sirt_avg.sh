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
    mkdir "$root/rec$k""_sirt"
    ran=$RANDOM;
    for m in `seq 1 $number`
    do
        mpirun -n 30 volrec_sglm -i "$root/group$k/$m/proj.mrc" -o "$root/rec$k""_sirt/rec$m.mrc" -a "$root/group$k/$m/simu.tlt" -g 0,0,0,200  -m SIRT,120,0.2
    done
    mrcavg -d "$root/rec$k""_sirt" -o "$root/group$k""_avg_sirt.mrc"
done
 
