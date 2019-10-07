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
    mkdir "$root/group$k"
    ran=$RANDOM;
    for m in `seq 1 $number`
    do
        mkdir "$root/group$k/$m"
        ./bin/simulator -i ./emd_3489_bin180_cut.map -o "$root/group$k/$m/simu.mrc" -s 200 -m 65 -r -1 -b 1
        start_angle=$(( ( $ran + $m * 360 / $number ) % 360 - 180 ));
        echo `seq $start_angle $interval $(( $start_angle + 2 * $proj_abs_range ))` > "$root/group$k/$m/simu.tlt"
        mpirun -n 30 volrec_sglm -i "$root/group$k/$m/simu.mrc" -o "$root/group$k/$m/proj.mrc" -a "$root/group$k/$m/simu.tlt" -g 0,0,0,200  -m RP
#         mpirun -n 30 volrec_sglm -i "$root/group$k/$m/proj.mrc" -o "$root/group$k/$m/rec.mrc" -a "$root/group$k/$m/simu.tlt" -g 0,0,0,200  -m SART,20,0.2
    done
done

for k in `seq 1 $groups`
do
    mkdir "$root/rec$k"
    ran=$RANDOM;
    for m in `seq 1 $number`
    do
        mpirun -n 30 volrec_sglm -i "$root/group$k/$m/proj.mrc" -o "$root/rec$k/rec$m.mrc" -a "$root/group$k/$m/simu.tlt" -g 0,0,0,200  -m SART,20,0.2
    done
    mrcavg -d "$root/rec$k" -o "$root/group$k""_avg_sart.mrc"
done
 
