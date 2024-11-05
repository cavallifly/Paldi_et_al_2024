#!/bin/sh

tag=TSA
resolution=5000 
chrlength=20000000
ncopies=$1
size=$(echo ${chrlength} ${ncopies} ${resolution} | awk '{print int($1/$3)}')
echo "Length of the chain ${chrlength} bp - ${size} ${resolution}bp-beads"
sigma=1.0

for phi in 0.03 ;
do

    b=$(grep -w bCG ../coarse_graining/coarse_graining_FS_20nm_CG_nu_5000bp_Mouse_mESC_phi_${phi}_TSA_${ncopies}_copies.txt | awk '{print $(NF-1)}')

    radius=$(grep radiusCG ../coarse_graining/coarse_graining_FS_20nm_CG_nu_5000bp_Mouse_mESC_phi_${phi}_TSA_${ncopies}_copies.txt | awk '{print $NF}')
    echo $sigma $phi $radius $b

    # Generate the initial rod-like conformation using TADphys
    for replica in $(seq $2 1 $2);
    do
	if [[ -e Initial_rosette_conformation_sphere_replica_${replica}.dat ]];
	then
	    continue
	fi
	touch Initial_rosette_conformation_sphere_replica_${replica}.dat
	seed=$(od -An -N3 -i /dev/random | awk '{print $1}')
	
	echo $replica $seed >> used_seeds.txt
	
	sed -e "s/XXXsigmaXXX/${sigma}/g" -e "s/XXXresolutionXXX/${resolution}/g" -e "s/XXXsizeXXX/${size}/g" -e "s/XXXchrlengthXXX/${chrlength}/g" -e "s/XXXseedXXX/${seed}/g" -e "s/XXXreplicaXXX/${replica}/g" -e "s/XXXncopiesXXX/${ncopies}/g" -e "s/XXXradiusXXX/${radius}/g" -e "s/XXXbXXX/${b}/g" ../scripts/00_Initial_conformation.py | python
	
    done
done
