#!/bin/bash

#SBATCH --job-name TR
##SBATCH -n 1                    # Number of cores. For now 56 is the number max of core available
#SBATCH --mem 15Gb
#SBATCH -t 4-00:00              # Runtime in D-HH:MM
#SBATCH -o 01_compute_contact_map_TR.out # File to which STDOUT will be written
#SBATCH -e 01_compute_contact_map_TR.out # File to which STDERR will be written 

ncopies=$1
replica=$2

for phi in 0.10 ;
do
    for phiDir in $(ls -1 | grep _${phi}_);
    do
	echo $phiDir
	cd  $phiDir
	
	for dir in $(ls -1 | grep EAA | grep _${ncopies}_copies | grep -v tau) ; #| grep "lkA_14.00_lkB_3.00\|lkA_1.00_lkB_0.00");		   
	do
	    
	    if [[ ! -d ${dir} ]];
	    then
		echo "${dir} is not a directory!"
		continue
	    fi
	    
	    #for minReplicates in $(seq 30 10 30);
	    for minReplicates in 4 ;
	    do
		if [[ ! -d ${dir}/replica_${minReplicates} ]];
		then
		    echo "${minReplicates} directories not available!"
		    continue
		fi

		tdelta=900000
		
		for tmax in $(seq 5400000 5400000 129600000);
		#for tmax in 21600000 ;			    
		do
		    tmin=$((${tmax}-5400000));
		    echo $tmin $tdelta $tmax

		    tag="OK"
		    check=$(ls -1 ${dir}/rep*/*_${tmax}.lammpstrj ${dir}/rep*/*_${tmax}.XYZ 2> /dev/null | wc -l)
		    if [[ $check -lt ${minReplicates} ]];
		    then
			echo "${dir} has ${check} replicate, less than ${minReplicates} done. Remember to include all replicates!"
			continue
			tag="tmp"
		    fi
		    
		    echo $dir
		    pwd
		    bash ../../scripts/01_compute_contact_map_genotoul_TR_and_CS.sh ${dir} ${phi} ${minReplicates} ${tmin} ${tmax} ${tdelta} ${tag} $replica

		done # Close cycle over $tmin
	    done # Close cycle over $minReplicates

	done # Close cycle over $dir
	
	cd .. # Exit phiDir
	
    done # Close cycle over $phiDir	
done # Close cycle over $phi
wait
