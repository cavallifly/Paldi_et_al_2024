#!/bin/bash

#SBATCH --job-name contacts
##SBATCH -n 1                    # Number of cores. For now 56 is the number max of core available
#SBATCH --mem 15Gb
#SBATCH -t 4-00:00              # Runtime in D-HH:MM
#SBATCH -o 02_compute_contact_map.out # File to which STDOUT will be written
#SBATCH -e 02_compute_contact_map.out # File to which STDERR will be written 

for phi in 0.10 ; #0.03 0.10 0.20 0.30 ;
do
    for phiDir in $(ls -1 | grep _${phi}_);
    do
	echo $phiDir
	cd  $phiDir
	
	for dir in $(ls -1 | grep _20_copies | grep "lkA_14.00_lkB_3.00\|lkA_1.00_lkB_0.00");
	do
	    
	    if [[ ! -d ${dir} ]];
	    then
		echo "${dir} is not a directory!"
		continue
	    fi
	    
	    #for minReplicates in $(seq 30 10 30);
	    for minReplicates in 5 ;
	    do
		if [[ ! -d ${dir}/replica_${minReplicates} ]];
		then
		    echo "${minReplicates} directories not available!"
		    continue
		fi

		tdelta=900000
		
		#for tmax in $(seq 5400000 5400000 129600000);
		#for tmax in 21600000 ;
		for tmax in 54000000 ;		
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
		    bash ../../scripts/02_compute_contact_map.sh ${dir} ${phi} ${minReplicates} ${tmin} ${tmax} ${tdelta} ${tag}

		done # Close cycle over $tmin
	    done # Close cycle over $minReplicates

	done # Close cycle over $dir
	
	cd .. # Exit phiDir
	
    done # Close cycle over $phiDir	
done # Close cycle over $phi
wait
