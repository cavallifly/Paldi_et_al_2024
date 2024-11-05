#!/bin/bash

#SBATCH --job-name radial
#SBATCH --mem 15Gb
#SBATCH -t 4-00:00              # Runtime in D-HH:MM
#SBATCH -o 04_compute_radialPositioning.out # File to which STDOUT will be written
#SBATCH -e 04_compute_radialPositioning.out # File to which STDERR will be written 

ncopies=$1
maxReplica=30
PI=3.141592653589793238462643383279502884197

phi=0.10
resolution=5000
CGfile=../../Modelling_setup/coarse_graining/coarse_graining_FS_20nm_CG_nu_${resolution}bp_Mouse_*_phi_${phi}_*_${ncopies}_copies.txt
 
b=$(cat ${CGfile} | grep -w bCG | awk '{print $9}')
lk=$(cat ${CGfile} | grep -w lkCG | awk '{print $10}')
radius=$(cat ${CGfile} | grep -w radiusCG | awk '{print $16}')
echo $b $lk $radius

compartmentFile=/zssd/scratch/mdistefano/2023_10_22_Project_hyperacetylation/Modelling_setup/scripts/_compartments_${resolution}bp_${ncopies}_copies
nBins=10

for phi in 0.10 ;
do
    for phiDir in $(ls -1 | grep _${phi}_);
    do
	echo $phiDir
	cd  $phiDir
	
	for dir in $(ls -1 | grep EAA | grep _${ncopies}_copies);		   
	do
	    
	    if [[ ! -d ${dir} ]];
	    then
		echo "${dir} is not a directory!"
		continue
	    fi
	    echo $dir
	    cd $dir
	    pwd
	    
	    tdelta=900000
	    
	    for tmax in $(seq 5400000 5400000 129600000);
	    #for tmax in $(seq 5400000 5400000 21600000);	    
	    do
		tmin=$(echo $tmax | awk '{print $1-5400000}')
		if [[ -e histogram_radial_position_from_${tmin}_to_${tmax}_${dir}_A.tab ]];
		then
		    ls -lrtha histogram_radial_position_from_${tmin}_to_${tmax}_${dir}_?.tab
		    continue
		fi
		

		for t in $(seq $tmin $tdelta $tmax) ;
		do
		    if [[ ! -e replica_1/compartmentalization_${t}.XYZ ]];
		    then
			continue
		    fi
		    
		    for replica in $(seq 1 1 ${maxReplica});
		    do
			if [[ ! -e replica_${replica}/compartmentalization_${t}.XYZ ]];
			then
			    continue
			fi
			
			echo $tmin $t ${tmax} $replica
			outFile=${PWD}/radial_position_from_${tmin}_to_${tmax}_${dir}_${replica}.tab				    
			ls -lrtha ./replica_${replica}/*_${t}.XYZ
			awk '{if(NR>9){print $1,sqrt($2*$2+$3*$3+$4+$4)}}' ./replica_${replica}/*_${t}.XYZ > _radial_position_${tmax}_${replica}

			paste _radial_position_${tmax}_${replica} ${compartmentFile} | awk -v t=$t -v r=${replica} '{if($1==$3) print t,r,$1,$NF,$2}' > _radial_position_${tmax}_${replica}_compartment
			cat _radial_position_${tmax}_${replica}_compartment >> ${outFile}
			for compartment in T A B;
			do
		    
			    #grep -w ${compartment} ${outFile} | awk -v R=${radius} -v nBins=${nBins} -v PI=${PI} 'BEGIN{for(i=0;i<=nBins;i++){r[i]=(i/nBins)^(1./3.)*R; print i,r[i],(4./3.*PI*r[i]*r[i]*r[i])/(4./3.*PI*R*R*R),R}}{for(i=1;i<=nBins;i++){if(r[i-1]<=$NF && $NF<r[i]){h[i]++}}}END{for(i=1;i<=nBins;i++){print i/nBins,h[i]}}' | head -15
			    grep -w ${compartment} _radial_position_${tmax}_${replica}_compartment | awk -v c=${compartment} -v R=${radius} -v nBins=${nBins} -v PI=${PI} 'BEGIN{for(i=0;i<=nBins;i++){r[i]=(i/nBins)^(1./3.)*R}}{for(i=1;i<=nBins;i++){if(r[i-1]<=$NF && $NF<r[i]){h[i]++}}}END{for(i=1;i<=nBins;i++){print i/nBins,h[i],c}}' >> histogram_radial_position_from_${tmin}_to_${tmax}_${dir}_${replica}_${compartment}.tab
			done
			rm -fvr _radial_position_${tmax}_${replica}_compartment _radial_position_${tmax}_${replica}
						
		    done # Close cycle over $replica
		done # Close cycle over $t

		for compartment in T A B;
		do

		    #ls -lrtha histogram_radial_position_from_${tmin}_to_${tmax}_${dir}_*_${compartment}.tab
		    awk -v c=${compartment} '{s[$1]+=$2; s2[$1]+=$2*$2; cnt[$1]++}END{for(i in s){avg=s[i]/cnt[i]; avg2=s2[i]/cnt[i]; stdDev=sqrt(avg2-avg*avg); print i,avg,stdDev,cnt[i],c}}' histogram_radial_position_from_${tmin}_to_${tmax}_${dir}_*_${compartment}.tab | sort -k 1,1n > histogram_radial_position_from_${tmin}_to_${tmax}_${dir}_${compartment}.tab
		    
		done
		rm -fvr _radial_position_*
		
	    done # Close cycle over $tmax    

	    
	    cd ..	    
	done # Close cycle over $dir
	
	cd .. # Exit phiDir
    done # Close cycle over $phiDir
    

done # Close cycle over $phi
wait
