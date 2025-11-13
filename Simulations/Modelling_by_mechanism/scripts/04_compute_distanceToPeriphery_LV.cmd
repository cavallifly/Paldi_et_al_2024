#!/bin/bash

#SBATCH --job-name radial
#SBATCH --mem 15Gb
#SBATCH -t 4-00:00              # Runtime in D-HH:MM
#SBATCH -o 04_compute_radialPositioning.out # File to which STDOUT will be written
#SBATCH -e 04_compute_radialPositioning.out # File to which STDERR will be written 

ncopies=$1
maxReplica=5
PI=3.141592653589793238462643383279502884197
b=54.1961
radius=46.4159

phi=0.10
resolution=5000
CGfile=../../Modelling_setup/coarse_graining/coarse_graining_FS_20nm_CG_nu_${resolution}bp_Mouse_*_phi_${phi}LV_*_${ncopies}_copies.txt
 
b=$(cat ${CGfile} | grep -w bCG | awk '{print $9}')
lk=$(cat ${CGfile} | grep -w lkCG | awk '{print $10}')
radius=$(cat ${CGfile} | grep -w radiusCG | awk '{print $16}')
echo $b $lk $radius

compartmentFile=../../Modelling_setup/scripts/_compartments_${resolution}bp_${ncopies}_copies
nBins=10

for phi in 0.10 ;
do
    for phiDir in $(ls -1 | grep _${phi}_ | grep -v nikita | grep Rn_49.0);
    do
	echo $phiDir

	#radius=$(echo ${phiDir} | sed "s/_/ /g" | awk '{print $6}')
	#echo $radius
	cd  $phiDir

	pwd
	for dir in $(ls -1 | grep EAA | grep _${ncopies}_copies | grep -v tau | grep "lkA_1.00_lkB_0.00"); # | grep "lkA_18.00_lkB_3.00\|lkA_1.00_lkB_0.00");
	do
	    echo $dir

	    if [[ ! -d ${dir} ]];
	    then
		echo "${dir} is not a directory!"
		continue
	    fi
	    echo $dir
	    cd $dir
	    pwd
	    radius=$(grep -i radius= replica_1/*py | sed "s/=/ /g" | awk '{print $NF}' | uniq)
	    echo "Radius $radius"

	    tdelta=900000
	    
	    for tmax in $(seq 5400000 5400000 129600000);
	    do
		tmin=$(echo $tmax | awk '{print $1-5400000}')
		if [[ ! -e replica_1/compartmentalization_${tmin}.XYZ ]];
                then
                    continue
                fi
		if [[ ! -e replica_1/compartmentalization_${tmax}.XYZ ]];
                then
                    continue
                fi

		if [[ ! -e histogram_distance_to_periphery_from_${tmin}_to_${tmax}_${dir}_A.tab ]];
		then
		    
		    for t in $(seq $tmin $tdelta $tmax) ;
		    do
			echo $t
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
			    outFile=${PWD}/distance_to_periphery_from_${tmin}_to_${tmax}_${dir}_${replica}.tab

			    ls -lrtha ./replica_${replica}/*_${t}.XYZ
			    awk -v r=$radius '{if(NR>9){print $1,r-sqrt($2*$2+$3*$3+$4*$4)}}' ./replica_${replica}/*_${t}.XYZ > _distance_to_periphery_${tmax}_${replica}
			    
			    paste _distance_to_periphery_${tmax}_${replica} ${compartmentFile} | awk -v t=$t -v r=${replica} '{if($1==$3) print t,r,$1,$NF,$2}' > _distance_to_periphery_${tmax}_${replica}_compartment
			    cat _distance_to_periphery_${tmax}_${replica}_compartment >> ${outFile}
			    for compartment in T A B;
			    do
				
				#grep -w ${compartment} ${outFile} | awk -v R=${radius} -v nBins=${nBins} -v PI=${PI} 'BEGIN{for(i=0;i<=nBins;i++){r[i]=(i/nBins)^(1./3.)*R; print i,r[i],(4./3.*PI*r[i]*r[i]*r[i])/(4./3.*PI*R*R*R),R}}{for(i=1;i<=nBins;i++){if(r[i-1]<=$NF && $NF<r[i]){h[i]++}}}END{for(i=1;i<=nBins;i++){print i/nBins,h[i]}}' | head -15
				grep -w ${compartment} _distance_to_periphery_${tmax}_${replica}_compartment | awk -v c=${compartment} -v R=${radius} -v nBins=${nBins} -v PI=${PI} 'BEGIN{for(i=0;i<=nBins;i++){r[i]=(i/nBins)^(1./3.)*R}}{for(i=1;i<=nBins;i++){if(r[i-1]<=$NF && $NF<r[i]){h[i]++}}}END{for(i=1;i<=nBins;i++){print i/nBins,r[i-1],r[i],h[i],c}}' >> histogram_distance_to_periphery_from_${tmin}_to_${tmax}_${dir}_${replica}_${compartment}.tab
				#head histogram_distance_to_periphery_from_${tmin}_to_${tmax}_${dir}_${replica}_${compartment}.tab

			    done
			    #rm -fvr _distance_to_periphery_${tmax}_${replica}_compartment _distance_to_periphery_${tmax}_${replica}
			    
			done # Close cycle over $replica
		    done # Close cycle over $t
		    
		    for compartment in T A B;
		    do
			
			#ls -lrtha histogram_distance_to_periphery_from_${tmin}_to_${tmax}_${dir}_*_${compartment}.tab
			awk -v c=${compartment} '{s[$1]+=$4; s2[$1]+=$4*$4; st[$1]=$2; end[$1]=$3; cnt[$1]++}END{for(i in s){avg=s[i]/cnt[i]; avg2=s2[i]/cnt[i]; stdDev=sqrt(avg2-avg*avg); if(cnt[i]==0){avg=0; stdDev=0}; print i,st[i],end[i],avg,stdDev,cnt[i],c}}' histogram_distance_to_periphery_from_${tmin}_to_${tmax}_${dir}_*_${compartment}.tab | sort -k 1,1n > histogram_distance_to_periphery_from_${tmin}_to_${tmax}_${dir}_${compartment}.tab
			cat histogram_distance_to_periphery_from_${tmin}_to_${tmax}_${dir}_${compartment}.tab
		    done
		    #rm -fvr _distance_to_periphery_*
		else
		    ls -lrtha histogram_distance_to_periphery_from_${tmin}_to_${tmax}_${dir}_?.tab
		fi
		if [[ ! -e histogram_distance_to_periphery_from_${tmin}_to_${tmax}_${dir}_A.tab ]];
		then
		    continue
		fi

		if [[ ! -e histogram_distance_to_periphery_from_${tmin}_to_${tmax}_${dir}_A.pdf ]];
		then
		    binwidth=0.1
		    ymin=0
		    ymax=0.4		    
		    
		    for compartment in T A B;
		    do
			nparticles=$(echo 300 | awk '{print $1*6}')
			color="blue"
			if [[ ${compartment} == "A" ]];
			then
			    color="red"
			fi
			if [[ ${compartment} == "T" ]];
			then
			    color="gray"
			    nparticles=$(echo 200 | awk '{print $1*2}')
			fi

			echo $compartment
			head histogram_distance_to_periphery_from_${tmin}_to_${tmax}_${dir}_${compartment}.tab 
			awk -v np=${nparticles} -v nc=${ncopies} -v b=$b -v r=$radius '{print $1,($3+$2)/2/r,$4/(np*nc),($3-$2)/r}' histogram_distance_to_periphery_from_${tmin}_to_${tmax}_${dir}_${compartment}.tab > data.dat
			cat data.dat
			sort -k 2,2n data.dat | tail -1
						
			sed -e "s/XXXyminXXX/${ymin}/g" -e "s/XXXymaxXXX/${ymax}/g" -e "s/XXXcolorXXX/${color}/g" -e "s/XXXbinwidthXXX/${binwidth}/g" /zssd/scratch/mdistefano/2023_10_22_Project_hyperacetylation/Modelling_by_mechanism/scripts/04_plot_radialPositioning_genotoul.gp | gnuplot

			ps2pdf data.ps			
			mv data.pdf histogram_distance_to_periphery_from_${tmin}_to_${tmax}_${dir}_${compartment}.pdf
			ls -lrtha histogram_distance_to_periphery_from_${tmin}_to_${tmax}_${dir}_${compartment}.pdf
			rm data.*
		    done
		else
		    ls -lrtha histogram_distance_to_periphery_from_${tmin}_to_${tmax}_${dir}_?.pdf
		fi
	    done # Close cycle over $tmax    
	    pwd

	    cd ..	    
	done # Close cycle over $dir
	
	cd .. # Exit phiDir
    done # Close cycle over $phiDir
    

done # Close cycle over $phi
wait
