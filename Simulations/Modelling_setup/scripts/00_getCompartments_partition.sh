#!/bin/sh

#SBATCH --job-name rewrite
#SBATCH -n 1                   # Number of cores. For now 56 is the number max of core available
#SBATCH --mem 15Gb
#SBATCH -t 4-00:00              # Runtime in D-HH:MM
#SBATCH -o rewrite.out # File to which STDOUT will be written
#SBATCH -e rewrite.out # File to which STDERR will be written

for nCopies in 1 ; #2 4 6 8 12 20 40 ;
do
    outFile5000=_compartments_5000bp_${nCopies}_copies
    if [[ ! -e ${outFile5000} ]];
    then
	grep "set atom" ../Modelling_by_mechanism/05_ABcompartmentalization_increased_stiffness/mESC_phi_0.10_rosette/EAA_0.010_EBB_0.050_fA_1.00_fB_1.00_lkA_8.00_lkB_4.00_${nCopies}_copies/replica_1/log.lammps | grep "type 2" | awk '{print $3,"A"}' >  scripts/${outFile5000}
	grep "set atom" ../Modelling_by_mechanism/05_ABcompartmentalization_increased_stiffness/mESC_phi_0.10_rosette/EAA_0.010_EBB_0.050_fA_1.00_fB_1.00_lkA_8.00_lkB_4.00_${nCopies}_copies/replica_1/log.lammps | grep "type 3" | awk '{print $3,"B"}' >> scripts/${outFile5000}
	awk -v nc=${nCopies} 'BEGIN{for(i=1;i<=4000*nc;i++){h[i]="T"}}{h[$1]=$2}END{for(i in h){print i,h[i]}}' scripts/${outFile5000} | sort -k 1,1n > _tmp ; mv _tmp scripts/${outFile5000}
    
	rsync -avz scripts/${outFile5000} 03_estimate_time_conversion_mESC_phi_0.10_rosette_${nCopies}_copies/scripts/
	echo T $(grep T scripts/${outFile5000} | wc -l)
	echo A $(grep A scripts/${outFile5000} | wc -l)
	echo B $(grep B scripts/${outFile5000} | wc -l)
	echo
	awk '{if(NR==1){cp=$2;pp=$1;fp=$1}else{if(cp!=$2 || (pp)%4000==0){lp=pp; print fp"-"lp,cp; fp=$1}; cp=$2;pp=$1}}END{lp=pp;print fp"-"lp,cp}' scripts/${outFile5000}
	echo
    fi

    outFile10000=_compartments_10000bp_${nCopies}_copies
    if [[ ! -e ${outFile10000} ]];
    then
	awk '{print int(($1-1)/2),$2}' scripts/${outFile5000} | uniq > scripts/${outFile10000}
	rsync -avz scripts/${outFile10000} 03_estimate_time_conversion_mESC_phi_0.10_rosette_${nCopies}_copies/scripts/
	echo T $(grep T scripts/${outFile10000} | wc -l)
	echo A $(grep A scripts/${outFile10000} | wc -l)
	echo B $(grep B scripts/${outFile10000} | wc -l)   
	echo
	awk '{if(NR==1){cp=$2;pp=$1;fp=$1}else{if(cp!=$2 || (pp+1)%2000==0){lp=pp; print fp"-"lp,cp; fp=$1}; cp=$2;pp=$1}}END{lp=pp;print fp"-"lp,cp}' scripts/${outFile10000}
	echo
    fi
    rsync -avz scripts/${outFile10000} scripts/${outFile5000} ../Modelling_by_mechanism/scripts/ --progress
    
done # Close cycle over ${nCopies}
