#!/bin/bash

#SBATCH --job-name p0.10
#SBATCH -n 1                    # Number of cores. For now 56 is the number max of core available
#SBATCH --mem-per-cpu 15Gb
#SBATCH -t 2-00:00              # Runtime in D-HH:MM
#SBATCH -o phi_0_10.out # File to which STDOUT will be written
#SBATCH -e phi_0_10.out # File to which STDERR will be written 

factor=1
ncopies=$1
phi=0.10
prevDir=02_generate_relaxed_mESC_phi_${phi}_rosette_${ncopies}_copies

R=0.5
Vol=$(awk -v R=$R -v PI=$PI 'BEGIN{print 4./3.*PI*R*R*R}')
nparticles=8000
echo $R $Vol $nparticles $(echo $nparticles $R | awk '{print "Lc = "2*$2*$1" sigma"}')

Rnew=$(echo $factor $R | awk '{printf("%.10f",$2/(($1)^(1/3)))}')
Vol=$(awk -v f=$factor -v R=$Rnew -v PI=$PI 'BEGIN{print f*4./3.*PI*R*R*R}')
nparticles=$(echo $nparticles $factor | awk '{print $1*$2}')
echo $Rnew $Vol $nparticles $(echo $nparticles $Rnew | awk '{print "Lc = "2*$2*$1}')
sigma=$(echo $Rnew | awk '{print 2*$1}')

tag=TSA
resolution=5000 
FSmodel=20nm
cellType=mESC

CGfile=../coarse_graining/coarse_graining_FS_${FSmodel}_CG_nu_${resolution}bp_Mouse_${cellType}_phi_${phi}_${tag}_${ncopies}_copies.txt

chrlength=$(cat ${CGfile} | grep "DNA content to simulate" | awk '{print $5}' | sed "s/bp//g")
size=$(echo ${chrlength} ${ncopies} ${resolution} | awk '{print int($2*$1/$3)}')
echo "Length of the chain ${chrlength} bp - ${size} ${resolution}bp-beads"
echo $b $lk $Tradius $chrlength ${size}

chrlength=20000000
size=4000

b=$(cat ${CGfile} | grep -w bCG | awk -v s=${sigma} '{print $9*s}')
lk=$(cat ${CGfile} | grep -w lkCG | awk '{print $10}')
radius=$(cat ${CGfile} | grep -w radiusCG | awk '{print $16}')

timestep=0.012

echo $b $lk

echo "Length of the chain ${chrlength} bp - ${size} ${resolution}bp-beads"

for replica in $(seq $2 1 $2) ;
do
    # Decompaction
    replicadir=replica_${replica}
    if [[ ! -e ../${prevDir}/replica_${replica}/final_conformation.txt ]];
    then
	ls -lrtha ../${prevDir}/replica_${replica}/final_conformation.txt
	exit
    fi
    
    if [[ ! -d ${replicadir} ]];
    then
	
	echo ${replicadir}
	mkdir -p ${replicadir}
	cd ${replicadir}

	pyScript=03_estimate_time_conversion.py
	
	# 50,000,000 relaxation to collect MSD
	run=$(echo 5000000 ${timestep} | awk '{print int($1*(0.012/$2))}')

	rsync -avz ../../${prevDir}/replica_${replica}/final_conformation.txt .
	
	sed -e "s/XXXsigmaXXX/${sigma}/g" -e "s/XXXbXXX/${b}/g" -e "s/XXXradiusXXX/${radius}/g" -e "s/XXXlkXXX/${lk}/g" -e "s/XXXrunXXX/${run}/g" -e "s/XXXreplicaXXX/${replica}/g" -e "s/XXXnparticlesXXX/${nparticles}/g" -e "s/XXXtimestepXXX/${timestep}/g" ../../scripts/${pyScript} > ${pyScript%.py}_${replica}.py

	python 03_estimate_time_conversion_${replica}.py &>> replica_${replica}.out
	rm -fvr mini*XYZ
	
	cd .. # Exit ${replicadir}
    fi
done # Close cycle over replica
