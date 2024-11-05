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
prevDir=01_compression_to_mESC_phi_${phi}_rosette_${ncopies}_copies

R=0.5
Vol=$(awk -v R=$R -v PI=$PI 'BEGIN{print 4./3.*PI*R*R*R}')
nparticles=8000
echo $R $Vol $nparticles $(echo $nparticles $R | awk '{print "Lc = "2*$2*$1" sigma"}')

Rnew=$(echo $factor $R | awk '{printf("%.10f",$2/(($1)^(1/3)))}')
Vol=$(awk -v f=$factor -v R=$Rnew -v PI=$PI 'BEGIN{print f*4./3.*PI*R*R*R}')
nparticles=$(echo $nparticles $factor | awk '{print $1*$2}')
echo $Rnew $Vol $nparticles $(echo $nparticles $Rnew | awk '{print "Lc = "2*$2*$1}')
sigma=$(echo $Rnew | awk '{print 2*$1}')

chrlength=$(cat ${CGfile} | grep "DNA content to simulate" | awk '{print $5}' | sed "s/bp//g")
size=$(echo ${chrlength} ${ncopies} ${resolution} | awk '{print int($2*$1/$3)}')
echo "Length of the chain ${chrlength} bp - ${size} ${resolution}bp-beads"
echo $b $lk $Tradius $chrlength ${size}

tag=TSA
resolution=5000 
FSmodel=20nm
cellType=mESC

CGfile=../coarse_graining/coarse_graining_FS_${FSmodel}_CG_nu_${resolution}bp_Mouse_${cellType}_phi_${phi}_${tag}_${ncopies}_copies.txt

b=$(cat ${CGfile} | grep -w bCG | awk '{print $9}')
lk=$(cat ${CGfile} | grep -w lkCG | awk '{print $10}')
radius=$(cat ${CGfile} | grep -w radiusCG | awk '{print $16}')

nparticles=4000

chrlength=$(cat ${CGfile} | grep "DNA content to simulate" | awk '{print $5}' | sed "s/bp//g")
size=$(echo ${chrlength} ${ncopies} ${resolution} | awk '{print int($2*$1/$3)}')
echo "Length of the chain ${chrlength} bp - ${size} ${resolution}bp-beads"
echo $b $lk $radius $chrlength ${size}

for replica in $(seq $2 1 $2);
do

    if [[ ! -e ../${prevDir}/replica_${replica}/langevin_dynamics_100000.XYZ ]];
    then
	echo "The file ../${prevDir}/replica_${replica}/langevin_dynamics_100000.XYZ is not available!"
	continue
    fi
    
    replicadir=replica_${replica}
    if [[ -d ${replicadir} ]];
    then
	continue
    fi

    echo $replicadir
    mkdir -p ${replicadir}
    cd ${replicadir}
    
    seed=$(od -An -N3 -i /dev/urandom)
    pyfile=02_generate_relaxed.py
    initial_conformation=../../${prevDir}/replica_${replica}/final_conformation.txt    
    extrusionTime=2500   # Use 1kb / s
    runtime=$(echo $extrusionTime | awk '{print int($1*40)}')      # Check that you can extrude at least the size of the target loops
    n=0

    sed -e "s/XXXreplicaXXX/${replica}/g" -e "s/XXXruntimeXXX/${runtime}/g" -e "s/XXXextrusionTimeXXX/${extrusionTime}/g" -e "s,XXXwdirXXX,${PWD},g" -e "s,XXXreset_timestepXXX,${reset_timestep},g" -e "s,XXXinitial_conformationXXX,${initial_conformation},g" -e "s/XXXseedXXX/${seed}/g" -e "s/XXXpressureXXX/${pressure}/g" -e "s/XXXlkXXX/${lk}/g" -e "s/XXXbXXX/${b}/g" -e "s/XXXnparticlesXXX/${nparticles}/g" -e "s/XXXradiusXXX/${radius}/g" -e "s/XXXsigmaXXX/${sigma}/g" -e "s/XXXncopiesXXX/${ncopies}/g" ../../scripts/${pyfile} > ${pyfile%.py}_${replica}.py
	
    n=$((${n}+1))
    python ${pyfile%.py}_${replica}.py &>> replica_${replica}.out
    
    cd .. # Exit ${replicadir}
    
done # Close cycle over $replica
