#!/bin/bash

#SBATCH --job-name p0.10
#SBATCH -n 1
#SBATCH --mem-per-cpu 15Gb
#SBATCH -t 2-00:00              # Runtime in D-HH:MM
#SBATCH -o phi_0_10.out # File to which STDOUT will be written
#SBATCH -e phi_0_10.out # File to which STDERR will be written 


factor=1

resolution=$1
nparticles=$2
ncopies=$3
replica=$4

tag=TSA
FSmodel=20nm
phi=0.10
cellType=mESC

CGfile=../coarse_graining/coarse_graining_FS_${FSmodel}_CG_nu_${resolution}bp_Mouse_${cellType}_phi_${phi}LV_${tag}_${ncopies}_copies.txt

b=$(cat ${CGfile} | grep -w bCG | awk '{print $9}')
lk=$(cat ${CGfile} | grep -w lkCG | awk '{print $10}')
Tradius=$(cat ${CGfile} | grep -w radiusCG | awk '{print $16}')

factor=1
R=0.5
Vol=$(awk -v R=$R -v PI=$PI 'BEGIN{print 4./3.*PI*R*R*R}')
echo $R $Vol $nparticles $(echo $nparticles $R | awk '{print "Lc = "2*$2*$1" sigma"}')

Rnew=$(echo $factor $R | awk '{printf("%.10f",$2/(($1)^(1/3)))}')
Vol=$(awk -v f=$factor -v R=$Rnew -v PI=$PI 'BEGIN{print f*4./3.*PI*R*R*R}')
nparticles=$(echo $nparticles $factor | awk '{print $1*$2}')
echo $Rnew $Vol $nparticles $(echo $nparticles $Rnew | awk '{print "Lc = "2*$2*$1}')
sigma=$(echo $Rnew | awk '{print 2*$1}')

chrlength=$(cat ${CGfile} | grep "DNA content to simulate" | awk '{print $5}' | sed "s/bp//g")
size=$(echo ${chrlength} ${ncopies} ${resolution} | awk '{print int($2*$1/$3)}')
echo "Length of the chain ${chrlength} bp - ${size} ${resolution}bp-beads"
echo "Bead diameter $b nm - Kuhn-legth $lk nm - Target radius $Tradius b - Ncopies ${ncopies} - Nparticles ${size}"

prevDir=00_generate_initial_conformation_rosette_phi_0.03_resolution_${resolution}bp_${ncopies}_copies

for replica in $(seq $replica 1 $replica);
do
    if [[ ! -e ../${prevDir}/Initial_rosette_conformation_sphere_replica_${replica}.dat ]];
    then
	echo "The file ../${prevDir}/Initial_rosette_conformation_sphere_replica_${replica}.dat is not generated!"
	continue
    fi   

    #runtime=100
    runtime=100000    
    replicadir=replica_${replica}
    echo ${Iradius} ${Tradius}
    if [[ -d ${replicadir} ]];
    then
	continue
    fi

    echo $replicadir
    mkdir -p ${replicadir}
    cd ${replicadir}
    
    seed=$(od -An -N3 -i /dev/urandom)
    pyfile=01_compression_to_the_target_density.py

    initial_conformation=../../${prevDir}/Initial_rosette_conformation_sphere_replica_${replica}.dat
    rsync -avz ${initial_conformation} final_conformation.txt

    Iradius=$(awk -v s=${sigma} '{if(NF==6){dx2=$4*$4;dy2=$5*$5;dz2=$6*$6;d=sqrt(dx2+dy2+dz2);if(d>max){max=d}}}END{print max+s*1.2}' final_conformation.txt)
    awk -v r=${Iradius} '{if($3=="xlo" || $3=="ylo" || $3=="zlo"){print -r,r,$3,$4}else{print $0}}' final_conformation.txt > _tmp ; mv _tmp final_conformation.txt
    echo $Iradius

    #for radius in $(seq ${Iradius} -0.1 ${Tradius})
    for radius in ${Iradius} ;
    do
	reset_timestep=$(echo $n | awk -v rt=${runtime} '{print $1*rt}')
	echo ${reset_timestep}
	sed -e "s/XXXsigmaXXX/${sigma}/g" -e "s/XXXreplicaXXX/${replica}/g" -e "s/XXXruntimeXXX/${runtime}/g" -e "s,XXXwdirXXX,${PWD},g" -e "s,XXXreset_timestepXXX,${reset_timestep},g" -e "s,XXXinitial_conformationXXX,${initial_conformation},g" -e "s/XXXseedXXX/${seed}/g" -e "s/XXXlkXXX/${lk}/g" -e "s/XXXbXXX/${b}/g" -e "s/XXXnparticlesXXX/${nparticles}/g" -e "s/XXXradiusXXX/${radius}/g" -e "s/XXXIradiusXXX/${Iradius}/g" -e "s/XXXTradiusXXX/${Tradius}/g" ../../scripts/${pyfile} > ${pyfile%.py}_${replica}.py
	
	#n=$((${n}+1))
	python ${pyfile%.py}_${replica}.py &>> replica_${replica}.out
	initial_conformation=final_conformation.txt
	bradius=${Tradius}
	awk -v r=${bradius} '{if($3=="xlo" || $3=="ylo" || $3=="zlo"){print -(r+1),(r+1),$3,$4}else{print $0}}' final_conformation.txt > _tmp ; mv _tmp final_conformation.txt	
	rm -fvr minim*
    done
    
    cd .. # Exit ${replicadir}
    
done # Close cycle over $replica
