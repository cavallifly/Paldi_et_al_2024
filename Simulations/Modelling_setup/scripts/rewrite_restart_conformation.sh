#!/bin/sh

#SBATCH --job-name rewrite
#SBATCH -n 1                   # Number of cores. For now 56 is the number max of core available
#SBATCH --mem 15Gb
#SBATCH -t 4-00:00              # Runtime in D-HH:MM
#SBATCH -o rewrite.out # File to which STDOUT will be written
#SBATCH -e rewrite.out # File to which STDERR will be written

mkdir -p relaxed_conformations

nCopies=$1
if [[ $nCopies == "" ]];
then
    echo "Ncopies ${nCopies} not provided!"
    exit
fi

awk '{print $0,"Compartment"}' scripts/_compartments_5000bp_$1_copies > _comp

for replica in $(seq 1 1 100);
do
    echo $replica

    outFile=./relaxed_conformations/relaxed_conformation_replica_${replica}.txt
    if [[ ! -e ${outFile} ]];
    then
	cat _comp ./replica_${replica}/relaxed_conformation.txt | awk 'BEGIN{f=0}{if($NF=="Compartment"){comp[$1]=$2}else{if($2=="angle" && $3=="types"){print 2,$2,$3; next}; if($2=="atom" && $3=="types"){print 3,$2,$3; next}; if(f==0 || NF==0){print $0}else{t=1;if(comp[$3]=="A" && comp[$4]=="A" && comp[$5]=="A"){t=2}; print $1,t,$3,$4,$5}; if($1=="Angles"){f=1}}}' > ${outFile}
	
	#diff replica_${replica}/relaxed_conformation.txt ${outFile}
	
	awk '{if(NF==5 && $2==2){print $0}}' ./replica_${replica}/relaxed_conformation.txt | wc -l
	awk '{if(NF==5 && $2==1){print $0}}' ./replica_${replica}/relaxed_conformation.txt | wc -l
	
	awk '{if(NF==5 && $2==2){print $0}}' ${outFile} | wc -l
	awk '{if(NF==5 && $2==1){print $0}}' ${outFile} | wc -l
    fi
done
rm _comp    
