#!/bin/bash
outfile=trajectory.lammpstrj

tag=$1
echo "Selected tag: ${tag}"
mode=$2
echo "Selected mode: ${mode}"
delta=$3
echo "Delta beteew snapshots: ${delta}"

rm -f -r  ${outfile}
rm -f -r _tmp _tmp2

#ls -lrth *${tag}*.lammpstrj | awk '{print $NF}' > _tmp
ls -lrth *${tag}*.XYZ *${tag}*.lammpstrj 2> /dev/null | awk '{print $NF}' > _tmp

cat _tmp

file=$(head -1 _tmp)
L=$(head -4 ${file} | tail -1 | awk '{print $1}')
head -4 $file
echo $L
for file in $(cat _tmp); 
do
    echo ${file}
    N=$(echo $L | awk '{print $1+9}')
    echo $N $L

    if [[ $mode == "initial" ]];
    then
	echo "Mode: initial"
	cat $file | head -$N | awk -v N=$N -v L=$L -v f=$file '{if((NR % N) == 9){print L; print f}; if((NR % N) >= 10 || (NR % N) == 0){if(NF==5){printf(" CA  %15.6f %15.6f %15.6f\n",$3,$4,$5)};if(NF==4){printf(" CA  %15.6f %15.6f %15.6f\n",$2,$3,$4)}}}' >> ${outfile}
	continue
    fi
    if [[ $mode == "final" ]];
    then
	echo "Mode: final"
	cat $file | tail -$N | awk -v N=$N -v L=$L -v f=$file '{if((NR % N) == 9){print L; print f}; if((NR % N) >= 10 || (NR % N) == 0){if(NF==5){printf(" CA  %15.6f %15.6f %15.6f\n",$3,$4,$5)};if(NF==4){printf(" CA  %15.6f %15.6f %15.6f\n",$2,$3,$4)}}}' >> ${outfile}
	continue
    fi
    if [[ $mode == "all" ]];
    then
	cat $file | awk -v N=$N -v L=$L -v fi=$file -v d=${delta} '{if((NR % N) == 2){f=1;if($1%d==0){f=0;t=$1}}; if((NR % N) == 9){if(NR>9 && f==0){print L; print t-d; for(i=1;i<=L;i++){printf(" CA  %15.6f %15.6f %15.6f\n",x[i],y[i],z[i])}}}; if((NR % N) >= 10 || (NR % N) == 0){if(NF==5 && f==0){x[$1]=$3;y[$1]=$4;z[$1]=$5};if(NF==4 && f==0){x[$1]=$2;y[$1]=$3;z[$1]=$4}}}END{print L; print f; for(i=1;i<=L;i++){printf(" CA  %15.6f %15.6f %15.6f\n",x[i],y[i],z[i])}}' >> ${outfile}
    fi
done	
rm -f -r _tmp 
