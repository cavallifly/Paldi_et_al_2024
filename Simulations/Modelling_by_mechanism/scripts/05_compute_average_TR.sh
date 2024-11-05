mainDir=$(pwd)

outFileDist=${mainDir}/mESC_phi_0.10_rosette/distances_to_TR.txt
outFileDistA=${mainDir}/mESC_phi_0.10_rosette/distances_to_TRA.txt
outFileDistB=${mainDir}/mESC_phi_0.10_rosette/distances_to_TRB.txt
rm -fvr ${outFileDist} ${outFileDistA} ${outFileDistB}

# Get median values per each experimental condition for comparisons
cd mESC_phi_0.10_rosette

nc=$1

for dir in $(ls -1 | grep lkB | grep _${nc}_copies | grep 0.080) ; # | grep -v "lkB_0.50\|lkB_0.25\|lkB_0.75\|lkA_16.00")
do
    if [[ ! -d ${dir} ]];
    then
	continue
    fi
    cd ${dir}
    for t in $(seq 5400000 5400000 129600000);
    do
	if [[ ! -e $(ls -1 rc_150nm_res_10000bp_from_*_to_${t}_every_900000_replica_2_TR_perChrom.tab 2> /dev/null) ]];
	then
	    continue
	fi
	th=$(echo $t | awk '{print int($1/5400000)"h"}')
	cat rc_150nm_res_10000bp_from_*_to_${t}_*_TR_perChrom.tab | grep -v transChromContactsFraction | awk '{print $NF}' | sort -k 1,1n > _tmp
	nlines=$(wc -l _tmp | awk '{print $1}')
	med=$(awk -v nl=${nlines} 'BEGIN{s=0;n1=int(nl/2)+1; if(nl%2==1){n2=n1}else{n2=n1-1}}{if(NR==n1){s+=$1}; if(NR==n2){s+=$1}}END{print "Median TR of ",nl,s/2}' _tmp)
	#echo $th $med
	medAll=$(echo $med | awk '{print $NF}')
	
	cat rc_150nm_res_10000bp_from_*_to_${t}_*_TRA_perChrom.tab | grep -v transChromContactsFraction | awk '{print $NF}' | sort -k 1,1n > _tmp
	nlines=$(wc -l _tmp | awk '{print $1}')
	med=$(awk -v nl=${nlines} 'BEGIN{s=0;n1=int(nl/2)+1; if(nl%2==1){n2=n1}else{n2=n1-1}}{if(NR==n1){s+=$1}; if(NR==n2){s+=$1}}END{print "Median TR of ",nl,s/2}' _tmp)
	#echo A $th $med
	medA=$(echo $med | awk '{print $NF}')
	
	cat rc_150nm_res_10000bp_from_*_to_${t}_*_TRB_perChrom.tab | grep -v transChromContactsFraction | awk '{print $NF}' | sort -k 1,1n > _tmp
	nlines=$(wc -l _tmp | awk '{print $1}')
	med=$(awk -v nl=${nlines} 'BEGIN{s=0;n1=int(nl/2)+1; if(nl%2==1){n2=n1}else{n2=n1-1}}{if(NR==n1){s+=$1}; if(NR==n2){s+=$1}}END{print "Median TR of ",nl,s/2}' _tmp)
	#echo B $th $med
	medB=$(echo $med | awk '{print $NF}')
	analyzedChains=$(echo $med | awk '{print $(NF-1)}')
	
	echo $dir $t ${analyzedChains}
	
	for condition in DMSO TSA TSA24hREC ;
	do
	    medianAll=$(cat ${mainDir}/../microC_TR_perChromosome/median_TR_perChrom.txt  | grep -w ${condition} | awk '{print $NF}')
	    medianA=$(  cat ${mainDir}/../microC_TR_perChromosome/median_TRA_perChrom.txt | grep -w ${condition} | awk '{print $NF}')
	    medianB=$(  cat ${mainDir}/../microC_TR_perChromosome/median_TRB_perChrom.txt | grep -w ${condition} | awk '{print $NF}')
	    #echo ${medianAll} ${medianA} ${medianB}

	    distanceAll=$(echo ${medianAll} ${medAll} | awk '{dp=$1-$2; print sqrt(dp*dp)}')
	    distanceA=$(  echo ${medianA}   ${medA}   | awk '{dp=$1-$2; print sqrt(dp*dp)}')
	    distanceB=$(  echo ${medianB}   ${medB}   | awk '{dp=$1-$2; print sqrt(dp*dp)}')
	    #echo ${distanceAll} ${distanceA} ${distanceB}
	    
	    echo $dir $t ${condition} ${distanceAll} >> ${outFileDist}
	    echo $dir $t ${condition} ${distanceA}   >> ${outFileDistA}
	    echo $dir $t ${condition} ${distanceB}   >> ${outFileDistB}	    
	done
	
	
	rm -fr _tmp
    done
    cd ..
done
cd ..

sort -k 3,3n -k 4,4n ${outFileDist} > _tmp ; mv _tmp ${outFileDist}
sort -k 3,3n -k 4,4n ${outFileDistA} > _tmp ; mv _tmp ${outFileDistA}
sort -k 3,3n -k 4,4n ${outFileDistB} > _tmp ; mv _tmp ${outFileDistB}
#head ${outFileDist} ${outFileDistA} ${outFileDistB}

echo
for condition in DMSO TSA TSA24hREC ;
do
    echo $condition
    outFileRankAll=${mainDir}/mESC_phi_0.10_rosette/ranking_distances_to_TR_${condition}.txt
    outFileRankA=${mainDir}/mESC_phi_0.10_rosette/ranking_distances_to_TRA_${condition}.txt
    outFileRankB=${mainDir}/mESC_phi_0.10_rosette/ranking_distances_to_TRB_${condition}.txt
    outFileRankMerged=${mainDir}/mESC_phi_0.10_rosette/ranking_distances_to_TR_${condition}_merged.txt
    
    cat ${outFileDist}  | grep -w ${condition} | sort -k 4,4n | grep -w 21600000 | awk '{h[$1]=$4"-"NR}END{for(i in h){print i,h[i],"All"}}' | sort -k 1,1n | sed "s/-/ /g" > ${outFileRankAll}
    cat ${outFileDistA} | grep -w ${condition} | sort -k 4,4n | grep -w 21600000 | awk '{h[$1]=$4"-"NR}END{for(i in h){print i,h[i],"A"}}'   | sort -k 1,1n | sed "s/-/ /g" > ${outFileRankA}
    cat ${outFileDistB} | grep -w ${condition} | sort -k 4,4n | grep -w 21600000 | awk '{h[$1]=$4"-"NR}END{for(i in h){print i,h[i],"B"}}'   | sort -k 1,1n | sed "s/-/ /g" > ${outFileRankB}
    
    paste ${outFileRankAll} ${outFileRankA} ${outFileRankB} | awk '{if($1==$5 && $1==$9){printf("%s\t%.3f\t%d\t%.3f\t%d\t%.3f\t%d\t%.2f\n",$1,$2,$3,$6,$7,$10,$11,($3+$7+$11)/3)}}' > ${outFileRankMerged}
    sort -k 8,8n ${outFileRankMerged} > _tmp ; mv _tmp ${outFileRankMerged}
    head -5 ${outFileRankMerged}
    echo
done
echo
