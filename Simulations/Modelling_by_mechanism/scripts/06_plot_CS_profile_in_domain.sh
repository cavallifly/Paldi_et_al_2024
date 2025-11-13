mainDir=$(pwd)

checkSymmetry=$(echo $mainDir | grep "05_AB\|06_AB" | wc -l | awk '{print $1}')
echo $checkSymmetry

nc=20
minReplicates=5
if [[ $nc -eq 1 ]];
then
    minReplicates=30
fi

cd mESC_phi_0.10_rosette

for condition in DMSO TSA ;
do
    tag=HomoOverAll
	    
    inMicrocFileA=$(ls -1 ../../../microC_CS_perChromosome/CS_${tag}_within_domains_rescaled_to_150_${condition}_*_A.txt)
    inMicrocFileB=$(ls -1 ../../../microC_CS_perChromosome/CS_${tag}_within_domains_rescaled_to_150_${condition}_*_B.txt)
    outFileDist=${mainDir}/mESC_phi_0.10_rosette/distances_to_CS_${tag}_within_domains_rescaled_to_150_${condition}_${nc}_copies.txt
    outFileDistA=${mainDir}/mESC_phi_0.10_rosette/distances_to_CS_${tag}_within_domains_rescaled_to_150_${condition}_A_${nc}_copies.txt
    outFileDistB=${mainDir}/mESC_phi_0.10_rosette/distances_to_CS_${tag}_within_domains_rescaled_to_150_${condition}_B_${nc}_copies.txt	    
    outFileRankA=${mainDir}/mESC_phi_0.10_rosette/ranking_distances_to_CS_${tag}_within_domains_rescaled_to_150_${condition}_A_${nc}_copies.txt
    outFileRankB=${mainDir}/mESC_phi_0.10_rosette/ranking_distances_to_CS_${tag}_within_domains_rescaled_to_150_${condition}_B_${nc}_copies.txt
    outFileRank=${mainDir}/mESC_phi_0.10_rosette/ranking_distances_to_CS_${tag}_within_domains_rescaled_to_150_${condition}_${nc}_copies.txt
    
    ls -lrtha ${inMicrocFileA} ${inMicrocFileB}

    rm -fvr ${outFileDist} ${outFileDistA} ${outFileDistB}
    
    for dir in $(ls -1 | grep _${nc}_copies | grep -v tau);
    do
	echo $dir
	if [[ ! -d ${dir} ]];
	then
	    continue
	fi
	
	cd $dir
	pwd
	
	
	for t in $(seq 5400000 5400000 129600000);
	#for t in 21600000 ;
	do
	    rm -fvr _?_${t}
	    from=$(echo $t | awk '{printf("%d\n",$1-5400000)}')
	    
	    nFiles=$(ls -1 rc_150nm_res_10000bp_from_${from}_to_${t}_every_900000_replica_?_0_obsOverExp_CS${tag}.tab 2> /dev/null | wc -l | awk '{print $1}')
	    if [[ ${nFiles} -lt ${minReplicates} ]];
	    then
		continue
	    fi
	    for replica in $(seq 1 1 ${minReplicates});
	    do
		for file in $(ls -1 rc_150nm_res_10000bp*to_${t}_every_900000_replica_${replica}_*_obsOverExp_CS${tag}.tab 2> /dev/null);
		do
		    echo $file
		    grep -w A ${file} | awk '{if(NR==1){pp=$1;n=1}else{if($1-pp!=1){n=1}}; print n++,$1,$2,$3,$4,$5,$6,$7; pp=$1}' >> _A_${t}
		    grep -w B ${file} | awk '{if(NR==1){pp=$1;n=1}else{if($1-pp!=1){n=1}}; print n++,$1,$2,$3,$4,$5,$6,$7; pp=$1}' >> _B_${t}
		done
	    done
	    tail -150 ${inMicrocFileA} > _A_${condition}
	    tail -150 ${inMicrocFileB} > _B_${condition}	
	    
	    dirName=$(echo $dir | sed "s/_/ /g" | awk '{printf("%s_%s_%s_%s",$1,$2,$3,$4); for(i=5;i<=NF;i++){printf("_%s",$i)}; printf("\n")}')
	    outFile=rc_150nm_res_10000bp_from_${from}_to_${t}_every_900000_CS${tag}_within_domains_${dirName}_${condition}.png
	    if [[ ! -e ${outFile} ]];
	    then
		cat _A_${t} | awk '{h[$1]+=$NF;h2[$1]+=$NF*$NF;c[$1]++}END{for(i in h){avg=h[i]/c[i]; avg2=h2[i]/c[i]; stddev=sqrt(avg2-avg*avg); print i,avg,stddev,c[i]}}' | sort -k 1,1n > _A
		cat _B_${t} | awk '{h[$1]+=$NF;h2[$1]+=$NF*$NF;c[$1]++}END{for(i in h){avg=h[i]/c[i]; avg2=h2[i]/c[i]; stddev=sqrt(avg2-avg*avg); print i,avg,stddev,c[i]}}' | sort -k 1,1n > _B
		
		sed "s/XXXconditionXXX/${condition}/g" ../../../scripts/XX_plot_CS_profile_in_domain.gp	| gnuplot
		mv data.ps ${outFile%.png}.ps
		ps2pdf ${outFile%.png}.ps
		rm ${outFile%.png}.ps
		#~/pdfcrop.pl --margin 0 ${outFile%.png}.pdf _tmp_${outFile%.png}.pdf
		#mv _tmp_${outFile%.png}.pdf ${outFile%.png}.pdf
		convert -background white -alpha remove -density 600 ${outFile%png}pdf ${outFile%png}png
		#conda run -n misha bash ../../../scripts/ps2pdf.sh #&> /dev/null	    
	    fi
	    cat _A_${t} | awk '{h[$1]+=$NF;h2[$1]+=$NF*$NF;c[$1]++}END{for(i in h){avg=h[i]/c[i]; avg2=h2[i]/c[i]; stddev=sqrt(avg2-avg*avg); print i,avg,stddev,c[i]}}' | sort -k 1,1n > _A
	    cat _B_${t} | awk '{h[$1]+=$NF;h2[$1]+=$NF*$NF;c[$1]++}END{for(i in h){avg=h[i]/c[i]; avg2=h2[i]/c[i]; stddev=sqrt(avg2-avg*avg); print i,avg,stddev,c[i]}}' | sort -k 1,1n > _B
	    head _?
	    
	    mv _A rc_150nm_res_10000bp_from_${from}_to_${t}_every_900000_CS${tag}_within_domains_A.txt
	    mv _B rc_150nm_res_10000bp_from_${from}_to_${t}_every_900000_CS${tag}_within_domains_B.txt	
	    
	    dirName=$(echo $dir | sed "s/_/ /g" | awk '{printf("%s_%s_%s_%s",$1,$2,$3,$4); for(i=5;i<=NF;i++){printf("_%s",$i)}; printf("\n")}')
	    paste rc_150nm_res_10000bp_from_${from}_to_${t}_every_900000_CS${tag}_within_domains_A.txt <(tail -150 ${inMicrocFileA}) | awk -v dir=${dir} -v t=${t} 'BEGIN{d=0}{if($1==$5){d+=sqrt(($2-$6)*($2-$6)); cnt++}}END{print "A",dir,t,d/cnt,d,cnt}' >> ${outFileDist}
	    paste rc_150nm_res_10000bp_from_${from}_to_${t}_every_900000_CS${tag}_within_domains_B.txt <(tail -150 ${inMicrocFileB}) | awk -v dir=${dir} -v t=${t} 'BEGIN{d=0}{if($1==$5){d+=sqrt(($2-$6)*($2-$6)); cnt++}}END{print "B",dir,t,d/cnt,d,cnt}' >> ${outFileDist}
	    head ${outFileDist}
	    
	    if [[ ${checkSymmetry} -eq 0 ]];
	    then
		dirName=$(echo $dir | sed "s/_/ /g" | awk '{printf("%s_%s_%s_%s",$1,$4,$3,$2); for(i=5;i<=NF;i++){printf("_%s",$i)}; printf("\n")}')
		paste rc_150nm_res_10000bp_from_${from}_to_${t}_every_900000_CS${tag}_within_domains_B.txt <(tail -150 ${inMicrocFileA}) | awk -v dirName=${dirName} -v t=${t} 'BEGIN{d=0}{if($1==$5){d+=sqrt(($2-$6)*($2-$6)); cnt++}}END{print "A",dirName,t,d/cnt,d,cnt}' >> ${outFileDist}
		paste rc_150nm_res_10000bp_from_${from}_to_${t}_every_900000_CS${tag}_within_domains_A.txt <(tail -150 ${inMicrocFileB}) | awk -v dirName=${dirName} -v t=${t} 'BEGIN{d=0}{if($1==$5){d+=sqrt(($2-$6)*($2-$6)); cnt++}}END{print "B",dirName,t,d/cnt,d,cnt}' >> ${outFileDist}

		outFile=rc_150nm_res_10000bp_from_${from}_to_${t}_every_900000_CS${tag}_within_domains_${dirName}_${condition}.png
		if [[ ! -e ${outFile} ]];
		then
		    cat _B_${t} | awk '{h[$1]+=$NF;h2[$1]+=$NF*$NF;c[$1]++}END{for(i in h){avg=h[i]/c[i]; avg2=h2[i]/c[i]; stddev=sqrt(avg2-avg*avg); print i,avg,stddev,c[i]}}' | sort -k 1,1n > _A
		    cat _A_${t} | awk '{h[$1]+=$NF;h2[$1]+=$NF*$NF;c[$1]++}END{for(i in h){avg=h[i]/c[i]; avg2=h2[i]/c[i]; stddev=sqrt(avg2-avg*avg); print i,avg,stddev,c[i]}}' | sort -k 1,1n > _B
		    
		    sed "s/XXXconditionXXX/${condition}/g" ../../../scripts/XX_plot_CS_profile_in_domain.gp	| gnuplot
		    mv data.ps ${outFile%.png}.ps
		    ps2pdf ${outFile%.png}.ps
		    rm ${outFile%.png}.ps
		    #~/pdfcrop.pl --margin 0 ${outFile%.png}.pdf _tmp_${outFile%.png}.pdf
		    #mv _tmp_${outFile%.png}.pdf ${outFile%.png}.pdf
		    convert -background white -alpha remove -density 600 ${outFile%png}pdf ${outFile%png}png
		#conda run -n misha bash ../../../scripts/ps2pdf.sh #&> /dev/null	    
		fi
	    fi																			    
	    rm -fvr _?_${t} _? _?_${condition}
	done
	cd ..
	pwd
    done
    
    sort -k 4,4n ${outFileDist} > _tmp ; mv _tmp ${outFileDist}
    head ${outFileDist}
    sort -k 4,4n ${outFileDist} | grep -w A > _tmp ; mv _tmp ${outFileDistA}
    sort -k 4,4n ${outFileDist} | grep -w B > _tmp ; mv _tmp ${outFileDistB}	    
    head ${outFileDistA}  head ${outFileDistB}

    sort -k 4,4n ${outFileDistA} | awk '{print $4}' | awk '{if(NR==1){pr=NR;pv=$1;print $1,pr}else{if($1!=pv){pr=NR}; print $1,pr; pv=$1}}' > _order
    #awk '{if(NF==2){r[$1]=$2}else{print $2,r[$4],"A"}}' _order <( grep -w A ${outFileDistA} | grep -w 21600000) | sort -k 2,2n > ${outFileRankA}
    awk '{if(NF==2){r[$1]=$2}else{print $2,r[$4],"A"}}' _order <( grep -w A ${outFileDistA}) | sort -k 2,2n > ${outFileRankA}
    head ${outFileRankA}
    rm -fvr _order

    sort -k 4,4n ${outFileDistB} | awk '{print $4}' | awk '{if(NR==1){pr=NR;pv=$1;print $1,pr}else{if($1!=pv){pr=NR}; print $1,pr; pv=$1}}' > _order
    #awk '{if(NF==2){r[$1]=$2}else{print $2,r[$4],"B"}}' _order <( grep -w B ${outFileDistB} | grep -w 21600000) | sort -k 2,2n > ${outFileRankB}
    awk '{if(NF==2){r[$1]=$2}else{print $2,r[$4],"B"}}' _order <( grep -w B ${outFileDistB}) | sort -k 2,2n > ${outFileRankB}
    head ${outFileRankB}
    rm -fvr _order

    cat ${outFileRankA} ${outFileRankB} | awk '{h[$1]+=$2; cnt[$1]; if($NF=="A"){rA[$1]=$2}; if($NF=="B"){rB[$1]=$2}}END{for(i in h) print i,21600000,rA[i],rB[i],h[i]/2}' | sort -k 5,5n > ${outFileRank}
    echo ${outFileRank}
    head ${outFileRank}
done

cd ..
