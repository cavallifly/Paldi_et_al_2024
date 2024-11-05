nc=$1

for phi in 0.10 ;
do
    if [[ ! -d mESC_phi_${phi}_rosette ]];
    then
	continue
    fi
    cd mESC_phi_${phi}_rosette
    
    #for dir in $(ls -1 | grep _${nc}_copies | grep "_16.00\|32.00" | grep 3.00);
    #for dir in $(ls -1 | grep _${nc}_copies | grep EAA_0.080_EBB_0.000_fA_1.00_fB_0.00_lkA_8.00_lkB_3.00);
    for dir in $(ls -1 | grep _${nc}_copies | grep 16.00);        
    do
	if [[ ! -d ${dir} ]];
	then
	    continue
	fi
	cd ${dir}
	pwd
	for rc in 150 ;
	do
	    for t in $(seq 5400000 5400000 129600000);
	    do	    
		if [[ ! -e $(ls -1 rc_${rc}nm_res_10000bp_from_*_to_${t}_every_900000_replica_2_CS*_perChrom.tab 2> /dev/null) ]];
		then
		    continue
		fi
		th=$(echo $t | awk '{print int($1/5400000)"h"}')
		#cat rc_150nm_res_10000bp_from_*_to_${t}_*_CS*_perChrom.tab | grep -v transChromContactsFraction | awk '{print $NF}' | head

		for compartment in A B;
		do
		    cat rc_150nm_res_10000bp_from_*_to_${t}_*_CS*_perChrom.tab | grep -v "#" | awk -v c=${compartment} '{if($1==c)print $NF}' | sort -k 1,1n > _tmp
		    nlines=$(wc -l _tmp | awk '{print $1}')
		    
		    med=$(awk -v nl=${nlines} -v c=${compartment} 'BEGIN{s=0;n1=int(nl/2); if(nl%2==1){n2=n1}else{n2=n1+1}}{if(NR==n1){s+=$1}; if(NR==n2){s+=$1}}END{print c,"Median CS of ",nl,s/2}' _tmp)
		    echo $th $med
		done
	    done
	done

	#rm -fvr _tmp
	cd ..
    done
    cd ..
done

echo
nlines=19
for compartment in A B;
do
    cat /zdata/data/mdistefano/2023_10_22_Project_hyperacetylation/Modelling_by_mechanism/microC_CS_perChromosome/compartment_strength_HomoOverAll_perChrom.tab | grep -v "chrX\|TSA\|None" | grep _${compartment} | awk '{print $4}' | sort -k 1,1n > _tmp
    nlines=$(wc -l _tmp | awk '{print $1}')
    awk -v nl=${nlines} -v c=${compartment} 'BEGIN{s=0;n1=int(nl/2); if(nl%2==1){n2=n1}else{n2=n1+1}}{if(NR==n1){s+=$1}; if(NR==n2){s+=$1}}END{print "DMSO",c,"Median CS of ",nl,s/2}' _tmp
done
for compartment in A B;
do
    cat /zdata/data/mdistefano/2023_10_22_Project_hyperacetylation/Modelling_by_mechanism/microC_CS_perChromosome/compartment_strength_HomoOverAll_perChrom.tab | grep -v "chrX\|DMSO\|None\|REC" | grep _${compartment} | awk '{print $4}' | sort -k 1,1n > _tmp
    nlines=$(wc -l _tmp | awk '{print $1}')	    
    awk -v nl=${nlines} -v c=${compartment} 'BEGIN{s=0;n1=int(nl/2); if(nl%2==1){n2=n1}else{n2=n1+1}}{if(NR==n1){s+=$1}; if(NR==n2){s+=$1}}END{print "TSA",c,"Median CS of ",nl,s/2}' _tmp
done
for compartment in A B;
do
    cat /zdata/data/mdistefano/2023_10_22_Project_hyperacetylation/Modelling_by_mechanism/microC_CS_perChromosome/compartment_strength_HomoOverAll_perChrom.tab | grep -v "chrX\|DMSO\|None" | grep REC | grep _${compartment} | awk '{print $4}' | sort -k 1,1n > _tmp
    nlines=$(wc -l _tmp | awk '{print $1}')	    
    awk -v nl=${nlines} -v c=${compartment} 'BEGIN{s=0;n1=int(nl/2); if(nl%2==1){n2=n1}else{n2=n1+1}}{if(NR==n1){s+=$1}; if(NR==n2){s+=$1}}END{print "24hREC",c,"Median CS of ",nl,s/2}' _tmp
done	
echo
echo
echo
rm _tmp
