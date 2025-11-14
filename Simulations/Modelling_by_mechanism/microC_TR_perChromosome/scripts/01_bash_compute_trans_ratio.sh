inDir=../01_cool_files/balanced_text_matrices/

mapResolution=10000 # in bp

condition=$1

for file in $(ls -1 ${inDir}/full${condition}_*${mapResolution}bp* 2> /dev/null | grep -v TR);
do
    echo $file
    head ${file}

    name=$(echo $file | sed "s,/, ,g" | awk '{print $NF}')

    #TR definition is based on the observed-over-expected (OoE) map. This map is computed as the entries of the contact or Hi-C map
    #divided by the P(s) for the corresponding genomic distance s. For each bin (b) of a given epigenomic state, the TR was then
    #quantified as the ratio between the average OoE value of b with bins of the same epigenomic state and the average OoE value of
    #b with any other bin in the same chromosome. TR scores of all the bins of a given epigenomic state are then pulled together to
    #form the TR distribution of that state, that were shown in the figures of this work. In this metric, no compartmentalization
    #corresponds to TR = 1, whereas any pattern of compartmentalization yields to TR > 1.

    for mode in ICE;
    do
    
	# Compute TR per bin
	outFile=${name%.tab}_TR${mode}.tab    
	if [[ ! -e ${outFile} ]];
	then
	    touch ${outFile}
	    echo ${outFile}
	
	    echo "#Bin cisChromContacts nCisChromBins transChromContacts nTransChromBins totContacts nTotBins cisChromContactsFraction transChromContactsFraction" | awk '{for(i=1;i<NF;i++){printf("%s\t",$i)}; printf("%s\n",$NF)}' > ${outFile}
	    head -10 ${file}
	    # ICE-normalized scores
	    if [[ ${mode} == "ICE" ]];
	    then
		# ICE-normalized contacts not-normalized by the number of valid-pairs -> Preferred
		awk '{if(NF!=8){next;}; ind1=$1"_"$2"_"$3; ind2=$4"_"$5"_"$6; if($1==$4){cis[ind1]+=$8; Ccis[ind1]++; if(ind1!=ind2){cis[ind2]+=$8; Ccis[ind2]++}}else{trans[ind1]+=$8; Ctrans[ind1]++; if(ind1!=ind2){trans[ind2]+=$8; Ctrans[ind2]++}}; tot[ind1]+=$8; Ctot[ind1]++; if(ind1!=ind2){tot[ind2]+=$8; Ctot[ind2]++}}END{for(i in Ctot){print i,cis[i],Ccis[i],trans[i],Ctrans[i],tot[i],Ctot[i],cis[i]/tot[i],trans[i]/tot[i]}}' ${file} | awk '{for(i=1;i<NF;i++){printf("%s\t",$i)}; printf("%s\n",$NF)}' >> ${outFile}
	    fi
	fi
	
	# Compute TR per chromosome
	inFile=${outFile}
	outFile=${outFile%.tab}_perChrom.tab
	
	echo "#chrom cisChromContactsFraction transChromContactsFraction" | awk '{for(i=1;i<NF;i++){printf("%s\t",$i)}; printf("%s\n",$NF)}' > ${outFile}
	for chrom in $(cat chrom_sizes_mm10_higlass.txt | awk '{printf("%s ",$1)}')
	do
	    echo $chrom
	    
	    sed "s/_/ /g" ${inFile} | grep -w $chrom | awk -v c=${chrom} -v n=${name} '{cis+=$10; trans+=$11; cnt++}END{print n,c,cis/cnt,trans/cnt}' >> ${outFile}	    
	done # Close cycle over ${chrom}
	cat ${outFile}

    done # Close cycle over ${mode}
done # Close cycle over ${file}

# Assign A or B compartment to the TR
for file in $(ls -1 *${mapResolution}bp*TRICE.tab | grep TR);
do
   
    outFile=${file%.tab}_compartment.tab
    echo $file ${outFile}
    if [[ -e ${outFile} ]];
    then
	continue
    fi

    condition=DMSO
    sort -k 1,1 -k 2,2n -k 3,3n <(cat A_compartments_cis_s${condition}.bed | awk '{print $1,$2,$3,"A"}') <(cat B_compartments_cis_s${condition}.bed | awk '{print $1,$2,$3,"B"}') > _tmp
    
    awk '{if(NF==4){n++; chrom[n]=$1;start[n]=$2;end[n]=$3;compartment[n]=$4;}else{for(i=0;i<=n;i++){if(chrom[i]==$1 && (start[i]<=($2+1) && ($2+1)<=end[i])){print $0,chrom[i],start[i],end[i],compartment[i]; next}}}}' _tmp <( sed "s/_/ /g" ${file} | sort -k 1,1 -k 2,2n -k 3,3n) | awk '{printf("%s_%s_%s\t",$1,$2,$3); for(i=4;i<NF;i++){printf("%s\t",$i)}; printf("%s\n",$NF)}' > ${outFile}
    rm _tmp *~ ./*/*~

    # Compute TR per chromosome
    inFile=${outFile}
    outFile=${outFile%.tab}_perChrom.tab

    rm -fvr ${outFile}
    echo "#chrom cisChromContactsFraction transChromContactsFraction" | awk '{for(i=1;i<NF;i++){printf("%s\t",$i)}; printf("%s\n",$NF)}' > ${outFile}
    for compartment in A B;
    do
	for chrom in $(cat /home/Programs/chrom_sizes_mm10_higlass.txt | grep -v "X\|Y\|M" | awk '{printf("%s ",$1)}')
	do
	    echo $chrom

	    name=${file%.tab}
	    sed "s/_/ /g" ${inFile} | grep -w ${compartment} | grep -w $chrom | awk -v c=${chrom} -v n=${name} '{cis+=$10; trans+=$11; cnt++}END{print n,c,cis/cnt,trans/cnt,$NF}' >> ${outFile} 
	done # Close cycle over ${chrom}
	cat ${outFile}
    done
    
done

rm -fvr median_TR_perChrom.txt
for file in $(ls -1 *ICE_perChrom.tab);
do
    echo $file

    # Compute median values
    grep -v "#" $file | grep -v "chrX\|nan\|#" | awk '{print $4}' | sort -k 1,1n > _tmp
    nlines=$(wc -l _tmp | awk '{print $1}')
    med=$(awk -v nl=${nlines} 'BEGIN{s=0;n1=int(nl/2)+1; if(nl%2==1){n2=n1}else{n2=n1-1}}{if(NR==n1){s+=$1}; if(NR==n2){s+=$1}}END{print "Median TR of ",nl,s/2}' _tmp)
    echo $file $med | sed -e "s/_perChrom//g" -e "s/\.tab//g" -e "s/_TRObs//g" -e "s/_TRICE//g" -e "s/full//g" -e "s/_10000bp//g" >> median_TR_perChrom.txt
    rm _tmp
done

for compartment in A B ;
do
    rm -fvr median_TR${compartment}_perChrom.txt
    for file in $(ls -1 *ICE_compartment_perChrom.tab);
    do
	echo $file

	# Compute median values
	grep -v "chrX\|nan\|#" $file | grep -w ${compartment} | awk '{print $4}' | sort -k 1,1n > _tmp
	nlines=$(wc -l _tmp | awk '{print $1}')
	med=$(awk -v nl=${nlines} 'BEGIN{s=0;n1=int(nl/2)+1; if(nl%2==1){n2=n1}else{n2=n1-1}}{if(NR==n1){s+=$1}; if(NR==n2){s+=$1}}END{print "Median TR of ",nl,s/2}' _tmp)
	echo $file $med | sed -e "s/_compartment//g" -e "s/_perChrom//g" -e "s/\.tab//g" -e "s/_TRObs//g" -e "s/_TRICE//g" -e "s/full//g" -e "s/_10000bp//g" >> median_TR${compartment}_perChrom.txt
	rm _tmp
    done
done





