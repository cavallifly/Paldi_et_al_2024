# Prepare micro-C datasets
fraction=1.00

nc=$1
conditions="DMSO TSA TSA24hREC"
if [[ $nc -eq 1 ]];
then
    conditions="DMSO"
fi

outFile=CS_perChrom.tab
outFileAll=CS_perChrom_All.tab
rm -fvr ${outFile}
mkdir -p CS_distributions

echo "sample CS" | awk '{print $1,$2}' > CS_perChrom_exp.tab
for condition in DMSO TSA TSA24hREC ;
do
    awk -v c=${condition} '{print c"_A",$2}' ./microC_CS_perChromosome/CS_HomoOverAll_within_domains_rescaled_to_150_${condition}_*A.txt >> CS_perChrom_exp.tab
    awk -v c=${condition} '{print c"_B",$2}' ./microC_CS_perChromosome/CS_HomoOverAll_within_domains_rescaled_to_150_${condition}_*B.txt >> CS_perChrom_exp.tab
    #cat microC_CS_perChromosome/${outFile} | grep -v "chrX\|nan\|#\|mC\|F-DMSO\|G-\|None" | sed -e "s/mC_E14_mESC_//g" -e "s/_NO_chrY_chrM//g" -e "s/bins100_//g" -e "s/NO_chrM_chrY_//g" -e "s/\.tab//g" -e "s/_mm10//g" -e "s/_at_10000bp/_${tag}_10kb/g" | awk '{print $1,$4}' >> CS_perChrom_exp.tab
done
cat CS_perChrom_exp.tab
cat CS_perChrom_exp.tab > ${outFileAll}

for phi in 0.10 ;
do
    
    for mainDir in 06_ABcompartmentalization_increased_stiffness_recovery 04_ABcompartmentalization 05_ABcompartmentalization_increased_stiffness ;
    do
	outFileMain=CS_perChrom_mainDir.tab
	cat CS_perChrom_exp.tab > ${outFileMain}
	
	for dir in $(ls -1 ./${mainDir}/mESC_phi_${phi}_rosette | grep _${nc}_copies | grep -v tau); #"lkA_14.00_lkB_3.00\|lkA_16.00_lkB_3.00"); #grep "lkA_14.00_lkB_3.00\|lkA_1.00_lkB_0.00");
	do
	    if [[ ! -d ./${mainDir}/mESC_phi_0.10_rosette/${dir} ]];
	    then
		continue
	    fi
	    echo $dir
	    
	    cat CS_perChrom_exp.tab > ${outFile}
	    
	    for th in $(seq 1 1 24);
	    do
		t=$(echo $th | awk '{print int($1*5400000)}') #; echo $t
		for rc in 150 ;
		do
		    inFile=rc_150nm_res_10000bp_from_*_to_${t}_every_900000_CSHomoOverAll_within_domains_A.txt
		    cat ./${mainDir}/mESC_phi_0.10_rosette/${dir}/${inFile} 2> /dev/null | grep -v "#" | awk -v d=$dir -v rc=${rc} -v t=$th '{print d"_"rc"nm_"t"h_A",$2}' >> ${outFile}
		    cat ./${mainDir}/mESC_phi_0.10_rosette/${dir}/${inFile} 2> /dev/null | grep -v "#" | awk -v d=$dir -v rc=${rc} -v t=$th '{print d"_"rc"nm_"t"h_A",$2}' >> ${outFileAll}
		    cat ./${mainDir}/mESC_phi_0.10_rosette/${dir}/${inFile} 2> /dev/null | grep -v "#" | awk -v d=$dir -v rc=${rc} -v t=$th '{print d"_"rc"nm_"t"h_A",$2}' >> ${outFileMain}
		    inFile=rc_150nm_res_10000bp_from_*_to_${t}_every_900000_CSHomoOverAll_within_domains_B.txt
		    cat ./${mainDir}/mESC_phi_0.10_rosette/${dir}/${inFile} 2> /dev/null | grep -v "#" | awk -v d=$dir -v rc=${rc} -v t=$th '{print d"_"rc"nm_"t"h_B",$2}' >> ${outFile}
		    cat ./${mainDir}/mESC_phi_0.10_rosette/${dir}/${inFile} 2> /dev/null | grep -v "#" | awk -v d=$dir -v rc=${rc} -v t=$th '{print d"_"rc"nm_"t"h_B",$2}' >> ${outFileAll}
		    cat ./${mainDir}/mESC_phi_0.10_rosette/${dir}/${inFile} 2> /dev/null | grep -v "#" | awk -v d=$dir -v rc=${rc} -v t=$th '{print d"_"rc"nm_"t"h_B",$2}' >> ${outFileMain}
		    #cat ./${mainDir}/mESC_phi_0.10_rosette/${dir}/rc_${rc}nm*_to_${t}_*_CSH*_perChr* | grep -v "#" | awk -v d=$dir -v rc=${rc} -v t=$th '{if($(NF-2)=="A" || $(NF-2)=="B") print d"_"rc"nm_"t"h_"$(NF-2),$NF}' >> ${outFile}
		done # Close cycle over $rc
	    done # Close cycle over $t
	    
	    awk '{print $1}' ${outFile} | sort | uniq | wc -l	    
	    
	    # Make the plots
	    Rscript scripts/07_make_CS_violinPlots.R
	    mv CS_perChrom.pdf CS_distributions/CS_within_domains_${mainDir}_${dir}.pdf
	    cp CS_perChrom.tab CS_distributions/CS_within_domains_${mainDir}_${dir}.tab
	    
	    for condition in ${conditions} ; 
	    do
		Rscript scripts/XX_make_CS_violinPlots_condition.R ${condition}
		mv CS_perChrom.pdf CS_distributions/CS_within_domains_${mainDir}_${dir}_${condition}.pdf
	    done
	    #exit
	done # Close cycle over $dir

	cp ${outFileMain} ${outFile}
	
	# Make the plot
	Rscript scripts/XX_make_CS_violinPlots.R
	mv CS_perChrom.pdf CS_distributions/CS_within_domains_all_${nc}_copies_${mainDir}.pdf
	mv CS_perChrom.tab CS_distributions/CS_within_domains_all_${nc}_copies_${mainDir}.pdf  
	
	rm -fvr CS_perChrom.tab
	exit
    done # Close cycle over $mainDir
done
mv ${outFileAll} ${outFile}
cat ${outFile} | awk '{print $1}' | uniq | sort | uniq

# Make the plot
Rscript scripts/XX_make_CS_violinPlots.R
mv CS_perChrom.pdf CS_distributions/CS_within_domains_all_${nc}_copies.pdf
mv CS_perChrom.tab CS_distributions/CS_within_domains_all_${nc}_copies.pdf    

rm -fvr CS_perChrom.tab CS_perChrom_exp.tab
