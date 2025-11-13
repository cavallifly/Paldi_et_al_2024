# Prepare micro-C datasets
echo "sample chrom avg_cis_ratio avg_trans_ratio" | awk '{for(i=1;i<NF;i++){printf("%s\t",$i)}; printf("%s\n",$NF)}' > AvgTransRatio_perChrom_exp.tab

for tag in ICE ;
do
    cat microC_TR_perChromosome/*_TR${tag}_perChrom.tab | grep -v "chrX\|nan\|#" | sed -e "s/\.tab//g" -e "s/_TRObs//g" -e "s/_TRICE//g" -e "s/_10000bp//g" -e "s/full//g" | awk '{for(i=1;i<NF;i++){printf("%s\t",$i)}; printf("%s\n",$NF)}' | grep -v 1h >> AvgTransRatio_perChrom_exp.tab

    cat microC_TR_perChromosome/*_TR${tag}_compartment_perChrom.tab | grep -v "chrX\|nan\|#" | sed -e "s/\.tab//g" -e "s/_TRObs//g" -e "s/_TRICE//g" -e "s/full//g" -e "s/_10000bp//g" | awk '{printf("%s_%s",$1,$NF); for(i=2;i<NF;i++){printf("\t%s",$i)}; printf("\n")}' | grep -v 1h >> AvgTransRatio_perChrom_exp.tab

done
#awk '{print $1}' AvgTransRatio_perChrom_exp.tab | uniq 

nc=$1
mainDir=05_ABcompartmentalization_increased_stiffness

for dir in $(ls -1 ./${mainDir}/mESC_phi_0.10_rosette | grep -v ":" | grep ${nc}_copies);
do
    if [[ ! -d ./${mainDir}/mESC_phi_0.10_rosette/${dir} ]];
    then
	continue
    fi
    echo $dir

    cp AvgTransRatio_perChrom_exp.tab AvgTransRatio_perChrom.tab
    awk '{print $1}' AvgTransRatio_perChrom.tab | uniq   

    for t in 21600000 ;    
    do
	th=$(echo $t | awk '{print int($1/5400000)}') ; echo $th
	for rc in 150 ;
	do
	    cat ./${mainDir}/mESC_phi_0.10_rosette/${dir}/rc_${rc}nm_res_10000bp_from_*_to_${t}_every_900000_replica_*_TR_perChrom.tab  | grep -v "#" | awk -v d=$dir -v rc=${rc} -v t=${th} '{print d"_"t"h",$(NF-2),$(NF-1),$NF}' | sed -e "s/fA_1.00_//g" -e "s/fB_1.00_//g" -e "s/fA_0.00_//g" -e "s/fB_0.00_//g" -e "s/EBB_0.000_//g" -e "s/EAA_0.080_//g" >> AvgTransRatio_perChrom.tab
	    cat ./${mainDir}/mESC_phi_0.10_rosette/${dir}/rc_${rc}nm_res_10000bp_from_*_to_${t}_every_900000_replica_*_TRA_perChrom.tab | grep -v "#" | awk -v d=$dir -v rc=${rc} -v t=${th} '{print d"_A_"t"h",$(NF-2),$(NF-1),$NF}' | sed -e "s/fA_1.00_//g" -e "s/fB_1.00_//g" -e "s/fA_0.00_//g" -e "s/fB_0.00_//g" -e "s/EBB_0.000_//g" -e "s/EAA_0.080_//g" >> AvgTransRatio_perChrom.tab
	    cat ./${mainDir}/mESC_phi_0.10_rosette/${dir}/rc_${rc}nm_res_10000bp_from_*_to_${t}_every_900000_replica_*_TRB_perChrom.tab | grep -v "#" | awk -v d=$dir -v rc=${rc} -v t=${th} '{print d"_B_"t"h",$(NF-2),$(NF-1),$NF}' | sed -e "s/fA_1.00_//g" -e "s/fB_1.00_//g" -e "s/fA_0.00_//g" -e "s/fB_0.00_//g" -e "s/EBB_0.000_//g" -e "s/EAA_0.080_//g" >> AvgTransRatio_perChrom.tab
	    
	    #tail AvgTransRatio_perChrom.tab
	done # Close cycle over $rc
    done # Close cycle over $t

    awk '{print $1}' AvgTransRatio_perChrom.tab | uniq
    wc -l AvgTransRatio_perChrom.tab

    # Make the plot
    Rscript scripts/08_make_TR_violinPlots.R
    
    mv AvgTransRatio_perChrom.pdf AvgTransRatio_perChrom_${mainDir}_${dir}.pdf
    mkdir -p TR_distributions
    mv AvgTransRatio_perChrom_${mainDir}_${dir}.pdf TR_distributions

done # Close cycle over $di

