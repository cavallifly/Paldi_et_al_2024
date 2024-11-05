scriptsdir=/work/user/mdistefano/2023_10_22_Project_hyperacetylation/Modelling_by_mechanism/scripts/
mapResolution=10000

condition=TSA
#cd $1

for modMatrix in $(ls -1 rc*900000.tab 2> /dev/null) ;
do
    echo $modMatrix
    size=$(tail -1 ${modMatrix} | awk '{print $1}')
    echo "Size $size"
    
    modMatrix0=${modMatrix}
    
    awk -v size=$size '{if($1!=$2 && ($1!=size && $2!=size)) print $1,$2,$3}' ${modMatrix} | awk -v np=${nparticles} '{print $0}' > _tmp_mod
    modMatrix=_tmp_mod
    head ${modMatrix}
    
    for exp in ${condition};
    do
	exp=$(echo $condition)

	name=$(echo $modMatrix0 | sed -e "s/\.tab//g")
	outFile=${name}_${exp}.png	
	if [[ -e ${outFile} ]];
	then
	    ls -lrtha ${outFile}
	    rm -fvr _tmp_mod
	    continue
	fi
	
	#cHiCMatrix="/work/user/mdistefano/2021_04_13_Project_with_Ivana/cHiCMatrices/chic_${exp}_merge_mm10_IJ_modelledRegion_chr18_53700000_56700000_at_5000bp_comparison.tab"
	#awk -v size=$size '{if($1!=$2 && ($1!=size && $2!=size)) print $0}' ${cHiCMatrix} > _tmp_cHiC
	#awk -v r=${mapResolution} -v start=${start} '{if($2<$5){print ($2-start)/r,($5-start)/r,$7}}' ${cHiCMatrix} > _tmp_cHiC
	#expMatrix=_tmp_cHiC
	#expFactor=$(awk '{d=sqrt(($1-$2)*($1-$2)); if(d==10){sum+=$3;cnt++}}END{print sum/cnt}' ${expMatrix})
	expFactor=1
	
	modFactor=$(awk '{d=sqrt(($1-$2)*($1-$2)); if(d==10){sum+=$3;cnt++}}END{print sum/cnt}' ${modMatrix})
	#head ${modMatrix}
	echo $size $modFactor
	#maxRank=$(wc -l $expMatrix | awk '{print $1}')
	#wc -l $expMatrix $modMatrix

	diagOff=5
	#paste ${modMatrix} ${expMatrix} | head
	#paste ${modMatrix} ${expMatrix} | tail

	# Model vs Hi-C
	#paste ${modMatrix} ${expMatrix} | awk -v fM=${modFactor} -v fE=${expFactor} -v diagOFF=${diagOFF} '{if($1!=$4 || $2!=$5){next}; v=$3/fM; if(sqrt(($1-$2)*($1-$2))<diagOFF){v=0/fM}; if($1<$2){print $1,$2,v}}' >   contact_map.tab
	#paste ${modMatrix} ${expMatrix} | awk -v fM=${modFactor} -v fE=${expFactor} -v diagOFF=${diagOFF} '{if($1!=$4 || $2!=$5){next}; v=$6/fE; if(sqrt(($1-$2)*($1-$2))<diagOFF || $6=="NaN"){v=0/fE}; if($1<$2){print $2,$1,v}}' >>  contact_map.tab

	# Model
	#paste ${modMatrix} ${expMatrix} | awk -v fM=${modFactor} -v fE=${expFactor} -v diagOFF=${diagOFF} '{if($1==$4 && $2==$5){v=$3/fM; if($7=="MASKED"){v=0/fM}; if(sqrt(($1-$2)*($1-$2))<diagOFF){v=0/fM}; if($1<=$2) print $2,$1,v}}' >  contact_map.tab
	#paste ${modMatrix} ${expMatrix} | awk -v fM=${modFactor} -v fE=${expFactor} -v diagOFF=${diagOFF} '{if($1==$4 && $2==$5){v=$3/fM; if($7=="MASKED"){v=0/fM}; if(sqrt(($1-$2)*($1-$2))<diagOFF){v=0/fM}; if($1<=$2) print $1,$2,v}}' >>  contact_map.tab
	cat ${modMatrix} | awk -v fM=${modFactor} -v diagOFF=${diagOFF} '{v=$3/fM; if(sqrt(($1-$2)*($1-$2))<diagOFF){v=0/fM}; if($1<=$2) print $2,$1,v}' >  contact_map.tab
	cat ${modMatrix} | awk -v fM=${modFactor} -v diagOFF=${diagOFF} '{v=$3/fM; if(sqrt(($1-$2)*($1-$2))<diagOFF){v=0/fM}; if($1<=$2) print $1,$2,v}' >>  contact_map.tab	
	head contact_map.tab

	# Hi-C
	#paste ${modMatrix} ${expMatrix} | awk -v fM=${modFactor} -v fE=${expFactor} -v diagOFF=${diagOFF} '{if($1==$4 && $2==$5){v=$6/fE; if($7=="MASKED"){v=0/fE}; if(sqrt(($1-$2)*($1-$2))<diagOFF){v=0/fE}; if($1< $2) print $2,$1,v}}' >  contact_map.tab
	#paste ${modMatrix} ${expMatrix} | awk -v fM=${modFactor} -v fE=${expFactor} -v diagOFF=${diagOFF} '{if($1==$4 && $2==$5){v=$6/fE; if($7=="MASKED"){v=0/fE}; if(sqrt(($1-$2)*($1-$2))<diagOFF){v=0/fE}; if($1< $2) print $1,$2,v}}' >> contact_map.tab 

	### Compute ranks
	# Model vs Hi-C
	#Observed -> general ranking -> check which ranks have the same expected and assign to all of those the same ranking!
	#paste ${modMatrix} ${expMatrix} | sort -k 3,3n | awk '{print $1,$2,$3,NR,$4,$5,$6}' | awk '{if(r[$3]==""){r[$3]=$4}; print $1,$2,r[$3],$5,$6,$7}' | awk -v m=$maxRank '{if($6=="NaN" || $6==0){print $1,$2,-100,$4,$5,$6}else{print $1,$2,($3-1)/(m-1)*200-100,$4,$5,$6}}' | sort -k 1,1n -k 2,2n | awk -v fM=1 -v fE=1 -v diagOFF=${diagOFF} '{if($1!=$4 || $2!=$5){next}; v=$3/fM; if(sqrt(($1-$2)*($1-$2))<diagOFF || $6=="NaN"){v=0/fM}; if($1<$2){print $1,$2,v}}' >   contact_map.tab
	#paste ${modMatrix} ${expMatrix} | sort -k 3,3n | awk '{print $1,$2,$3,NR,$4,$5,$6}' | awk '{if(r[$3]==""){r[$3]=$4}; print $1,$2,r[$3],$5,$6,$7}' | awk -v m=$maxRank '{if($6=="NaN" || $6==0){print $1,$2,-100,$4,$5,$6}else{print $1,$2,($3-1)/(m-1)*200-100,$4,$5,$6}}' | sort -k 1,1n -k 2,2n | awk -v fM=1 -v fE=1 -v diagOFF=${diagOFF} '{if($1!=$4 || $2!=$5){next}; v=$3/fM; if(sqrt(($1-$2)*($1-$2))<diagOFF || $6=="NaN"){v=0/fM}; if($1<$2){print $1,$2,0}}' >   contact_map.tab
	#paste ${modMatrix} ${expMatrix} | sort -k 6,6n | awk '{print $1,$2,$3,$4,$5,$6,NR}' | awk '{if(r[$6]==""){r[$6]=$7}; print $1,$2,$3,$4,$5,r[$6]}' | awk -v m=$maxRank '{if($6=="NaN" || $6==0){print $1,$2,-100,$4,$5,$6}else{print $1,$2,$3,$4,$5,($6-1)/(m-1)*200-100}}' | sort -k 1,1n -k 2,2n | awk -v fM=1 -v fE=1 -v diagOFF=${diagOFF} '{if($1!=$4 || $2!=$5){next}; v=$6/fM; if(sqrt(($1-$2)*($1-$2))<diagOFF || $6=="NaN"){v=0/fM}; if($1<$2){print $2,$1,v}}' >>   contact_map.tab	

	# Model
	#paste ${modMatrix} ${expMatrix} | sort -k 3,3n | awk -v m=$maxRank '{if($6=="NaN"){print $1,$2,0,$4,$5,$6}else{print $1,$2,NR/m*200-100,$4,$5,$6}}' | sort -k 1,1n -k 2,2n | awk -v fM=1 -v fE=1 -v diagOFF=${diagOFF} '{if($1!=$4 || $2!=$5){next}; v=$3/fM; if(sqrt(($1-$2)*($1-$2))<diagOFF || $6=="NaN"){v=0/fM}; if($1<$2){print $1,$2,v}}' >   contact_map.tab
	#paste ${modMatrix} ${expMatrix} | sort -k 3,3n | awk -v m=$maxRank '{if($6=="NaN"){print $1,$2,0,$4,$5,$6}else{print $1,$2,NR/m*200-100,$4,$5,$6}}' | sort -k 1,1n -k 2,2n | awk -v fM=1 -v fE=1 -v diagOFF=${diagOFF} '{if($1!=$4 || $2!=$5){next}; v=$3/fM; if(sqrt(($1-$2)*($1-$2))<diagOFF || $6=="NaN"){v=0/fM}; if($1<$2){print $2,$1,v}}' >>  contact_map.tab

	# Hi-C
	#paste ${modMatrix} ${expMatrix} | sort -k 6,6n | awk -v m=$maxRank '{if($6=="NaN"){print $1,$2,$3,$4,$5,0}else{print $1,$2,$3,$4,$5,NR/m*200-100}}' | sort -k 1,1n -k 2,2n | awk -v fM=1 -v fE=1 -v diagOFF=${diagOFF} '{if($1!=$4 || $2!=$5){next}; v=$6/fE; if(sqrt(($1-$2)*($1-$2))<diagOFF || $6=="NaN"){v=0/fE}; if($1<$2){print $1,$2,v}}' >  contact_map.tab	
	#paste ${modMatrix} ${expMatrix} | sort -k 6,6n | awk -v m=$maxRank '{if($6=="NaN"){print $1,$2,$3,$4,$5,0}else{print $1,$2,$3,$4,$5,NR/m*200-100}}' | sort -k 1,1n -k 2,2n | awk -v fM=1 -v fE=1 -v diagOFF=${diagOFF} '{if($1!=$4 || $2!=$5){next}; v=$6/fE; if(sqrt(($1-$2)*($1-$2))<diagOFF || $6=="NaN"){v=0/fE}; if($1<$2){print $2,$1,v}}' >>  contact_map.tab	

	### Compute ranks differences
	#paste ${modMatrix} ${expMatrix} | sort -k 3,3n | awk -v m=$maxRank '{if($6=="NaN"){print $1,$2,0,$4,$5,$6}else{print $1,$2,NR/m*200-100,$4,$5,$6}}' | sort -k 1,1n -k 2,2n | awk -v fM=1 -v fE=1 -v diagOFF=${diagOFF} '{if($1!=$4 || $2!=$5){next}; vM=$3/fM; vE=$6/fE; if(sqrt(($1-$2)*($1-$2))<diagOFF){vM=0/fM}; if(sqrt(($1-$2)*($1-$2))<diagOFF){vE=0/fE}; if($1<$2){print $1,$2,vM-VE}}' >   contact_map_diff.tab
	#paste ${modMatrix} ${expMatrix} | sort -k 3,3n | awk -v m=$maxRank '{if($6=="NaN"){print $1,$2,0}}'
	#paste ${modMatrix} ${expMatrix} | sort -k 3,3n | awk -v m=$maxRank '{if($6=="NaN"){print $4,$5,0}}'
	#paste ${modMatrix} ${expMatrix} | sort -k 3,3n | awk -v m=$maxRank '{if($6=="NaN"){print $1,$2,0}else{print $1,$2,NR/m*200-100}}' | sort -k 1,1n -k 2,2n > _mod_rank
	#paste ${modMatrix} ${expMatrix} | sort -k 6,6n | awk -v m=$maxRank '{if($6=="NaN"){print $4,$5,0}else{print $4,$5,NR/m*200-100}}' | sort -k 1,1n -k 2,2n > _exp_rank
	#paste _mod_rank _exp_rank | awk -v fM=1 -v fE=1 -v diagOFF=${diagOFF} '{if($1!=$4 || $2!=$5){next}; vE=$6/fE; vM=$3/fM; if(sqrt(($1-$2)*($1-$2))<diagOFF){vM=0/fM; vE=0/fE}; if($1<$2){print $1,$2,vM-vE}}' >   contact_map_diff.tab
	#paste _mod_rank _exp_rank | awk -v fM=1 -v fE=1 -v diagOFF=${diagOFF} '{if($1!=$4 || $2!=$5){next}; vE=$6/fE; vM=$3/fM; if(sqrt(($1-$2)*($1-$2))<diagOFF){vM=0/fM; vE=0/fE}; if($1<$2){print $2,$1,vM-vE}}' >>   contact_map_diff.tab
	#wc -l _mod_rank _exp_rank
	#awk -v size=${size} 'BEGIN{for(i=0;i<(size-80);i++){print i,i,0}}' >> contact_map_diff.tab
	#sort -k 1,1n -k 2,2n contact_map_diff.tab > _a.tab ; mv _a.tab contact_map_diff.tab


	#if [[ ! -e contacts_vs_L_${name}_${exp}.txt ]];
	#then
	#    paste ${modMatrix} ${expMatrix} | sort -k 3,3n | awk '{if($6!="NaN"){print $1,$2,$3}}' | grep -v NaN | awk '{d=sqrt(($1-$2)*($1-$2)); h[d]+=$3; cnt[d]++}END{f=h[10]/cnt[10]; for(i in h){print i,h[i]/cnt[i]/f}}' | sort -k 1n > contacts_vs_L_${name}_${exp}.txt
	#fi
	#if [[ ! -e contacts_vs_L_${name}_${exp}.png ]];
	#then
	#    #cp /work/user/mdistefano/2021_04_13_Project_with_Ivana/cHiCMatrices/contacts_vs_L_chic_${exp}_merge_mm10_IJ_modelledRegion_chr18_53700000_56700000_at_5000bp_comparison.tab exp.txt
	#    cp /work/user/mdistefano/2021_04_13_Project_with_Ivana/cHiCMatrices/average_contacts_chic_${exp}_merge_mm10_IJ_modelledRegion_chr18_53700000_56700000_at_5000bp_comparison.tab exp.txt
	#    cp contacts_vs_L_${name}_${exp}.txt data.txt
	#    gnuplot ${scriptsdir}/03_plot_contacts_vs_L_genotoul.gp
	#    mv data.ps contacts_vs_L_${name}_${exp}.ps
	#    head data.txt exp.txt
	#    rm -fvr data.txt exp.txt
	#fi

	awk -v size=${size} 'BEGIN{for(i=0;i<size;i++){print i,i,0}}' >> contact_map.tab
	sort -k 1,1n -k 2,2n contact_map.tab > _a.tab ; mv _a.tab contact_map.tab
	#sort -k 1,1n -k 2,2n contact_map.tab | awk '{if((50<=$1 && $1<=550) && (50<=$2 && $2<=550)) print $0}' > _a.tab ; mv _a.tab contact_map.tab
	#head contact_map.tab
	#tail contact_map.tab
	#wc -l contact_map.tab
	#cat contact_map.tab
	
	cbmin=$(awk '{if(NR==1){min=$3}; if($3<min) min=$3}END{print min}' contact_map.tab)
	cbmax=$(awk '{if($3>max) max=$3}END{print max}' contact_map.tab)
	echo $cbmin $cbmax
	for maxFactor in 0.1 #1 2 3 4 5 0.001 0.01 0.1 1 ;
	do
	    for minFactor in 100 #1 2 3 4 5 0.001 0.01 0.1 1 10 100 ;
	    do
		for scale in 0.75 ; #0.80 0.90 0.95 ;
		do
		    for factor in 1 ;
		    do
			if [[ ! -e ${outFile} ]];
			then
			    wc -l contact_map.tab
			    sed -e "s/XXXmaxFactorXXX/${maxFactor}/g" -e "s/XXXminFactorXXX/${minFactor}/g" -e "s/XXXcbmaxXXX/${cbmax}/g" -e "s/XXXcbminXXX/${cbmin}/g" -e "s/XXXsizeXXX/$((${size}-80))/g" -e "s/XXXfactorXXX/${factor}/g" -e "s/XXXscaleXXX/${scale}/g" -e "s/XXXconditionXXX/${exp}/g" ${scriptsdir}/03_plot_contact_matrix.gp | gnuplot #2> /dev/null
			    mv contact_map.ps ${outFile%.png}.ps
			    bash ${scriptsdir}/ps2pdf.sh #&> /dev/null

			    #wc -l contact_map_diff.tab
			    #mv contact_map_diff.tab contact_map.tab
			    #sed -e "s/XXXmaxFactorXXX/${maxFactor}/g" -e "s/XXXminFactorXXX/${minFactor}/g" -e "s/XXXcbmaxXXX/${cbmax}/g" -e "s/XXXcbminXXX/${cbmin}/g" -e "s/XXXsizeXXX/$((${size}-80))/g" -e "s/XXXfactorXXX/${factor}/g" -e "s/XXXscaleXXX/${scale}/g" -e "s/XXXconditionXXX/${exp}/g" ${scriptsdir}/03_plot_contact_matrix_genotoul.gp | gnuplot #2> /dev/null
			    #mv contact_map.ps ${outFile%.png}_diff.ps
			    #bash ${scriptsdir}/ps2pdf.sh #&> /dev/null
			fi
			ls -lrtha ${outFile}
		    done
		done
	    done
	done
	
    done # Close cycle over $modMatrix
done # Close cycle over $exp
rm -fr contact_map.tab contact_map_diff.tab _tmp_mod _tmp_cHiC _exp_rank _mod_rank
	
cd .. # Exit ${condition}
