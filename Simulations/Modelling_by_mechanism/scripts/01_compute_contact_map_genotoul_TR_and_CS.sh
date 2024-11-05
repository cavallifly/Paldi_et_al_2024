# Arguments

condition=${1}
cd ${condition}
echo "wDir: $PWD"

modelResolution=5000
mapResolution=10000

res=$(awk -v cres=${modelResolution} -v mres=${mapResolution} 'BEGIN{print int(mres/cres)}')

ncopies=1 # In this way the 2 chains won't be considered as two copies of the same chromosome

phi=${2}

nCopies=2
tmpCopies=$(echo $PWD | grep copies | sed "s/_/ /g" | awk '{print $(NF-1)}')
if [[ ${tmpCopies} != "" ]];
then
    nCopies=${tmpCopies}
fi

CGfile=/zdata/data/mdistefano/2023_10_22_Project_hyperacetylation/Modelling_setup/coarse_graining/coarse_graining_FS_20nm_CG_nu_${modelResolution}bp_Mouse_mESC_phi_${phi}_TSA_${nCopies}_copies.txt
b=$(cat ${CGfile} | grep -w bCG | awk '{print $9}')


dir=$(echo ${PWD} | sed "s,/, ,g" | awk '{print $NF}')

for rc in 150 ; #100 150 200 ; # From tests on loop-extrusion only 50 100 150 ; # in nm
do
    echo "rc (nm) = ${rc}"
    for nreplicas in ${3} ;
    do
	echo "nreplicas = ${nreplicas}"	
	for tmin in ${4} ;
	do
	    for tmax in ${5} ;
	    do
		if [[ ${tmin} -ge ${tmax} ]];
		then
		    echo "t-min ${tmin} is larger than t-max ${tmax}"
		    continue
		fi	    
		for tdelta in ${6} ;		      
		do
		    echo "timestep: From ${tmin} to ${tmax} by ${tdelta}"	
		    
		    tag=${7}
		    
		    for copies in intra ;
		    do
			dirName=$(echo $dir | sed -e "s/\./_/g")
			echo "dirname : $dirName"

			for replica in ${8}; #$(seq 1 1 ${nreplicas})
			#for replica in $(seq 1 1 ${nreplicas})			
			do
			    nparticles=$(echo 4000 $nCopies | awk '{print $1*$2}') # All the particles in the system
			    echo ${res} ${nparticles} ${rc} ${b} ${nCopies}		    
			    if [[ ! -d replica_${replica} ]];
			    then
				echo "replica_${replica} doesn't exist!"
				continue
			    fi
			    
			    outMatrix=rc_${rc}nm_res_${mapResolution}bp_from_${tmin}_to_${tmax}_every_${tdelta}_replica_${replica}.tab
			    stringCopy="Considering only cis-copy contacts per each chain"
			    i=0
			    if [[ ${copies} == "inter" ]];
			    then
				outMatrix=rc_${rc}nm_res_${mapResolution}bp_${dirName}_from_${tmin}_to_${tmax}_every_${tdelta}_replica_${replica}_inter.tab			    
				stringCopy="Considering cis-copy and trans-copy contacts per each chain"
				i=0
			    fi
			    if [[ $tag == "tmp" ]];
			    then
				outMatrix=${outMatrix%.tab}_${tag}.tab
			    fi
			    echo "outMatrix : ${outMatrix}"
			    echo "Copies $i ${stringCopy}"


			    tmpDir=_tmp_rc_${replica}_TR
			    if [[ -d ${tmpDir} ]];
			    then
				echo "${tmpDir} exists!"
				continue
			    fi
			    
			    outFile=${outMatrix%.tab}_TR.tab
			    ls -lrtha ${outFile}
			    if [[ ! -e ${outFile} ]];
			    then
				mkdir ${tmpDir}
							    
				cd ${tmpDir}
				touch ${outMatrix}
				
				rm -fvr DATA_FILES_INPUT.txt			
				echo ${tmin} ${tdelta} ${tmax}
				for t in $(seq ${tmin} ${tdelta} ${tmax});
				do
				    echo $t
				    for r in ${replica} ; #$(seq 1 1 ${nreplicas});
				    do
					ls -1 ${PWD}/../replica_${r}/*_${t}.lammpstrj ${PWD}/../replica_${r}/*_${t}.XYZ 2> /dev/null >> DATA_FILES_INPUT.txt
				    done # Close cycle over ${r}
				done # Close cycle over ${t}		       
				cat DATA_FILES_INPUT.txt
				
				nlines=$(cat DATA_FILES_INPUT.txt | wc -l | awk '{print $1}')
				if [[ $nlines -lt 7 ]];
				then
				    rm -fr ${outMatrix} ${outMatrix%.tab}.png output.log distances.txt DATA_FILES_INPUT.txt	    
				    continue
				fi
				echo "Computing the model matrix using ${nlines} snapshots"
				
				rsync -avz /zdata/data/mdistefano/2023_10_22_Project_hyperacetylation/Modelling_by_mechanism/scripts/_compartments_10000bp_${nCopies}_copies _compartments
				sed -e "s/T/0/g" -e "s/A/1/g" -e "s/B/2/g" _compartments > _tmp ; mv _tmp _compartments
				if [[ ! -e _chromosomes ]];
                                then
                                    rm -fvr _chromosomes
                                    for nc in $(seq 1 1 ${nCopies});
                                    do
					offset=$(echo $nc 2000 | awk '{print int((($1-1)*$2))}')
                                        awk -v nc=$nc 'BEGIN{for(i=0;i<=40;i++){chr[i]=i}; s1=(nc-1)*2000; s2=nc*2000; for(i=s1;i<s2;i++){print i,chr[nc-1]}}' >> _chromosomes
                                    done    
				fi
				
				echo $outMatrix #$outProb
				touch ${outMatrix}
				#touch ${outFile}
				
				# -r:		 Resolution of the map in monomers
				# -p:		 Number of particles in the system
				# -d:		 Contact distance cutoff between particles (nm)
				# -b:		 Monomer diameter (nm)
				# -k:		 Ploidy of the system
				# -i:                Consider (0) or not (1) inter copies contacts
				echo ${res} ${nparticles} ${rc} ${b} ${ncopies} 0
				np=$(echo ${nparticles} | awk '{print $1}')
				../../../scripts/compute_contact_map_TR_and_CS -r ${res} -p ${np} -d ${rc} -b ${b} -k ${ncopies} -i 0
				mv TR_per_bin.txt ${outFile}
				
				mv contacts.tab ${outMatrix}
				rm -fr output.log distances.txt DATA_FILES_INPUT.txt contacts.tab

				#rm ${outMatrix}
				rsync -avz ${outMatrix%.tab}_TR.tab ../
				cd ..
			    fi

			    # Compute TR per chromosome
			    inFile=${outFile}
			    outFile=${outMatrix%.tab}_TR_perChrom.tab
			    if [[ ! -e ${outFile} ]];
			    then
				cd ${tmpDir}
				rsync -avz ../${inFile} .
				
				# TR all
				echo "#chrom cisChromContactsFraction transChromContactsFraction" | awk '{for(i=1;i<NF;i++){printf("%s\t",$i)}; printf("%s\n",$NF)}' > ${outFile}
				for chrom in $(seq 0 1 $((${nCopies}-1)));
				do
				    echo $chrom
				    sed "s/_/ /g" ${inFile} | grep -w $chrom | awk -v c=${chrom} '{if($2==c){cis+=$(NF-1); trans+=$NF; cnt++}}END{if(cnt!=0){print c,cis/cnt,trans/cnt}}' >> ${outFile}
				    
				    echo ${1}
				    cat ${outFile}
				done # Close cycle over ${chrom}
				rsync -avz ${outMatrix%.tab}_TR_perChrom.tab ../
				cd ..
			    fi

			    outFileComp=${outMatrix%.tab}_TRA_perChrom.tab
			    if [[ ! -e ${outFileComp} ]];
			    then
				cd ${tmpDir}
				rsync -avz ../${inFile} .
				
				# TR per compartment
				for compartment in A B ;
				do
				    outFileComp=${outMatrix%.tab}_TR${compartment}_perChrom.tab
				    echo "#chrom cisChromContactsFraction transChromContactsFraction" | awk '{for(i=1;i<NF;i++){printf("%s\t",$i)}; printf("%s\n",$NF)}' > ${outFileComp}
				    for chrom in $(seq 0 1 ${nCopies});
				    do
					echo $chrom
					sed "s/_/ /g" ${inFile} | grep -w $chrom | grep -w ${compartment} | awk -v c=${chrom} '{if($2==c){cis+=$(NF-1); trans+=$NF; cnt++}}END{if(cnt!=0){print c,cis/cnt,trans/cnt}}' >> ${outFileComp}
					
					echo ${1}
					cat ${outFileComp}
				    done # Close cycle over ${chrom}
				done
				rsync -avz ${outMatrix%.tab}_TR?_perChrom.tab ..				
				cd ..
			    fi

			    inFile=${outMatrix}
			    outFile=${outMatrix%.tab}_CSHomoOverAll_perChrom.tab
			    if [[ ! -e ${outFile} ]];
			    then
				cd ${tmpDir}

				echo "#compartment chrom CS" | awk '{for(i=1;i<NF;i++){printf("%s\t",$i)}; printf("%s\n",$NF)}' > ${outFile}
				
				for chrom in $(awk '{print $2}' _chromosomes | uniq)
				do
				    CSPerbin=${inFile%.tab}_${chrom}_obsOverExp_CSHomoOverAll.tab
				    
				    mv CS_per_bin_${chrom}.txt ${CSPerbin}
				    
				    echo "# Compute CS per chromosome ${chrom}"
				    sed "s/_/ /g" ${CSPerbin} | awk -v c=${chrom} '{CS[$2]+=$NF; cnt[$2]++}END{for(i in CS){print i,c,CS[i]/cnt[i]}}' | awk '{for(i=1;i<NF;i++){printf("%s\t",$i)}; printf("%s\n",$NF)}' >> ${outFile}				
				done
				rsync -avz *CS*.tab ..
				cd ..
			    fi
			    
			    #rsync -avz ${outMatrix%.tab}_TR_perChrom.tab ${outMatrix%.tab}_TR?_perChrom.tab ${outMatrix%.tab}_TR.tab ../
			    #cd ../
			    rm -fvr ${tmpDir}
			    if [[ $nCopies -eq 1 ]];
			    then
				rm -fvr rc*TR*
			    fi
			    #exit
			done # Close cycle over ${replica}
		    done # Close cycle over ${copies}			
		done # Close cycle over ${nreplicas}		    		
	    done # Close cycle over ${tmin}
	done # Close cycle over ${tmax}
    done # Close cycle over ${tdelta}		    	
done # Close cycle over ${rc}
