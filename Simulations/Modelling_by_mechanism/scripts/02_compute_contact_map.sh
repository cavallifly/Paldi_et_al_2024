# Arguments

condition=${1}
cd ${condition}
echo "wDir: $PWD"

modelResolution=5000
mapResolution=10000

res=$(awk -v cres=${modelResolution} -v mres=${mapResolution} 'BEGIN{print int(mres/cres)}')

ncopies=1 # In this way the 2 chains won't be considered as two copies of the same chromosome

phi=${2}

tmpCopies=$(echo $PWD | grep copies | sed "s/_/ /g" | awk '{print $(NF-1)}')
nCopies=1
if [[ ${tmpCopies} != "" ]];
then
    nCopies=${tmpCopies}
fi

CGfile=/zssd/scratch/mdistefano/2023_10_22_Project_hyperacetylation/Modelling_setup/coarse_graining/coarse_graining_FS_20nm_CG_nu_${modelResolution}bp_Mouse_mESC_phi_${phi}_TSA_${nCopies}_copies.txt
b=$(cat ${CGfile} | grep -w bCG | awk '{print $9}')

dir=$(echo ${PWD} | sed "s,/, ,g" | awk '{print $NF}')

for rc in 150 ;# 150 200 ; # From tests on loop-extrusion only 50 100 150 ; # in nm
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
		    continue
		fi	    
		for tdelta in ${6} ;		      
		do
		    echo "timestep: From ${tmin} to ${tmax} by ${tdelta}"	

		    if [[ -e rc_${rc}nm_res_${mapResolution}bp_from_${tmin}_to_${tmax}_every_${tdelta}_TSA.png ]];
		    then
			continue
		    fi

		    #if [[ -e rc_${rc}nm_res_${mapResolution}bp_from_${tmin}_to_${tmax}_every_${tdelta}.tab ]];
		    #then
		    #continue
		    #fi
		    
		    tag=${7}
		    
		    for copies in intra ;
		    do
			dirName=$(echo $dir | sed -e "s/\./_/g")
			echo "dirname : $dirName"
			
			confDir=${PWD}
			
			nparticles=$(echo ${modelResolution} ${mapResolution} ${nCopies} | awk '{print int(4000*$3)}') # All the particles in the system
			echo ${res} ${nparticles} ${rc} ${b} ${nCopies}
			
			outMatrix=rc_${rc}nm_res_${mapResolution}bp_from_${tmin}_to_${tmax}_every_${tdelta}.tab
			stringCopy="Considering only cis-copy contacts per each chain"
			i=0
			if [[ ${copies} == "inter" ]];
			then
			    outMatrix=rc_${rc}nm_res_${mapResolution}bp_${dirName}_from_${tmin}_to_${tmax}_every_${tdelta}_inter.tab			    
			    stringCopy="Considering cis-copy and trans-copy contacts per each chain"
			    i=0
			fi
			if [[ $tag == "tmp" ]];
			then
			    outMatrix=${outMatrix%.tab}_${tag}.tab
			fi
			echo "outMatrix : ${outMatrix}"
			echo "Copies $i ${stringCopy}"
			
			#touch ${outMatrix}
			
			rm -fvr DATA_FILES_INPUT.txt			
			echo ${tmin} ${tdelta} ${tmax}
			for r in $(seq 1 1 ${nreplicas})
			do
			    for t in $(seq ${tmin} ${tdelta} ${tmax});
			    do
				echo $t
				ls -1 ${confDir}/replica_${r}/*_${t}.lammpstrj ${confDir}/replica_${r}/*_${t}.XYZ 2> /dev/null >> DATA_FILES_INPUT.txt
			    done # Close cycle over ${t}
			done
			cat DATA_FILES_INPUT.txt
			
			nlines=$(cat DATA_FILES_INPUT.txt | wc -l | awk '{print $1}')
			if [[ $nlines -eq 0 ]];
			then
			    rm -fr ${outMatrix} ${outMatrix%.tab}.png output.log distances.txt DATA_FILES_INPUT.txt	    
			    continue
			fi
			echo "Computing the model matrix using ${nlines} snapshots"
			
			echo $outMatrix #$outProb
			if [[ ! -e ${outMatrix} ]];
			then
			    touch ${outMatrix}
			    
			    # -r:		 Resolution of the map in monomers
			    # -p:		 Number of particles in the system
			    # -d:		 Contact distance cutoff between particles (nm)
			    # -b:		 Monomer diameter (nm)
			    # -k:		 Ploidy of the system
			    # -i:                Consider (0) or not (1) inter copies contacts
			    echo ${res} ${nparticles} ${rc} ${b} ${nCopies} 0
			    np=$(echo ${nparticles} | awk '{print $1}')
			    ./scripts/compute_contact_map -r ${res} -p ${np} -d ${rc} -b ${b} -k ${nCopies} -i 0

			    mv contacts.tab ${outMatrix}
			    rm -fr output.log distances.txt DATA_FILES_INPUT.txt contacts.tab			    
			fi
			size=$(awk -v np=${nparticles} -v r=${res} 'BEGIN{print int(np/r)}')
			echo "Size of the matrix ${size}"

			bash ../../../scripts/03_plot_contact_analysis.sh

			# Create .mcool file
			head -2001 ${outMatrix} | tail
			#cat ${outMatrix} | awk -v s=$size -v nc=${nCopies} -v res=${mapResolution} 'BEGIN{size=int(s/nc)}{if($3==0){next}; chrom1="chr"int($1/size)+1; start1=$1%size*res; end1=start1+res; chrom2="chr"int($2/size)+1 ; start2=$2%size*res; end2=start2+res; printf("%s\t%d\t%d\t%s\t%d\t%d\t%d\n",chrom1,start1,end1,chrom2,start2,end2,$3)}' > ${outMatrix%tab}ints
			cat ${outMatrix} | awk -v s=$size -v nc=1 -v res=${mapResolution} 'BEGIN{size=int(s/nc)}{if($3==0){next}; chrom1="chr"int($1/size)+1; start1=$1%size*res; end1=start1+res; chrom2="chr"int($2/size)+1 ; start2=$2%size*res; end2=start2+res; printf("%s\t%d\t%d\t%s\t%d\t%d\t%d\n",chrom1,start1,end1,chrom2,start2,end2,$3)}' > ${outMatrix%tab}ints
			
			#rm -fr ${outMatrix}
			
		    done # Close cycle over ${copies}			
		done # Close cycle over ${nreplicas}		    		
	    done # Close cycle over ${tmin}
	done # Close cycle over ${tmax}
    done # Close cycle over ${tdelta}		    	
done # Close cycle over ${rc}
