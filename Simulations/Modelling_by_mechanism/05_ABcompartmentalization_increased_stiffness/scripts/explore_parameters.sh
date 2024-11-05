scriptsDir=${PWD}/scripts/

condition=05_ABcompartmentalization_increased_stiffness

# 20_copies
lkFAs="2.00 4.00 6.00 8.00 10.00 12.00 13.00 14.00 15.00"
lkFBs="0.00 1.00 2.00 3.00 4.00"

AAattrEnergies="0.080" # Attraction energies
BBattrEnergies="0.000" # Attraction energies

pythonScript=${scriptsDir}/01_runs_template.py

tag=TSA
resolution=5000 

ncopies=20
FSmodel=20nm
phi=0.10
cellType=mESC
t=$(echo $condition | awk '{print substr($1,1,1)}')
echo ""$t

CGfile=../../Modelling_setup/coarse_graining/coarse_graining_FS_${FSmodel}_CG_nu_${resolution}bp_Mouse_${cellType}_phi_${phi}_${tag}_${ncopies}_copies.txt
 
b=$(cat ${CGfile} | grep -w bCG | awk '{print $9}')
lk=$(cat ${CGfile} | grep -w lkCG | awk '{print $10}')
radius=$(cat ${CGfile} | grep -w radiusCG | awk '{print $16}')

dir=mESC_phi_${phi}_rosette
mkdir -p ${dir}
cd ${dir}

for lkFA in ${lkFAs} ;
do
    for lkFB in ${lkFBs} ;
    do
	if (( $(echo "${lkFA} <= ${lkFB}" |bc -l) ));
        then
	    continue
        fi

	for AAattrEnergy in ${AAattrEnergies} ;
	do
	    for BBattrEnergy in ${BBattrEnergies} ;
	    do
		lkSA=$(awk -v lk=${lk} -v lkFA=${lkFA} 'BEGIN{print lkFA*lk}')
		lkSB=$(awk -v lk=${lk} -v lkFB=${lkFB} 'BEGIN{print lkFB*lk}')		
		
		wDir=EAA_${AAattrEnergy}_EBB_${BBattrEnergy}_lkA_${lkFA}_lkB_${lkFB}_${ncopies}_copies		
		
		echo $wDir
		mkdir -p ${wDir}
		cd ${wDir}
		
		if [[ ! -e 01_runs_${condition}.py ]];
		then
		    awk '{print $0}' ${pythonScript} | sed -e "s/BBattrEnergy/${BBattrEnergy}/g" -e "s/AAattrEnergy/${AAattrEnergy}/g" > 01_runs_${condition}.py
		    #awk '{print $0}' ${pythonScript} > 01_runs_${condition}.py
		fi
		
		for replica in $(seq 1 1 5);
		do
		    replicaDir=replica_${replica}
		    if [[ -d ${replicaDir} ]];
		    then
			echo "${replicaDir} exists in ${condition}. Remove/rename it, if you need to re-run it!"
			continue
		    fi
		    echo $replicaDir
		    
		    prevDir=../../../..//Modelling_setup/03_estimate_time_conversion_mESC_phi_${phi}_rosette_${ncopies}_copies/relaxed_conformations/
		    if [[ ! -e ${prevDir}/relaxed_conformation_replica_${replica}.txt ]];
		    then
			ls -lrtha ${prevDir}
			continue
		    fi		
		    
		    mkdir -p ${replicaDir}		
		    cd ${replicaDir}
		    
		    rsync -avz ../${prevDir}/relaxed_conformation_replica_${replica}.txt initial_conformation.txt &> /dev/null
		    
		    sed -e "s/XXXreplicaXXX/${replica}/g" -e "s/XXXradiusXXX/${radius}/g" -e "s/XXXlkXXX/${lk}/g" -e "s/XXXlkFAXXX/${lkFA}/g" -e "s/XXXlkFBXXX/${lkFB}/g" -e "s/XXXbXXX/${b}/g" -e "s/XXXlkXXX/${lk}/g" -e "s/XXXphiXXX/${phi}/g" -e "s/XXXncopiesXXX/${ncopies}/g"	../01_runs_${condition}.py > 01_runs_${condition}_replica_${replica}.py
		    sed -e "s/XXXreplicaXXX/${replica}/g" -e "s/XXXconditionXXX/${condition}/g" -e "s,XXXdirXXX,${PWD},g" -e "s,XXXncXXX,${ncopies},g" ${scriptsDir}/01_runs_template.cmd > jobscript_${condition}_replica_${replica}.cmd
		    
		    cd .. # Exit $replicaDir
		done # Close cycle over $replica
		cd .. # Exit ${wDir}
	    done    
	done
    done
done
cd .. # Exit ${dir}
	
