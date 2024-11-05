scriptsDir=${PWD}/scripts/

condition=AB_compartmentalization
# Parameters' exploration
AAattrEnergies="0.100 0.090 0.080 0.070 0.060 0.050 0.040 0.030 0.020 0.010 0.000" # Attraction energies
BBattrEnergies="0.100 0.090 0.080 0.070 0.060 0.050 0.040 0.030 0.020 0.010 0.000" # Attraction energies

# DMSO case
#EAA_0.080_EBB_0.000_fA_1.00_fB_1.00_1_copies
#AAattrEnergies="0.080" # Attraction energies
#BBattrEnergies="0.000" # Attraction energies

# TSA case
#EAA_0.070_EBB_0.010_fA_1.00_fB_1.00_1_copies
#AAattrEnergies="0.050" # Attraction energies
#BBattrEnergies="0.010" # Attraction energies

pythonScript=${scriptsDir}/01_runs_template.py

tag=TSA
resolution=5000 

ncopies=1
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

for AAattrEnergy in ${AAattrEnergies} ;
do
    for BBattrEnergy in ${BBattrEnergies} ;
    do
	if (( $(echo "${BBattrEnergy} > ${AAattrEnergy}" | bc -l) ));
	then
	    continue
	fi
	
	wDir=EAA_${AAattrEnergy}_EBB_${BBattrEnergy}_${ncopies}_copies
	
	echo $wDir
	mkdir -p ${wDir}
	cd ${wDir}
	
	if [[ ! -e 01_runs_${condition}.py ]];
	then
	    awk '{print $0}' ${pythonScript} | sed -e "s/BBattrEnergy/${BBattrEnergy}/g" -e "s/AAattrEnergy/${AAattrEnergy}/g" > 01_runs_${condition}.py
	fi
	
	for replica in $(seq 1 1 100);
	do
	    replicaDir=replica_${replica}
	    if [[ -d ${replicaDir} ]];
	    then
		#echo "${replicaDir} exists in ${condition}. Remove/rename it, if you need to re-run it!"
		continue
	    fi
	    echo $replicaDir
	    
	    prevDir=/zssd/scratch/mdistefano/2023_10_22_Project_hyperacetylation/Modelling_setup/03_estimate_time_conversion_mESC_phi_${phi}_rosette_${ncopies}_copies/replica_${replica}/
	    if [[ ! -e ${prevDir}/relaxed_conformation.txt ]];
	    then
		ls -lrtha ${prevDir}
		continue
	    fi		
	    
	    mkdir -p ${replicaDir}		
	    cd ${replicaDir}
	    
	    rsync -avz ${prevDir}/relaxed_conformation.txt initial_conformation.txt &> /dev/null
	    
	    sed -e "s/XXXreplicaXXX/${replica}/g" -e "s/XXXradiusXXX/${radius}/g" -e "s/XXXlkXXX/${lk}/g" -e "s/XXXbXXX/${b}/g" -e "s/XXXlkXXX/${lk}/g" -e "s/XXXphiXXX/${phi}/g" -e "s/XXXncopiesXXX/${ncopies}/g"	../01_runs_${condition}.py > 01_runs_${condition}_replica_${replica}.py
	    sed -e "s/XXXreplicaXXX/${replica}/g" -e "s/XXXconditionXXX/${condition}/g" -e "s,XXXdirXXX,${PWD},g" ${scriptsDir}/01_runs_template.cmd > jobscript_${condition}_replica_${replica}.cmd 
	    
	    cd .. # Exit $replicaDir
	done # Close cycle over $replica
	cd .. # Exit ${wDir}
    done    
done
cd .. # Exit ${dir}
