ncpus=40
nReplicas=5
resolution=5000


for ncopies in 20 ; #1 2 4 6 8 10 12 14 16 18 20 30 40 ;
do
    #cd coarse_graining
    #bash scripts/coarse_graining.sh ${ncopies} ${resolution} &> coarse_graining_${ncopies}_copies_at_resolution_${resolution}.out
    #cd ..

    nparticles=$(grep NCG coarse_graining/coarse_graining_FS_20nm_CG_nu_${resolution}bp_Mouse_mESC_phi_0.10LV_TSA_${ncopies}_copies.txt | awk '{print $NF}')
    echo $nparticles

    ### NOT USED !!! ###
    #ncpu=0
    #outDir=00_generate_initial_conformation_rosette_phi_0.01_resolution_${resolution}bp_${ncopies}_copies    

    #mkdir -p ${outDir}
    #cd ${outDir}
    #echo "Generating initial conformation with ${ncopiess} rosettes each of ${nparticles} at ${resolution} bp resolution and at volumic density phi=0.01"
    #for replica in $(seq 1 1 ${nReplicas});
    #do
    #echo "Generating replica ${replica}..."
    #resolution=$1
    #nparticles=$2
    #ncopies=$3
    #replica=$4
    #bash ../scripts/00_Initial_conformation.sh $resolution $nparticles $ncopies $replica &> 00_Initial_conformation_replica_${replica}.out &
	
    #ncpu=$((${ncpu}+1))
	
    #if [[ ${ncpu} -eq ${ncpus} ]];
    #then
    #wait
    #ncpu=0
    #fi
    #done
    #cd ..
    ### NOT USED !!! ###
    
    ncpu=0
    outDir=01_compression_to_mESC_rosette_phi_0.10LV_resolution_${resolution}bp_${ncopies}_copies
    mkdir -p ${outDir}
    cd ${outDir}
    echo "Compressing system to the target volumic density phi=0.10"
    for replica in $(seq 1 1 ${nReplicas});
    do
	echo "Compressing replica ${replica}..."
	bash ../scripts/01_compression_to_the_target_density_LV.sh $resolution $nparticles $ncopies $replica &> 01_compression_to_the_target_density_${replica}.out &
	ncpu=$((${ncpu}+1))
	
	if [[ ${ncpu} -eq ${ncpus} ]];
	then
	    wait
	    ncpu=0
	fi
    done
    cd ..
    wait

    ncpu=0
    outDir=02_generate_relaxed_mESC_rosette_phi_0.10LV_resolution_${resolution}bp_${ncopies}_copies
    mkdir -p ${outDir}
    cd ${outDir}
    echo "Relaxing system after compression"
    for replica in $(seq 1 1 ${nReplicas});
    do
	echo "Relaxing replica ${replica}..."
	bash ../scripts/02_generate_relaxed_LV.sh $resolution $nparticles $ncopies $replica &> 02_generate_relaxed_${replica}.out &
	ncpu=$((${ncpu}+1))
    
	if [[ ${ncpu} -eq ${ncpus} ]];
	then
	    wait
	    ncpu=0
	fi
    done
    cd ..
    wait

    ncpu=0
    outDir=03_estimate_time_conversion_mESC_rosette_phi_0.10LV_resolution_${resolution}bp_${ncopies}_copies    
    if [[ -e ${outDir}.tar.gz ]];
    then
	exit
    fi
    mkdir -p ${outDir}
    cd ${outDir}
    echo "Relaxing system for time conversion"
    for replica in $(seq 1 1 ${nReplicas});
    do
	if [[ ! -e replica_${replica} ]];
	then
	    echo "Relaxing for time conversion replica ${replica}..."
	    bash ../scripts/03_estimate_time_conversion_LV.sh $resolution $nparticles $ncopies $replica &> 03_estimate_time_conversion_${replica}.out &
	    ncpu=$((${ncpu}+1))
	    
	    if [[ ${ncpu} -eq ${ncpus} ]];
	    then
		wait
		ncpu=0
	    fi
	fi
    done
    cd ..
    wait
    
done
