ncpus=40

for ncopies in 1 20 ;
do
    if [[ $ncopies -eq 1 ]];
    then
	nReplicas=30
    fi
    if [[ $ncopies -eq 20 ]];
    then
	nReplicas=5
    fi    
    
    cd coarse_graining
    bash scripts/coarse_graining.sh ${ncopies} &> coarse_graining_${ncopies}_copies.out
    cd ..

    ncpu=0
    outDir=00_generate_initial_conformation_rosette_phi_0.03_${ncopies}_copies    
    #if [[ ! -d ${outDir} ]];
    #then
    #if [[ ! -e ${outDir}.tar.gz ]];
    #then
    mkdir -p ${outDir}
    cd ${outDir}
    echo "Generating initial conformation with ${ncopies} rosettes of 20Mb at 5000bp resolution and at volumic density phi=0.03"
    for replica in $(seq 1 1 ${nReplicas});
    do
	echo "Generating replica ${replica}..."
	bash ../scripts/00_Initial_conformation.sh $ncopies $replica &> 00_Initial_conformation_replica_${replica}.out &
	
	ncpu=$((${ncpu}+1))
	
	if [[ ${ncpu} -eq ${ncpus} ]];
	then
	    wait
	    ncpu=0
	fi
    done
    cd ..
    #fi
    #fi
    wait
    
    ncpu=0
    outDir=01_compression_to_mESC_phi_0.10_rosette_${ncopies}_copies
    #if [[ ! -d ${outDir} ]];
    #then
    #if [[ ! -e ${outDir}.tar.gz ]];
    #then
    mkdir -p ${outDir}
    cd ${outDir}
    echo "Compressing system to the target volumic density phi=0.10"
    for replica in $(seq 1 1 ${nReplicas});
    do
	echo "Compressing replica ${replica}..."
	bash ../scripts/01_compression_to_the_target_density.sh $ncopies $replica &> 01_compression_to_the_target_density_${replica}.out &
	ncpu=$((${ncpu}+1))
	
	if [[ ${ncpu} -eq ${ncpus} ]];
	then
	    wait
	    ncpu=0
	fi
    done
    cd ..
    #fi
    #fi
    wait
    
    ncpu=0
    outDir=02_generate_relaxed_mESC_phi_0.10_rosette_${ncopies}_copies
    #if [[ ! -d ${outDir} ]];
    #then
    #if [[ ! -e ${outDir}.tar.gz ]];
    #then
    mkdir -p ${outDir}
    cd ${outDir}
    echo "Relaxing system after compression"
    for replica in $(seq 1 1 ${nReplicas});
    do
	echo "Relaxing replica ${replica}..."
	bash ../scripts/02_generate_relaxed.sh $ncopies $replica &> 02_generate_relaxed_${replica}.out &
	ncpu=$((${ncpu}+1))
	
	if [[ ${ncpu} -eq ${ncpus} ]];
	then
	    wait
		    ncpu=0
	fi
    done
    cd ..
    #fi
    #fi
    wait

    ncpu=0
    outDir=03_estimate_time_conversion_mESC_phi_0.10_rosette_${ncopies}_copies
    #if [[ ! -d ${outDir} ]];
    #then
#	if [[ ! -e ${outDir}.tar.gz ]];
#	then
    mkdir -p ${outDir}
    cd ${outDir}
    echo "Relaxing system for time conversion"
    for replica in $(seq 1 1 ${nReplicas});
    do
	if [[ ! -e replica_${replica} ]];
	then
	    echo "Relaxing for time conversion replica ${replica}..."
	    bash ../scripts/03_estimate_time_conversion.sh $ncopies $replica &> 03_estimate_time_conversion_${replica}.out &
	    ncpu=$((${ncpu}+1))
	    
	    if [[ ${ncpu} -eq ${ncpus} ]];
	    then
		wait
		ncpu=0
	    fi
	fi
    done
    cd ..
    #	fi
    #    fi
    wait
    
done
