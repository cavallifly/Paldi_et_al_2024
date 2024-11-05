tMax=$( echo 129600000 | awk '{print $1/2}')
echo $tMax

for dir in $(ls -1 | grep EAA_0.080_EBB_0.000_fA_1.00_fB_0.00_lkA_18.00_lkB_3.00_lkArec_1.00_lkBrec_0.25_20_copies_tau_0.012);
do
    echo $dir

    cd $dir

    for replicaDir in $(ls -1 | grep replica_);
    do
	if [[ ! -d ${replicaDir} ]];
	then
	    continue
	fi
	
	cd ${replicaDir}
	pwd
	if [[ -e compartmentalization_${tMax}.XYZ ]];
	then
	    ls -lrtha compartmentalization_${tMax}.XYZ
	    cd ..
	    continue
	fi
	ls -lrtha | tail

	t=$(ls -lrtha compartmentalization_*.XYZ | tail -1 | sed -e "s/_/ /g" -e "s/\./ /g" | awk '{print $(NF-1)}')
	echo "Last timestep $t"
	
	tail -80000 compartmentalization_${t}.XYZ | awk '{print $1,1,1,$2,$3,$4,0,0,0}' > _coords
	
	for file in $(ls -1 replica_?.out);
	do
	    mv ${file} ${file%.out}_part1.out
	done
	
	for file in $(ls -1 log.lammps);
	do
	    mv ${file} log_part1.lammps
	done

	pyFile=$(ls -1 *py)
	cat ${pyFile} | sed -e "s/int(4\*60\*60\*second) #/int(4\*60\*60\*second) - ${t} #/g" -e "s/\"runtime\"   : runtime,/\"runtime\"   : runtime, \"reset_timestep\" : ${t},/g" > _tmp
	mv _tmp ${pyFile}

	head -15 initial_conformation.txt >  _tmp
	cat _coords                       >> _tmp
	tail -159946 initial_conformation.txt >> _tmp
	mv _tmp initial_conformation.txt

	cmdFile=$(ls -1 *cmd*)
	sed -e "s/k 4/k 8/g" -e "s/p 4/p 8/g" ${cmdFile} > _tmp
	mv _tmp ${cmdFile}
	mv -v ${cmdFile} ${cmdFile%_launched}
	ls -lrtha

	cd ..
    done

    cd ..    
done
