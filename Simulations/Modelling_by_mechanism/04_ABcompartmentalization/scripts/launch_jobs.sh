n=0

for job in $(ls -1 ./*/*/rep*/*cmd 2> /dev/null);
do 
    
    if [[ -e $job ]];
    then
	echo $job
	#continue
	mv ${job} ${job}_launched
    	bash ${job}_launched

    	n=$(($n+1))
    	echo "Launched job number $n"
    	echo
	if [[ $n -eq 40 ]];
	then
	    wait
	    n=0
	fi
    fi

done # Close cycle over $dir
