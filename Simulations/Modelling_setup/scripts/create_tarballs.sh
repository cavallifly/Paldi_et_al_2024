for dir in $(ls -1 | grep 03 | grep _12_copies | grep -v tar); 
do 
    echo $dir 
    if [[ ! -e ${dir}.tar.gz ]];
    then
        tar -czvf ${dir}.tar.gz ${dir} 
    fi
    rm -fvr $dir
done
