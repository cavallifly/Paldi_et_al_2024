
outDir=balanced_text_matrices
mkdir -p ${outDir} 

for resolution in 10000;
do
   for coolFile in $(ls -1 full*.mcool );
   do
       outFile=${outDir}/${coolFile%.mcool}_at_${resolution}bp.tab
       if [[ -e ${outFile} ]];
       then
	   continue
       fi
       touch ${outFile}
       conda run -n cooler cooler dump --balanced --join ${coolFile}::/resolutions/${resolution} > ${outFile}
       
   done # Close cycle over $coolFile
done # Close cycle over $resolution
