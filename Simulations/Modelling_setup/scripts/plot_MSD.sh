PERL5LIB=""
b=$(bash ./*/0*.sh | head -1 | awk '{print $1}')
t=$(grep "timestep=" scripts/0*py | sed -e "s/timestep=//g" -e "s/,//g" | awk '{print $1}')
echo "Monomer diameter = $b nm" 
echo "Tau_LJ = $t dt"

gpfile=../scripts/plot_MSD.gp
ylabel="<MSD> ({/Symbol m}m^2)" 
xlabel="Time lag ({/Symbol t}_{LJ})"
xmin=$(head -3 replica_1_*/MSD.txt | tail -1 | awk -v t=${t} '{print $1*t}')
xmax=$(tail -1 replica_1_*/MSD.txt | awk -v t=${t} '{print $1*t}')
ymin=$(awk -v b=$b 'BEGIN{print 0.00005}')
ymax=$(awk -v b=$b 'BEGIN{print 10}')
psoutfile=plot_MSD.ps
title=
datatitle="All replicas"
colour=green

rsync -avz ~/TOOLS/*_Gu_* .
rsync -avz ~/TOOLS/MSD_Mach_et_al_mESC_Figure1eTetO-exCTCFNOCre.txt .

rm -fvr output_plot_MSD.log
for dir in MSD_analysis ; #tauLJ_0_012 ;
do
    txtinfile=${dir}/MSD.txt
    #set label 2 "XXXthresholdXXX {/Symbol t}_{LJ} {/Symbol \273} ${second}s"
    #1.0285920697311175 0.013094617074388578

    # From Mach et al
    sec=10    
    MSDthreshold=$(awk '{if($1==10) printf("%.4f",$2)}' MSD_Mach_et_al_mESC_*.txt | head -1 | awk '{print $1}')
    echo "MSD at 10s in experiment $MSDthreshold"    
    threshold=$(awk -v m=${MSDthreshold} '{if(NR==1){p=$1;Mp=$2}; if(Mp<m && m<$2){print (p+$1)*0.5}; p=$1;Mp=$2}' ${txtinfile})
    echo "Mach et al tau_LJ for ${sec} second $threshold"

    sed -e "s/XXXsecondXXX/${sec}/g" -e "s/XXXthresholdXXX/${threshold}/g" -e "s/XXXMSDthresholdXXX/${MSDthreshold}/g" -e "s/XXXcolourXXX/${colour}/g" -e "s,XXXylabelXXX,${ylabel},g" -e "s,XXXxlabelXXX,${xlabel},g" -e "s,XXXxminXXX,${xmin},g" -e "s,XXXxmaxXXX,${xmax},g" -e "s,XXXyminXXX,${ymin},g" -e "s,XXXymaxXXX,${ymax},g" -e "s,XXXpsoutfileXXX,${psoutfile},g" -e "s,XXXtitleXXX,${title},g" -e "s,XXXtxtinfileXXX,${txtinfile},g" -e "s,XXXdatatitleXXX,${datatitle},g" ${gpfile} | gnuplot &>> output_plot_MSD.log
    bash ~/TOOLS/ps2pdf.sh
    mv plot_MSD.png plot_MSD_mapping_on_10s_Mach_et_al.png

    # From Gu et al
    sec=1        
    MSDthreshold=$(awk '{if($1==1.0285920697311175) printf("%.4f",$2*3./2.)}' MSD_Gu_et_al_mESC_noTranscription.txt | head -1 | awk '{print $1}')
    echo "MSD at ${sec}s in experiment $MSDthreshold"    
    threshold=$(awk -v m=${MSDthreshold} '{if(NR==1){p=$1;Mp=$2}; if(Mp<m && m<$2){print (p+$1)*0.5}; p=$1;Mp=$2}' ${txtinfile})
    echo "Gu et al tau_LJ for ${sec} second $threshold"

    sed -e "s/XXXsecondXXX/${sec}/g" -e "s/XXXthresholdXXX/${threshold}/g" -e "s/XXXMSDthresholdXXX/${MSDthreshold}/g" -e "s/XXXcolourXXX/${colour}/g" -e "s,XXXylabelXXX,${ylabel},g" -e "s,XXXxlabelXXX,${xlabel},g" -e "s,XXXxminXXX,${xmin},g" -e "s,XXXxmaxXXX,${xmax},g" -e "s,XXXyminXXX,${ymin},g" -e "s,XXXymaxXXX,${ymax},g" -e "s,XXXpsoutfileXXX,${psoutfile},g" -e "s,XXXtitleXXX,${title},g" -e "s,XXXtxtinfileXXX,${txtinfile},g" -e "s,XXXdatatitleXXX,${datatitle},g" ${gpfile} | gnuplot &>> output_plot_MSD.log
    bash ~/TOOLS/ps2pdf.sh
    mv plot_MSD.png plot_MSD_mapping_on_1s_Gu_et_al.png

    mv plot_MSD* fit.log ${dir}


    
done
