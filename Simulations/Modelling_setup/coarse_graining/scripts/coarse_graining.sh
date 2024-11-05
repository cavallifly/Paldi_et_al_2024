# Ou et al 2017 between 2 and 5 nucleosomes per clutch ~ 1 kb per clutch

#Compute the parameters of the coarse-graining (CG) model applying the ‘yeast’ fine-scale model with nu=1kbp, sigma=20nm, and lk=100nm
#and imposing volumic densities of Φ=0.10 and DNA density of the mESC cells, ρ=0.006 bp/nm3.

#nuCGtarget=10000
nuCGtarget=5000

# General parameters
pi=3.14159265359
c=19 # Numerical constant

# DNA content to simulate
chr=TSA

ncopies=20
DNAcontent=20000000
#ncopies=1
#DNAcontent=5451075338 # diploid autosomes + a single copy of X and Y
#DNAcontent=$(awk -v nu=$nuCGtarget '{if($1=="X" || $1=="Y"){s+=(int($2/nu+1)*nu); next}; s+=(int($2/nu+1)*nu)*2;}END{print s}' <(cat /work/user/mdistefano/mishaDB/trackdb/mm10/chrom_sizes.txt | grep -v "M"))
echo $chr $DNAcontent

# FS model    
#"yeast" fiber
nuFS=1000
bFS=20
lkFS=100

# Data for mESC to derive the volumetric density of the CG model
assembly=mm10
#GenomeLength=$(awk -v nu=$nuCGtarget '{s+=$2}END{print (int(s/nu)*nu+1)*2}' <( grep -v M /media/data/mishaDB/trackdb/${assembly}/chrom_sizes.txt))

### mESC ###
cellType=mESC
rho=0.006

for phiCG in 0.03 0.10 0.20 0.30;
do

    outfile=coarse_graining_FS_${bFS}nm_CG_nu_${nuCGtarget}bp_Mouse_${cellType}_phi_${phiCG}_${chr}_${ncopies}_copies.txt
    
    echo "Defining the quantities for the fine-scale system:"
    NFS=$(awk -v nc=$ncopies -v nu=$nuFS -v DNA=${DNAcontent} 'BEGIN{print int(DNA/nu)*nc}')
    NkFS=$(awk -v N=$NFS -v b=$bFS -v lk=$lkFS 'BEGIN{print N*b/lk}')
    rhoFS=$rho
    rhokFS=$(awk -v rho=$rhoFS -v nu=$nuFS -v b=$bFS -v lk=$lkFS 'BEGIN{print rho/nu*b/lk}')
    LFS=$(awk -v b=$bFS -v N=$NFS 'BEGIN{print N*b}')
    LeFS=$(awk -v lk=$lkFS -v c=$c -v rhok=$rhokFS 'BEGIN{print lk*(c/rhok/lk^3)^2}')
    lkFSbFS_nuFS=$(awk -v lk=$lkFS -v b=$bFS -v nu=$nuFS 'BEGIN{print lk*b/nu}')
    lkFSnuFS_bFS=$(awk -v lk=$lkFS -v b=$bFS -v nu=$nuFS 'BEGIN{print lk*nu/b}')
    echo "Fine scale model is the ${bFS}-nm fiber"
    echo "DNA content to simulate $((${DNAcontent}*${ncopies}))bp"
    
    echo
    
    echo "nuFS (DNA content of a monomer in b) = ${nuFS} bp"
    echo "bFS (Diameter of a bead in nm) = ${bFS} nm"
    echo "lkFS (Kuhn length of the chain in nm) = ${lkFS} nm"
    echo "NFS (Number of monomers to represent the chromosome) = $NFS"
    echo "NkFS (Number if Kuhn lengths of the chain) = ${NkFS}"
    echo "rhoFS (Genome density in bp/nm³) = ${rhoFS}"
    echo "rhokFS (Genome density in Kuhn lengths bp/nm³) = ${rhokFS}"   
    echo "LFS (Polymer contour length) = $LFS nm"
    echo "LeFS (Entanglement length of the chain in nm) = ${LeFS} nm"
    echo "Number of monomers in a Kuhn length FS = "$(awk -v lk=$lkFS -v b=$bFS 'BEGIN{print lk/b}')
    echo "lkFSnuFS_bFS (DNA content of a Kuhn length FS) = ${lkFSnuFS_bFS} bp"
    echo
    
    
    for phi in ${phiCG} ;
    do
	nuI=100 #$nuCGtarget
	nuF=20000000
	flag=0
	
	for dnu in 20000 10000 1000 100 10 1 0.5 ;
	do
	    echo $dnu
	    for DNAperlk in $(seq ${nuI} ${dnu} ${nuF})
	    do
		if [[ $flag -eq 1 ]];
		then
		    continue
		fi
	   	
		# CG model
		lkCGnuCG_bCG=${DNAperlk}
		lkCG=$(awk -v lkCGnuCG_bCG=${lkCGnuCG_bCG} -v lkFSbFS_nuFS=${lkFSbFS_nuFS} 'BEGIN{print sqrt(lkCGnuCG_bCG*lkFSbFS_nuFS)}')
		#bCG=$(awk -v phi=${phiCG} -v rho=${rhoFS} -v lkCGnuCG_bCG=${lkCGnuCG_bCG} -v lkFSbFS_nuFS=${lkFSbFS_nuFS} -v pi=$pi 'BEGIN{print sqrt(sqrt(lkCGnuCG_bCG/lkFSbFS_nuFS)/rho*6/pi*phi)}')
		bCG=$(awk -v phi=${phiCG} -v rho=${rhoFS} -v lkCGnuCG_bCG=${lkCGnuCG_bCG} -v lkFSbFS_nuFS=${lkFSbFS_nuFS} -v pi=$pi 'BEGIN{print sqrt(sqrt(lkCGnuCG_bCG/lkFSbFS_nuFS)/rho*6/pi*phi)}')
		echo ${bCG} ${lkCGnuCG_bCG} ${lkFSbFS_nuFS}
		nuCG=$(awk -v b=${bCG} -v lkCGnuCG_bCG=${lkCGnuCG_bCG} -v lkFSbFS_nuFS=${lkFSbFS_nuFS} 'BEGIN{print int(sqrt(lkCGnuCG_bCG/lkFSbFS_nuFS)*b)}')
		NCG=$(awk -v nc=$ncopies -v nu=$nuCG -v DNA=$DNAcontent 'BEGIN{print int(DNA/nu)*nc}')
		radiusCG=$(awk -v N=$NCG -v b=$bCG -v nu=$nuCG -v rho=$rhoFS -v pi=$pi -v phi=${phiCG} 'BEGIN{print (N/phi)^(1/3)/2}')
		LCG=$(awk -v b=$bCG -v N=$NCG 'BEGIN{print int(N*b)}')
		NkCG=$(awk -v N=$NCG -v b=$bCG -v lk=$lkCG 'BEGIN{print N*b/lk}')
		rhoCG=$rhoFS
		rhokCG=$(awk -v rho=$rhoCG -v nu=$nuCG -v b=$bCG -v lk=$lkCG 'BEGIN{print rho/nu*b/lk}')
		LeCG=$(awk -v lk=$lkCG -v c=$c -v rhok=$rhokCG 'BEGIN{print lk*(c/rhok/lk^3)^2}')  
		
		echo $nuCGtarget $phi $DNAperlk ${nuCG} ${nuI} ${nuF} ${dnu}
		
		if [[ $nuCG -eq ${nuCGtarget} ]];
		then
		    
		    (
			echo "Fine scale model is the ${bFS}-nm fiber"
			echo "DNA content to simulate ${DNAcontent}bp"
			
			echo
			
			echo "nuFS (DNA content of a monomer in b) = ${nuFS} bp"
			echo "bFS (Diameter of a bead in nm) = ${bFS} nm"
			echo "lkFS (Kuhn length of the chain in nm) = ${lkFS} nm"
			echo "NFS (Number of monomers to represent the chromosome) = $NFS"
			echo "NkFS (Number if Kuhn lengths of the chain) = ${NkFS}"
			echo "rhoFS (Genome density in bp/nm³) = ${rhoFS}"
			echo "rhokFS (Genome density in Kuhn lengths bp/nm³) = ${rhokFS}"   
			echo "LFS (Polymer contour length) = $LFS nm"
			echo "LeFS (Entanglement length of the chain in nm) = ${LeFS} nm"
			echo "Number of monomers in a Kuhn length FS = "$(awk -v lk=$lkFS -v b=$bFS 'BEGIN{print lk/b}')
			echo "lkFSnuFS_bFS (DNA content of a Kuhn length FS) = ${lkFSnuFS_bFS} bp"
			echo
			
			echo "Coarse-grained model at ${nuCGtarget} bp"
			echo "Volume of the ${celltype} nucleus ${Vnucleus}"
			echo "Genome length in ${celltype} ${GenomeLength}"
			echo "phiCG (Volumetric density of the chain in the CG model ${celltype} parameters) = ${phiCG}"
			echo
			
			echo "nuCG (DNA content of a monomer in b) = ${nuCG} bp"
			echo "bCG (Diameter of a bead in nm) = ${bCG} nm"
			echo "lkCG (Kuhn length of the chain in nm) = ${lkCG} nm"
			echo "NCG (Number of monomers to represent the chromosome) = $NCG"
			echo "radiusCG (Radius of the spherical simulation environment for the CG model in monomer diameters) = ${radiusCG}"
			echo "NkCG (Number if Kuhn lengths of the chain) = ${NkCG}"
			echo "rhoCG (Genome density in bp/nm³) = ${rhoCG} nm"
			echo "rhokCG (Genome density in Kuhn lengths bp/nm³) = ${rhokCG} nm"
			echo "LCG (Polymer contour length) = $LCG nm"
			echo "LeCG (Entanglement length of the chain in nm) = ${LeCG} nm"
			echo "Number of monomers in a Kuhn length CG = "$(awk -v lk=$lkCG -v b=$bCG 'BEGIN{print lk/b}')
			echo "Number of monomers in a persistence length CG = "$(awk -v lk=$lkCG -v b=$bCG 'BEGIN{print lk/b/2.}')
			echo "lkCGnuCG_bCG (DNA content of a Kuhn length CG) = ${lkCGnuCG_bCG} bp"
			
			echo
			
			echo "These quantities have to be the same in the FS and CG models:"
			echo "rho (The Genome densities of FS and CG are imposed to be the same) = ${rhoFS}"
			
			echo "Ratio of L and Le FS = "$(awk -v L=$LFS -v Le=$LeFS 'BEGIN{print L/Le}')
			echo "Ratio of L and Le CG = "$(awk -v L=$LCG -v Le=$LeCG 'BEGIN{print L/Le}')
			
			echo "Volume of the simulation box FS "$(awk -v N=$NFS -v nu=$nuFS -v rho=$rhoFS 'BEGIN{print N*nu/rho}')
			echo "Volume of the simulation box CG "$(awk -v N=$NCG -v nu=$nuCG -v rho=$rhoCG 'BEGIN{print N*nu/rho}')
			
			flag=1
			
		    ) > ${outfile}
		fi
		if [[ $nuCG -gt ${nuCGtarget} ]];
		then
		    nuI=$(awk -v D=${DNAperlk} -v d=${dnu} 'BEGIN{print D-d}')
		    nuF=${DNAperlk}
		    break
		fi
	    done
	done
    done
done
