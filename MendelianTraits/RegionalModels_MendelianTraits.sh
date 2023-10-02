#!/bin/sh
#$ -cwd
#$ -l h_rt=47:30:00
#$ -pe sharedmem 4
#$ -l h_vmem=8G
#$ -o _o_files/
#$ -e _e_files/
. /etc/profile.d/modules.sh

module load igmm/apps/dissect/1.15.2c
module load igmm/apps/plink/1.90b4
module load igmm/apps/gcta/1.91.2beta

i=$SGE_TASK_ID

traitinfo=$(sed "${i}q;d" ../traits.txt)

IFS=','; trait=($traitinfo); unset IFS;

mkdir -p ${trait[0]}

cd ${trait[0]}



##0.1
maf=0.1
output=01
mkdir -p Rare_${output}

##0.1 Chrom
plink --bfile ../../InputData/sheep_geno_imputed_Plates1to97_20220627 --chr ${trait[3]} --max-maf ${maf} --nonfounders --write-snplist --out Rare_${output}/RareSNPs${output}_Chrom_${trait[3]}
gcta64 --bfile ../../InputData/sheep_geno_imputed_Plates1to97_20220627 --extract Rare_${output}/RareSNPs${output}_Chrom_${trait[3]}.snplist --make-grm --out Rare_${output}/RareSNPs${output}_Chrom_${trait[3]} 
gcta64 --grm Rare_${output}/RareSNPs${output}_Chrom_${trait[3]} --make-grm-gz --out Rare_${output}/RareSNPs${output}_Chrom_${trait[3]}   

gcta64 --bfile ../../InputData/sheep_geno_imputed_Plates1to97_20220627 --chr ${trait[3]} --make-grm --out Rare_${output}/Chrom_${trait[3]}   
gcta64 --grm Rare_${output}/Chrom_${trait[3]}    --make-grm-gz --out Rare_${output}/Chrom_${trait[3]}      
 
echo "MAF_${maf} Rare_${output}/RareSNPs${output}_Chrom_${trait[3]}
Full Rare_${output}/Chrom_${trait[3]}" > Rare_${output}/GRM_File.txt
mpirun -np 8 dissect.mpich --reml --gcta-grms-gz Rare_${output}/GRM_File.txt --pheno ${trait[0]}Phenotypes.txt --pheno-col 1 --blue --reml-maxit 100  --out Rare_${output}/${trait[0]}_MAF${output}_H2_Chrom_${trait[3]}

## 0.1 Chrom region
plink --bfile ../../InputData/sheep_geno_imputed_Plates1to97_20220627 --chr ${trait[3]} --from-kb ${trait[4]} --to-kb ${trait[5]} --max-maf ${maf} --nonfounders --write-snplist --out Rare_${output}/RareSNPs${output}_Chrom_${trait[3]}_region
gcta64 --bfile ../../InputData/sheep_geno_imputed_Plates1to97_20220627 --extract Rare_${output}/RareSNPs${output}_Chrom_${trait[3]}_region.snplist --out Rare_${output}/RareSNPs${output}_Chrom_${trait[3]}_region 
gcta64 --grm Rare_${output}/RareSNPs${output}_Chrom_${trait[3]}_region --make-grm-gz --out Rare_${output}/RareSNPs${output}_Chrom_${trait[3]}_region   

plink --bfile ../../InputData/sheep_geno_imputed_Plates1to97_20220627 --chr ${trait[3]} --from-kb ${trait[4]} --to-kb ${trait[5]} --nonfounders --write-snplist --out Rare_${output}/Chrom_${trait[3]}_region
gcta64 --bfile ../../InputData/sheep_geno_imputed_Plates1to97_20220627 --extract Rare_${output}/Chrom_${trait[3]}_region.snplist --make-grm --out Rare_${output}/Chrom_${trait[3]}_region   
gcta64 --grm Rare_${output}/Chrom_${trait[3]}_region --make-grm-gz --out Rare_${output}/Chrom_${trait[3]}_region  
 
echo "MAF_${maf} Rare_${output}/RareSNPs${output}_Chrom_${trait[3]}_region
Full Rare_${output}/Chrom_${trait[3]}_region" > Rare_${output}/GRM_File.txt
mpirun -np 8 dissect.mpich --reml --gcta-grms-gz Rare_${output}/GRM_File.txt --pheno ${trait[0]}Phenotypes.txt --pheno-col 1 --blue --reml-maxit 100  --out Rare_${output}/${trait[0]}_MAF${output}_H2_Chrom_${trait[3]}_region


##0.05
maf=0.05
output=005
mkdir -p Rare_${output}

##0.05 Chrom
plink --bfile ../../InputData/sheep_geno_imputed_Plates1to97_20220627 --chr ${trait[3]} --max-maf ${maf} --nonfounders --write-snplist --out Rare_${output}/RareSNPs${output}_Chrom_${trait[3]}
gcta64 --bfile ../../InputData/sheep_geno_imputed_Plates1to97_20220627 --extract Rare_${output}/RareSNPs${output}_Chrom_${trait[3]}.snplist --make-grm --out Rare_${output}/RareSNPs${output}_Chrom_${trait[3]} 
gcta64 --grm Rare_${output}/RareSNPs${output}_Chrom_${trait[3]} --make-grm-gz --out Rare_${output}/RareSNPs${output}_Chrom_${trait[3]}   

gcta64 --bfile ../../InputData/sheep_geno_imputed_Plates1to97_20220627 --chr ${trait[3]} --make-grm --out Rare_${output}/Chrom_${trait[3]}   
gcta64 --grm Rare_${output}/Chrom_${trait[3]}    --make-grm-gz --out Rare_${output}/Chrom_${trait[3]}      
 
echo "MAF_${maf} Rare_${output}/RareSNPs${output}_Chrom_${trait[3]}
Full Rare_${output}/Chrom_${trait[3]}" > Rare_${output}/GRM_File.txt
mpirun -np 8 dissect.mpich --reml --gcta-grms-gz Rare_${output}/GRM_File.txt --pheno ${trait[0]}Phenotypes.txt --pheno-col 1 --blue --reml-maxit 100  --out Rare_${output}/${trait[0]}_MAF${output}_H2_Chrom_${trait[3]}

## 0.05 Chrom region
plink --bfile ../../InputData/sheep_geno_imputed_Plates1to97_20220627 --chr ${trait[3]} --from-kb ${trait[4]} --to-kb ${trait[5]} --max-maf ${maf} --nonfounders --write-snplist --out Rare_${output}/RareSNPs${output}_Chrom_${trait[3]}_region
gcta64 --bfile ../../InputData/sheep_geno_imputed_Plates1to97_20220627 --extract Rare_${output}/RareSNPs${output}_Chrom_${trait[3]}_region.snplist --out Rare_${output}/RareSNPs${output}_Chrom_${trait[3]}_region 
gcta64 --grm Rare_${output}/RareSNPs${output}_Chrom_${trait[3]}_region --make-grm-gz --out Rare_${output}/RareSNPs${output}_Chrom_${trait[3]}_region   

plink --bfile ../../InputData/sheep_geno_imputed_Plates1to97_20220627 --chr ${trait[3]} --from-kb ${trait[4]} --to-kb ${trait[5]} --nonfounders --write-snplist --out Rare_${output}/Chrom_${trait[3]}_region
gcta64 --bfile ../../InputData/sheep_geno_imputed_Plates1to97_20220627 --extract Rare_${output}/Chrom_${trait[3]}_region.snplist --make-grm --out Rare_${output}/Chrom_${trait[3]}_region   
gcta64 --grm Rare_${output}/Chrom_${trait[3]}_region --make-grm-gz --out Rare_${output}/Chrom_${trait[3]}_region  
 
echo "MAF_${maf} Rare_${output}/RareSNPs${output}_Chrom_${trait[3]}_region
Full Rare_${output}/Chrom_${trait[3]}_region" > Rare_${output}/GRM_File.txt
mpirun -np 8 dissect.mpich --reml --gcta-grms-gz Rare_${output}/GRM_File.txt --pheno ${trait[0]}Phenotypes.txt --pheno-col 1 --blue --reml-maxit 100  --out Rare_${output}/${trait[0]}_MAF${output}_H2_Chrom_${trait[3]}_region



##0.01
maf=0.01
output=001
mkdir -p Rare_${output}

##0.01 Chrom
plink --bfile ../../InputData/sheep_geno_imputed_Plates1to97_20220627 --chr ${trait[3]} --max-maf ${maf} --nonfounders --write-snplist --out Rare_${output}/RareSNPs${output}_Chrom_${trait[3]}
gcta64 --bfile ../../InputData/sheep_geno_imputed_Plates1to97_20220627 --extract Rare_${output}/RareSNPs${output}_Chrom_${trait[3]}.snplist --make-grm --out Rare_${output}/RareSNPs${output}_Chrom_${trait[3]} 
gcta64 --grm Rare_${output}/RareSNPs${output}_Chrom_${trait[3]} --make-grm-gz --out Rare_${output}/RareSNPs${output}_Chrom_${trait[3]}   

gcta64 --bfile ../../InputData/sheep_geno_imputed_Plates1to97_20220627 --chr ${trait[3]} --make-grm --out Rare_${output}/Chrom_${trait[3]}   
gcta64 --grm Rare_${output}/Chrom_${trait[3]}    --make-grm-gz --out Rare_${output}/Chrom_${trait[3]}      
 
echo "MAF_${maf} Rare_${output}/RareSNPs${output}_Chrom_${trait[3]}
Full Rare_${output}/Chrom_${trait[3]}" > Rare_${output}/GRM_File.txt
mpirun -np 8 dissect.mpich --reml --gcta-grms-gz Rare_${output}/GRM_File.txt --pheno ${trait[0]}Phenotypes.txt --pheno-col 1 --blue --reml-maxit 100  --out Rare_${output}/${trait[0]}_MAF${output}_H2_Chrom_${trait[3]}

## 0.01 Chrom region
plink --bfile ../../InputData/sheep_geno_imputed_Plates1to97_20220627 --chr ${trait[3]} --from-kb ${trait[4]} --to-kb ${trait[5]} --max-maf ${maf} --nonfounders --write-snplist --out Rare_${output}/RareSNPs${output}_Chrom_${trait[3]}_region
gcta64 --bfile ../../InputData/sheep_geno_imputed_Plates1to97_20220627 --extract Rare_${output}/RareSNPs${output}_Chrom_${trait[3]}_region.snplist --out Rare_${output}/RareSNPs${output}_Chrom_${trait[3]}_region 
gcta64 --grm Rare_${output}/RareSNPs${output}_Chrom_${trait[3]}_region --make-grm-gz --out Rare_${output}/RareSNPs${output}_Chrom_${trait[3]}_region   

plink --bfile ../../InputData/sheep_geno_imputed_Plates1to97_20220627 --chr ${trait[3]} --from-kb ${trait[4]} --to-kb ${trait[5]} --nonfounders --write-snplist --out Rare_${output}/Chrom_${trait[3]}_region
gcta64 --bfile ../../InputData/sheep_geno_imputed_Plates1to97_20220627 --extract Rare_${output}/Chrom_${trait[3]}_region.snplist --make-grm --out Rare_${output}/Chrom_${trait[3]}_region   
gcta64 --grm Rare_${output}/Chrom_${trait[3]}_region --make-grm-gz --out Rare_${output}/Chrom_${trait[3]}_region  
 
echo "MAF_${maf} Rare_${output}/RareSNPs${output}_Chrom_${trait[3]}_region
Full Rare_${output}/Chrom_${trait[3]}_region" > Rare_${output}/GRM_File.txt
mpirun -np 8 dissect.mpich --reml --gcta-grms-gz Rare_${output}/GRM_File.txt --pheno ${trait[0]}Phenotypes.txt --pheno-col 1 --blue --reml-maxit 100  --out Rare_${output}/${trait[0]}_MAF${output}_H2_Chrom_${trait[3]}_region


##0.005
maf=0.005
output=0005
mkdir -p Rare_${output}

##0.005 Chrom
plink --bfile ../../InputData/sheep_geno_imputed_Plates1to97_20220627 --chr ${trait[3]} --max-maf ${maf} --nonfounders --write-snplist --out Rare_${output}/RareSNPs${output}_Chrom_${trait[3]}
gcta64 --bfile ../../InputData/sheep_geno_imputed_Plates1to97_20220627 --extract Rare_${output}/RareSNPs${output}_Chrom_${trait[3]}.snplist --make-grm --out Rare_${output}/RareSNPs${output}_Chrom_${trait[3]} 
gcta64 --grm Rare_${output}/RareSNPs${output}_Chrom_${trait[3]} --make-grm-gz --out Rare_${output}/RareSNPs${output}_Chrom_${trait[3]}   

gcta64 --bfile ../../InputData/sheep_geno_imputed_Plates1to97_20220627 --chr ${trait[3]} --make-grm --out Rare_${output}/Chrom_${trait[3]}   
gcta64 --grm Rare_${output}/Chrom_${trait[3]}    --make-grm-gz --out Rare_${output}/Chrom_${trait[3]}      
 
echo "MAF_${maf} Rare_${output}/RareSNPs${output}_Chrom_${trait[3]}
Full Rare_${output}/Chrom_${trait[3]}" > Rare_${output}/GRM_File.txt
mpirun -np 8 dissect.mpich --reml --gcta-grms-gz Rare_${output}/GRM_File.txt --pheno ${trait[0]}Phenotypes.txt --pheno-col 1 --blue --reml-maxit 100  --out Rare_${output}/${trait[0]}_MAF${output}_H2_Chrom_${trait[3]}

## 0.005 Chrom region
plink --bfile ../../InputData/sheep_geno_imputed_Plates1to97_20220627 --chr ${trait[3]} --from-kb ${trait[4]} --to-kb ${trait[5]} --max-maf ${maf} --nonfounders --write-snplist --out Rare_${output}/RareSNPs${output}_Chrom_${trait[3]}_region
gcta64 --bfile ../../InputData/sheep_geno_imputed_Plates1to97_20220627 --extract Rare_${output}/RareSNPs${output}_Chrom_${trait[3]}_region.snplist --out Rare_${output}/RareSNPs${output}_Chrom_${trait[3]}_region 
gcta64 --grm Rare_${output}/RareSNPs${output}_Chrom_${trait[3]}_region --make-grm-gz --out Rare_${output}/RareSNPs${output}_Chrom_${trait[3]}_region   

plink --bfile ../../InputData/sheep_geno_imputed_Plates1to97_20220627 --chr ${trait[3]} --from-kb ${trait[4]} --to-kb ${trait[5]} --nonfounders --write-snplist --out Rare_${output}/Chrom_${trait[3]}_region
gcta64 --bfile ../../InputData/sheep_geno_imputed_Plates1to97_20220627 --extract Rare_${output}/Chrom_${trait[3]}_region.snplist --make-grm --out Rare_${output}/Chrom_${trait[3]}_region   
gcta64 --grm Rare_${output}/Chrom_${trait[3]}_region --make-grm-gz --out Rare_${output}/Chrom_${trait[3]}_region  
 
echo "MAF_${maf} Rare_${output}/RareSNPs${output}_Chrom_${trait[3]}_region
Full Rare_${output}/Chrom_${trait[3]}_region" > Rare_${output}/GRM_File.txt
mpirun -np 8 dissect.mpich --reml --gcta-grms-gz Rare_${output}/GRM_File.txt --pheno ${trait[0]}Phenotypes.txt --pheno-col 1 --blue --reml-maxit 100  --out Rare_${output}/${trait[0]}_MAF${output}_H2_Chrom_${trait[3]}_region


##0.001
maf=0.001
output=0001
mkdir -p Rare_${output}

##0.001 Chrom
plink --bfile ../../InputData/sheep_geno_imputed_Plates1to97_20220627 --chr ${trait[3]} --max-maf ${maf} --nonfounders --write-snplist --out Rare_${output}/RareSNPs${output}_Chrom_${trait[3]}
gcta64 --bfile ../../InputData/sheep_geno_imputed_Plates1to97_20220627 --extract Rare_${output}/RareSNPs${output}_Chrom_${trait[3]}.snplist --make-grm --out Rare_${output}/RareSNPs${output}_Chrom_${trait[3]} 
gcta64 --grm Rare_${output}/RareSNPs${output}_Chrom_${trait[3]} --make-grm-gz --out Rare_${output}/RareSNPs${output}_Chrom_${trait[3]}   

gcta64 --bfile ../../InputData/sheep_geno_imputed_Plates1to97_20220627 --chr ${trait[3]} --make-grm --out Rare_${output}/Chrom_${trait[3]}   
gcta64 --grm Rare_${output}/Chrom_${trait[3]}    --make-grm-gz --out Rare_${output}/Chrom_${trait[3]}      
 
echo "MAF_${maf} Rare_${output}/RareSNPs${output}_Chrom_${trait[3]}
Full Rare_${output}/Chrom_${trait[3]}" > Rare_${output}/GRM_File.txt
mpirun -np 8 dissect.mpich --reml --gcta-grms-gz Rare_${output}/GRM_File.txt --pheno ${trait[0]}Phenotypes.txt --pheno-col 1 --blue --reml-maxit 100  --out Rare_${output}/${trait[0]}_MAF${output}_H2_Chrom_${trait[3]}

## 0.001 Chrom region
plink --bfile ../../InputData/sheep_geno_imputed_Plates1to97_20220627 --chr ${trait[3]} --from-kb ${trait[4]} --to-kb ${trait[5]} --max-maf ${maf} --nonfounders --write-snplist --out Rare_${output}/RareSNPs${output}_Chrom_${trait[3]}_region
gcta64 --bfile ../../InputData/sheep_geno_imputed_Plates1to97_20220627 --extract Rare_${output}/RareSNPs${output}_Chrom_${trait[3]}_region.snplist --out Rare_${output}/RareSNPs${output}_Chrom_${trait[3]}_region 
gcta64 --grm Rare_${output}/RareSNPs${output}_Chrom_${trait[3]}_region --make-grm-gz --out Rare_${output}/RareSNPs${output}_Chrom_${trait[3]}_region   

plink --bfile ../../InputData/sheep_geno_imputed_Plates1to97_20220627 --chr ${trait[3]} --from-kb ${trait[4]} --to-kb ${trait[5]} --nonfounders --write-snplist --out Rare_${output}/Chrom_${trait[3]}_region
gcta64 --bfile ../../InputData/sheep_geno_imputed_Plates1to97_20220627 --extract Rare_${output}/Chrom_${trait[3]}_region.snplist --make-grm --out Rare_${output}/Chrom_${trait[3]}_region   
gcta64 --grm Rare_${output}/Chrom_${trait[3]}_region --make-grm-gz --out Rare_${output}/Chrom_${trait[3]}_region  
 
echo "MAF_${maf} Rare_${output}/RareSNPs${output}_Chrom_${trait[3]}_region
Full Rare_${output}/Chrom_${trait[3]}_region" > Rare_${output}/GRM_File.txt
mpirun -np 8 dissect.mpich --reml --gcta-grms-gz Rare_${output}/GRM_File.txt --pheno ${trait[0]}Phenotypes.txt --pheno-col 1 --blue --reml-maxit 100  --out Rare_${output}/${trait[0]}_MAF${output}_H2_Chrom_${trait[3]}_region
