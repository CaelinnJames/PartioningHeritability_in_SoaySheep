#!/bin/sh
#$ -cwd
#$ -l h_rt=23:00:00
#$ -pe sharedmem 1
#$ -l h_vmem=75G
#$ -o o_files/
#$ -e e_files/
. /etc/profile.d/modules.sh

module load igmm/apps/dissect/1.15.2c
module load igmm/apps/plink/1.90b4
module load igmm/apps/gcta/1.91.2beta

i=$SGE_TASK_ID

traitinfo=$(sed "${i}q;d" ../traits.txt)

IFS=','; trait=($traitinfo); unset IFS;

mkdir -p ${trait[0]}

cd ${trait[0]}

##Preparing genotype and phenotype data
sort -k 3,3 ${trait[0]}Phenotypes.txt > Pheno1HD.txt
awk '($2 >0)' ../../InputData/sheep_geno_imputed_Plates1to97_20220627.ped > Ped1HD.ped
sort -k 2,2 Ped1HD.ped > Ped2HD.ped
join -1 3 -2 2 Pheno1HD.txt Ped2HD.ped > Ped3HD.ped
cut -d ' ' -f 1,4,5 --complement Ped3HD.ped > ${trait[0]}HD.ped
rm Ped*HD.ped
rm Pheno1HD.txt
plink --file ../../InputData/sheep_geno_imputed_Plates1to97_20220627 --ped ${trait[0]}HD.ped --sheep --chr 1-26 --make-bed --out ${trait[0]}BFile
gcta64 --bfile ${trait[0]}BFile --make-grm --autosome --autosome-num 26 --out GRM_GCTA
gcta64 --grm GRM_GCTA --make-grm-gz --out GRM_GCTA
awk '{$2=$4="";print $0}' ${trait[0]}Phenotypes.txt  > ${trait[0]}_IDs.txt
sort ${trait[0]}_IDs.txt > tmp.txt
uniq tmp.txt > ${trait[0]}_IDs.txt

##No Thresholded GRM (model 1)
mkdir -p GRM_Null
echo "GRM GRM_GCTA" > GRM_Null/GRM_File.txt
mpirun -np 8 dissect.mpich --reml --gcta-grms-gz GRM_Null/GRM_File.txt --pheno ${trait[0]}Phenotypes.txt --pheno-col 2 --blue --reml-maxit 100 --covar ${trait[0]}Assoc.txt --qcovar ${trait[0]}QAssoc.txt --random-effects ${trait[0]}Random.txt --random-effects-cols ${trait[1]} --out GRM_Null/${trait[0]}_Full_H2

##Relatedness thresholded model (model 2) t=0.05
mkdir -p GRM_005
gcta64 --grm GRM_GCTA --make-bK 0.050 --out GRM_005/GRM_005_GCTA
gcta64 --grm GRM_005/GRM_005_GCTA --make-grm-gz --out GRM_005/GRM_005_GCTA
echo "GRM_0.05 GRM_005/GRM_005_GCTA
GRM GRM_GCTA" > GRM_005/GRM_File.txt
mpirun -np 8 dissect.mpich --reml --gcta-grms-gz GRM_005/GRM_File.txt --pheno ${trait[0]}Phenotypes.txt --pheno-col 2 --blue --reml-maxit 100 --covar ${trait[0]}Assoc.txt --qcovar ${trait[0]}QAssoc.txt --random-effects ${trait[0]}Random.txt --random-effects-cols ${trait[1]} --out GRM_005/${trait[0]}_005_H2

##Relatedness thresholded model (model 2) t=0.1
mkdir -p GRM_01
gcta64 --grm GRM_GCTA --make-bK 0.1 --out GRM_01/GRM_01_GCTA
gcta64 --grm GRM_01/GRM_01_GCTA --make-grm-gz --out GRM_01/GRM_01_GCTA
echo "GRM_0.1 GRM_01/GRM_01_GCTA
GRM GRM_GCTA" > GRM_01/GRM_File.txt
mpirun -np 8 dissect.mpich --reml --gcta-grms-gz GRM_01/GRM_File.txt --pheno ${trait[0]}Phenotypes.txt --pheno-col 2 --blue --reml-maxit 100 --covar ${trait[0]}Assoc.txt --qcovar ${trait[0]}QAssoc.txt --random-effects ${trait[0]}Random.txt --random-effects-cols ${trait[1]} --out GRM_01/${trait[0]}_01_H2

##MAF thresholded model (model 3) MAF=0.1
mkdir -p Rare_01
gcta64 --bfile ${trait[0]}BFile --autosome-num 26 --autosome --extract ../../GRMs/RareSNPs01.snplist --make-grm-gz --out  Rare_01/RareSNPs01  
echo "MAF_0.1  Rare_01/RareSNPs01
GRM GRM_GCTA" > Rare_01/GRM_File.txt
mpirun -np 8 dissect.mpich --reml --gcta-grms-gz Rare_01/GRM_File.txt --pheno ${trait[0]}Phenotypes.txt --pheno-col 2 --blue --reml-maxit 100 --covar ${trait[0]}Assoc.txt --qcovar ${trait[0]}QAssoc.txt --random-effects ${trait[0]}Random.txt --random-effects-cols ${trait[1]} --out Rare_01/${trait[0]}_MAF01_H2

##MAF thresholded model (model 3) MAF=0.05
mkdir -p Rare_005
gcta64 --bfile ${trait[0]}BFile --autosome-num 26 --autosome --extract ../../GRMs/RareSNPs005.snplist --make-grm-gz --out  Rare_005/RareSNPs005  
echo "MAF_0.05 Rare_005/RareSNPs005
GRM GRM_GCTA" > Rare_005/GRM_File.txt
mpirun -np 8 dissect.mpich --reml --gcta-grms-gz Rare_005/GRM_File.txt --pheno ${trait[0]}Phenotypes.txt --pheno-col 2 --blue --reml-maxit 100 --covar ${trait[0]}Assoc.txt --qcovar ${trait[0]}QAssoc.txt --random-effects ${trait[0]}Random.txt --random-effects-cols ${trait[1]} --out Rare_005/${trait[0]}_MAF005_H2

#MAF thresholded model (model 3) MAF=0.01
mkdir -p Rare_001
gcta64 --bfile ${trait[0]}BFile --autosome-num 26 --autosome --extract ../../GRMs/RareSNPs001.snplist --make-grm-gz --out  Rare_001/RareSNPs001  
echo "MAF_0.01  Rare_001/RareSNPs001
GRM GRM_GCTA" > Rare_001/GRM_File.txt
mpirun -np 8 dissect.mpich --reml --gcta-grms-gz Rare_001/GRM_File.txt --pheno ${trait[0]}Phenotypes.txt --pheno-col 2 --blue --reml-maxit 100 --covar ${trait[0]}Assoc.txt --qcovar ${trait[0]}QAssoc.txt --random-effects ${trait[0]}Random.txt --random-effects-cols ${trait[1]} --out Rare_001/${trait[0]}_MAF001_H2

##MAF thresholded model (model 3) MAF=0.005
mkdir -p Rare_0005
gcta64 --bfile ${trait[0]}BFile --autosome-num 26 --autosome --extract ../../GRMs/RareSNPs0005.snplist --make-grm-gz --out  Rare_0005/RareSNPs0005  
echo "MAF_0.005 Rare_0005/RareSNPs0005
GRM GRM_GCTA" > Rare_0005/GRM_File.txt
mpirun -np 8 dissect.mpich --reml --gcta-grms-gz Rare_0005/GRM_File.txt --pheno ${trait[0]}Phenotypes.txt --pheno-col 2 --blue --reml-maxit 100 --covar ${trait[0]}Assoc.txt --qcovar ${trait[0]}QAssoc.txt --random-effects ${trait[0]}Random.txt --random-effects-cols ${trait[1]} --out Rare_0005/${trait[0]}_MAF0005_H2

##MAF thresholded model (model 3) MAF=0.001
mkdir -p Rare_0001
gcta64 --bfile ${trait[0]}BFile --autosome-num 26 --autosome --extract ../../GRMs/RareSNPs0001.snplist --make-grm-gz --out  Rare_0001/RareSNPs0001  
echo "MAF_0.001  Rare_0001/RareSNPs0001
GRM GRM_GCTA" > Rare_0001/GRM_File.txt
mpirun -np 8 dissect.mpich --reml --gcta-grms-gz Rare_0001/GRM_File.txt --pheno ${trait[0]}Phenotypes.txt --pheno-col 2 --blue --reml-maxit 100 --covar ${trait[0]}Assoc.txt --qcovar ${trait[0]}QAssoc.txt --random-effects ${trait[0]}Random.txt --random-effects-cols ${trait[1]} --out Rare_0001/${trait[0]}_MAF0001_H2

##Pedigree model (model 4)
mkdir -p Pedigree
echo "Pedigree Ped_2_4" > Pedigree/Pedigree_File.txt
mpirun -np 8 dissect.mpich --reml --gcta-grms-gz Pedigree/Pedigree_File.txt --pheno ${trait[0]}Phenotypes.txt --pheno-col 2 --blue --reml-maxit 100 --covar ${trait[0]}Assoc.txt --qcovar ${trait[0]}QAssoc.txt --random-effects ${trait[0]}Random.txt --random-effects-cols ${trait[1]} --out Pedigree/${trait[0]}_Full_H2
