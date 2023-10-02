#!/bin/sh
#$ -cwd
#$ -l h_rt=47:30:00
#$ -pe sharedmem 4
#$ -l h_vmem=8G
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

##No Thresholded GRM (model 1)
mkdir -p GRM_Null
echo "GRM ../../GRMs/GRM_GCTA" > GRM_Null/GRM_File.txt
mpirun -np 8 dissect.mpich --reml --gcta-grms-gz GRM_Null/GRM_File.txt --pheno ${trait[0]}Phenotypes.txt --pheno-col 1 --blue --reml-maxit 100 --covar ${trait[0]}Assoc.txt --qcovar ${trait[0]}QAssoc.txt --random-effects ${trait[0]}Random.txt --random-effects-cols ${trait[1]} --out GRM_Null/${trait[0]}_Full_H2

##Relatedness thresholded model (model 2) t=0.05
mkdir -p GRM_005
echo "GRM_0.05 ../../GRMs/GRM_005
GRM ../../GRMs/GRM_GCTA" > GRM_005/GRM_File.txt
mpirun -np 8 dissect.mpich --reml --gcta-grms-gz GRM_005/GRM_File.txt --pheno ${trait[0]}Phenotypes.txt --pheno-col 1 --blue --reml-maxit 100 --covar ${trait[0]}Assoc.txt --qcovar ${trait[0]}QAssoc.txt --random-effects ${trait[0]}Random.txt --random-effects-cols ${trait[1]} --out GRM_005/${trait[0]}_005_H2

##Relatedness thresholded model (model 2) t=0.1
mkdir -p GRM_01
echo "GRM_0.1 ../../GRMs/GRM_01
GRM ../../GRMs/GRM_GCTA" > GRM_01/GRM_File.txt
mpirun -np 8 dissect.mpich --reml --gcta-grms-gz GRM_01/GRM_File.txt --pheno ${trait[0]}Phenotypes.txt --pheno-col 1 --blue --reml-maxit 100 --covar ${trait[0]}Assoc.txt --qcovar ${trait[0]}QAssoc.txt --random-effects ${trait[0]}Random.txt --random-effects-cols ${trait[1]} --out GRM_100/${trait[0]}_01_H2

##MAF thresholded model (model 3) MAF=0.1
mkdir -p Rare_01
echo "MAF_0.1 ../../GRMs/GRM_RareSNPs01
GRM ../../GRMs/GRM_GCTA" > Rare_01/GRM_File.txt
mpirun -np 8 dissect.mpich --reml --gcta-grms-gz Rare_01/GRM_File.txt --pheno ${trait[0]}Phenotypes.txt --pheno-col 1 --blue --reml-maxit 100 --covar ${trait[0]}Assoc.txt --qcovar ${trait[0]}QAssoc.txt --random-effects ${trait[0]}Random.txt --random-effects-cols ${trait[1]} --out Rare_01/${trait[0]}_MAF01_H2

##MAF thresholded model (model 3) MAF=0.05
mkdir -p Rare_005
echo "MAF_0.05 ../../GRMs/GRM_RareSNPs005
GRM ../../GRMs/GRM_GCTA" > Rare_005/GRM_File.txt
mpirun -np 8 dissect.mpich --reml --gcta-grms-gz Rare_005/GRM_File.txt --pheno ${trait[0]}Phenotypes.txt --pheno-col 1 --blue --reml-maxit 100 --covar ${trait[0]}Assoc.txt --qcovar ${trait[0]}QAssoc.txt --random-effects ${trait[0]}Random.txt --random-effects-cols ${trait[1]} --out Rare_005/${trait[0]}_MAF005_H2

##MAF thresholded model (model 3) MAF=0.01
mkdir -p Rare_001
echo "MAF_0.01 ../../GRMs/GRM_RareSNPs001
GRM ../../GRMs/GRM_GCTA" > Rare_001/GRM_File.txt
mpirun -np 8 dissect.mpich --reml --gcta-grms-gz Rare_001/GRM_File.txt --pheno ${trait[0]}Phenotypes.txt --pheno-col 1 --blue --reml-maxit 100 --covar ${trait[0]}Assoc.txt --qcovar ${trait[0]}QAssoc.txt --random-effects ${trait[0]}Random.txt --random-effects-cols ${trait[1]} --out Rare_001/${trait[0]}_MAF001_H2

##MAF thresholded model (model 3) MAF=0.005
mkdir -p Rare_0005
echo "MAF_0.005 ../../GRMs/GRM_RareSNPs0005
GRM ../../GRMs/GRM_GCTA" > Rare_0005/GRM_File.txt
mpirun -np 8 dissect.mpich --reml --gcta-grms-gz Rare_0005/GRM_File.txt --pheno ${trait[0]}Phenotypes.txt --pheno-col 1 --blue --reml-maxit 100 --covar ${trait[0]}Assoc.txt --qcovar ${trait[0]}QAssoc.txt --random-effects ${trait[0]}Random.txt --random-effects-cols ${trait[1]} --out Rare_0005/${trait[0]}_MAF0005_H2

##MAF thresholded model (model 3) MAF=0.001
mkdir -p Rare_0001
echo "MAF_0.001 ../../GRMs/GRM_RareSNPs0001
GRM ../../GRMs/GRM_GCTA" > Rare_0001/GRM_File.txt
mpirun -np 8 dissect.mpich --reml --gcta-grms-gz Rare_0001/GRM_File.txt --pheno ${trait[0]}Phenotypes.txt --pheno-col 1 --blue --reml-maxit 100 --covar ${trait[0]}Assoc.txt --qcovar ${trait[0]}QAssoc.txt --random-effects ${trait[0]}Random.txt --random-effects-cols ${trait[1]} --out Rare_0001/${trait[0]}_MAF0001_H2

##Pedigree model (model 4)
mkdir -p Pedigree
echo "Pedigree ../../GRMs/AMatrix_2_4/Ped_2_4" > Pedigree/Pedigree_File.txt
mpirun -np 8 dissect.mpich --reml --gcta-grms-gz Pedigree/Pedigree_File.txt --pheno ${trait[0]}Phenotypes.txt --pheno-col 1 --blue --reml-maxit 100 --covar ${trait[0]}Assoc.txt --qcovar ${trait[0]}QAssoc.txt --random-effects ${trait[0]}Random.txt --random-effects-cols ${trait[1]} --out Pedigree/${trait[0]}_Full_H2
