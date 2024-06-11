#This script performs conditional analysis using GCTA tool. 
#GCTA-COJO: multi-SNP-based conditional & joint association analysis using GWAS summary data
#https://yanglab.westlake.edu.cn/software/gcta/#COJO

#Step 1
#Create a file including SNPs that are conditioned: *.cond.snplist

#Step 2
#Prepare input file: *.ma

#Step 3
#Prepare LD reference file 
#Here, we use 1000G reference. 
#UKB reference in discussion.

ref=/mnt/zfs/P001_annots/1000G/data/ph3v5a/vol1/ftp/release/20130502/plink
in=/mnt/zfs/P201_UKB_pQTL/results/TREM2/cojo/input
out=/mnt/zfs/P201_UKB_pQTL/results/TREM2/cojo/output

plink --bfile $ref/ALL.chr.phase3v5a.genotypes --chr 6 --make-bed --out $in/TREM2_Region_1KG

#Stpe 4 
#Run conditional analysis 

assay="ASGR2 CCL20 NEFL TREM2 TREML2"

for i in $assay; 

do 
gcta64  --bfile $in/TREM2_Region_1KG  --chr 6 --cojo-file $in/"$i"_pQTL_in_TREM2_Region.ma --cojo-cond $in/"$i".cond.snplist --out $out/"$i"_pQTL_in_TREM2_Region_cond

done

#Note: --cojo-cond is the option relevant to our analysis. --cojo-slct, --cojo-top-SNPs and --cojo-join can be considered to understand independent signal in the gene region of interest.
#1. Select multiple associated SNPs through a stepwise selection procedure
#gcta64  --bfile $in/TREM2_Region_1KG  --chr 6 --maf 0.01 --cojo-p 5e-8 --cojo-file $in/TREML2_pQTL_in_TREM2_Region.ma --cojo-slct --out #$out/TREML2_pQTL_in_TREM2_Region_slct

#2. Select a fixed number of of top associated SNPs through a stepwise selection procedure
#gcta64  --bfile $in/TREM2_Region_1KG  --chr 6 --maf 0.01 --cojo-p 5e-8 --cojo-file $in/TREML2_pQTL_in_TREM2_Region.ma --cojo-top-SNPs 10 --out #$out/TREML2_pQTL_in_TREM2_Region_topSNPs

#3. Estimate the joint effects of a subset of SNPs (given in the file test.snplist) without model selection
#gcta64  --bfile $in/TREM2_Region_1KG  --chr 6 --extract test.snplist --cojo-file $in/TREML2_pQTL_in_TREM2_Region.ma --cojo-joint --out TREML2_pQTL_in_TREM2_Region_joint

#4. Perform single-SNP association analyses conditional on a set of SNPs (given in the file cond.snplist) without model selection
#gcta64  --bfile $in/TREM2_Region_1KG  --chr 6 --maf 0.01 --cojo-file $in/TREML2_pQTL_in_TREM2_Region.ma --cojo-cond $in/cond.snplist --out $out/TREML2_pQTL_in_TREM2_Region_cond

