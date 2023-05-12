#snp pos files

#!/bin/bash
#SBATCH --account=def-rieseber
#SBATCH --time=1-5
#SBATCH --mem=8G
module load StdEnv/2020 vcftools/0.1.16

for i in $(cat Chromosomes.txt);do  vcftools --gzvcf Annuus.ann_env.tranche90_snps_bi_AN50_AF99_MAF05_AN50_biallelic_DP5_12_NoTE_noMap_PASS_beagle_25_texas.vcf.gz --ldhelmet --chr $i --phase
d --out "$i""_Ldhelmet_SNPs"; done

#more Chromosomes.txt
#Ha412HOChr01
#Ha412HOChr02
#Ha412HOChr03
#...
#Ha412HOChr17
