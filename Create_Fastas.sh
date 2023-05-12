#!/bin/bash
#SBATCH --account=def-rieseber
#SBATCH --time=1-5
#SBATCH --mem=8G
module load StdEnv/2020 vcftools/0.1.16 gcc/9.3.0 bcftools/1.11 samtools/1.12

tabix -p vcf Annuus.ann_env.tranche90_snps_bi_AN50_AF99_MAF05_AN50_biallelic_DP5_12_NoTE_noMap_PASS.vcf.gz

bcftools view -S tex_selected.txt  Annuus.ann_env.tranche90_snps_bi_AN50_AF99_MAF05_AN50_biallelic_DP5_12_NoTE_noMap_PASS.vcf.gz -Oz > Annuus.ann_env.tranche90_snps_bi_AN50_AF99_MAF05_AN50_b
iallelic_DP5_12_NoTE_noMap_PASS_25.vcf.gz

zcat Annuus.ann_env.tranche90_snps_bi_AN50_AF99_MAF05_AN50_biallelic_DP5_12_NoTE_noMap_PASS_25.vcf.gz | grep -v "#" | grep "PASS" |  awk -v OFS="\t" '{print $1, $2}' > pos_to_keept.txt

grep -wFf pos_to_keept.txt <(zcat Annuus.ann_env.tranche90_snps_bi_AN50_beagle_AF99.vcf.gz) > Annuus.ann_env.tranche90_snps_bi_AN50_beagle_AF99_filtered.vcf

zcat Annuus.ann_env.tranche90_snps_bi_AN50_beagle_AF99.vcf.gz | grep "#" > head.txt

cat head.txt Annuus.ann_env.tranche90_snps_bi_AN50_beagle_AF99_filtered.vcf | bgzip -c > Annuus.ann_env.tranche90_snps_bi_AN50_AF99_MAF05_AN50_biallelic_DP5_12_beagle_PASS_25.vcf.gz

rm head.txt
m Annuus.ann_env.tranche90_snps_bi_AN50_beagle_AF99_filtered.vcf

tabix -p vcf Annuus.ann_env.tranche90_snps_bi_AN50_AF99_MAF05_AN50_biallelic_DP5_12_beagle_PASS_25.vcf.gz

bcftools view -S tex_selected.txt Annuus.ann_env.tranche90_snps_bi_AN50_AF99_MAF05_AN50_biallelic_DP5_12_beagle_PASS_25.vcf.gz -Oz > Annuus.ann_env.tranche90_snps_bi_AN50_AF99_MAF05_AN50_bia
llelic_DP5_12_NoTE_noMap_PASS_beagle_25_texas.vcf.gz


tabix -p vcf Annuus.ann_env.tranche90_snps_bi_AN50_AF99_MAF05_AN50_biallelic_DP5_12_NoTE_noMap_PASS_beagle_25_texas.vcf.gz


for i in $(cat tex_selected.txt)
do
   for j in $(cat Chromosomes.txt)
   do
   samtools faidx  Ha412HOv2.0-20181130.fasta $j | vcf-consensus -H 1 -s $i Annuus.ann_env.tranche90_snps_bi_AN50_AF99_MAF05_AN50_biallelic_DP5_12_NoTE_noMap_PASS_beagle_25_texas.vcf.gz > "$
i""_1_""$j"".fasta"
   samtools faidx  Ha412HOv2.0-20181130.fasta $j | vcf-consensus -H 2 -s $i Annuus.ann_env.tranche90_snps_bi_AN50_AF99_MAF05_AN50_biallelic_DP5_12_NoTE_noMap_PASS_beagle_25_texas.vcf.gz > "$
i""_2_""$j"".fasta"
   done
done

for i in $(cat Chromosomes.txt)
do
cat *"$i".fasta > "$i""_ALL.fasta"
done

for i in $(ls Ha412HOChr*_ALL.fasta)
do
name=$(echo $i | cut -d "." -f 1)
awk 'BEGIN{FS=" "}{if(!/>/){print toupper($0)}else{print $1}}' $i > "$name""_Cap.fasta"
done

Number=1
#change names
for j in $(ls -1v *_Cap.fasta)
do
awk '/^>/{print $0"_"++i; next}{print}' $j > "$Number""_""$j"
let Number=Number+1
done
