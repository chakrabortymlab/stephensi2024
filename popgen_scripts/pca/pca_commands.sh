#commands for performing principal component analysis of genome-wide SNPs with maf>0.05

#filter snps 
module load GCC/11.2.0 VCFtools/0.1.16
vcftools --vcf fivepops115_X23_snps.vcf.gz --maf 0.05  --recode --out fivepops115_X23_snps_filteredmaf


module load GCC/12.3.0 PLINK/2.00a3.7
# perform linkage pruning -  identify prune sites
plink --vcf fivepops115_X23_snps_filteredmaf.recode.vcf --double-id --allow-extra-chr \
--set-missing-var-ids @:# \
--indep-pairwise 50 10 0.1 --out fivepops115_pca

#perform PCA
plink --vcf fivepops115_X23_snps_filteredmaf.recode.vcf --double-id --allow-extra-chr --set-missing-var-ids @:# \
--extract fivepops115_pca.prune.in \
--make-bed --pca --out fivepops115_pca
