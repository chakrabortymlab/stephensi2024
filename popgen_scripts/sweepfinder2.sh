#script used to convert population snp data in vcf format to 
#input files for SweepFinder2, then run SW2

#remove indels and keep only chr2,3,X
vcftools --remove-indels --minQ 20 --recode --chr chr2 --chr chr3 --chr chrX --vcf tvm_wild_04012023.vcf --out tvm_X23_onlysnps_filtered
#get snp allele counts
vcftools --vcf tvm_X23_onlysnps_filtered.vcf --counts2 --out tvm_X23_snps
#reformat to create sweepfinder input files
tail -n+2 tvm_X23_snps.frq.count | awk -v OFS='\t' '{print $1,$2,$6,$4,1}' > tvm_X23_snps.input

awk '$1 =="chrX"{print $2"\t" $3"\t" $4"\t" $5}' tvm_X23_snps.input  > tvm_chromX_sf_input.txt
awk '$1 =="chr2"{print $2"\t" $3"\t" $4"\t" $5}' tvm_X23_snps.input  > tvm_chrom2_sf_input.txt
awk '$1 =="chr3"{print $2"\t" $3"\t" $4"\t" $5}' tvm_X23_snps.input  > tvm_chrom3_sf_input.txt

echo -e "position\tx\tn\tfolded"| cat - tvm_chromX_sf_input.txt > temp.txt && mv temp.txt tvm_chromX_sf_input.txt
echo -e "position\tx\tn\tfolded"| cat - tvm_chrom2_sf_input.txt > temp.txt && mv temp.txt tvm_chrom2_sf_input.txt
echo -e "position\tx\tn\tfolded"| cat - tvm_chrom3_sf_input.txt > temp.txt && mv temp.txt tvm_chrom3_sf_input.txt

#run sweepfinder2 on each chromosome, test site every 5kb
SweepFinder2 -sg 5000 tvm_chromX_sf_input.txt tvm_chrX_sf.sweeps
SweepFinder2 -sg 5000 tvm_chrom2_sf_input.txt tvm_chr2_sf.sweeps
SweepFinder2 -sg 5000 tvm_chrom3_sf_input.txt tvm_chr3_sf.sweeps


