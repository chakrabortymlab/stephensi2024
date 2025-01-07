#commands for filtering and merging SNP and SV calls
#Alex Samano 2024

##### SNPs #####
# merge population VCFs 
bcftools merge -0 -l fivepops115_vcf_list.txt > ansteph_fivepops115_merged_freebayes.vcf
#keep only SNPs, chromosomes 2,3,X, non missing data in at least 75% of individuals, and minQ=20
vcftools --remove-indels --max-missing 0.75 --minQ 20 --chr chr2 --chr chr3 --chr chrX --vcf ansteph_fivepops115_merged_freebayes.vcf --recode --recode-INFO-all --out asteph_fivepops115_snps_filt
#annotate SNPs using SnpEff
snpEff -v -stats -nodownload fivepops115_stats.html Anopheles_stephensi_build asteph_fivepops115_snps_filt.vcf > asteph_fivepops115_snps_filt_annot.vcf



##### SVs #####
# filter individual VCFs by size, SV type, and quality filter
# then add strain and chromosome to SV ID so that this information is retained when merged by Jasmine
while read strain 
do
bcftools view -i '(INFO/END-POS) >= 100 && (INFO/END-POS) <= 100000 && FILTER="PASS" && ALT="<DUP>"' $strain.dellycalls.vcf > $strain.DUPS_filtered.vcf
bcftools view -i '(INFO/END-POS) >= 100 && (INFO/END-POS) <= 100000 && FILTER="PASS" && ALT="<DEL>"' $strain.dellycalls.vcf > $strain.DELS_filtered.vcf
awk '$1 ~ /^#/ {OFS="\t"; print $0;next}{OFS="\t"; print $1,$2,$1"_"$3"_'"$strain"',$4,$5,$6,$7,$8,$9,$10}' $strain.DUPS_filtered.vcf > $strain.DUPS_forjasmine_filtered.vcf
awk '$1 ~ /^#/ {OFS="\t"; print $0;next}{OFS="\t"; print $1,$2,$1"_"$3"_'"$strain"',$4,$5,$6,$7,$8,$9,$10}' $strain.DELS_filtered.vcf > $strain.DELS_forjasmine_filtered.vcf
done < strain_id_list.txt 

# merge with Jasmine (conda package)
ls *.DUPS_forjasmine_filtered.vcf > dups_list.txt
ls *.DELS_forjasmine_filtered.vcf > dels_list.txt

jasmine file_list=dups_list.txt out_file=stephensi_fivepops115_DUPS.vcf --ignore_strand --output_genotypes
jasmine file_list=dels_list.txt out_file=stephensi_fivepops115_DELS.vcf --ignore_strand --output_genotypes

#convert VCFs to bed
bcftools query -f '%CHROM\t%POS\t%END\n' stephensi_fivepops115_DUPS.vcf | sort -k1,1 -k2,2n > stephensi_fivepops115_DUPS.bed
bcftools query -f '%CHROM\t%POS\t%END\n' stephensi_fivepops115_DELS.vcf | sort -k1,1 -k2,2n > stephensi_fivepops115_DELS.bed

#identify deletions with over 90% overlap with annotated reference TEs
te_annot=asteph.fasta.mod.EDTA.TEanno.gff
bedtools intersect -f 0.9 -wa -a stephensi_fivepops115_DELS.bed -b $te_annot> stephensi_fivepops115_DELS_TEoverlap.bed
bedtools intersect -v -f 0.9 -wa -a stephensi_fivepops115_DELS.bed -b $te_annot> stephensi_fivepops115_DELS_noTE.bed

