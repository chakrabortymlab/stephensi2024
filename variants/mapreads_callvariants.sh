###### Read quality checking and mapping commands ######
# sample id in this code is 'B4'

# Quality checking of raw reads
fastqc B4_R1.fastq.gz -o B4_R1_pre_FastQC.txt
fastqc B4_R2.fastq.gz -o B4_R2_pre_FastQC.txt

# Trimming of adapters
java -jar /softwares/trimmomatic/bin/trimmomatic.jar PE -phred33 \
B4_R1.fastq.gz B4_R2.fastq.gz \
B4_R1_paired_trim.fq.gz B4_R1_unpaired.fq.gz \
B4_R2_paired_trim.fq.gz B4_R2_unpaired.fq.gz \
ILLUMINACLIP:/softwares/trimmomatic/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

# check the quality of reads after removing the adapters
fastqc B4_R1_paired_trim.fastq.gz -o B4_R1_post_FastQC.txt
fastqc B4_R2_paired_trim.fastq.gz -o B4_R2_post_FastQC.txt

# Indexing of reference genome, bwa-0.7.17
bwa index anstep_genomev2.fasta

# Mapping of trimmed reads against a genome, bwa-0.7.17
bwa mem -t 16 -M anstep_genomev2.fasta B4_R1_paired_trim.fastq.gz B4_R2_paired_trim.fastq.gz > B4.sam

# SAM to BAM format, samtools-1.9
samtools view -bS -@12 B4.sam > B4.bam

# BAM to Sorted.bam, samtools-1.9
samtools sort B4.bam > B4.sorted.bam

# Marking the duplicates (Picard tools)
java -jar /home/softwares/picard/build/libs/picard.jar MarkDuplicates \
REMOVE_DUPLICATES=true I=B4.sorted.bam O=B4.sorted.nodup.bam \
M=B4_duplicatedata.txt

#Add read group name to the samples
java -jar /home/softwares/picard/build/libs/picard.jar AddOrReplaceReadGroups I=B4.sorted.nodup.bam \
O=B4.sorted.RGnodup.bam RGID=B4 RGLB=WGS RGPL=illumina RGPU=unit1 RGSM=B4

# Indexing the final bam file, samtools-1.9
samtools index B4.sorted.RGnodup.bam 

#Checking of Bam statistics using BAMQC(Qualimap)
bamqc B4.sorted.RGnodup.bam

#bam validate to check the bam file status
bam validate B4.sorted.RGnodup.bam




###### Variant Calling (SNPs, indels, SVs)
# call SNPs and indels by population bam file list (ex. bangalore_bamfiles_list.txt)
freebayes -f anstep_genomev2.fasta -0 -C 2 -L bangalore_bamfiles_list.txt > bangalore_wild_04012023.vcf
vcftools --gzvcf bangalore_wild_04012023.vcf.gz --max-missing 0.5 --minQ 20 --recode --recode-INFO-all --out bangalore_wild_04012023_filt

# run delly to call SVs by individual bam (ex. B4)
delly call -g anstep_genomev2.fasta B4.sorted.RGnodup.bam > B4.dellycalls.vcf


