# alignment of An. gambiae and An. stephensi genomes to identify SVs in syntenic sequences to infer ancestral state
# Alex Samano 2024

# mask repeat sequences
RepeatMasker -nolow -norna -lib AGAP_repeat_lib.fa AnoGambNW_F1_1_genomic_X23.fna
RepeatMasker -nolow -norna -lib AGAP_repeat_lib.fa AnoGambNW_F1_1_genomic_X23.fna

asteph=UCI_ANSTEP_V1.0_X23_masked.fna
agam=AnoGambNW_F1_1_genomic_X23_masked.fna

# alignment using LASTZ and Chain/Net pipeline 
module load GCC/10.3.0  OpenMPI/4.1.1  LASTZ/1.04.15
lastz $asteph[multiple] $agam K=4000 L=5000 Y=3400 E=30 H=3000 O=400 T=1 --format=axt --out=asteph_agamb.axt 

# make chains
axtChain -linearGap=medium asteph_agamb.axt -faT $asteph -faQ $agam asteph_agamb.chain 

# make nets and get syntenic regions
chainPreNet asteph_agamb.chain asteph.sizes agam.sizes asteph_agamb_preNet.chain 
chainNet asteph_agamb_preNet.chain -minSpace=1 asteph.sizes agam.sizes asteph.net agam.net 
netSyntenic asteph.net asteph_syn.net 

# convert fastas to 2bit 
faToTwoBit $asteph asteph.2bit
faToTwoBit $agam agamb.2bit

# convert net to maf
netToAxt asteph_syn.net asteph_agamb_preNet.chain asteph.2bit agamb.2bit asteph_agamb_syn.axt 
axtSort asteph_agamb_syn.axt asteph_agamb_syn_sorted.axt 
axtToMaf asteph_agamb_syn_sorted.axt asteph.sizes agam.sizes -tPrefix=asteph_ -qPrefix=agamb_ asteph_agamb_syn.maf

# convert maf to bed
grep 'asteph' asteph_agamb_syn.maf | awk '{OFS="\t"; print $2,$3,$3+$4,$4}' | tail -n+2 > asteph_agamb_syn.bed
grep 'agamb' asteph_agamb_syn.maf | awk '{OFS="\t"; print $2,$3,$3+$4,$4}' | tail -n+2 > agamb_syn.bed

# overlap syntenic regions with sv breakpoints
bedtools intersect -wa -wb -f 1.0 -a fivepops115_dups.bed -b asteph_agamb_syn.bed > dups_syntenic_overlap.txt
bedtools intersect -wa -wb -f 1.0 -a fivepops115_dels.bed -b asteph_agamb_syn.bed > dels_syntenic_overlap.txt