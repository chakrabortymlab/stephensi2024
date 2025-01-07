
#shuffle coordinates of SVs with allele frequency 0.25, count how many are within 5kb on either side of a sweep window
for i in {1..100000}
do
bedtools shuffle -chrom -g chromSizes.txt -i lkd_svs_over25.bed > lkd_svs_over25_shuffled.bed
count=$(bedtools window -w 5000 -a lkd_sweep_windows.bed -b lkd_svs_over25_shuffled.bed | cut -f8 | sort | uniq | wc -l)
echo "rep_$i $count"
done > lkd_100Kshuff_svs_nearpeak.txt.txt
