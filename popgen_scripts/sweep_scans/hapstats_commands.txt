#these commands were used to convert a multi-sample VCF with population SNP data into the proper input file for the SelectionHapStats pipeline (Garud 2015, Harris 2018)
#to calculate G123 and G2/G1 statistics

#Author: Alex Samano 2024


module load GCC/11.2.0 BCFtools/1.14  
#restructure snp vcf data
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' kochisnpsonly_X23.annot.vcf > kochi_X23_onlysnps_filt.reformat

#use custom python script (make_mlg.py) to make multi-locus genotype format 

#then reformat the MLG for G12 input
cut -f1,2,5-22 kochi_X23_onlysnps.mlg | awk -v OFS="," '$1=="X"{$1=$1; print}' | cut -d ","  -f2-21 > kochi_X.input
cut -f1,2,5-22 kochi_X23_onlysnps.mlg | awk -v OFS="," '$1=="2"{$1=$1; print}' | cut -d ","  -f2-21 > kochi_2.input
cut -f1,2,5-22 kochi_X23_onlysnps.mlg | awk -v OFS="," '$1=="3"{$1=$1; print}' | cut -d ","  -f2-21 > kochi_3.input


#run hapstats
H12_H2H1=/scratch/user/asamano/tools/SelectionHapStats/scripts/H12_H2H1.py

python2 $H12_H2H1 kochi_X.input 18 -o kochi_chrX_G12.output.txt -w 100 -j 50 -d 0
python2 $H12_H2H1 kochi_2.input 18 -o kochi_chr2_G12.output.txt -w 100 -j 50 -d 0
python2 $H12_H2H1 kochi_3.input 18 -o kochi_chr3_G12.output.txt -w 100 -j 50 -d 0


#calculate threshold by simulating neutral evolution with msms
#assumming length 1 Gb, mu 1.36e-9, and bottlegrowth parameters inferred by dadi
#ms command used: msms 36 1 -t 9562.062 -r 703092.8 1000000 -en 0.0004881431 1 1.30724133e-03 -eg 0.0004881431 1 14941.06 -I 1 36
#highest G123 from 1000 neutral sims used as threshold
