# variants

mapreads_callvariants.sh: Commands to perform quality checking on Illumina reads, mapping to An. stephensi reference genome, and calling SNPs and SVs.

filter_merge_variants.sh: Commands for filtering SNP and SV calls and merging VCFs.

VCFs: Full VCF output from merging individual SV VCFs with Jasmine. To download the annotated SNP VCF, follow this [link](https://zenodo.org/records/14542541)

sv_enrichment_analysis: Commands and data used to perform enrichment analysis of high-frequency SVs associated with CLR peaks. 

agam_asteph_synteny.sh: Alignment An. gambiae and An. stephensi genomes (chrX,2,3) to identify syntenic regions, then identify SVs within these syntenic blocks. 

simulate_missing_genotypes: Script for simulating missing genotypes in SNP data to evaluate effect on non-synonymous site frequency spectrum.
