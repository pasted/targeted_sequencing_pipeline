mkdir ../assembly
mkdir ../cnv_analysis
mkdir ../coverage
mkdir ../duplicates
mkdir ../intervals
mkdir ../logs
mkdir ../metrics
mkdir ../raw_reads
mkdir ../results
mkdir ../summary
mkdir ../variants
mkdir ../variants/haplotyper
mkdir ../variants/freebayes
mkdir ../variants/platypus
mkdir ../variants_6q24
mkdir ../variants_t1d
#Update the regions to dismiss / keep from the main repo
cp /mnt/Data1/targeted_sequencing/intervals/variants_regions_to_dismiss/unwanted_variants.csv configuration/unwanted_variants.csv
cp /mnt/Data1/targeted_sequencing/intervals/variants_regions_to_process/wanted_regions.csv configuration/wanted_regions.csv
