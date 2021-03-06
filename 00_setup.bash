mkdir ../assembly
mkdir ../cnv_analysis
mkdir ../cnv_analysis/female
mkdir ../cnv_analysis/female/v501
mkdir ../cnv_analysis/female/v501/results
mkdir ../cnv_analysis/female/v501/results/plots
mkdir ../cnv_analysis/female/v603
mkdir ../cnv_analysis/female/v603/results
mkdir ../cnv_analysis/female/v603/results/plots
mkdir ../cnv_analysis/male
mkdir ../cnv_analysis/male/v501
mkdir ../cnv_analysis/male/v501/results
mkdir ../cnv_analysis/male/v501/results/plots
mkdir ../cnv_analysis/male/v603
mkdir ../cnv_analysis/male/v603/results
mkdir ../cnv_analysis/male/v603/results/plots
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
mkdir ../clean_call
mkdir r_scripts
#Update the regions to dismiss / keep from the main repo
if [ "$(hostname)" == "dna-prd-app01.exe.nhs.uk" ]; then
	cp /mnt/data1/resources/diagnostic_intervals/variants_regions_to_dismiss/unwanted_variants.csv configuration/unwanted_variants.csv
	cp /mnt/data1/resources/diagnostic_intervals/variants_regions_to_process/wanted_regions.csv configuration/wanted_regions.csv
else
	cp /mnt/Data1/targeted_sequencing/intervals/variants_regions_to_dismiss/unwanted_variants.csv configuration/unwanted_variants.csv
	cp /mnt/Data1/targeted_sequencing/intervals/variants_regions_to_process/wanted_regions.csv configuration/wanted_regions.csv
fi

bundle install
