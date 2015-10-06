class Pipeline
		require 'yaml'
		require 'logger'
		require 'parallel'
		require 'awesome_print'
		
		require_relative '../shared/batch'
		require_relative '../shared/panel'
		require_relative 'sample_parser'
		require_relative 'wget'
		require_relative 'rename'
		require_relative 'bwa_mem'
		require_relative 'sam_to_bam'
		require_relative 'fix_mate_pair'
		require_relative 'sort_bam'
		require_relative 'mark_duplicates'
		require_relative 'realign_indels'
		require_relative 'calculate_metrics'
		require_relative 'parse_metrics'
		require_relative 'variant_caller'
		require_relative 'filter_variants'
		require_relative 'alamut_annotation'
		require_relative 'variants_to_table'
		
		def error_check(out, sample, stage, logger)
			if out && (out[0] > 0)
				puts "ERROR :: #{stage} :: Halting pipeline at Sample #{sample.panel_version}_#{sample.sample_id}\n"
				ap "#{out[1].inspect}"
				logger.info('error') { "SAMPLE :: #{sample.sample_id} #{out[1].inspect}" }
				
				raise Parallel::Break
			end
		end
		
		logger = Logger.new('../../logs/pipeline.log')
		
		#Load the config YAML file and pass the settings to local variables
		this_batch = YAML.load_file('../configuration/config.yaml')
		
		this_parser = SampleParser.new
		this_pipeline = Pipeline.new
		
		samples = this_parser.parse_sample_list(this_batch.sample_list_path)
		
		samples_first = []
		
		samples.each do |this_sample|
			if [ "v501_1146", "v501_1147", "v501_1148", "v501_1149", "v501_1150", "v501_1151", "v501_1152" ].include? this_sample.capture_number
				samples_first.push this_sample	
			else
				#samples_first.push this_sample
			end
		end
		
		##ERRORS to check
		#sample.gender being blank - any of the require fields being set to blank
		#add sample validation code
		
		#sequencing service panel version having variable formats (or just the wrong one)
		#add translation step as early as possible
		
		#Add INTERMEDIATE file ext to all assembly files except final bam
		
		###ERRORS
		
		#Seperate loop for the WGET cmd due to a throttling issue
#		results = Parallel.map(samples, :in_processes=>1 ) do |this_sample|
#			this_wget = Wget.new
#			out = this_wget.fetch_fastq(this_sample, this_batch, logger)
#			this_pipeline.error_check(out, this_sample, "Wget FastQ", logger)
#			
#			this_rename = Rename.new
#			out = this_rename.remove_fastq_adaptor_string(this_sample, this_batch, logger)
#			this_pipeline.error_check(out, this_sample, "FastQ file rename", logger)
#		end      
		
		results = Parallel.map(samples_first, :in_processes=>32 ) do |this_sample|
			puts this_sample.inspect
			this_bwa_mem = BwaMem.new
			out = this_bwa_mem.run_bwa_mem(this_sample, this_batch, logger)
			this_pipeline.error_check(out, this_sample, "BWA-MEM alignment", logger)
			
			this_sam_to_bam = SamToBam.new
			out = this_sam_to_bam.convert_sam_to_bam(this_sample, this_batch, logger)
			this_pipeline.error_check(out, this_sample, "SAM to BAM", logger)
			
			this_fix_mate_pair = FixMatePair.new
			out = this_fix_mate_pair.fix_information(this_sample, this_batch, logger)
			this_pipeline.error_check(out, this_sample, "Fix MatePair info", logger)
			
			this_sort_bam = SortBam.new
			out = this_sort_bam.sort(this_sample, this_batch, logger)
			this_pipeline.error_check(out, this_sample, "Sort BAM", logger)
			
			this_mark_dups = MarkDuplicates.new
			out = this_mark_dups.remove_dups(this_sample, this_batch, logger)
			this_pipeline.error_check(out, this_sample, "Remove Duplicates", logger)
			
			this_realign_indels = RealignIndels.new
			out = this_realign_indels.create_targets(this_sample, this_batch, logger)
			this_pipeline.error_check(out, this_sample, "Create targets", logger)
			
			out = this_realign_indels.realign(this_sample, this_batch, logger)
			this_pipeline.error_check(out, this_sample, "Indel Realignment", logger)

	## Metrics overall
      
	#puts this_sample.inspect
			this_metric = CalculateMetrics.new
			out = this_metric.overall_metrics(this_sample, this_batch, logger)
			this_pipeline.error_check(out, this_sample, "Overall metrics", logger)
			
			out = this_metric.phenotype_metrics(this_sample, this_batch, logger)
			this_pipeline.error_check(out, this_sample, "Phenotype metrics", logger)
			
			out = this_metric.phenotype_coverage(this_sample, this_batch, logger)
			this_pipeline.error_check(out, this_sample, "Phenotype coverage", logger)
		
			### Phenotype specifc section
			### Haplotype Caller
        
			this_caller = VariantCaller.new
			out = this_caller.haplotype_caller(this_sample, this_batch, logger)
			
			this_pipeline.error_check(out, this_sample, "Variant Calling", logger)
			
			### Variant filtration
        
			this_selection = FilterVariants.new
			out = this_selection.annotate_filters(this_sample, this_batch, logger)
			this_pipeline.error_check(out, this_sample, "Annotate filters", logger)
			
			## Select Variants - discordance with No Known Clinical Significance 
			
			out = this_selection.select_discordant_variants(this_sample, this_batch, logger, "filtered", "nkmi", this_batch.common_variants_nkmi_path)
			this_pipeline.error_check(out, this_sample, "NKMI variant discordance", logger)
			
			## Select Variants - discordance with Common Artefacts
        
			out = this_selection.select_discordant_variants(this_sample, this_batch, logger, "nkmi", "ca", this_batch.common_artefacts_path)
			this_pipeline.error_check(out, this_sample, "CA variant discordance", logger)
			
			## Select Variants - phenotype specific intervals
			
			out = this_selection.select_phenotype_variants(this_sample, this_batch, logger, "ca", "phenotype")
			this_pipeline.error_check(out, this_sample, "Phenotype specific variant selection", logger)
					
		end
		
		results = Parallel.map(samples_first, :in_processes=>1 ) do |this_sample|
			# Alamut Annotation
			this_annotation = AlamutAnnotation.new
			out = this_annotation.annotate(this_sample, this_batch, logger)
			puts "#{out.inspect}"
			this_pipeline.error_check(out, this_sample, "Alamut annotation", logger)
		end
	
	#Parse batch metrics in order
		this_metric = ParseMetrics.new
		phenotype_metric_line = this_batch.phenotype_metric_line
		overall_metric_line = this_batch.overall_metric_line
		duplicate_metric_line = this_batch.duplicate_metric_line
    
		phenotype_metric_hash = Hash.new
		overall_metric_array = Array.new
		duplicate_array = Array.new
		
		samples.each do |this_sample|
			phenotype_metric_file_path = "#{this_batch.base_path}/#{this_batch.batch_id}/metrics/#{this_sample.panel_version.downcase}_#{this_sample.sample_id}_#{this_sample.gender.upcase}_#{this_sample.phenotype.upcase}.phenotype.bait_capture_metrics"
			overall_metric_file_path = "#{this_batch.base_path}/#{this_batch.batch_id}/metrics/#{this_sample.panel_version.downcase}_#{this_sample.sample_id}_#{this_sample.gender.upcase}.overall.bait_capture_metrics"
			duplicates_file_path = "#{this_batch.base_path}/#{this_batch.batch_id}/duplicates/#{this_sample.panel_version.downcase}_#{this_sample.sample_id}.duplicates"
			
			this_phenotype_metric = this_metric.parse_metrics(phenotype_metric_file_path, phenotype_metric_line)
			if phenotype_metric_hash.has_key?(this_sample.phenotype.upcase)
				this_phenotype_array = phenotype_metric_hash[this_sample.phenotype.upcase]
				this_phenotype_array.push(this_phenotype_metric)
				phenotype_metric_hash[this_sample.phenotype.upcase] = this_phenotype_array
			else
				this_phenotype_array = Array.new
				this_phenotype_array.push(this_phenotype_metric)
				phenotype_metric_hash[this_sample.phenotype.upcase] = this_phenotype_array
			end
			
			this_overall_metric = this_metric.parse_metrics(overall_metric_file_path, overall_metric_line)
			overall_metric_array.push(this_overall_metric)
			
			this_duplicate = this_metric.parse_duplicates(duplicates_file_path, duplicate_metric_line)
			#ap this_duplicate.inspect
			duplicate_array.push(this_duplicate)
		end
    
		#ap phenotype_metric_hash
		
		#write out phenotype metrics
		phenotype_metric_hash.each_pair do |phenotype, this_metric_array|
			csv_data = CSV.generate(col_sep: "\t") do |csv|
				csv << this_metric_array.first.variable_order
				this_metric_array.each do |this_metric|
					csv << this_metric.print_attributes
				end
			end
			# Writing the csv data back to the same file, (also specify UTF-8 format)
			File.open("#{this_batch.base_path}/#{this_batch.batch_id}/metrics/#{this_batch.batch_id}_#{phenotype}.phenotype.metrics", 'w:UTF-8') { |file| file.write(csv_data)}
		end
		
		csv_data = CSV.generate(col_sep: "\t") do |csv|
				csv << overall_metric_array.first.variable_order
				overall_metric_array.each do |this_metric|
					csv << this_metric.print_attributes
				end
		end
		
		File.open("#{this_batch.base_path}/#{this_batch.batch_id}/metrics/#{this_batch.batch_id}.overall.metrics", 'w:UTF-8') { |file| file.write(csv_data)}
		
		csv_data = CSV.generate(col_sep: "\t") do |csv|
				csv << duplicate_array.first.variable_order
				duplicate_array.each do |this_duplicate|
					csv << this_duplicate.print_attributes
				end
		end
			
		File.open("#{this_batch.base_path}/#{this_batch.batch_id}/duplicates/#{this_batch.batch_id}.overall.duplicates", 'w:UTF-8') { |file| file.write(csv_data)}



	#SNP typing
		allowed_genesets = ["v5","v501"]
		input_file_string = ""
		samples.each do |this_sample|
			if allowed_genesets.include?(this_sample.panel_version)
				input_file_string += "-I #{this_batch.base_path}/#{this_batch.batch_id}/assembly/#{this_sample.panel_version}_#{this_sample.sample_id}_#{this_sample.gender.upcase}.realigned.bam "
			end
		end

	#	6q24 SNPs
	#	Only parse samples with the 6q24 region targeted   
  
		output_file_string = "#{this_batch.base_path}/#{this_batch.batch_id}/variants/haplotyper/#{this_batch.batch_id}.6q24.ug_call.vcf"	
		annotated_file_string = "#{this_batch.base_path}/#{this_batch.batch_id}/variants/haplotyper/#{this_batch.batch_id}.filtered.6q24.ug_call.vcf"	
		genotype_file_string = "#{this_batch.base_path}/#{this_batch.batch_id}/variants_6q24/#{this_batch.batch_id}.6q24snps.GATK-#{this_batch.gatk_version}.GT.table"
		allele_depth_file_string = "#{this_batch.base_path}/#{this_batch.batch_id}/variants_6q24/#{this_batch.batch_id}.6q24snps.GATK-#{this_batch.gatk_version}.AD.table"
		
		this_caller = VariantCaller.new
		out = this_caller.genotype_caller_snps(input_file_string, output_file_string, "6q24", this_batch, logger)
			
		this_annotation = FilterVariants.new
		out = this_annotation.annotate_snp_filters(output_file_string, annotated_file_string, this_batch, logger)
		
		this_table = VariantsToTable.new
		out = this_table.genotype_table(annotated_file_string, genotype_file_string, this_batch, logger)
		out = this_table.allelic_depth_table(annotated_file_string, allele_depth_file_string, this_batch, logger)

	#	type_one_snps

		output_file_string = "#{this_batch.base_path}/#{this_batch.batch_id}/variants/haplotyper/#{this_batch.batch_id}.t1d.hc_call.vcf"	
		annotated_file_string = "#{this_batch.base_path}/#{this_batch.batch_id}/variants/haplotyper/#{this_batch.batch_id}.filtered.t1d.hc_call.vcf"	
		genotype_file_string = "#{this_batch.base_path}/#{this_batch.batch_id}/variants_t1d/#{this_batch.batch_id}.t1d.GATK-#{this_batch.gatk_version}.GT.table"
		allele_depth_file_string = "#{this_batch.base_path}/#{this_batch.batch_id}/variants_t1d/#{this_batch.batch_id}.t1d.GATK-#{this_batch.gatk_version}.AD.table"
		
		this_caller = VariantCaller.new
		out = this_caller.haplotype_caller_snps(input_file_string, output_file_string, "t1d", "#{this_batch.type_one_snps_path}", this_batch, logger)
			
		this_annotation = FilterVariants.new
		out = this_annotation.annotate_snp_filters(output_file_string, annotated_file_string, this_batch, logger)
		
		this_table = VariantsToTable.new
		out = this_table.genotype_table(annotated_file_string, genotype_file_string, this_batch, logger)
		out = this_table.allelic_depth_table(annotated_file_string, allele_depth_file_string, this_batch, logger)
			
end
