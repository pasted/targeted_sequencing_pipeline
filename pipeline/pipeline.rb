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
				return false
			else
				return true
			end
		end
			
	 def generate_input_file_string(samples, allowed_genesets)
			#allowed_genesets = ["v5","v501"]
			input_file_string = ""
			samples.each do |this_sample|
				if allowed_genesets.include?(this_sample.panel_version)
					input_file_string += "-I ../../assembly/#{this_sample.panel_version}_#{this_sample.sample_id}_#{this_sample.gender.upcase}.realigned.bam "
				end
			end
			return input_file_string
	 end
	 
	 def annotate_variants(samples, this_batch, logger, this_pipeline)
			results = Parallel.map(samples, :in_processes=>1 ) do |this_sample|
				# Alamut Annotation
				this_annotation = AlamutAnnotation.new
				out = this_annotation.annotate(this_sample, this_batch, logger)
				puts "#{out.inspect}"
				this_pipeline.error_check(out, this_sample, "Alamut annotation", logger)
			end
			return true
	end

	
	
		logger = Logger.new('../../logs/pipeline.log')
		
		#Load the config YAML file and pass the settings to local variables
		this_batch = YAML.load_file('../configuration/config.yaml')
		
		this_parser = SampleParser.new
		this_pipeline = Pipeline.new
		
		samples = this_parser.parse_sample_list(this_batch.sample_list_path)
		
		samples_first = []
		
		samples.each do |this_sample|
			if [ "v5_0510", "v5_0511", "v5_0552", "v5_0557", "v5_0558", "v5_0562", "v5_0564", "v5_0569", "v5_0572", "v5_0574", "v5_0580", "v5_0581", "v5_0583", "v5_0584", "v5_0599", "v5_0605", "v5_0609", "v5_0623", "v5_0624", "v5_0632", "v5_0753", "v5_0760", "v5_0762", "v5_0786", "v5_0788", "v5_0791", "v5_0797", "v5_0799", "v5_0800", "v5_0801", "v5_0802", "v5_0822", "v5_0827", "v5_0828", "v5_0829", "v5_0832", "v5_0834", "v5_0839", "v5_0846", "v5_0854", "v5_0863", "v5_0883", "v5_0884", "v5_0888", "v5_0894", "v5_0896", "v5_0897", "v5_0929", "v5_0946", "v5_0979", "v501_0952", "v501_0971", "v501_1004", "v501_1006", "v501_1008", "v501_1012", "v501_1014", "v501_1015", "v501_1016", "v501_1018" ].include? this_sample.capture_number
				#samples_first.push this_sample	
			else
				samples_first.push this_sample
			end
		end
		
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
		
#		results = Parallel.map(samples_first, :in_processes=>20 ) do |this_sample|
#			puts this_sample.inspect
#			this_bwa_mem = BwaMem.new
#			out = this_bwa_mem.run_bwa_mem(this_sample, this_batch, logger)
#			this_pipeline.error_check(out, this_sample, "BWA-MEM alignment", logger)
#			
#			this_sam_to_bam = SamToBam.new
#			out = this_sam_to_bam.convert_sam_to_bam(this_sample, this_batch, logger)
#			this_pipeline.error_check(out, this_sample, "SAM to BAM", logger)
#			
#			this_fix_mate_pair = FixMatePair.new
#			out = this_fix_mate_pair.fix_information(this_sample, this_batch, logger)
#			this_pipeline.error_check(out, this_sample, "Fix MatePair info", logger)
#			
#			this_sort_bam = SortBam.new
#			out = this_sort_bam.sort(this_sample, this_batch, logger)
#			this_pipeline.error_check(out, this_sample, "Sort BAM", logger)
#			
#			this_mark_dups = MarkDuplicates.new
#			out = this_mark_dups.remove_dups(this_sample, this_batch, logger)
#			this_pipeline.error_check(out, this_sample, "Remove Duplicates", logger)
#			
#			this_realign_indels = RealignIndels.new
#			out = this_realign_indels.create_targets(this_sample, this_batch, logger)
#			this_pipeline.error_check(out, this_sample, "Create targets", logger)
#			
#			out = this_realign_indels.realign(this_sample, this_batch, logger)
#			this_pipeline.error_check(out, this_sample, "Indel Realignment", logger)
#
#			## Metrics overall
#     
#			this_metric = CalculateMetrics.new
#			out = this_metric.overall_metrics(this_sample, this_batch, logger)
#			this_pipeline.error_check(out, this_sample, "Overall metrics", logger)
#			
#			out = this_metric.phenotype_metrics(this_sample, this_batch, logger)
#			this_pipeline.error_check(out, this_sample, "Phenotype metrics", logger)
#			
#			out = this_metric.phenotype_coverage(this_sample, this_batch, logger)
#			this_pipeline.error_check(out, this_sample, "Phenotype coverage", logger)
#		
#			### Phenotype specifc section
#			### Haplotype Caller
#      
#			this_caller = VariantCaller.new
#			out = this_caller.haplotype_caller(this_sample, this_batch, logger, false)
#			
#			this_pipeline.error_check(out, this_sample, "Variant Calling", logger)
#
#			### Variant filtration
#       
#			this_selection = FilterVariants.new
#			out = this_selection.annotate_filters(this_sample, this_batch, logger)
#			this_pipeline.error_check(out, this_sample, "Annotate filters", logger)
#		
#			## Select Variants - discordance with No Known Clinical Significance 
#		
#			out = this_selection.select_discordant_variants(this_sample, this_batch, logger, "filtered", "nkmi", this_batch.common_variants_nkmi_path)
#			this_pipeline.error_check(out, this_sample, "NKMI variant discordance", logger)
#		
#			## Select Variants - discordance with Common Artefacts
#       
#			out = this_selection.select_discordant_variants(this_sample, this_batch, logger, "nkmi", "ca", "#{this_batch.common_artefacts_path}")
#			this_pipeline.error_check(out, this_sample, "CA variant discordance", logger)
#		
#			## Select Variants - phenotype specific intervals
#		
#			out = this_selection.select_phenotype_variants(this_sample, this_batch, logger, "ca", "phenotype")
#			this_pipeline.error_check(out, this_sample, "Phenotype specific variant selection", logger)
#					
#		end
		

		#annotate variants
		this_pipeline.annotate_variants(samples, this_batch, logger, this_pipeline)
		#Parse batch metrics in order
		this_metric = ParseMetrics.new
		this_metric.parse_batch_metrics(this_batch, samples)

#		#SNP typing
#
#		input_file_string = this_pipeline.generate_input_file_string(samples, ["v5","v501"])
#		if input_file_string != ""
#		#	6q24 SNPs
#		#	Only parse samples with the 6q24 region targeted   
#			this_pipeline.call_6q24_snps(this_batch, logger, input_file_string)
#		#	type_one_snps
#			this_pipeline.call_type_one_snps(this_batch, logger, input_file_string)
#		else
#			puts "No v5 or v501 samples to run through snp typing"
#			logger.info('stage') { "Variant caller - SNP Typing :: No v5 or v501 samples present." }
#		end
		
end
