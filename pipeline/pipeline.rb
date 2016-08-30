# @author Garan Jones
# @note Pipeline class: main Pipeline class to run the samples through the stages
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
		require_relative 'cnv_caller'
		require_relative 'file_compressor'
		
	# @author Garan Jones
	# Check the IO.pipe STDOUT returned from Parallel from a non-zero value (error)
	# If non-zero send break signal to Parallel, add error log to logger
  	#
  	# @param out [Array<String, Object>] An array with command exit status and an IO.pipe object
  	# @param sample [Object] A Sample object
  	# @param stage [String] Text description of which stage of the pipeline is being run
  	# @param logger [Object] A Ruby Logger object
  	# @return [Boolean] If IO.pipe STDOUT has passed the error check (TRUE, Parallel::Break not called) or not (FALSE, Parallel::Break called)
	 def error_check(out, sample, stage, logger)
			if out && (out[0] > 0)
				puts "ERROR :: #{stage} :: Halting pipeline at Sample #{sample.panel_version}_#{sample.sample_id}\n"
				ap "#{out[1].inspect}"
				#Write error message out to logger
				logger.info('error') { "SAMPLE :: #{sample.sample_id} #{out[1].inspect}" }
				#Raise Parallel::Break, stop all current processes run by Parallel
				raise Parallel::Break
				return false
			else
				return true
			end
	 end

	# @author Garan Jones
	# Use an array of Sample objects to generate a concatenated list of BAM filepaths
	# based on the allowed genesets
  	#
  	# @param samples [Array<Object>] An array containing the sample objects
  	# @param allowed_genesets [Array<String>] A list of allowed genesets
  	# @return [String] A string containing a concatenated list of BAM filepaths
	 def generate_input_file_string(batch, samples, allowed_genesets)
			#allowed_genesets = ["v5","v501"]
			input_file_string = ""
			samples.each do |this_sample|
				if allowed_genesets.include?(this_sample.panel_version)
					input_file_string += "-I #{batch.base_path}/#{batch.batch_id}/assembly/#{this_sample.panel_version}_#{this_sample.sample_id}_#{this_sample.gender.upcase}.realigned.bam "
				end
			end
			return input_file_string
	 end
	 
	# @author Garan Jones
	# Annotate variants using AlamutAnnotation class
  	#
  	# @param samples [Array<Object>] An array containing the sample objects
  	# @param this_batch [Object] A Batch object containing details on the batch
  	# @param logger [Object] A Logger object
  	# @param this_pipeline [Object] An instance of the Pipeline class
  	# @return [true]
	 def annotate_variants(samples, this_batch, logger, this_pipeline)
			results = Parallel.map(samples, :in_processes=>1 ) do |this_sample|
				# Alamut Annotation
				this_annotation = AlamutAnnotation.new
				out = this_annotation.annotate(this_sample, this_batch, logger)
				puts "#{out.inspect}"
				error_check(out, this_sample, "Alamut annotation", logger)
			end
			return true
	 end

	 def run_assembly(this_sample, this_batch, logger)
	 	 	## Align forward and reverse Fastq files to produce BAM file
 			this_bwa_mem = BwaMem.new
 			out = this_bwa_mem.run_bwa_mem(this_sample, this_batch, logger)
 			error_check(out, this_sample, "BWA-MEM alignment", logger)
 			
 			this_sam_to_bam = SamToBam.new
 			out = this_sam_to_bam.convert_sam_to_bam(this_sample, this_batch, logger)
 			error_check(out, this_sample, "SAM to BAM", logger)
 			
 			this_fix_mate_pair = FixMatePair.new
 			out = this_fix_mate_pair.fix_information(this_sample, this_batch, logger)
 			error_check(out, this_sample, "Fix MatePair info", logger)
 			
 			this_sort_bam = SortBam.new
 			out = this_sort_bam.sort(this_sample, this_batch, logger)
 			error_check(out, this_sample, "Sort BAM", logger)
 			
 			this_mark_dups = MarkDuplicates.new
 			out = this_mark_dups.remove_dups(this_sample, this_batch, logger)
 			error_check(out, this_sample, "Remove Duplicates", logger)
 			
 			this_realign_indels = RealignIndels.new
 			out = this_realign_indels.create_targets(this_sample, this_batch, logger)
 			error_check(out, this_sample, "Create targets", logger)
 			
 			out = this_realign_indels.realign(this_sample, this_batch, logger)
 			error_check(out, this_sample, "Indel Realignment", logger)
 			
	 end
	 
	 def run_metrics(this_sample, this_batch, logger)
	 	 	## Metrics overall
      
 			this_metric = CalculateMetrics.new
 			out = this_metric.overall_metrics(this_sample, this_batch, logger)
 			error_check(out, this_sample, "Overall metrics", logger)
 			
 			out = this_metric.phenotype_metrics(this_sample, this_batch, logger)
 			error_check(out, this_sample, "Phenotype metrics", logger)
 			
 			out = this_metric.phenotype_coverage(this_sample, this_batch, logger)
 			error_check(out, this_sample, "Phenotype coverage", logger)
 			
	 end
	 
	 def run_variant_caller(this_sample, this_batch, logger)		
 			## Phenotype specifc section
 			## Haplotype Caller
       
 			this_caller = VariantCaller.new
 			#bam_out: outputs the HaplotypeCaller representation of the BAM file, should include the ArtificialHaplotypes 
 			#set to false to reduce run-time
 			bam_out = false
 			out = this_caller.haplotype_caller(this_sample, this_batch, logger, bam_out)
 			
 			error_check(out, this_sample, "Variant Calling", logger)
	 end
	 
	 def run_select_variants(this_sample, this_batch, logger) 
			## Variant filtration
       
			this_selection = FilterVariants.new
			out = this_selection.annotate_filters(this_sample, this_batch, logger)
			error_check(out, this_sample, "Annotate filters", logger)
		
			## Select Variants - discordance with No Known Clinical Significance 
		
			out = this_selection.select_discordant_variants(this_sample, this_batch, logger, "filtered", "nkmi", this_batch.common_variants_nkmi_path)
			error_check(out, this_sample, "NKMI variant discordance", logger)
		
			## Select Variants - discordance with Common Artefacts
       
			out = this_selection.select_discordant_variants(this_sample, this_batch, logger, "nkmi", "ca", "#{this_batch.common_artefacts_path}")
			error_check(out, this_sample, "CA variant discordance", logger)
		
			## Select Variants - phenotype specific intervals
		
			out = this_selection.select_phenotype_variants(this_sample, this_batch, logger, "ca", "phenotype")
 			error_check(out, this_sample, "Phenotype specific variant selection", logger)
	 end
	 
	 def run_exome_depth(samples, this_batch, logger)
 	  	["FEMALE", "MALE"].each do |this_gender|
 	  			panels = Array.new
 	  		  samples.each do |this_sample|
 	  		  	panels.push(this_sample.sequencing_panel_version)
 	  			end
 	  			panels.uniq!

	 		  	["v501", "v603"].each do |this_panel_version|
						run_type = "interbatch"
	 		  		this_cnv_caller = CnvCaller.new
	 		  		out = this_cnv_caller.exome_depth(this_panel_version, this_gender, run_type, this_batch, logger)
	 		  		puts out
	 		    end
 	  	end
 	 end
	
	 def run_pipeline(this_pipeline)
	 	 	path = File.expand_path(__FILE__)
			base_path = path.split("/scripts/").first
  		
			logger = Logger.new("#{base_path}/logs/pipeline.log")
			
			#Load the config YAML file and pass the settings to local variables
			this_batch = YAML.load_file("#{base_path}/scripts/configuration/config.yaml")
			
			this_parser = SampleParser.new
			
			samples = this_parser.parse_sample_list("#{base_path}/#{this_batch.sample_list_path}")
			
			samples_first = []
			
 			samples.each do |this_sample|
 				if [ "v603_EX1508144" ].include? this_sample.capture_number
 					samples_first.push this_sample	
 				else
 					#samples_first.push this_sample
 				end
 			end
 	
 	  #Seperate loop for the WGET cmd due to a throttling issue
 		results = Parallel.map(samples, :in_processes=>1 ) do |this_sample|
 			this_wget = Wget.new
 			out = this_wget.fetch_fastq(this_sample, this_batch, logger)
 			this_pipeline.error_check(out, this_sample, "Wget FastQ", logger)
 		
 			this_rename = Rename.new
 			out = this_rename.remove_fastq_adaptor_string(this_sample, this_batch, logger)
 			this_pipeline.error_check(out, this_sample, "FastQ file rename", logger)

 			out = this_rename.rename_symlink(this_sample, this_batch, logger)
 			this_pipeline.error_check(out, this_sample, "Symlink rename", logger)
		end


#Only required if the Fastq files are not gzipped
# 		results = Parallel.map(samples, :in_processes=>20 ) do |this_sample|
#				this_compressor = FileCompressor.new
#				out = this_compressor.gzip_fastq("R1", this_sample, this_batch, logger)
#				if out
#					this_pipeline.error_check(out, this_sample, "Gzip FastQ", logger)
#				else
#					puts "SYMLINK :: #{this_sample.inspect}"
#				end
#
#				out = this_compressor.gzip_fastq("R2", this_sample, this_batch, logger)
#				if out 
#					this_pipeline.error_check(out, this_sample, "Gzip FastQ", logger)
#			  else
#					puts "SYMLINK :: #{this_sample.inspect}"
#				end
# 		end  
 	
#		Main pipeline loop, set to the number of concurrent processes to reflect server load
 		results = Parallel.map(samples, :in_processes=>20 ) do |this_sample|
 			puts this_sample.inspect
 			
 			run_assembly(this_sample, this_batch, logger)
 			
 			run_metrics(this_sample, this_batch, logger)
 			
 			run_variant_caller(this_sample, this_batch, logger)
 			
 			run_select_variants(this_sample, this_batch, logger)
					
 		end
 	
 
	 	#Run ExomeDepth over gender specific batches
	 	run_exome_depth(samples, this_batch, logger)
 	  
	 	#annotate variants
		this_pipeline.annotate_variants(samples, this_batch, logger, this_pipeline)
		
		#double tab in HSmetrics columns is throwing the metrics out of alignment, use sed to remove
 		`sed -i $'s/\t\t/\t/g' #{base_path}/metrics/*`
 		
 		#Parse batch metrics in order
 		this_metric = ParseMetrics.new
 		this_metric.parse_batch_metrics(this_batch, samples)


		#SNP typing

 		input_file_string = this_pipeline.generate_input_file_string(this_batch, samples, ["v5","v501"])
 		if input_file_string != ""
 			this_caller = VariantCaller.new
 		#	6q24 SNPs
 		#	Only parse samples with the 6q24 region targeted   
 			this_caller.call_6q24_snps(this_batch, logger, input_file_string)
 		#	type_one_snps
 			this_caller.call_type_one_snps(this_batch, logger, input_file_string)
 		else
 			puts "No v5 or v501 samples to run through snp typing"
 			logger.info('stage') { "Variant caller - SNP Typing :: No v5 or v501 samples present." }
 		end
	end#run_pipeline method
	
	
	if __FILE__ == $0
		this_pipeline = Pipeline.new
		this_pipeline.run_pipeline(this_pipeline)
	end
end
