class VariantCaller
	require_relative 'wrapper'
	require_relative 'filter_variants'
	require_relative 'variants_to_table'
	
	# Runs GATK HaplotypeCaller (version specified by the config file)
	#	@param sample [Sample] The Sample instance
	# @param batch [Batch] The Batch instance, derived from the Config YAML file
	# @param logger [Logger] The Logger instance
	#	@param bam_out_required [Boolean] A Boolean to select whether or not to output the BAM file
	# @return output [Array] An array containing the ExitStatus of the Wrapper class and it's STDOUT	
	def haplotype_caller(sample, batch, logger, bam_out_required)
		
			this_wrapper = Wrapper.new
			puts "Running GATK HaplotypeCaller on sample...#{sample.sequencing_panel_version.downcase}_#{sample.sample_id}"
			
			input_file_string = "../../assembly/#{sample.panel_version.downcase}_#{sample.sample_id}_#{sample.gender.upcase}.realigned.bam"
			
			output_file_string = "../../variants/haplotyper/#{sample.panel_version.downcase}_#{sample.sample_id}_#{sample.gender.upcase}_#{sample.phenotype.upcase}.hap_call.vcf"
			this_panel = batch.select_panel(sample.panel_version)
			
			phenotype_intervals_path="#{this_panel.intervals_directory}/#{sample.panel_version.downcase}_#{sample.phenotype}_variant_calling.bed"
			
			cmd = "#{batch.java_path} -Xmx4g -jar #{batch.gatk_path} -T HaplotypeCaller -R #{batch.reference_path} -rf BadCigar -stand_call_conf 50.0 -stand_emit_conf 30.0 "
			
			cmd += "-L #{phenotype_intervals_path} -I #{input_file_string} -o #{output_file_string} -log ../../logs/#{sample.panel_version}_#{sample.sample_id}_haplotype_caller.log"

			if bam_out_required
				#Output BAM file
				output_bam_string = "../../assembly/#{sample.panel_version}_#{sample.sample_id}_#{sample.gender.upcase}_#{sample.phenotype.upcase}.haplotype_caller.bam"
				#OPTIONAL Output BAM file to check informative reads
				cmd += "--bamOutput #{output_bam_string} "
			end
			
			logger.info('stage') { "Variant caller :: Sample #{sample.panel_version.downcase}_#{sample.sample_id}" }
			output = this_wrapper.run_command(cmd, logger)

			return output
	end
	
	# Runs GATK HaplotypeCaller over the T1D snps, converts the genotypes and allele depths to table format
	# @param this_batch [Batch] The Batch instance, derived from the Config YAML file
	# @param logger [Logger] The Logger instance
	#	@param bam_out_required [Boolean] A Boolean to select whether or not to output the BAM file
	# @return output [Array] An array containing the ExitStatus of the Wrapper class and it's STDOUT	
	def call_type_one_snps(this_batch, logger, input_file_string)
	
		output_file_string = "../../variants/haplotyper/#{this_batch.batch_id}.t1d.hc_call.vcf"	
		annotated_file_string = "../../variants/haplotyper/#{this_batch.batch_id}.filtered.t1d.hc_call.vcf"	
		genotype_file_string = "../../variants_t1d/#{this_batch.batch_id}.t1d.GATK-#{this_batch.gatk_version}.GT.table"
		allele_depth_file_string = "../../variants_t1d/#{this_batch.batch_id}.t1d.GATK-#{this_batch.gatk_version}.AD.table"
		
		this_caller = VariantCaller.new
		out = this_caller.haplotype_caller_snps(input_file_string, output_file_string, "t1d", "#{this_batch.type_one_snps_path}", this_batch, logger)
			
		this_annotation = FilterVariants.new
		out = this_annotation.annotate_snp_filters(output_file_string, annotated_file_string, this_batch, logger)
		
		this_table = VariantsToTable.new
		out = this_table.genotype_table(annotated_file_string, genotype_file_string, this_batch, logger)
		out = this_table.allelic_depth_table(annotated_file_string, allele_depth_file_string, this_batch, logger)
		
		return true
	end
	
	# Runs GATK HaplotypeCaller over the 6q24 snps, converts the genotypes and allele depths to table format
	# @param this_batch [Batch] The Batch instance, derived from the Config YAML file
	# @param logger [Logger] The Logger instance
	#	@param bam_out_required [Boolean] A Boolean to select whether or not to output the BAM file
	# @return output [Array] An array containing the ExitStatus of the Wrapper class and it's STDOUT	
	def call_6q24_snps(this_batch, input_file_string) 
		
		output_file_string = "../../variants/haplotyper/#{this_batch.batch_id}.6q24.ug_call.vcf"	
		annotated_file_string = "../../variants/haplotyper/#{this_batch.batch_id}.filtered.6q24.ug_call.vcf"	
		genotype_file_string = "../../variants_6q24/#{this_batch.batch_id}.6q24snps.GATK-#{this_batch.gatk_version}.GT.table"
		allele_depth_file_string = "../../variants_6q24/#{this_batch.batch_id}.6q24snps.GATK-#{this_batch.gatk_version}.AD.table"
		
		this_caller = VariantCaller.new
		out = this_caller.genotype_caller_snps(input_file_string, output_file_string, "6q24", this_batch, logger)
			
		this_annotation = FilterVariants.new
		out = this_annotation.annotate_snp_filters(output_file_string, annotated_file_string, this_batch, logger)
		
		this_table = VariantsToTable.new
		out = this_table.genotype_table(annotated_file_string, genotype_file_string, this_batch, logger)
		out = this_table.allelic_depth_table(annotated_file_string, allele_depth_file_string, this_batch, logger)
		
		return true
	end
	
	def platypus(sample, batch, logger)
		
			this_wrapper = Wrapper.new
			puts "Running Platypus on sample...#{sample.sequencing_panel_version.downcase}_#{sample.sample_id}"
			
			input_file_string = "../../assembly/assembly/#{sample.panel_version.downcase}_#{sample.sample_id}_#{sample.gender.upcase}.realigned.bam"
			
			output_file_string = "../../assembly/variants/platypus/#{sample.panel_version.downcase}_#{sample.sample_id}_#{sample.gender.upcase}_#{sample.phenotype.upcase}.platypus.vcf"
			this_panel = batch.select_panel(sample.panel_version)
			
			cmd = "python /usr/share/platypus/Platypus_0.8.1/Platypus.py callVariants --bamFiles=#{input_file_string} --refFile=#{batch.reference_path} --output=#{output_file_string}"
	

			logger.info('stage') { "Variant caller Platypus :: Sample #{sample.panel_version.downcase}_#{sample.sample_id}" }
			output = this_wrapper.run_command(cmd, logger)

			return output
	end
	
	#DO NOT USE until new version of Freebayes downloaded and checked
	#Current version uses too much system resources
	def freebayes(sample, batch, logger)
		
			this_wrapper = Wrapper.new
			puts "Running Freebayes on sample...#{sample.sequencing_panel_version.downcase}_#{sample.sample_id}"
			
			input_file_string = "../../assembly/assembly/#{sample.panel_version.downcase}_#{sample.sample_id}_#{sample.gender.upcase}.realigned.bam"
			
			output_file_string = "../../assembly/variants/freebayes/#{sample.panel_version.downcase}_#{sample.sample_id}_#{sample.gender.upcase}_#{sample.phenotype.upcase}.freebayes.vcf"
			this_panel = batch.select_panel(sample.panel_version)
			
						
			cmd = "freebayes -f #{batch.reference_path} #{input_file_string} > #{output_file_string}"
			#cmd += ""

			logger.info('stage') { "Variant caller Freebayes :: Sample #{sample.panel_version.downcase}_#{sample.sample_id}" }
			output = this_wrapper.run_command(cmd, logger)

			return output
	end
	
	def genotype_caller_snps(input_file_path, output_file_path, snp_type, batch, logger)
		
			this_wrapper = Wrapper.new
			output = ""
			
				puts "Running GATK Unified Genotyper #{snp_type} on batch."
				
				cmd = "#{batch.java_path} -Xmx4g -jar #{batch.gatk_path} -T UnifiedGenotyper -R #{batch.reference_path} -D #{batch.dbsnp_path} -L #{batch.ndm_snps_path} -log #{batch.base_path}/logs/#{batch.batch_id}_#{snp_type}_ug.log "
				cmd += "#{input_file_path} -o #{output_file_path} --output_mode EMIT_ALL_SITES "
      	
				logger.info('stage') { "Variant caller UG snps :: Batch #{batch.batch_id}" }
				output = this_wrapper.run_command(cmd, logger)
			
			return output
	end
	
	def haplotype_caller_snps(input_file_path, output_file_path, snp_type, snps_path, batch, logger)
		
			this_wrapper = Wrapper.new
			output = ""
			
				puts "Running GATK Haplotype Caller #{snp_type} on batch."
				
				cmd = "#{batch.java_path} -Xmx4g -jar #{batch.gatk_path} -T HaplotypeCaller -R #{batch.reference_path} -D #{batch.dbsnp_path} -L #{snps_path} -log #{batch.base_path}/logs/#{batch.batch_id}_#{snp_type}_hc.log "
				cmd += "#{input_file_path} -o #{output_file_path} --output_mode EMIT_ALL_SITES "
      	
				logger.info('stage') { "Variant caller HC snps :: Batch #{batch.batch_id}" }
				output = this_wrapper.run_command(cmd, logger)
			
			return output
	end
	
	

end
