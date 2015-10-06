class VariantCaller
	require_relative 'wrapper'			
			
	def haplotype_caller(sample, batch, logger)
		
			this_wrapper = Wrapper.new
			puts "Running GATK HaplotypeCaller on sample...#{sample.sequencing_panel_version.downcase}_#{sample.sample_id}"
			
			input_file_string = "#{batch.base_path}/#{batch.batch_id}/assembly/#{sample.panel_version.downcase}_#{sample.sample_id}_#{sample.gender.upcase}.realigned.bam"
			#Output BAM file
			#output_bam_string = "#{batch.base_path}/#{batch.batch_id}/assembly/#{sample.panel_version}_#{sample.sample_id}_#{sample.gender.upcase}_#{sample.phenotype.upcase}.haplotype_caller.bam"
			
			output_file_string = "#{batch.base_path}/#{batch.batch_id}/variants/haplotyper/#{sample.panel_version.downcase}_#{sample.sample_id}_#{sample.gender.upcase}_#{sample.phenotype.upcase}.hap_call.vcf"
			this_panel = batch.select_panel(sample.panel_version)
			
			phenotype_intervals_path="#{this_panel.intervals_directory}/#{sample.panel_version.downcase}_#{sample.phenotype}_variant_calling.bed"
			
			cmd = "#{batch.java_path} -Xmx4g -jar #{batch.gatk_path} -T HaplotypeCaller -R #{batch.reference_path} -rf BadCigar -stand_call_conf 50.0 -stand_emit_conf 30.0 "
			#OPTIONAL Output BAM file to check informative reads
			#cmd += "--bamOutput #{output_bam_string} "
			cmd += "-L #{phenotype_intervals_path} -I #{input_file_string} -o #{output_file_string} -log #{batch.base_path}/logs/#{sample.panel_version}_#{sample.sample_id}_haplotype_caller.log"

			logger.info('stage') { "Variant caller :: Sample #{sample.panel_version.downcase}_#{sample.sample_id}" }
			output = this_wrapper.run_command(cmd, logger)

			return output
	end
	
 def platypus(sample, batch, logger)
		
			this_wrapper = Wrapper.new
			puts "Running Platypus on sample...#{sample.sequencing_panel_version.downcase}_#{sample.sample_id}"
			
			input_file_string = "#{batch.base_path}/#{batch.batch_id}/assembly/#{sample.panel_version.downcase}_#{sample.sample_id}_#{sample.gender.upcase}.realigned.bam"
			
			output_file_string = "#{batch.base_path}/#{batch.batch_id}/variants/platypus/#{sample.panel_version.downcase}_#{sample.sample_id}_#{sample.gender.upcase}_#{sample.phenotype.upcase}.platypus.vcf"
			this_panel = batch.select_panel(sample.panel_version)
			
			phenotype_intervals_path="#{this_panel.intervals_directory}/#{sample.panel_version.downcase}_#{sample.phenotype}_variant_calling.bed"
			
			cmd = "python /usr/share/platypus/Platypus_0.8.1/Platypus.py callVariants --bamFiles=#{input_file_string} --refFile=#{batch.reference_path} --output=#{output_file_string}"
			#cmd += ""

			logger.info('stage') { "Variant caller Platypus :: Sample #{sample.panel_version.downcase}_#{sample.sample_id}" }
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
