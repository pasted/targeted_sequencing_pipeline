class FilterVariants
	require_relative 'wrapper'			
			
	def annotate_filters(sample, batch, logger)
		
			this_wrapper = Wrapper.new
			puts "Adding VariantFiltration annotation on sample...#{sample.sequencing_panel_version}_#{sample.sample_id}"
			
			input_file_string="#{batch.base_path}/#{batch.batch_id}/variants/haplotyper/#{sample.panel_version.downcase}_#{sample.sample_id}_#{sample.gender.upcase}_#{sample.phenotype.upcase}.hap_call.vcf"
			output_file_string="#{batch.base_path}/#{batch.batch_id}/variants/haplotyper/#{sample.panel_version.downcase}_#{sample.sample_id}_#{sample.gender.upcase}_#{sample.phenotype.upcase}.filtered.hap_call.vcf"
			this_panel = batch.select_panel(sample.panel_version)
			
			
			cmd = "#{batch.java_path} -Xmx4g -jar #{batch.gatk_path} -T VariantFiltration -R #{batch.reference_path} --filterExpression 'QD < 2.0' --filterName 'QD2' --filterExpression 'MQ < 40.0' --filterName 'MQ40' "
			cmd += "--filterExpression 'ReadPosRankSum < -8.0' --filterName 'RPRS-8' --filterExpression 'FS > 60.0' --filterName 'FS60' --filterExpression 'MQRankSum < -12.5' --filterName 'MQRankSum-12.5' "
			cmd += "-o #{output_file_string} --variant #{input_file_string} -log #{batch.base_path}/logs/#{sample.panel_version}_#{sample.sample_id}_filtered_hap_call.log"
			
			logger.info('stage') { "Annotate Filters :: Sample #{sample.panel_version}_#{sample.sample_id}" }
			output = this_wrapper.run_command(cmd, logger)
		  
			return output
	end
	
	def annotate_snp_filters(input_file_string, output_file_string, batch, logger)
		
			this_wrapper = Wrapper.new
			puts "Adding VariantFiltration annotation on sample...#{batch.batch_id}"
						
			
			cmd = "#{batch.java_path} -Xmx4g -jar #{batch.gatk_path} -T VariantFiltration -R #{batch.reference_path} --filterExpression 'QD < 2.0' --filterName 'QD2' --filterExpression 'MQ < 40.0' --filterName 'MQ40' "
			cmd += "--filterExpression 'ReadPosRankSum < -8.0' --filterName 'RPRS-8' --filterExpression 'FS > 60.0' --filterName 'FS60' --filterExpression 'MQRankSum < -12.5' --filterName 'MQRankSum-12.5' "
			cmd += "-o #{output_file_string} --variant #{input_file_string} -log #{batch.base_path}/logs/#{batch.batch_id}_filtered_6q24_ug_call.log"
			
			logger.info('stage') { "Annotate UG Batch Filters :: #{batch.batch_id}" }
			output = this_wrapper.run_command(cmd, logger)
		  
			return output
	end
	
	def select_discordant_variants(sample, batch, logger, prefix_in, prefix_out, discordant_file_path)
				
			this_wrapper = Wrapper.new
			puts "Running SelectVariants for #{prefix_out} on sample...#{sample.sequencing_panel_version}_#{sample.sample_id}"
			
			input_file_string="#{batch.base_path}/#{batch.batch_id}/variants/haplotyper/#{sample.panel_version.downcase}_#{sample.sample_id}_#{sample.gender.upcase}_#{sample.phenotype.upcase}.#{prefix_in}.hap_call.vcf"
			output_file_string="#{batch.base_path}/#{batch.batch_id}/variants/haplotyper/#{sample.panel_version.downcase}_#{sample.sample_id}_#{sample.gender.upcase}_#{sample.phenotype.upcase}.#{prefix_out}.hap_call.vcf"
			
						
			cmd = "#{batch.java_path} -Xmx4g -jar #{batch.gatk_path} -T SelectVariants -R #{batch.reference_path} "
			cmd += "-o #{output_file_string} --variant #{input_file_string} --discordance #{discordant_file_path} -U LENIENT_VCF_PROCESSING"
			
			logger.info('stage') { "Select Variants #{prefix_out} :: Sample #{sample.panel_version}_#{sample.sample_id}" }
			output = this_wrapper.run_command(cmd, logger)
		  
			return output
		
	end
	
	def select_phenotype_variants(sample, batch, logger, prefix_in, prefix_out)
				
			this_wrapper = Wrapper.new
			puts "Running SelectVariants for phenotype specific regions of interest #{sample.phenotype} on sample...#{sample.sequencing_panel_version}_#{sample.sample_id}"
			
			input_file_string="#{batch.base_path}/#{batch.batch_id}/variants/haplotyper/#{sample.panel_version.downcase}_#{sample.sample_id}_#{sample.gender.upcase}_#{sample.phenotype.upcase}.#{prefix_in}.hap_call.vcf"
			output_file_string="#{batch.base_path}/#{batch.batch_id}/variants/haplotyper/#{sample.panel_version.downcase}_#{sample.sample_id}_#{sample.gender.upcase}_#{sample.phenotype.upcase}.#{prefix_out}.hap_call.vcf"
			
			this_panel = batch.select_panel(sample.panel_version)
			phenotype_intervals_path="#{this_panel.intervals_directory}/#{sample.panel_version.downcase}_#{sample.phenotype}_variant_calling.bed"
						
			cmd = "#{batch.java_path} -Xmx4g -jar #{batch.gatk_path} -T SelectVariants -R #{batch.reference_path} "
			cmd += "-o #{output_file_string} --variant #{input_file_string} -L #{phenotype_intervals_path} -log #{batch.base_path}/logs/#{sample.panel_version}_#{sample.sample_id}_#{sample.phenotype}.log"
			
			logger.info('stage') { "Select Variants #{sample.phenotype} :: Sample #{sample.panel_version}_#{sample.sample_id}" }
			output = this_wrapper.run_command(cmd, logger)
		  
			return output
		
	end
	
	
	

end
