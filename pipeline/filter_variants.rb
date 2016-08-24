# @author Garan Jones
# FilterVariants class: Various methods for filtering and adding annotation for filtering VCF files
class FilterVariants
	require_relative 'wrapper'			
	
		# @author Garan Jones
		# Add Filter field to VCF based on various criteria
		# 'PASS' variant is not annotated with any other expression (i.e it passes quality control)
		# 'QD2' variant falls below depth of quality of 2.0 (less than 2 read support this call)
		# 'MQ40' variant has a mapping quality of less than 40
		# 'RPRS-8' ReadPosRankSum; Z-score from Wilcoxon rank sum test of Alt vs. Ref read mapping qualities. If the alternate bases are more likely to be found on reads with lower MQ than reference bases then the site is likely mismapped
		# 'FS60' FisherStrand Bias; Phred-scaled p-value using Fisher's exact test to detect strand bias. If the reference carrying reads are balanced between forward and reverse strands then the alternate carrying reads should be as well
		# 'MQRankSum-12.5' 
  	# @param sample [Object] A Sample object
  	# @param batch [Object] A Batch object containing details on the batch
  	# @param logger [Object] A Ruby Logger object
  	# @return [Array<String, Object>] An array with command exit status and an IO.pipe object
	def annotate_filters(sample, batch, logger)
		
			this_wrapper = Wrapper.new
			puts "Adding VariantFiltration annotation on sample...#{sample.sequencing_panel_version}_#{sample.ex_number}"
			
			input_file_string="#{batch.base_path}/#{batch.batch_id}/variants/haplotyper/#{sample.panel_version.downcase}_#{sample.ex_number}_#{sample.gender.upcase}.hap_call.vcf"
			output_file_string="#{batch.base_path}/#{batch.batch_id}/variants/haplotyper/#{sample.panel_version.downcase}_#{sample.ex_number}_#{sample.gender.upcase}.filtered.hap_call.vcf"
			this_panel = batch.select_panel(sample.panel_version)
			
			#Build the command to be passed to the Wrapper object
			cmd = "#{batch.java_path} -Xmx4g -jar #{batch.gatk_path} -T VariantFiltration -R #{batch.reference_path} --filterExpression 'QD < 2.0' --filterName 'QD2' --filterExpression 'MQ < 40.0' --filterName 'MQ40' "
			cmd += "--filterExpression 'ReadPosRankSum < -8.0' --filterName 'RPRS-8' --filterExpression 'FS > 60.0' --filterName 'FS60' --filterExpression 'MQRankSum < -12.5' --filterName 'MQRankSum-12.5' "
			cmd += "-o #{output_file_string} --variant #{input_file_string} -log #{batch.base_path}/#{batch.batch_id}/logs/#{sample.panel_version}_#{sample.ex_number}_filtered_hap_call.log"
			
			#Add an entry to the Logger object, stating that the filters are being annotated on the sample vcf
			logger.info('stage') { "Annotate Filters :: Sample #{sample.panel_version}_#{sample.ex_number}" }
			
			#Send the command to the Wrapper object, return the IO.pipe output
			output = this_wrapper.run_command(cmd, logger)
		  
			return output
	end
	
		# @author Garan Jones
		# Add Filter field to SNP-typing specific VCF based on various criteria
		# 'PASS' variant is not annotated with any other expression (i.e it passes quality control)
		# 'QD2' variant falls below depth of quality of 2.0 (less than 2 read support this call)
		# 'MQ40' variant has a mapping quality of less than 40
		# 'RPRS-8' ReadPosRankSum; Z-score from Wilcoxon rank sum test of Alt vs. Ref read mapping qualities. If the alternate bases are more likely to be found on reads with lower MQ than reference bases then the site is likely mismapped
		# 'FS60' FisherStrand Bias; Phred-scaled p-value using Fisher's exact test to detect strand bias. If the reference carrying reads are balanced between forward and reverse strands then the alternate carrying reads should be as well
		# 'MQRankSum-12.5' 
  	# @param input_file_string [String] Filepath to input file
  	# @param output_file_string [String] Filepath to output file
  	# @param batch [Object] A Batch object containing details on the batch
  	# @param logger [Object] A Ruby Logger object
  	# @return [Array<String, Object>] An array with command exit status and an IO.pipe object
	def annotate_snp_filters(input_file_string, output_file_string, batch, logger)
		
			this_wrapper = Wrapper.new
			puts "Adding VariantFiltration annotation on sample...#{batch.batch_id}"
						
			#Build the command to be passed to the Wrapper object
			cmd = "#{batch.java_path} -Xmx4g -jar #{batch.gatk_path} -T VariantFiltration -R #{batch.reference_path} --filterExpression 'QD < 2.0' --filterName 'QD2' --filterExpression 'MQ < 40.0' --filterName 'MQ40' "
			cmd += "--filterExpression 'ReadPosRankSum < -8.0' --filterName 'RPRS-8' --filterExpression 'FS > 60.0' --filterName 'FS60' --filterExpression 'MQRankSum < -12.5' --filterName 'MQRankSum-12.5' "
			cmd += "-o #{output_file_string} --variant #{input_file_string} -log #{batch.base_path}/#{batch.batch_id}/logs/#{batch.batch_id}_filtered_6q24_ug_call.log"
			
			#Add an entry to the Logger object, stating that the filters are being annotated on the sample vcf
			logger.info('stage') { "Annotate UG Batch Filters :: #{batch.batch_id}" }
			
			#Send the command to the Wrapper object, return the IO.pipe output
			output = this_wrapper.run_command(cmd, logger)
		  
			return output
	end
	
		# @author Garan Jones
		# Select discordant variants from given VCF file of variants to remove (common artefacts / non known clinical significance)
		# @param sample [Object] A Sample oject
  	# @param batch [Object] A Batch object containing details on the batch
  	# @param logger [Object] A Ruby Logger object
  	# @param prefix_in [String] Prefix of input filename, reflects variants that are being selected against (ca / nkmi) at previous step
  	# @param prefix_out [String] Prefix of output file, reflects the variants that are being selected against (ca / nkmi) at this step
  	# @param discordant_file_path [String] Filepath of VCF file that contains the variants to be removed from sample VCF file
  	# @return [Array<String, Object>] An array with command exit status and an IO.pipe object
	def select_discordant_variants(sample, batch, logger, prefix_in, prefix_out, discordant_file_path)
				
			this_wrapper = Wrapper.new
			puts "Running SelectVariants for #{prefix_out} on sample...#{sample.sequencing_panel_version}_#{sample.ex_number}"
			
			input_file_string="#{batch.base_path}/#{batch.batch_id}/variants/haplotyper/#{sample.panel_version.downcase}_#{sample.ex_number}_#{sample.gender.upcase}.#{prefix_in}.hap_call.vcf"
			output_file_string="#{batch.base_path}/#{batch.batch_id}/variants/haplotyper/#{sample.panel_version.downcase}_#{sample.ex_number}_#{sample.gender.upcase}.#{prefix_out}.hap_call.vcf"
			
			#Build the command to be passed to the Wrapper object	
			cmd = "#{batch.java_path} -Xmx4g -jar #{batch.gatk_path} -T SelectVariants -R #{batch.reference_path} "
			cmd += "-o #{output_file_string} --variant #{input_file_string} --discordance #{discordant_file_path} -U LENIENT_VCF_PROCESSING"
			
			#Add an entry to the Logger object
			logger.info('stage') { "Select Variants #{prefix_out} :: Sample #{sample.panel_version}_#{sample.ex_number}" }
			
			#Send the command to the Wrapper object, return the IO.pipe output
			output = this_wrapper.run_command(cmd, logger)
		  
			return output
		
	end
	
		# @author Garan Jones
		# Select variants from sample VCF based on supplied phenotype intervals
		# @param sample [Object] A Sample oject
  	# @param batch [Object] A Batch object containing details on the batch
  	# @param logger [Object] A Ruby Logger object
  	# @param prefix_in [String] Prefix of input filename
  	# @param prefix_out [String] Prefix of output file
  	# @return [Array<String, Object>] An array with command exit status and an IO.pipe object
	def select_phenotype_variants(sample, batch, logger, prefix_in, prefix_out)
				
			this_wrapper = Wrapper.new
			puts "Running SelectVariants for phenotype specific regions of interest #{sample.phenotype} on sample...#{sample.sequencing_panel_version}_#{sample.ex_number}"
			
			input_file_string="#{batch.base_path}/#{batch.batch_id}/variants/haplotyper/#{sample.panel_version.downcase}_#{sample.ex_number}_#{sample.gender.upcase}.#{prefix_in}.hap_call.vcf"
			output_file_string="#{batch.base_path}/#{batch.batch_id}/variants/haplotyper/#{sample.panel_version.downcase}_#{sample.ex_number}_#{sample.gender.upcase}_#{sample.phenotype.upcase}.#{prefix_out}.hap_call.vcf"
			
			this_panel = batch.select_panel(sample.panel_version)
			phenotype_intervals_path="#{this_panel.intervals_directory}/#{sample.panel_version.downcase}_#{sample.phenotype}_variant_calling.bed"
			
			#Build the command to be passed to the Wrapper object	
			cmd = "#{batch.java_path} -Xmx4g -jar #{batch.gatk_path} -T SelectVariants -R #{batch.reference_path} "
			cmd += "-o #{output_file_string} --variant #{input_file_string} -L #{phenotype_intervals_path} -log #{batch.base_path}/#{batch.batch_id}/logs/#{sample.panel_version}_#{sample.ex_number}_#{sample.phenotype}.log"
			
			#Add an entry to the Logger object
			logger.info('stage') { "Select Variants #{sample.phenotype} :: Sample #{sample.panel_version}_#{sample.ex_number}" }
			
			#Send the command to the Wrapper object, return the IO.pipe output
			output = this_wrapper.run_command(cmd, logger)
		  
			return output
		
	end
	
	
	

end
