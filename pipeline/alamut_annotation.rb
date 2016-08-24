# @author Garan Jones
# @note AlamutAnnotation class: use Alamut (version and filepath specified by config.yaml)
class AlamutAnnotation
	require_relative 'wrapper'			
	
		# @author Garan Jones
		# Annotate each sample final VCF file with Alamut
  	# --alltrans; return annotation for all available transcripts for each variant
  	# --ssIntronicRange; set varLocation as 'splice site' if variant is intronic and within this range
  	# --outputVCFInfo; return additional fields directly from VCF input
  	# --outputVCFGenotypeData; return additional fields from VCF input for the sample genotype fields
  	# --outputVCFQuality; return VCF specific quality fields
  	# --outputVCFFilter; return additional VCF filter field
  	# --hgmdUser; taken from Batch object
  	# --hgmdPasswd; taken from Batch object
  	# @param sample [Object] A Sample object
  	# @param batch [Object] A Batch object containing details on the batch
  	# @param logger [Object] A Ruby Logger object
  	# @return [Array<String, Object>] An array with command exit status and an IO.pipe object
	def annotate(sample, batch, logger)

			this_wrapper = Wrapper.new
			puts "Running Alamut annotation on sample...#{sample.sequencing_panel_version}_#{sample.ex_number}"
			
			input_file_string = "#{batch.base_path}/#{batch.batch_id}/variants/haplotyper/#{sample.panel_version}_#{sample.ex_number}_#{sample.gender.upcase}_#{sample.phenotype.upcase}.phenotype.hap_call.vcf"
			annotated_file_string = "#{batch.base_path}/#{batch.batch_id}/variants/haplotyper/#{sample.panel_version}_#{sample.ex_number}_#{sample.gender.upcase}_#{sample.phenotype.upcase}.alamut.alltrans.txt"
			unannotated_file_string = "#{batch.base_path}/#{batch.batch_id}/variants/haplotyper/#{sample.panel_version}_#{sample.ex_number}_#{sample.gender.upcase}_#{sample.phenotype.upcase}.alamut.alltrans.unannotated.txt"
			
			#Build command to send via Wrapper object to the shell to be executed
			cmd = "#{batch.alamut_path} --hgmdUser #{batch.hgmd_user} --hgmdPasswd #{batch.hgmd_pass} --in #{input_file_string} --ann #{annotated_file_string} --unann #{unannotated_file_string} "
			cmd += "--assbly GRCh37 --alltrans --ssIntronicRange 2 --outputVCFInfo AC AF AN DP FS MQ MQ0 QD --outputVCFGenotypeData AD DP GQ GT PL --outputVCFQuality --outputVCFFilter"
			
			#Add info entry to log to say Wrapper object has been sent command
			logger.info('stage') { "Alamut annotation :: Sample #{sample.panel_version}_#{sample.ex_number}" }
			
			#Run command and return IO.pipe object from Wrapper class
			output = this_wrapper.run_command(cmd, logger)
		  
			return output
	end
	

end
