# @author Garan Jones
# VerifyBamId class: Container class for various variant calling software calls.
class VerifyBamId
	require_relative 'wrapper'
	
		# Runs VerifyBamID 
		#	@param sample [Object] The Sample instance
		# @param batch [Object] The Batch instance, derived from the Config YAML file
		# @param logger [Object] The Logger instance
		#	@param bam_out_required [Boolean] A Boolean to select whether or not to output the BAM file
		# @return [Array<String, Object>] An array with command exit status and an IO.pipe object
	def verifybam_check(sample, batch, logger)
		
			this_wrapper = Wrapper.new
			puts "Running VerifyBamID on sample...#{sample.sequencing_panel_version.downcase}_#{sample.ex_number}"
			
			input_file_string = "#{batch.base_path}/#{batch.batch_id}/assembly/#{sample.panel_version.downcase}_#{sample.ex_number}_#{sample.gender.upcase}.realigned.bam"
			id_label = "#{sample.panel_version.downcase}_#{sample.ex_number}"
			
			cmd = "verifyBamID --vcf #{batch.cleancall_resources_path}/#{sample.panel_version.downcase}_exac_subset.vcf --bam #{input_file_string} --out #{batch.base_path}/#{batch.batch_id}/clean_call/verify_bam_id.#{id_label}"			
						
			puts cmd
			
			logger.info('stage') { "VerifyBamID :: Sample #{sample.panel_version.downcase}_#{sample.ex_number}" }
			output = this_wrapper.run_command(cmd, logger)

			return output
	end
	
end
