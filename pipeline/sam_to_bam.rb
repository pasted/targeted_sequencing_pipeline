# @author Garan Jones
# SamToBam class: use Picard to convert the SAM input to a BAM output
class SamToBam
	require_relative 'wrapper'			
		
		# @author Garan Jones
		# Use Picard to convert the SAM input to a BAM output
  	# @param sample [Object] A Sample object
  	# @param batch [Object] A Batch object containing details on the batch
  	# @param logger [Object] A Ruby Logger object
  	# @return [Array<String, Object>] An array with command exit status and an IO.pipe object
	def convert_sam_to_bam(sample, batch, logger)
		
			this_wrapper = Wrapper.new
			puts "Converting SAM to BAM...#{sample.sequencing_panel_version}_#{sample.ex_number}"

			#Build the command to be passed to the Wrapper object
		 	cmd = "#{batch.java_path} -Xmx4g -jar #{batch.picard_path}/picard.jar SamFormatConverter I=#{batch.base_path}/#{batch.batch_id}/assembly/#{sample.panel_version}_#{sample.ex_number}_#{sample.gender.upcase}.sam "
		 	cmd += "O=#{batch.base_path}/#{batch.batch_id}/assembly/#{sample.panel_version}_#{sample.ex_number}_#{sample.gender.upcase}.bam TMP_DIR=#{batch.tmp_path} VALIDATION_STRINGENCY=SILENT"
			
		 	#Add an entry to the Logger object
		 	logger.info('stage') { "Converting SAM to BAM :: Sample #{sample.panel_version}_#{sample.ex_number}" }
		 	
		 	#Send the command to the Wrapper object, return the IO.pipe output
			output = this_wrapper.run_command(cmd, logger)
		  
			return output
	end
	

end
