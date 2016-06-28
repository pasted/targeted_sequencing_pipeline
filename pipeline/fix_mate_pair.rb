# @author Garan Jones
# FixMatePair class: method for using Picard to fix the mate pair information in the selected BAM file
class FixMatePair
	require_relative 'wrapper'			

		# @author Garan Jones
		# Use Picard to fix the mate pair information of sample BAM file
  	# @param sample [Object] A Sample object
  	# @param batch [Object] A Batch object containing details on the batch
  	# @param logger [Object] A Ruby Logger object
  	# @return [Array<String, Object>] An array with command exit status and an IO.pipe object
	def fix_information(sample, batch, logger)
		
			this_wrapper = Wrapper.new
			puts "Fixing Mate-Pair information...#{sample.sequencing_panel_version}_#{sample.ex_number}"

			#Build the command to be passed to the Wrapper object
		 	cmd = "#{batch.java_path} -Xmx4g -jar #{batch.picard_path}/picard.jar FixMateInformation INPUT=#{batch.base_path}/#{batch.batch_id}/assembly/#{sample.panel_version}_#{sample.ex_number}_#{sample.gender.upcase}.bam "
		 	cmd += "OUTPUT=#{batch.base_path}/#{batch.batch_id}/assembly/#{sample.panel_version}_#{sample.ex_number}_#{sample.gender.upcase}.fixmate.bam TMP_DIR=#{batch.tmp_path} SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true"
			
		 	#Add an entry to the Logger object
		 	logger.info('stage') { "Fixmate pair information :: Sample #{sample.panel_version}_#{sample.ex_number}" }
			
		 	#Send the command to the Wrapper object, return the IO.pipe output
		 	output = this_wrapper.run_command(cmd, logger)
		  
			return output
	end
	

end
