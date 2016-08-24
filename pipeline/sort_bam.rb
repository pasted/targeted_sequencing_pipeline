# @author Garan Jones
# SortBam class: Use Picard to order the BAM file by coordinates
class SortBam
	require_relative 'wrapper'			
			
	  # @author Garan Jones
		# Use Picard to order the BAM file by coordinate
  	# @param sample [Object] The Sample object
  	# @param batch [Object] The Batch object
  	# @param logger [Object] The Logger object
  	# @return [Array<String, Object>] An array with command exit status and an IO.pipe object
	def sort(sample, batch, logger)
		
			this_wrapper = Wrapper.new
			puts "Sorting BAM...#{sample.sequencing_panel_version}_#{sample.ex_number}"

			#Build the command to be passed to the Wrapper object
		 	cmd = "#{batch.java_path} -Xmx4g -Djava.io.tmpdir=#{batch.tmp_path} -jar #{batch.picard_path}/picard.jar SortSam INPUT=#{batch.base_path}/#{batch.batch_id}/assembly/#{sample.panel_version}_#{sample.ex_number}_#{sample.gender.upcase}.fixmate.bam "
		 	cmd += "OUTPUT=#{batch.base_path}/#{batch.batch_id}/assembly/#{sample.panel_version}_#{sample.ex_number}_#{sample.gender.upcase}.fixmate.sorted.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true"
			
		 	#Add an entry to the Logger object
		 	logger.info('stage') { "Sorting BAM :: Sample #{sample.panel_version}_#{sample.ex_number}" }
		 	
		 	#Send the command to the Wrapper object, return the IO.pipe output
			output = this_wrapper.run_command(cmd, logger)
		  
			return output
	end
	

end
