# @author Garan Jones
# MarkDuplicates class: method for using Picard to mark the duplicate reads in the selected BAM file, results of library prep artefacts
class MarkDuplicates
	require_relative 'wrapper'			
		
		# @author Garan Jones
		# Use Picard to mark the duplicate reads in the selected BAM file
  	# @param sample [Object] A Sample object
  	# @param batch [Object] A Batch object containing details on the batch
  	# @param logger [Object] A Ruby Logger object
  	# @return [Array<String, Object>] An array with command exit status and an IO.pipe object
	def remove_dups(sample, batch, logger)
		
			this_wrapper = Wrapper.new
			puts "Running MarkDuplicates on sample...#{sample.sequencing_panel_version}_#{sample.ex_number}"

			#Build the command to be passed to the Wrapper object
		 	cmd = "#{batch.java_path} -Xmx4g -Djava.io.tmpdir=#{batch.tmp_path} -jar #{batch.picard_path}/picard.jar MarkDuplicates INPUT=#{batch.base_path}/#{batch.batch_id}/assembly/#{sample.panel_version}_#{sample.ex_number}_#{sample.gender.upcase}.fixmate.sorted.bam "
		 	cmd += "OUTPUT=#{batch.base_path}/#{batch.batch_id}/assembly/#{sample.panel_version}_#{sample.ex_number}_#{sample.gender.upcase}.fixmate.rmdup.bam METRICS_FILE=#{batch.base_path}/#{batch.batch_id}/duplicates/#{sample.panel_version}_#{sample.ex_number}.duplicates REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true"
		 	
		 	#Add an entry to the Logger object
		 	logger.info('stage') { "Mark duplicates :: Sample #{sample.panel_version}_#{sample.ex_number}" }
			
		 	#Send the command to the Wrapper object, return the IO.pipe output
		 	output = this_wrapper.run_command(cmd, logger)
		  
			return output
	end
	

end
