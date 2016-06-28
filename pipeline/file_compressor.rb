# @author Garan Jones
# FileCompressor class: Compression methods
class FileCompressor

	require_relative 'wrapper'

		# @author Garan Jones
		# Gzip files
		# @param read_direction [String] String containing the read direction
  	# @param panel_version [String] String containing the panel version
  	# @param gender [String] Current gender of batched samples to the checked
  	# @param batch [Object] A Batch object containing details on the batch
  	# @param logger [Object] A Ruby Logger object
  	# @return [Array<String, Object>] An array with command exit status and an IO.pipe object
	def gzip_fastq(read_direction, sample, batch, logger)
		
			this_wrapper = Wrapper.new
			puts "Running Gzip on ..."

			#Find files matching the pattern
			file_array = Dir["#{batch.base_path}/#{batch.batch_id}/raw_reads/#{sample.sequencing_panel_version}_#{sample.ex_number}_*_#{read_direction}_001.fastq"]
		 	
		 	#Only use the first file name returned - should only be one 
		 	reads=file_array.first
		 	puts reads
		 	if File.symlink?("#{reads}") != nil
		 		cmd = "gzip #{reads}"
      	
				#Add an entry to the Logger object, stating that Gzip is being run on the sample fastqs
				logger.info('stage') { "Gzip Fastq :: Sample #{sample.panel_version}_#{sample.ex_number}" }
				
				#Send the command to the Wrapper object, return the IO.pipe output
				output = this_wrapper.run_command(cmd, logger)
		  	
				return output
			else
				return nil
			end

	end
	
		# @author Garan Jones
		# UnGzip files
		# @param read_direction [String] String containing the read direction
  	# @param panel_version [String] String containing the panel version
  	# @param gender [String] Current gender of batched samples to the checked
  	# @param batch [Object] A Batch object containing details on the batch
  	# @param logger [Object] A Ruby Logger object
  	# @return [Array<String, Object>] An array with command exit status and an IO.pipe object
	def gunzip_fastq(read_direction, sample, batch, logger)
		
			this_wrapper = Wrapper.new
			puts "Running Gzip on ..."

			#Find files matching the pattern
			file_array = Dir["#{batch.base_path}/#{batch.batch_id}/raw_reads/#{sample.sequencing_panel_version}_#{sample.sample_id}_*_#{read_direction}_001.fastq"]
		 	
		 	#Only use the first file name returned - should only be one 
		 	reads=file_array.first
		 	
		 	if File.symlink?("#{reads}") != nil
		 		cmd = "gunzip #{reads}"
				
				#Add an entry to the Logger object, stating that Gzip is being run on the sample fastqs
				logger.info('stage') { "Gunzip Fastq :: Sample #{sample.panel_version}_#{sample.sample_id}" }
				
				#Send the command to the Wrapper object, return the IO.pipe output
				output = this_wrapper.run_command(cmd, logger)
		  	
				return output
			else
				return nil
			end

	end

end
