# @author Garan Jones
# RealignIndels class: method for using GATK to preform Indel realignment
class RealignIndels
	require_relative 'wrapper'			
			
		# @author Garan Jones
		# Use GATK to select regions which require local realignment
  	# @param sample [Object] A Sample object
  	# @param batch [Object] A Batch object containing details on the batch
  	# @param logger [Object] A Ruby Logger object
  	# @return [Array<String, Object>] An array with command exit status and an IO.pipe object	
	def create_targets(sample, batch, logger)
		
			this_wrapper = Wrapper.new
			puts "Running RealignerTargetCreator on sample...#{sample.sequencing_panel_version}_#{sample.ex_number}"
			
			input_file_string="-I #{batch.base_path}/#{batch.batch_id}/assembly/#{sample.panel_version}_#{sample.ex_number}_#{sample.gender.upcase}.fixmate.rmdup.bam"
			
			#Build the command to be passed to the Wrapper object		
			cmd = "#{batch.java_path} -Djava.io.tmpdir=#{batch.tmp_path} -Xmx6g -jar #{batch.gatk_path} -T RealignerTargetCreator -known #{batch.known_path} -R #{batch.reference_path} "
			cmd += "-log #{batch.base_path}/#{batch.batch_id}/logs/RealignerTargetCreator.#{sample.panel_version}_#{sample.ex_number}.nodup.log #{input_file_string} -o #{batch.base_path}/#{batch.batch_id}/intervals/#{sample.panel_version}_#{sample.ex_number}.intervals"

			#Add an entry to the Logger object
			logger.info('stage') { "Realigner Target Creator :: Sample #{sample.panel_version}_#{sample.ex_number}" }
			
			#Send the command to the Wrapper object, return the IO.pipe output
			output = this_wrapper.run_command(cmd, logger)
		  
			return output
	end
	
		# @author Garan Jones
		# Use GATK to preform local realignment to correct Indels
  	# @param sample [Object] A Sample object
  	# @param batch [Object] A Batch object containing details on the batch
  	# @param logger [Object] A Ruby Logger object
  	# @return [Array<String, Object>] An array with command exit status and an IO.pipe object
	def realign(sample, batch, logger)
		
			this_wrapper = Wrapper.new
			puts "Running IndelRealigner on sample...#{sample.sequencing_panel_version}_#{sample.ex_number}"
						
			input_file_string="#{batch.base_path}/#{batch.batch_id}/assembly/#{sample.panel_version}_#{sample.ex_number}_#{sample.gender.upcase}.fixmate.rmdup.bam"
			output_file_string="#{batch.base_path}/#{batch.batch_id}/assembly/#{sample.panel_version}_#{sample.ex_number}_#{sample.gender.upcase}.realigned.bam "
			
			#Build the command to be passed to the Wrapper object
		  cmd = "#{batch.java_path} -Djava.io.tmpdir=#{batch.tmp_path} -Xmx6g -jar #{batch.gatk_path} -T IndelRealigner -known #{batch.known_path} -R #{batch.reference_path} "
			cmd += "-log #{batch.base_path}/#{batch.batch_id}/logs/IndelRealigner.#{sample.panel_version}_#{sample.ex_number}.nodup.log -compress 0 --maxReadsInMemory 1000000000 --maxReadsForRealignment 6000000 "
			cmd += "-I #{input_file_string} -o #{output_file_string} -targetIntervals #{batch.base_path}/#{batch.batch_id}/intervals/#{sample.panel_version}_#{sample.ex_number}.intervals"
	
			#Add an entry to the Logger object
			logger.info('stage') { "IndelRealigner :: Sample #{sample.panel_version}_#{sample.ex_number}" }
			
			#Send the command to the Wrapper object, return the IO.pipe output
			output = this_wrapper.run_command(cmd, logger)
		  
			return output
	end
	

end
