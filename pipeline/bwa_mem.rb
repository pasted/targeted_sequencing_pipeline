# @author Garan Jones
# BwaMem class: use BWA MEM to align provided forward and reverse fastq to produce SAM file (version and filepath specified by config.yaml)
class BwaMem
	require_relative 'wrapper'			

		# @author Garan Jones
		# Use BWA MEM to align provided forward and reverse fastq to produce SAM file
  	# @param sample [Object] A Sample object
  	# @param batch [Object] A Batch object containing details on the batch
  	# @param logger [Object] A Ruby Logger object
  	# @return [Array<String, Object>] An array with command exit status and an IO.pipe object
	def run_bwa_mem(sample, batch, logger)
		
			this_wrapper = Wrapper.new
			puts "Renaming FASTQ raw reads...#{sample.sequencing_panel_version}_#{sample.ex_number}"
			
			#Find files matching the pattern
			forward_file_array = Dir["#{batch.base_path}/#{batch.batch_id}/raw_reads/#{sample.sequencing_panel_version}_#{sample.ex_number}_*_R1_001.fastq.gz"]
			reverse_file_array = Dir["#{batch.base_path}/#{batch.batch_id}/raw_reads/#{sample.sequencing_panel_version}_#{sample.ex_number}_*_R2_001.fastq.gz"]
		 	
		 	#Only use the first file name returned - should only be one 
		 	forward_reads=forward_file_array.first
		 	reverse_reads=reverse_file_array.first
		 	
		 	#Obtain fastq header information from first line with matching '@',
		 	#use zgrep due to fastq files being compressed
		 	forward_header = `zgrep -m 1 '@' #{forward_reads}`
		 	
		 	#Split header string on ':' character return the third element which should contain the flowcell id
		 	forward_flowcell = forward_header.split(":")[2]
		 	
		 	#Override the collected flowcell id if the config file contains an alternative
		 	if batch.flowcell != ""
		 		forward_flowcell = "#{batch.flowcell}"
		 		puts forward_flowcell
		 	end
		 	
		 	#Construct the read group that BWA MEM requires
		 	#Build the command to be passed to the Wrapper object
		 	read_group="@RG\tID:#{forward_flowcell}\tPL:#{batch.sequencer_name}\tSM:#{sample.panel_version}_#{sample.ex_number}\tLB:#{batch.batch_id}"
		 	cmd = "#{batch.bwa_path} mem -t 4 -M -v 1 -R '#{read_group}' #{batch.reference_path} #{forward_reads} #{reverse_reads} > #{batch.base_path}/#{batch.batch_id}/assembly/#{sample.panel_version}_#{sample.ex_number}_#{sample.gender.upcase}.sam"
			puts cmd
			
			#Add an entry to the Logger object, stating that BWA MEM is being run on the sample
			logger.info('stage') { "BWA MEM alignment :: Sample #{sample.panel_version}_#{sample.ex_number}" }
			
			#Send the command to the Wrapper object, return the IO.pipe output
			output = this_wrapper.run_command(cmd, logger)
		  
			return output
	end
	

end
