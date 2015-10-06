class BwaMem
	require_relative 'wrapper'			
			
	def run_bwa_mem(sample, batch, logger)
		
			this_wrapper = Wrapper.new
			puts "Renaming FASTQ raw reads...#{sample.sequencing_panel_version}_#{sample.sample_id}"
			
			#Find files matching the pattern
			forward_file_array = Dir["#{batch.base_path}/#{batch.batch_id}/raw_reads/#{sample.sequencing_panel_version}_#{sample.sample_id}_*_R1_001.fastq"]
			reverse_file_array = Dir["#{batch.base_path}/#{batch.batch_id}/raw_reads/#{sample.sequencing_panel_version}_#{sample.sample_id}_*_R2_001.fastq"]
		 	
		 	#Only use the first file name returned - should only be one 
		 	forward_reads=forward_file_array.first
		 	reverse_reads=reverse_file_array.first
		 	
		 	forward_header = `grep -m 1 '@' #{forward_reads}`
		 	
		 	forward_flowcell = forward_header.split(":")[2]
		 	puts forward_flowcell
		 	
		 	read_group="@RG\tID:#{forward_flowcell}\tPL:#{batch.sequencer_name}\tSM:#{sample.panel_version}_#{sample.sample_id}\tLB:#{batch.batch_id}"
		 	cmd = "#{batch.bwa_path} mem -t 4 -M -v 1 -R '#{read_group}' #{batch.reference_path} #{forward_reads} #{reverse_reads} > #{batch.base_path}/#{batch.batch_id}/assembly/#{sample.panel_version}_#{sample.sample_id}_#{sample.gender.upcase}.sam"
			puts cmd
			
			logger.info('stage') { "BWA MEM alignment :: Sample #{sample.panel_version}_#{sample.sample_id}" }
			output = this_wrapper.run_command(cmd, logger)
		  
			return output
	end
	

end
