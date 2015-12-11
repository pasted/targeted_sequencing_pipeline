class MarkDuplicates
	require_relative 'wrapper'			
			
	def remove_dups(sample, batch, logger)
		
			this_wrapper = Wrapper.new
			puts "Running MarkDuplicates on sample...#{sample.sequencing_panel_version}_#{sample.sample_id}"

		 	cmd = "#{batch.java_path} -Xmx4g -Djava.io.tmpdir=#{batch.tmp_path} -jar #{batch.picard_path}/picard.jar MarkDuplicates INPUT=../../assembly/#{sample.panel_version}_#{sample.sample_id}_#{sample.gender.upcase}.fixmate.sorted.bam "
		 	cmd += "OUTPUT=../../assembly/#{sample.panel_version}_#{sample.sample_id}_#{sample.gender.upcase}.fixmate.rmdup.bam METRICS_FILE=../../duplicates/#{sample.panel_version}_#{sample.sample_id}.duplicates REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true"
		 	
		 	logger.info('stage') { "Mark duplicates :: Sample #{sample.panel_version}_#{sample.sample_id}" }
			output = this_wrapper.run_command(cmd, logger)
		  
			return output
	end
	

end
