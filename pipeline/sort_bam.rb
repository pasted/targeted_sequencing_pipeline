class SortBam
	require_relative 'wrapper'			
			
	def sort(sample, batch, logger)
		
			this_wrapper = Wrapper.new
			puts "Sorting BAM...#{sample.sequencing_panel_version}_#{sample.sample_id}"

		 	cmd = "#{batch.java_path} -Xmx4g -Djava.io.tmpdir=#{batch.tmp_path} -jar #{batch.picard_path}/picard.jar SortSam INPUT=#{batch.base_path}/#{batch.batch_id}/assembly/#{sample.panel_version}_#{sample.sample_id}_#{sample.gender.upcase}.fixmate.bam "
		 	cmd += "OUTPUT=#{batch.base_path}/#{batch.batch_id}/assembly/#{sample.panel_version}_#{sample.sample_id}_#{sample.gender.upcase}.fixmate.sorted.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true"
			
		 	logger.info('stage') { "Sorting BAM :: Sample #{sample.panel_version}_#{sample.sample_id}" }
			output = this_wrapper.run_command(cmd, logger)
		  
			return output
	end
	

end
