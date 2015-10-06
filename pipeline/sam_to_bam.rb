class SamToBam
	require_relative 'wrapper'			
			
	def convert_sam_to_bam(sample, batch, logger)
		
			this_wrapper = Wrapper.new
			puts "Converting SAM to BAM...#{sample.sequencing_panel_version}_#{sample.sample_id}"

		 	cmd = "#{batch.java_path} -Xmx4g -jar #{batch.picard_path}/picard.jar SamFormatConverter I=#{batch.base_path}/#{batch.batch_id}/assembly/#{sample.panel_version}_#{sample.sample_id}_#{sample.gender.upcase}.sam "
		 	cmd += "O=#{batch.base_path}/#{batch.batch_id}/assembly/#{sample.panel_version}_#{sample.sample_id}_#{sample.gender.upcase}.bam TMP_DIR=#{batch.tmp_path} VALIDATION_STRINGENCY=SILENT"
			
		 	logger.info('stage') { "Converting SAM to BAM :: Sample #{sample.panel_version}_#{sample.sample_id}" }
			output = this_wrapper.run_command(cmd, logger)
		  
			return output
	end
	

end
