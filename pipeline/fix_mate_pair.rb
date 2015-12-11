class FixMatePair
	require_relative 'wrapper'			
			
	def fix_information(sample, batch, logger)
		
			this_wrapper = Wrapper.new
			puts "Fixing Mate-Pair information...#{sample.sequencing_panel_version}_#{sample.sample_id}"

		 	cmd = "#{batch.java_path} -Xmx4g -jar #{batch.picard_path}/picard.jar FixMateInformation INPUT=../../assembly/#{sample.panel_version}_#{sample.sample_id}_#{sample.gender.upcase}.bam "
		 	cmd += "OUTPUT=../../assembly/#{sample.panel_version}_#{sample.sample_id}_#{sample.gender.upcase}.fixmate.bam TMP_DIR=#{batch.tmp_path} SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true"
			
		 	logger.info('stage') { "Fixmate pair information :: Sample #{sample.panel_version}_#{sample.sample_id}" }
			output = this_wrapper.run_command(cmd, logger)
		  
			return output
	end
	

end
