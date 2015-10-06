class RealignIndels
	require_relative 'wrapper'			
			
	def create_targets(sample, batch, logger)
		
			this_wrapper = Wrapper.new
			puts "Running RealignerTargetCreator on sample...#{sample.sequencing_panel_version}_#{sample.sample_id}"
			
			input_file_string="-I #{batch.base_path}/#{batch.batch_id}/assembly/#{sample.panel_version}_#{sample.sample_id}_#{sample.gender.upcase}.fixmate.rmdup.bam"
						
			cmd = "#{batch.java_path} -Djava.io.tmpdir=#{batch.tmp_path} -Xmx6g -jar #{batch.gatk_path} -T RealignerTargetCreator -known #{batch.known_path} -R #{batch.reference_path} "
			cmd += "-log #{batch.base_path}/#{batch.batch_id}/logs/RealignerTargetCreator.#{sample.panel_version}_#{sample.sample_id}.nodup.log #{input_file_string} -o #{batch.base_path}/#{batch.batch_id}/intervals/#{sample.panel_version}_#{sample.sample_id}.intervals"

			logger.info('stage') { "Realigner Target Creator :: Sample #{sample.panel_version}_#{sample.sample_id}" }
			output = this_wrapper.run_command(cmd, logger)
		  
			return output
	end
	
	def realign(sample, batch, logger)
		
			this_wrapper = Wrapper.new
			puts "Running IndelRealigner on sample...#{sample.sequencing_panel_version}_#{sample.sample_id}"
						
			input_file_string="#{batch.base_path}/#{batch.batch_id}/assembly/#{sample.panel_version}_#{sample.sample_id}_#{sample.gender.upcase}.fixmate.rmdup.bam"
			output_file_string="#{batch.base_path}/#{batch.batch_id}/assembly/#{sample.panel_version}_#{sample.sample_id}_#{sample.gender.upcase}.realigned.bam "
			
		  cmd = "#{batch.java_path} -Djava.io.tmpdir=#{batch.tmp_path} -Xmx6g -jar #{batch.gatk_path} -T IndelRealigner -known #{batch.known_path} -R #{batch.reference_path} "
			cmd += "-log #{batch.base_path}/#{batch.batch_id}/logs/IndelRealigner.#{sample.panel_version}_#{sample.sample_id}.nodup.log -compress 0 --maxReadsInMemory 1000000000 --maxReadsForRealignment 6000000 "
			cmd += "-I #{input_file_string} -o #{output_file_string} -targetIntervals #{batch.base_path}/#{batch.batch_id}/intervals/#{sample.panel_version}_#{sample.sample_id}.intervals"
	
			logger.info('stage') { "IndelRealigner :: Sample #{sample.panel_version}_#{sample.sample_id}" }
			output = this_wrapper.run_command(cmd, logger)
		  
			return output
	end
	

end
