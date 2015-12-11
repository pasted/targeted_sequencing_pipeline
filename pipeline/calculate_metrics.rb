class CalculateMetrics

	require_relative 'wrapper'
	
	def	create_to_interval_list(sample, batch, logger)
		
			this_wrapper = Wrapper.new
		#	puts "Running Phenotype specific Calculate HS Metrics on sample...#{sample.sequencing_panel_version}_#{sample.sample_id}"
		#		
		#	input_file_string		=	"#{batch.base_path}/#{batch.batch_id}/assembly/#{sample.panel_version.downcase}_#{sample.sample_id}_#{sample.gender.upcase}.realigned.bam"
		#	output_file_string	=	"#{batch.base_path}/#{batch.batch_id}/metrics/#{sample.panel_version.downcase}_#{sample.sample_id}_#{sample.gender.upcase}.phenotype.bait_capture_metrics"
		#	this_panel = batch.select_panel(sample.panel_version.upcase)
	  #
		#	phenotype_intervals_path	=	"#{this_panel.intervals_directory}/#{sample.panel_version.downcase}_#{sample.phenotype}_metrics.bed"
		#	
		#	cmd = "( #{batch.java_path} -Xmx4g -jar #{batch.picard_path}/picard.jar BedToIntervalList INPUT=#{input_file_string} OUTPUT=#{output_file_string} SEQUENCE_DICTIONARY=/mnt/Data1/resources/human_g1k_v37.dict "
		#  
		#	logger.info('stage') { "Phenotype specific HS Metrics :: Sample #{sample.panel_version}_#{sample.sample_id}" }
			output = this_wrapper.run_command(cmd, logger)
			
			return output
	
	end
	
	def overall_metrics(sample, batch, logger)
		
			this_wrapper = Wrapper.new
			puts "Running Overall Calculate HS Metrics on sample...#{sample.sequencing_panel_version}_#{sample.sample_id}"
			
			input_file_string="../../assembly/#{sample.panel_version.downcase}_#{sample.sample_id}_#{sample.gender.upcase}.realigned.bam"
			output_file_string="../../metrics/#{sample.panel_version.downcase}_#{sample.sample_id}_#{sample.gender.upcase}.overall.bait_capture_metrics"
			
			this_panel = batch.select_panel(sample.panel_version.upcase)
			puts sample.panel_version
			puts this_panel.inspect
			bait_intervals_path="#{this_panel.intervals_directory}/#{sample.panel_version.downcase}_all_covered_bases.interval_list"
			
			cmd = "( #{batch.java_path} -Xmx4g -jar #{batch.picard_path}/picard.jar CalculateHsMetrics I=#{input_file_string} o=#{output_file_string} LEVEL=SAMPLE "
			cmd += "BAIT_INTERVALS=#{bait_intervals_path} TARGET_INTERVALS=#{bait_intervals_path} VALIDATION_STRINGENCY=SILENT )"

			logger.info('stage') { "Overall HS Metrics :: Sample #{sample.panel_version}_#{sample.sample_id}" }
			output = this_wrapper.run_command(cmd, logger)
		  
			return output
	end
	
	def phenotype_metrics(sample, batch, logger)
		
			this_wrapper = Wrapper.new
			puts "Running Phenotype specific Calculate HS Metrics on sample...#{sample.sequencing_panel_version}_#{sample.sample_id}"
				
			input_file_string		=	"../../assembly/#{sample.panel_version.downcase}_#{sample.sample_id}_#{sample.gender.upcase}.realigned.bam"
			output_file_string	=	"../../metrics/#{sample.panel_version.downcase}_#{sample.sample_id}_#{sample.gender.upcase}_#{sample.phenotype.upcase}.phenotype.bait_capture_metrics"
			this_panel = batch.select_panel(sample.panel_version.upcase)
	
			phenotype_intervals_path	=	"#{this_panel.intervals_directory}/#{sample.panel_version.downcase}_#{sample.phenotype}_metrics.interval_list"
			
			cmd = "( #{batch.java_path} -Xmx4g -jar #{batch.picard_path}/picard.jar CalculateHsMetrics I=#{input_file_string} o=#{output_file_string} LEVEL=SAMPLE "
			cmd += "BAIT_INTERVALS=#{phenotype_intervals_path} TARGET_INTERVALS=#{phenotype_intervals_path} VALIDATION_STRINGENCY=SILENT )"
		  
			logger.info('stage') { "Phenotype specific HS Metrics :: Sample #{sample.panel_version}_#{sample.sample_id}" }
			output = this_wrapper.run_command(cmd, logger)
			
			return output
	end
	
	def phenotype_coverage(sample, batch, logger)
		
			this_wrapper = Wrapper.new
			puts "Running Phenotype specific GATK DepthOfCoverage on sample...#{sample.sequencing_panel_version}_#{sample.sample_id}"
				
			input_file_string  = "../../assembly/#{sample.panel_version.downcase}_#{sample.sample_id}_#{sample.gender.upcase}.realigned.bam"
			output_file_string = "../../coverage/#{sample.panel_version.downcase}_#{sample.sample_id}_#{sample.gender.upcase}_#{sample.phenotype.upcase}.phenotype.by_base_coverage"
			this_panel = batch.select_panel(sample.panel_version.upcase)
			
			phenotype_intervals_path	=	"#{this_panel.intervals_directory}/#{sample.panel_version.downcase}_#{sample.phenotype}_metrics.interval_list"
		
			cmd = "( #{batch.java_path} -Xmx4g -jar #{batch.gatk_path} -T DepthOfCoverage -L #{phenotype_intervals_path} -R #{batch.reference_path} -I #{input_file_string} -o #{output_file_string} "
			cmd += "-mmq 30 -mbq 30 -dels -ct 1 -ct 10 -ct 20 -ct 30 -ct 40 -ct 50 -omitLocusTable )"
		  
			logger.info('stage') { "Phenotype specific coverage :: Sample #{sample.panel_version}_#{sample.sample_id}" }
			output = this_wrapper.run_command(cmd, logger)
			
			return output
	end
	
	def batch_phenotype_metrics(sample, batch, logger)
		batch_file_string = "../../metrics/#{sample.panel_version.downcase}_#{sample.sample_id}_#{sample.gender.upcase}_#{sample.phenotype.upcase}.phenotype.bait_capture_metrics"
		metrics_array = self.parse_metrics(batch_file_string)
		puts metrics_array.inspect
	end

end
