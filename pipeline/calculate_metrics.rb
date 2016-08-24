# @author Garan Jones
# CalculateMetrics class: Various methods to generate metrics from BAM files
class CalculateMetrics
	require_relative 'wrapper'
	
		# @author Garan Jones
		# STUBBED_METHOD: use Picard to generate interval lists from BED files
  	# @param sample [Object] A Sample object
  	# @param batch [Object] A Batch object containing details on the batch
  	# @param logger [Object] A Ruby Logger object
  	# @return [Array<String, Object>] An array with command exit status and an IO.pipe object
	def	create_to_interval_list(sample, batch, logger)
		
			this_wrapper = Wrapper.new
		#	input_file_string		=	""
		#	output_file_string	=	""
		#	
		#	cmd = "( #{batch.java_path} -Xmx4g -jar #{batch.picard_path}/picard.jar BedToIntervalList INPUT=#{input_file_string} OUTPUT=#{output_file_string} SEQUENCE_DICTIONARY=/mnt/Data1/resources/human_g1k_v37.dict "
		#  
		#	logger.info('stage') { "" }
			output = this_wrapper.run_command(cmd, logger)
			
			return output
	
	end
	
	  # @author Garan Jones
		# Calculate the overall metrics for each sample using Picard CalculateHSMetrics
  	# @param sample [Object] A Sample object
  	# @param batch [Object] A Batch object containing details on the batch
  	# @param logger [Object] A Ruby Logger object
  	# @return [Array<String, Object>] An array with command exit status and an IO.pipe object
	def overall_metrics(sample, batch, logger)
		
			this_wrapper = Wrapper.new
			puts "Running Overall Calculate HS Metrics on sample...#{sample.sequencing_panel_version}_#{sample.ex_number}"
			
			input_file_string="#{batch.base_path}/#{batch.batch_id}/assembly/#{sample.panel_version.downcase}_#{sample.ex_number}_#{sample.gender.upcase}.realigned.bam"
			output_file_string="#{batch.base_path}/#{batch.batch_id}/metrics/#{sample.panel_version.downcase}_#{sample.ex_number}_#{sample.gender.upcase}.overall.bait_capture_metrics"
			
			this_panel = batch.select_panel(sample.panel_version.upcase)
			puts sample.panel_version
			puts this_panel.inspect
			bait_intervals_path="#{this_panel.intervals_directory}/#{sample.panel_version.downcase}_all_covered_bases.interval_list"
			
			cmd = "( #{batch.java_path} -Xmx4g -jar #{batch.picard_path}/picard.jar CalculateHsMetrics I=#{input_file_string} o=#{output_file_string} LEVEL=SAMPLE "
			cmd += "BAIT_INTERVALS=#{bait_intervals_path} TARGET_INTERVALS=#{bait_intervals_path} VALIDATION_STRINGENCY=SILENT )"

			logger.info('stage') { "Overall HS Metrics :: Sample #{sample.panel_version}_#{sample.ex_number}" }
			output = this_wrapper.run_command(cmd, logger)
		  
			return output
	end
	
		# @author Garan Jones
		# Calculate the phenotype specific metrics for each sample using Picard CalculateHSMetrics
  	# @param sample [Object] A Sample object
  	# @param batch [Object] A Batch object containing details on the batch
  	# @param logger [Object] A Ruby Logger object
  	# @return [Array<String, Object>] An array with command exit status and an IO.pipe object	
	def phenotype_metrics(sample, batch, logger)
		
			this_wrapper = Wrapper.new
			puts "Running Phenotype specific Calculate HS Metrics on sample...#{sample.sequencing_panel_version}_#{sample.ex_number}"
				
			input_file_string		=	"#{batch.base_path}/#{batch.batch_id}/assembly/#{sample.panel_version.downcase}_#{sample.ex_number}_#{sample.gender.upcase}.realigned.bam"
			output_file_string	=	"#{batch.base_path}/#{batch.batch_id}/metrics/#{sample.panel_version.downcase}_#{sample.ex_number}_#{sample.gender.upcase}_#{sample.phenotype.upcase}.phenotype.bait_capture_metrics"
			this_panel = batch.select_panel(sample.panel_version.upcase)
	
			phenotype_intervals_path	=	"#{this_panel.intervals_directory}/#{sample.panel_version.downcase}_#{sample.phenotype}_metrics.interval_list"
			
			#Build command to send via Wrapper object to the shell to be executed
			cmd = "( #{batch.java_path} -Xmx4g -jar #{batch.picard_path}/picard.jar CalculateHsMetrics I=#{input_file_string} o=#{output_file_string} LEVEL=SAMPLE "
			cmd += "BAIT_INTERVALS=#{phenotype_intervals_path} TARGET_INTERVALS=#{phenotype_intervals_path} VALIDATION_STRINGENCY=SILENT )"
		  
			#Add info entry to log to say Wrapper object has been sent command
			logger.info('stage') { "Phenotype specific HS Metrics :: Sample #{sample.panel_version}_#{sample.ex_number}" }
			
			#Run command and return IO.pipe object from Wrapper class
			output = this_wrapper.run_command(cmd, logger)
			
			return output
	end
	
		# @author Garan Jones
		# Calculate the DepthOfCoverage using GATK, based on provided phenotype intervals
		# -mmq; Minimum mapping quality of reads to count towards depth
		# -mbq; Minimum quality of bases to count towards depth
		# -dels; Include information on deletions
		# -ct; Coverage threshold (in percent) for summarizing statistics
		# -omitLocusTable; Do not calculate per-sample per-depth counts of loci
  	# @param sample [Object] A Sample object
  	# @param batch [Object] A Batch object containing details on the batch
  	# @param logger [Object] A Ruby Logger object
  	# @return [Array<String, Object>] An array with command exit status and an IO.pipe object	
	def phenotype_coverage(sample, batch, logger)
		
			this_wrapper = Wrapper.new
			puts "Running Phenotype specific GATK DepthOfCoverage on sample...#{sample.sequencing_panel_version}_#{sample.ex_number}"
				
			input_file_string  = "#{batch.base_path}/#{batch.batch_id}/assembly/#{sample.panel_version.downcase}_#{sample.ex_number}_#{sample.gender.upcase}.realigned.bam"
			output_file_string = "#{batch.base_path}/#{batch.batch_id}/coverage/#{sample.panel_version.downcase}_#{sample.ex_number}_#{sample.gender.upcase}_#{sample.phenotype.upcase}.phenotype.by_base_coverage"
			this_panel = batch.select_panel(sample.panel_version.upcase)
			
			phenotype_intervals_path	=	"#{this_panel.intervals_directory}/#{sample.panel_version.downcase}_#{sample.phenotype}_metrics.interval_list"
		
			#Build command to send via Wrapper object to the shell to be executed
			cmd = "( #{batch.java_path} -Xmx4g -jar #{batch.gatk_path} -T DepthOfCoverage -L #{phenotype_intervals_path} -R #{batch.reference_path} -I #{input_file_string} -o #{output_file_string} "
			cmd += "-mmq 30 -mbq 30 -dels -ct 1 -ct 10 -ct 20 -ct 30 -ct 40 -ct 50 -omitLocusTable )"
		  
			#Add info entry to log to say Wrapper object has been sent command
			logger.info('stage') { "Phenotype specific coverage :: Sample #{sample.panel_version}_#{sample.ex_number}" }
			
			#Run command and return IO.pipe object from Wrapper class
			output = this_wrapper.run_command(cmd, logger)
			
			return output
	end
	
		# @author Garan Jones
		# STUBBED_METHOD: Batched phenotype metrics
  	# @param sample [Object] A Sample object
  	# @param batch [Object] A Batch object containing details on the batch
  	# @param logger [Object] A Ruby Logger object
  	# @return [Array<String, Object>] An array with command exit status and an IO.pipe object
	def batch_phenotype_metrics(sample, batch, logger)
		batch_file_string = "#{batch.base_path}/#{batch.batch_id}/metrics/#{sample.panel_version.downcase}_#{sample.ex_number}_#{sample.gender.upcase}_#{sample.phenotype.upcase}.phenotype.bait_capture_metrics"
		metrics_array = self.parse_metrics(batch_file_string)
		puts metrics_array.inspect
	end

end
