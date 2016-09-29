# @author Garan Jones
# CleanCaller class: Container class for various variant calling software calls.
class CleanCaller
	require_relative 'wrapper'
	
		# Runs CleanCall 
		#	@param sample [Object] The Sample instance
		# @param batch [Object] The Batch instance, derived from the Config YAML file
		# @param logger [Object] The Logger instance
		#	@param bam_out_required [Boolean] A Boolean to select whether or not to output the BAM file
		# @return [Array<String, Object>] An array with command exit status and an IO.pipe object
	def cleancall_check(sample, batch, logger)
		
			this_wrapper = Wrapper.new
			puts "Running CleanCall on sample...#{sample.sequencing_panel_version.downcase}_#{sample.ex_number}"
			
			input_file_string = "#{batch.base_path}/#{batch.batch_id}/assembly/#{sample.panel_version.downcase}_#{sample.ex_number}_#{sample.gender.upcase}.realigned.bam"
			id_label = "#{sample.panel_version.downcase}_#{sample.ex_number}"
			
			cmd = "#{batch.cleancall_path}/bin/samtools view -q 20 -F 0x0704 -uh #{input_file_string} | #{batch.cleancall_path}/bin/samtools calmd -AEbr - #{batch.reference_path} | #{batch.cleancall_path}/bin/samtools mpileup"			
			cmd += " -s -O -f #{batch.reference_path} -d 255 -l #{batch.cleancall_resources_path}/#{sample.panel_version.downcase}_exac_subset.vcf - | #{batch.cleancall_path}/bin/bgzip -c > #{batch.base_path}/#{batch.batch_id}/clean_call/clean_call.pileup.#{id_label}.txt.gz \n "
			cmd += " #{batch.cleancall_path}/bin/tabix -f -s 1 -b 2 -e 2 #{batch.base_path}/#{batch.batch_id}/clean_call/clean_call.pileup.#{id_label}.txt.gz \n "
			cmd += " #{batch.cleancall_path}/bin/cleanCall verify --vcf #{batch.cleancall_resources_path}/#{sample.panel_version.downcase}_exac_subset.vcf --minAF 0.01 --minCallRate 0.95 --out #{batch.base_path}/#{batch.batch_id}/clean_call/clean_call.verify.#{id_label} --mpu #{batch.base_path}/#{batch.batch_id}/clean_call/clean_call.pileup.#{id_label}.txt.gz --smID #{id_label} --maxDepth 20 \n "
			cmd += " rm #{batch.base_path}/#{batch.batch_id}/clean_call/cleanCall.pileup.#{id_label}.txt.gz #{batch.base_path}/#{batch.batch_id}/clean_call/clean_call.pileup.#{id_label}.txt.gz.tbi #{batch.base_path}/#{batch.batch_id}/clean_call/clean_call.verify.#{id_label}.depthSM"
			
			puts cmd
			
			logger.info('stage') { "Clean caller :: Sample #{sample.panel_version.downcase}_#{sample.ex_number}" }
			output = this_wrapper.run_command(cmd, logger)

			return output
	end
	
	def cleancall_mpileup(sample, batch, logger)
		
			this_wrapper = Wrapper.new
			puts "Running CleanCall mpileup on sample...#{sample.sequencing_panel_version.downcase}_#{sample.ex_number}"
			
			input_file_string = "#{batch.base_path}/#{batch.batch_id}/assembly/#{sample.panel_version.downcase}_#{sample.ex_number}_#{sample.gender.upcase}.realigned.bam"
			id_label = "#{sample.panel_version.downcase}_#{sample.ex_number}"
			
			cmd = "#{batch.cleancall_path}/bin/samtools view -q 20 -F 0x0704 -uh #{input_file_string} | #{batch.cleancall_path}/bin/samtools calmd -AEbr - #{batch.reference_path} | #{batch.cleancall_path}/bin/samtools mpileup"			
			cmd += " -s -O -f #{batch.reference_path} -d 255 -l #{batch.cleancall_resources_path}/#{sample.panel_version.downcase}_exac_subset.vcf - | #{batch.cleancall_path}/bin/bgzip -c > #{batch.base_path}/#{batch.batch_id}/clean_call/clean_call.pileup.#{id_label}.txt.gz \n "
			
			puts cmd
			
			logger.info('stage') { "Clean caller mpileup :: Sample #{sample.panel_version.downcase}_#{sample.ex_number}" }
			output = this_wrapper.run_command(cmd, logger)

			return output
	end
	
	def cleancall_tabix(sample, batch, logger)
		
			this_wrapper = Wrapper.new
			puts "Running CleanCall tabix on sample...#{sample.sequencing_panel_version.downcase}_#{sample.ex_number}"
			
			input_file_string = "#{batch.base_path}/#{batch.batch_id}/assembly/#{sample.panel_version.downcase}_#{sample.ex_number}_#{sample.gender.upcase}.realigned.bam"
			id_label = "#{sample.panel_version.downcase}_#{sample.ex_number}"
			
			cmd = " #{batch.cleancall_path}/bin/tabix -f -s 1 -b 2 -e 2 #{batch.base_path}/#{batch.batch_id}/clean_call/clean_call.pileup.#{id_label}.txt.gz"			
			puts cmd
			
			logger.info('stage') { "Clean caller tabix :: Sample #{sample.panel_version.downcase}_#{sample.ex_number}" }
			output = this_wrapper.run_command(cmd, logger)

			return output
	end
	
	def cleancall_verify(sample, batch, logger)
		
			this_wrapper = Wrapper.new
			puts "Running CleanCall verify on sample...#{sample.sequencing_panel_version.downcase}_#{sample.ex_number}"
			
			input_file_string = "#{batch.base_path}/#{batch.batch_id}/assembly/#{sample.panel_version.downcase}_#{sample.ex_number}_#{sample.gender.upcase}.realigned.bam"
			id_label = "#{sample.panel_version.downcase}_#{sample.ex_number}"
			
			cmd = "#{batch.cleancall_path}/bin/cleanCall verify --vcf #{batch.cleancall_resources_path}/#{sample.panel_version.downcase}_exac_subset.vcf --minAF 0.01 --minCallRate 0.95 --out #{batch.base_path}/#{batch.batch_id}/clean_call/clean_call.verify.#{id_label} --mpu #{batch.base_path}/#{batch.batch_id}/clean_call/clean_call.pileup.#{id_label}.txt.gz --smID #{id_label} --maxDepth 20 "
			
			puts cmd
			
			logger.info('stage') { "Clean caller verify :: Sample #{sample.panel_version.downcase}_#{sample.ex_number}" }
			output = this_wrapper.run_command(cmd, logger)

			return output
	end
	
	def cleancall_cleanup(sample, batch, logger)
		
			this_wrapper = Wrapper.new
			puts "Running CleanCall cleanup on sample...#{sample.sequencing_panel_version.downcase}_#{sample.ex_number}"
			
			input_file_string = "#{batch.base_path}/#{batch.batch_id}/assembly/#{sample.panel_version.downcase}_#{sample.ex_number}_#{sample.gender.upcase}.realigned.bam"
			id_label = "#{sample.panel_version.downcase}_#{sample.ex_number}"
			
			cmd = " rm #{batch.base_path}/#{batch.batch_id}/clean_call/clean_call.pileup.#{id_label}.txt.gz #{batch.base_path}/#{batch.batch_id}/clean_call/clean_call.pileup.#{id_label}.txt.gz.tbi #{batch.base_path}/#{batch.batch_id}/clean_call/clean_call.verify.#{id_label}.depthSM"
		
			puts cmd
			
			logger.info('stage') { "Clean caller cleanup :: Sample #{sample.panel_version.downcase}_#{sample.ex_number}" }
			output = this_wrapper.run_command(cmd, logger)

			return output
	end
	
end
