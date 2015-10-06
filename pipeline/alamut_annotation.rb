class AlamutAnnotation
	require_relative 'wrapper'			
			
	def annotate(sample, batch, logger)
		
			this_wrapper = Wrapper.new
			puts "Running Alamut annotation on sample...#{sample.sequencing_panel_version}_#{sample.sample_id}"
			
			input_file_string = "#{batch.base_path}/#{batch.batch_id}/variants/haplotyper/#{sample.panel_version}_#{sample.sample_id}_#{sample.gender.upcase}_#{sample.phenotype.upcase}.phenotype.hap_call.vcf"
			annotated_file_string = "#{batch.base_path}/#{batch.batch_id}/variants/haplotyper/#{sample.panel_version}_#{sample.sample_id}_#{sample.gender.upcase}_#{sample.phenotype.upcase}.alamut.alltrans.txt"
			unannotated_file_string = "#{batch.base_path}/#{batch.batch_id}/variants/haplotyper/#{sample.panel_version}_#{sample.sample_id}_#{sample.gender.upcase}_#{sample.phenotype.upcase}.alamut.alltrans.unannotated.txt"

			puts input_file_string
			puts annotated_file_string
			puts unannotated_file_string
			
			cmd = "#{batch.alamut_path} --hgmdUser #{batch.hgmd_user} --hgmdPasswd #{batch.hgmd_pass} --in #{input_file_string} --ann #{annotated_file_string} --unann #{unannotated_file_string} "
			cmd += "--alltrans --ssIntronicRange 2 --outputVCFInfo AC AF AN DP FS MQ MQ0 QD --outputVCFGenotypeData AD DP GQ GT PL --outputVCFQuality --outputVCFFilter"

			logger.info('stage') { "Alamut annotation :: Sample #{sample.panel_version}_#{sample.sample_id}" }
			
			output = this_wrapper.run_command(cmd, logger)
		  
			return output
	end
	

end
