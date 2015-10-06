class VariantsToTable
	require_relative 'wrapper'
	
	def genotype_table(input_file_path, output_file_path, batch, logger)
			this_wrapper = Wrapper.new
			puts "Converting VCF to genotype table on batch #{batch.batch_id}..."

			cmd = "#{batch.java_path} -Xmx4g -jar #{batch.gatk_path} -T VariantsToTable -R #{batch.reference_path} --variant #{input_file_path} "
			cmd += "-o #{output_file_path} --allowMissingData -F CHROM -F POS -F ID -F REF -F ALT -F QUAL -F FILTER -F AC -F AF -GF GT "
			cmd += "-log #{batch.base_path}/logs/#{batch.batch_id}_GT_table_filtered_ug_call.log"
			
			logger.info('stage') { "Genotype table :: Batch #{batch.batch_id}" }
			output = this_wrapper.run_command(cmd, logger)
		  
			return output
	end
	
	def allelic_depth_table(input_file_path, output_file_path, batch, logger)
			this_wrapper = Wrapper.new
			puts "Converting VCF to allele depth table on batch #{batch.batch_id}..."

			cmd = "#{batch.java_path} -Xmx4g -jar #{batch.gatk_path} -T VariantsToTable -R #{batch.reference_path} --variant #{input_file_path} "
			cmd += "-o #{output_file_path} --allowMissingData -F CHROM -F POS -F ID -F REF -F ALT -F QUAL -F FILTER -F AC -F AF -GF AD "
			cmd += "-log #{batch.base_path}/logs/#{batch.batch_id}_GT_table_filtered_ug_call.log"
			
			logger.info('stage') { "Genotype table :: Batch #{batch.batch_id}" }
			output = this_wrapper.run_command(cmd, logger)
		  
			return output
	end
			
	
end
