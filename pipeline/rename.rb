class Rename
	require_relative 'wrapper'
	
	#	perl -w ${base_path}/scripts/rename.pl 's/(.*)_[ACGT]{6}_(.*).fastq/$1_$2.fastq/' ${base_path}/raw_reads/*.fastq
	#	perl -w ${base_path}/scripts/rename.pl 's/(.*)_[ACGT]{6}_(.*).fastq/$1_$2.fastq/' ${base_path}/raw_reads/md5sum/*.fastq.md5sum
			
			
	def remove_fastq_adaptor_string(sample, parser, logger)
		
			this_wrapper = Wrapper.new
			puts "Renaming FASTQ raw reads...#{sample.sequencing_panel_version}_#{sample.sample_id}"
			
			logger.info('stage') { "Renaming FASTQ raw reads :: Sample #{sample.panel_version}_#{sample.sample_id}" }
			output = this_wrapper.run_command("perl -w #{parser.base_path}/#{parser.batch_id}/scripts/perl_scripts/rename.pl 's/(.*)_[ACGT]{6}_(.*).fastq/$1_$2.fastq/' #{parser.base_path}/#{parser.batch_id}/raw_reads/*.fastq", logger)
		  
			return output
	end
	
	def remove_md5sum_adaptor_string(sample, parser, logger)
		
			this_wrapper = Wrapper.new
			puts "Renaming FASTQ MD5SUM ...#{sample.sequencing_panel_version}_#{sample.sample_id}"
			
			logger.info('stage') { "Renaming FASTQ MD5SUM :: Sample #{sample.panel_version}_#{sample.sample_id}" }
			output = this_wrapper.run_command("perl -w #{parser.base_path}/#{parser.batch_id}/scripts/perl_scripts/rename.pl 's/(.*)_[ACGT]{6}_(.*).fastq/$1_$2.fastq/' #{parser.base_path}/#{parser.batch_id}/raw_reads/md5sum/*.fastq.md5sum", logger)
		  
			return output
	end
end
