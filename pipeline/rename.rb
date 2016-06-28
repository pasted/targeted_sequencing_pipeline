# @author Garan Jones
# Rename class: call a perl script to rename fastq files according to a given regex
class Rename
	require_relative 'wrapper'

		# @author Garan Jones
		# Call a perl script to rename fastq files according to a given regex
		#	perl -w ${base_path}/scripts/rename.pl 's/(.*)_[ACGT]{6}_(.*).fastq/$1_$2.fastq/' ${base_path}/raw_reads/*.fastq
		#	perl -w ${base_path}/scripts/rename.pl 's/(.*)_[ACGT]{6}_(.*).fastq/$1_$2.fastq/' ${base_path}/raw_reads/md5sum/*.fastq.md5sum
  	# @param sample [Object] A Sample object
  	# @param parser [Object] A Parser object
  	# @param logger [Object] A Ruby Logger object
  	# @return [Array<String, Object>] An array with command exit status and an IO.pipe object
	def remove_fastq_adaptor_string(sample, parser, logger)
		
			this_wrapper = Wrapper.new
			puts "Renaming FASTQ raw reads...#{sample.sequencing_panel_version}_#{sample.ex_number}"
			
			#Add an entry to the Logger object
			logger.info('stage') { "Renaming FASTQ raw reads :: Sample #{sample.panel_version}_#{sample.ex_number}" }
			
			#Build the command to be passed to the Wrapper object and run it
			output = this_wrapper.run_command("perl -w #{parser.base_path}/#{parser.batch_id}/scripts/perl_scripts/rename.pl 's/(.*)_[ACGT]{6}_(.*).fastq/$1_$2.fastq/' #{parser.base_path}/#{parser.batch_id}/raw_reads/*.fastq", logger)
		  
			return output
	end
	
		# @author Garan Jones
	 	# Call a perl script to rename fastq files according to a given regex
   	# @param sample [Object] A Sample object
   	# @param parser [Object] A Parser object
   	# @param logger [Object] A Ruby Logger object
   	# @return [Array<String, Object>] An array with command exit status and an IO.pipe object
	def remove_md5sum_adaptor_string(sample, parser, logger)
		
			this_wrapper = Wrapper.new
			puts "Renaming FASTQ MD5SUM ...#{sample.sequencing_panel_version}_#{sample.ex_number}"
			
			#Add an entry to the Logger object
			logger.info('stage') { "Renaming FASTQ MD5SUM :: Sample #{sample.panel_version}_#{sample.ex_number}" }
			
			#Build the command to be passed to the Wrapper object and run it
			output = this_wrapper.run_command("perl -w #{parser.base_path}/#{parser.batch_id}/scripts/perl_scripts/rename.pl 's/(.*)_[ACGT]{6}_(.*).fastq/$1_$2.fastq/' #{parser.base_path}/#{parser.batch_id}/raw_reads/md5sum/*.fastq.md5sum", logger)
		  
			return output
	end
	
	def rename_symlink(sample, parser, logger)
			
			this_wrapper = Wrapper.new
			puts "Renaming Symlink ...#{sample.sequencing_panel_version}_#{sample.ex_number}"
			
			#Add an entry to the Logger object
			logger.info('stage') { "Renaming Symlink :: Sample #{sample.panel_version}_#{sample.ex_number}" }
			
			output = this_wrapper.run_command("perl -w #{parser.base_path}/#{parser.batch_id}/scripts/perl_scripts/rename.pl 's/#{sample.sequencing_panel_version}_#{sample.sample_id}_(.*).fastq.gz/#{sample.sequencing_panel_version}_#{sample.ex_number}_$1.fastq.gz/' #{parser.base_path}/#{parser.batch_id}/raw_reads/*.fastq.gz", logger)
			
			return output
	end
end
