class Wget
	require_relative 'wrapper'
	
	def fetch_files(sample, parser, logger)
		
			this_wrapper = Wrapper.new
			puts "Downloading FASTQ raw reads...#{sample.sequencing_panel_version}_#{sample.sample_id}"
			
			logger.info('stage') { "Downloading FASTQ raw reads :: Sample #{sample.panel_version}_#{sample.sample_id}" }
			this_folder = "Sample_#{sample.sequencing_panel_version}_#{sample.sample_id}/raw_illumina_reads/"
			output = this_wrapper.run_command("wget -r -l 1 --no-parent -nH --cut-dirs=3 -A \"*.fastq.md5sum\" -P \"../../raw_reads/md5sum/\" -w 1  #{parser.ftp_url}/#{parser.flowcell}/#{this_folder}", logger)
			output = this_wrapper.run_command("wget -r -l 1 --no-parent -nH --cut-dirs=3 -A \"*.fastq\" -P \"../../raw_reads/\" -w 1 #{parser.ftp_url}/#{parser.flowcell}/#{this_folder}", logger)
		  

			## Wget options
			## -r recursively download
			## -l maximum depth to descend to
			## --no-parent ignore upper directories
			## -nH dont save under hostname folder (zeus-galaxy.ex.ac.uk)
			## --cut-dirs number of directory levels to remove
			## -A select only files with this extension
			## -P directory prefix
			## -w wait time between requests
			#wget_output = `wget -r -l 1 --no-parent -nH --cut-dirs=3 -A "*.fastq" -P "#{base_path}/raw_reads/" -w 1 #{ftp_string}/#{this_folder}`				
			#wget_output = `wget -r -l 1 --no-parent -nH --cut-dirs=3 -A "*.fastq.md5sum" -P "#{parser.base_path}/#{parser.batch_id}/raw_reads/md5sum/" -w 1  #{parser.ftp_url}/#{this_folder}`
			
			return output
	end
	
	def fetch_md5sum(sample, parser, logger)
		
			this_wrapper = Wrapper.new
			puts "Downloading FASTQ MD5SUMs ...#{sample.sequencing_panel_version}_#{sample.sample_id}"
			
			logger.info('stage') { "Downloading FASTQ MD5SUMs :: Sample #{sample.panel_version}_#{sample.sample_id}" }
			this_folder = "Sample_#{sample.sequencing_panel_version}_#{sample.sample_id}/raw_illumina_reads/"
			output = this_wrapper.run_command("wget -r -l 1 --no-parent -nH --cut-dirs=3 -A \"*.fastq.md5sum\" -P \"../../raw_reads/md5sum/\" -w 1  #{parser.ftp_url}/#{parser.flowcell}/#{this_folder}", logger)
		  
			return output
	end
	
	def fetch_fastq(sample, parser, logger)
		
			this_wrapper = Wrapper.new
			puts "Downloading FASTQ raw reads...#{sample.sequencing_panel_version}_#{sample.sample_id}"
			
			logger.info('stage') { "Downloading FASTQ raw reads :: Sample #{sample.panel_version}_#{sample.sample_id}" }
			this_folder = "Sample_#{sample.sequencing_panel_version}_#{sample.sample_id}/raw_illumina_reads/"
			output = this_wrapper.run_command("wget -r -l 1 --no-parent -nH --cut-dirs=3 -A \"*.fastq\" -P \"../../raw_reads/\" -w 1 #{parser.ftp_url}/#{parser.flowcell}/#{this_folder}", logger)
		  
			return output
	end
	
end
