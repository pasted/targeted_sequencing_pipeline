# @author Garan Jones
# Wget class: Wrapper class to use wget to download fastq data from sequencing servers to processing server
class Wget
	require_relative 'wrapper'
	
		# Downloads Fastq and MD5sum data
		# @param sample [Object] The Sample instance
		# @param batch [Object] The Batch instance
		# @param logger [Object] The Logger instance
		# @return [Array<String, Object>] An array with command exit status and an IO.pipe object
	def fetch_files(sample, batch, logger)
		
			this_wrapper = Wrapper.new
			puts "Downloading FASTQ raw reads...#{sample.sequencing_panel_version}_#{sample.ex_number}"
			
			logger.info('stage') { "Downloading FASTQ raw reads :: Sample #{sample.panel_version}_#{sample.ex_number}" }
			this_folder = "Sample_#{sample.sequencing_panel_version}_#{sample.ex_number}/raw_illumina_reads/"
			output = this_wrapper.run_command("wget -r -l 1 --no-parent -nH --cut-dirs=3 -A \"*.fastq.md5sum\" -P \"#{batch.base_path}/#{batch.batch_id}/raw_reads/md5sum/\" -w 1  #{batch.ftp_url}/#{batch.flowcell}/#{this_folder}", logger)
			output = this_wrapper.run_command("wget -r -l 1 --no-parent -nH --cut-dirs=3 -A \"*.fastq.gz\" -P \"#{batch.base_path}/#{batch.batch_id}/raw_reads/\" -w 1 #{batch.ftp_url}/#{batch.flowcell}/#{this_folder}", logger)
		  

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
	
		# Downloads only the Fastq MD5sum files
		# @param sample [Object] The Sample instance
		# @param batch [Object] The Batch instance
		# @param logger [Object] The Logger instance
		# @return [Array<String, Object>] An array with command exit status and an IO.pipe object
	def fetch_md5sum(sample, batch, logger)
		
			this_wrapper = Wrapper.new
			puts "Downloading FASTQ MD5SUMs ...#{sample.sequencing_panel_version}_#{sample.ex_number}"
			
			logger.info('stage') { "Downloading FASTQ MD5SUMs :: Sample #{sample.panel_version}_#{sample.ex_number}" }
			this_folder = "Sample_#{sample.sequencing_panel_version}_#{sample.ex_number}/raw_illumina_reads/"
			output = this_wrapper.run_command("wget -r -l 1 --no-parent -nH --cut-dirs=3 -A \"*.fastq.md5sum\" -P \"#{batch.base_path}/#{batch.batch_id}/raw_reads/md5sum/\" -w 1  #{batch.ftp_url}/#{batch.flowcell}/#{this_folder}", logger)
		  
			return output
	end
	
		# Downloads only the Fastq files
		# @param sample [Object] The Sample instance
		# @param batch [Object] The Batch instance
		# @param logger [Object] The Logger instance
		# @return [Array<String, Object>] An array with command exit status and an IO.pipe object
	def fetch_fastq(sample, batch, logger)
		
			this_wrapper = Wrapper.new
			puts "Downloading FASTQ raw reads...#{sample.sequencing_panel_version}_#{sample.ex_number}"
			
			logger.info('stage') { "Downloading FASTQ raw reads :: Sample #{sample.panel_version}_#{sample.ex_number}" }
			this_folder = "Sample_#{sample.sequencing_panel_version}_#{sample.ex_number}/raw_illumina_reads/"
			output = this_wrapper.run_command("wget -r -l 1 --no-parent -nH --cut-dirs=3 -A \"*.fastq.gz\" -P \"#{batch.base_path}/#{batch.batch_id}/raw_reads/\" -w 1 #{batch.ftp_url}/#{batch.flowcell}/#{this_folder}", logger)
		  
			return output
	end
	
end
