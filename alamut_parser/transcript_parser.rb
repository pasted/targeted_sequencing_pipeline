#Script to retrieve Alamut annotation on a transcript by transcript basis when using the --alltrans flag in Alamut-HT
#transcript_list.txt contains the transcripts to be written out to file.
class TranscriptParser

		attr_accessor :input_file_name, :output_file_name, :sample_id, :transcripts_file_name, :batch_id, :base_path, :variants_dir		
		attr_accessor :phenotype
		
		def set_directories(base_path, variants_dir)
				self.base_path = base_path
				self.variants_dir = variants_dir
		end
		
		def set_input_file_name(batch_id, sample_id, panel_version, phenotype)
			self.batch_id ||= batch_id
			self.sample_id ||= sample_id
			self.phenotype ||= phenotype
			if (self.batch_id != nil) && (self.sample_id != nil) && (self.phenotype != nil)
				self.input_file_name = "#{self.base_path}/#{batch_id}/#{self.variants_dir}/#{panel_version}_#{sample_id}_#{phenotype}.alamut.alltrans.txt"
			else
				raise ArgumentError "Missing sample information for Input file name: batch_id, sample_id or phenotype"
			end
		end
		
		def set_output_file_name(batch_id, sample_id, panel_version, phenotype)
			self.batch_id ||= batch_id
			self.sample_id ||= sample_id
			self.phenotype ||= phenotype
			if (self.batch_id != nil) && (self.sample_id != nil) && (self.phenotype != nil)
				self.output_file_name = "#{self.base_path}/#{batch_id}/#{self.variants_dir}/#{panel_version}_#{sample_id}_#{phenotype}.alamut.selected.txt"
			else
				raise ArgumentError "Missing sample information for Output file name: batch_id, sample_id or phenotype"
			end
		end
		
		def set_file_names(batch_id, sample_id, phenotype, panel_version, transcripts_file_name)
			self.set_input_file_name(batch_id, sample_id, panel_version, phenotype)
			self.set_output_file_name(batch_id, sample_id, panel_version, phenotype)
			self.transcripts_file_name = transcripts_file_name
		end
		
		def select_transcripts
						
					output_file = File.open(self.output_file_name, "w+")
					
					File.readlines(self.transcripts_file_name).each do |this_transcript|
						this_transcript.strip!
					
						File.readlines(self.input_file_name).each do |line|
							if line =~ /#{this_transcript}/
									output_file.puts(line)
							end
						end
						
					end
				
					output_file.close
		end

end
