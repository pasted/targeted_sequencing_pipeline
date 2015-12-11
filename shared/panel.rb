class Panel
	attr_accessor :batch_id, :panel_id, :panel_version, :base_path, :variants_directory, :intervals_directory
	attr_accessor :transcript_file_path, :unwanted_file_path, :wanted_file_path, :sample_list_path
	attr_accessor :coverage_cutoff, :allowed_genesets
end
