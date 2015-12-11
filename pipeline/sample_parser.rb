class SampleParser
  require 'yaml'
  require 'smarter_csv'
  require_relative '../shared/sample'
  
  attr_accessor :batch_id, :ftp_url, :base_path, :variants_directory, :intervals_directory, :transcript_file_path 
  attr_accessor :unwanted_file_path, :wanted_file_path, :sample_list_path, :panel_id, :panel_version
  
  def parse_sample_list(sample_file_path)
  	options = { :col_sep => "," }
  	sample_array = Array.new
  	
  	SmarterCSV.process( sample_file_path, options ) do |csv|
  		this_sample = Sample.new
  		this_sample.sequencing_panel_version 	= csv.first[:sequencing_panel_version]
  		this_sample.capture_number 						= csv.first[:capture_number]
  		this_sample.mody_number								= csv.first[:mody_number]
  		this_sample.ex_number									= csv.first[:ex_number]
  		this_sample.gender										= csv.first[:gender]
  		this_sample.profile										= csv.first[:profile]
  		this_sample.phenotype									= csv.first[:phenotype]
  		this_sample.sample_type								= csv.first[:sample_type]
  		this_sample.comment										= csv.first[:comment]
  		this_sample.parse_panel_version
  		this_sample.parse_sample_id
  		if this_sample.phenotype == nil
  			this_sample.parse_phenotype
  		end
  		
  		sample_array.push(this_sample)

  	end
  	
  	return sample_array
  end
  
  
  
end
