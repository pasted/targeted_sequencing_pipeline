class SnpParser
  require 'yaml'
  require 'json'
  require 'smarter_csv'
  require 'spreadsheet'
  require 'awesome_print'
  
  require_relative '../shared/sample_store'
  require_relative '../shared/sample'
  require_relative '../shared/batch'
  require_relative '../shared/panel'
  require_relative 'snp'
  require_relative 'snp_store'
  require_relative 'allele'
  require_relative 'rank'

  
  attr_accessor :batch_id, :base_path, :variants_directory, :intervals_directory, :transcript_file_path 
  attr_accessor :unwanted_file_path, :wanted_file_path, :sample_list_path, :panel_id, :panel_version
  
  #Load the SNP data with genotypes from a TSV file
  #
  # @param file_name [String] A string of the filepath to the file
  # @capture_numbers [String] An array of capture ids
  # @return snp_array [Array] An array of Snp objects
  def parse_genotype_table_file(file_name, capture_numbers)
  	options = { :col_sep => "\t", :headers_in_file => true }
  	snp_array = Array.new 
  	
  	if File.exists?(file_name) && ( File.stat(file_name).size > 0 )
  		
  		SmarterCSV.process( file_name, options ) do |csv|
  			this_snp = Snp.new
  			this_snp.id 			= csv.first[:id]
  			this_snp.chrom 		= csv.first[:chrom]
  			this_snp.pos 			= csv.first[:pos]
  			this_snp.ref	 		= csv.first[:ref]
  			this_snp.alt			= csv.first[:alt]
  			this_snp.qual 		= csv.first[:qual]
  			this_snp.filter 	= csv.first[:filter]
  			this_snp.ac 			= csv.first[:ac]
  			this_snp.af 			= csv.first[:af]
  			
  			capture_numbers.each do |this_capture_number|
  				this_allele = Allele.new
  				this_allele.sample_id 	= this_capture_number
  				this_allele.genotype 	= csv.first[:"#{this_capture_number}.gt"]
  				this_snp.alleles.push(this_allele)
  			end
  			snp_array.push(this_snp)
  		end
  	else
  		
  		if File.exists?(file_name)
  			puts "File exists? #{File.exists?(file_name)}"
  			puts "File size? #{File.stat(file_name).size}"
  		else
  			puts "File exists? #{File.exists?(file_name)}"
  		end
  	end
  	return snp_array
  end
  
  #Load the sample list from a CSV file
  #
  # @param sample_file_path [String] A string of the filepath to the file
  # @return sample_array [Array] An array of Sample objects
  def parse_sample_list(sample_file_path)
  	options = { :col_sep => "," }
  	sample_array = Array.new
  	
  	SmarterCSV.process( sample_file_path, options ) do |csv|
  		this_sample = Sample.new
  		this_sample.capture_number 	= csv.first[:capture_number]
  		this_sample.mody_number		= csv.first[:mody_number]
  		this_sample.ex_number		= csv.first[:ex_number]
  		this_sample.gender			= csv.first[:gender]
  		this_sample.profile			= csv.first[:profile]
  		this_sample.phenotype		= csv.first[:phenotype]
  		this_sample.sample_type		= csv.first[:sample_type]
  		this_sample.comment			= csv.first[:comment]
  		this_sample.parse_panel_version
  		this_sample.parse_sample_id
  		if this_sample.phenotype == nil
  			this_sample.parse_phenotype
  		end
  		
  		sample_array.push(this_sample)

  	end
  	
  	return sample_array
  end
  
  #Load the ranked scores from a CSV file
  #
  # @param rank_file_path [String] A string of the filepath to the file
  # @return rank_array [Array] An array of Rank objects
  def load_ranked_scores(rank_file_path)
  	options = { :col_sep => "," }
  	rank_array = Array.new
  	this_order = 1
  	
  	SmarterCSV.process( rank_file_path, options ) do |csv|
  		this_rank = Rank.new
  		this_rank.order						= this_order
  		this_rank.score 					= csv.first[:score]
  		this_rank.percentage			= csv.first[:percentage]
  		
  		this_order += 1
  		
  		rank_array.push(this_rank)

  	end
  	
  	return rank_array
  end
  
  #Loop through results and write the sample_id, total score, T1D Percentage, Non-T1D Percentage
  #out to the Excel worksheet
  #
  # @param snps_hash [Hash] A hash containing the results
  # @param this_book [Object] A Spreadsheet::Workbook object
  # @param worksheet_name [String] A string with the name of the Workbook tab
  # @return [Object] A Spreadsheet::Workbook object with the results
  def write_snps_scores(snps_hash, this_book, worksheet_name)
		this_sheet = this_book.create_worksheet :name => "#{worksheet_name}"
		
		row_number = 0
		this_sheet.row(row_number).push "Sample Id", "Score", "T1D Percentage", "Non-T1D Percentage"
		row_number = row_number + 1
		snps_hash.each_pair do |this_key, this_score|		
			this_sheet.row(row_number).push "#{this_key}", "#{this_score[0]}", "#{this_score[1]}", "#{this_score[2]}"
			row_number = row_number + 1
		end

		return this_book
	end
	
	#Loop through the samples, hash key, this_key contains the sample_id
  #this_array contains the combined score for this sample in element 0 OR the string "ERROR"
  #this_array will contain the T1D percentage rank in element 1 OR the string "NA"
  #this_array will contain the Non-T1D percentage rank in element 2 OR the string "NA"
  #
  # @param sample_hash [Hash] A hash containing the Sample objects to be processed
  # @param t1d_ranked_scores [Hash] A hash containing the T1D ranked scores by percentage from the CSV file
  # @param nont1d_ranked_scores [Hash] A hash containing the Non-T1D ranked scores by percentage from the CSV file
  # @return [Hash] processed sample_hash with the updated T1D percentages and Non-T1D percentages
	def process_samples(sample_hash, t1d_ranked_scores, non_t1d_ranked_scores)
  	
  	sample_hash.each do |this_key, this_array|
  		this_sample_score = this_array[0]
  		current_highest_score = 0
  		
  		if this_sample_score != "ERROR"
  			t1d_ranked_scores.each do |this_rank|
  				
  				if this_rank.score < this_sample_score
  					current_highest_score = this_rank.score
  				else
  					this_array[1] = this_rank.percentage
  				end
  			end
  			non_t1d_ranked_scores.each do |this_rank|
  				if this_rank.score < this_sample_score
  					current_highest_score = this_rank.score
  				else
  					this_array[2] = this_rank.percentage
  				end
  			end
  		else
  			this_array[1] = "NA"
  			this_array[2] = "NA"
  		end
  		
  	end
  	
  	return sample_hash
  end
  
  #Load the config YAML file and pass the settings to local variables
  this_batch = YAML.load_file('../configuration/config.yaml')
  
  #Init AlamutParser class
  parser = SnpParser.new()
  
  parser.batch_id = this_batch.batch_id
  parser.base_path = this_batch.base_path
  parser.sample_list_path = this_batch.sample_list_path

  samples = parser.parse_sample_list(parser.sample_list_path)

  sample_store = SampleStore.new(samples)
  
  capture_numbers = sample_store.capture_numbers.collect {|capture_number| capture_number.match("v501") ? capture_number : nil}
  capture_numbers.compact!
  

  scores = YAML.load_file("t1d_snps_score.yaml")
  t1d_ranked_scores = parser.load_ranked_scores("t1d_ranked_scores.csv")
  non_t1d_ranked_scores = parser.load_ranked_scores("non-t1d_ranked_scores.csv")
  
  
  snp_array = parser.parse_genotype_table_file("#{this_batch.base_path}/#{this_batch.batch_id}/variants_t1d/#{this_batch.batch_id}.t1d.GATK-#{this_batch.gatk_version}.GT.table", capture_numbers)
  
  sample_hash = Hash[capture_numbers.map {|x| [:"#{x}", [0.0,0]]}]
  
  snp_store = SnpStore.new(snp_array)
  snp_store.process_merged_snps(["rs2187668","rs7454108"])
  
  snp_array.each do |this_snp|
  	#Collect the scores
  	#Find the relevant score row

  	scores.each_pair do |this_id, this_score|
  		if this_id == this_snp.id
  			
  			#Loop through the alleles and pull off the genotypes
  			this_snp.alleles.each do |this_allele|
  				begin
  					this_genotype_score = this_score.fetch(:"#{this_allele.genotype}")
  					current_score = sample_hash.fetch(:"#{this_allele.sample_id}")
  					if current_score[0] != "ERROR"
  						current_score[0] += this_genotype_score.to_f
  						sample_hash.store(:"#{this_allele.sample_id}", current_score)
  					else
  						puts "#{this_allele.sample_id} has an error."
  					end
  				rescue
  					puts "Sample #{this_allele.sample_id} has an unknown genotype #{this_allele.genotype}"
  					current_score = sample_hash.fetch(:"#{this_allele.sample_id}")
  					current_score[0] = "ERROR"
  					sample_hash.store(:"#{this_allele.sample_id}", current_score)
  					puts "#{sample_hash[:"#{this_allele.sample_id}"]}"
  				end#begin..rescue
  				
  			end#alleles loop
  			
  		end#if clause
  	end#scores loop
  	
  end#snp_array loop
  

  
  sample_hash = parser.process_samples(sample_hash, t1d_ranked_scores, non_t1d_ranked_scores)
  
  puts JSON.pretty_generate(sample_hash)
  
  this_book = Spreadsheet::Workbook.new
  		
  this_book = parser.write_snps_scores(sample_hash, this_book, "T1D scores")
  
  this_book.write "#{this_batch.base_path}/#{this_batch.batch_id}/results/#{this_batch.batch_id}_T1D_snps.xls"
  
end#end SnpParser class
