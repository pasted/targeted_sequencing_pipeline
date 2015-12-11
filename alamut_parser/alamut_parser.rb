class AlamutParser
  require 'yaml'
  require 'smarter_csv'
  require_relative 'interval'
  require_relative 'region'
  require_relative '../shared/gene'
  require_relative '../shared/phenotype'
  require_relative '../shared/transcript_store'
  require_relative 'sample_store'
  require_relative '../shared/variant'
  require_relative '../shared/sample'
  require_relative '../shared/batch'
  require_relative '../shared/panel'
  require_relative '../shared/metric'
  require 'spreadsheet'
  
  attr_accessor :batch_id, :base_path, :variants_directory, :intervals_directory, :transcript_file_path 
  attr_accessor :unwanted_file_path, :wanted_file_path, :sample_list_path, :panel_id, :panel_version
  
  def parse_file(file_name, sample_id_lc)
  	options = { :col_sep => "\t" }
  	variant_array = Array.new 
  	
  	if File.exists?(file_name) && ( File.stat(file_name).size > 0 )
  	
  		SmarterCSV.process( file_name, options ) do |csv|
  			this_variant = Variant.new
  			#puts csv.first.inspect
  			
  			this_variant.chromosome 					= csv.first[:chrom]
  			this_variant.position							= csv.first[:pos]
  			this_variant.genomic_dna_start		= csv.first[:gdnastart]
  			this_variant.genomic_dna_end			= csv.first[:gdnaend]
  			this_variant.gene 								= csv.first[:gene]
  			this_variant.genotype							=	csv.first[:"gt_(#{sample_id_lc})"]
  			this_variant.full_transcript 			= csv.first[:transcript]
  			this_variant.var_type							= csv.first[:vartype]
  			this_variant.coding_effect				= csv.first[:codingeffect]
  			this_variant.var_location 				= csv.first[:varlocation]
  			this_variant.genomic_nomen 				= csv.first[:gnomen]
  			this_variant.cdna_nomen	 					= csv.first[:cnomen]
  			this_variant.protein_nomen 				= csv.first[:pnomen]
  			this_variant.exon									= csv.first[:exon]
  			this_variant.intron								= csv.first[:intron]
  			this_variant.distance_nearest_splice_site 	= csv.first[:distnearestss]
  			this_variant.nearest_splice_site_type 			= csv.first[:nearestsstype]
  			this_variant.nearest_splice_site_change 		= csv.first[:nearestsschange]
  			this_variant.local_splice_effect						= csv.first[:localspliceeffect]
  			this_variant.rs_id 									= csv.first[:rsid]
  			this_variant.rs_validated 					= csv.first[:rsvalidated]
  			this_variant.rs_maf									= csv.first[:rsmaf]
  			this_variant.genomes_1000_freq			= csv.first[:"1000g_AF"]
  			this_variant.genomes_1000_afr_freq	= csv.first[:"1000g_AFR_AF"]
  			this_variant.genomes_1000_sas_freq  = csv.first[:"1000g_SAS_AF"]
  			this_variant.genomes_1000_eas_freq  = csv.first[:"1000g_EAS_AF"]
  			this_variant.genomes_1000_eur_freq  = csv.first[:"1000g_EUR_AF"]
  			this_variant.genomes_1000_amr_freq  = csv.first[:"1000g_AMR_AF"]
  			this_variant.exac_all_freq       		= csv.first[:exacallfreq]
  			this_variant.exac_afr_freq       	  = csv.first[:exacafrfreq]
  			this_variant.exac_amr_freq       	  = csv.first[:exacamrfreq]
  			this_variant.exac_eas_freq       	  = csv.first[:exaceasfreq]
  			this_variant.exac_sas_freq       	  = csv.first[:exacsasfreq]
  			this_variant.exac_nfe_freq       	  = csv.first[:exacnfefreq]
  			this_variant.exac_fin_freq       	  = csv.first[:exacfinfreq]
  			this_variant.exac_oth_freq       	  = csv.first[:exacothfreq]
  			this_variant.exac_afr_hmz       	  = csv.first[:exacafrhmz]
  			this_variant.exac_amr_hmz        	  = csv.first[:exacamrhmz]
  			this_variant.exac_eas_hmz           = csv.first[:exaceashmz]
  			this_variant.exac_sas_hmz        	  = csv.first[:exacsashmz]
  			this_variant.exac_nfe_hmz        	  = csv.first[:exacnfehmz]
  			this_variant.exac_fin_hmz        	  = csv.first[:exacfinhmz]
  			this_variant.exac_oth_hmz        	  = csv.first[:exacothhmz]
  			this_variant.exac_filter         	  = csv.first[:exacfilter]
  			this_variant.exac_read_depth        = csv.first[:exacreaddepth]
  			this_variant.clin_var_ids           = csv.first[:clinvarids]
  			this_variant.clin_var_origins    	  = csv.first[:clinvarorigins]
  			this_variant.clin_var_methods    	  = csv.first[:clinvarmethods]
  			this_variant.clin_var_clin_signifs	= csv.first[:clinvarclinsignifs]
  			this_variant.clin_var_review_status	= csv.first[:clinvarreviewstatus]
  			this_variant.clin_var_phenotypes	  = csv.first[:clinvarphenotypes]
  			this_variant.cosmic_ids         	  = csv.first[:cosmicids]
  			this_variant.cosmic_tissues      	  = csv.first[:cosmictissues]
  			this_variant.cosmic_freqs           = csv.first[:cosmicfreqs]
  			this_variant.cosmic_sample_counts	  = csv.first[:cosmicsamplecounts]
  			
  			this_variant.esp_all_maf						= csv.first[:espallmaf]
  			this_variant.hgmd_phenotype 				= csv.first[:hgmdphenotype]
  			this_variant.hgmd_pub_med_id 				= csv.first[:hgmdpubmedid]
  			this_variant.hgmd_sub_category 			= csv.first[:hgmdsubcategory]
  			this_variant.n_orthos								= csv.first[:northos]
  			this_variant.conserved_orthos				= csv.first[:conservedorthos]
  			this_variant.conserved_dist_species	= csv.first[:conserveddistspecies]
  			this_variant.grantham_dist					= csv.first[:granthamdist]
  			this_variant.agv_gd_class						= csv.first[:agvgdclass]
  			this_variant.sift_prediction				= csv.first[:siftprediction]
  			this_variant.wt_nuc									=	csv.first[:wtnuc]
  			this_variant.var_nuc								=	csv.first[:varnuc]
  			this_variant.ins_nucs								=	csv.first[:insnucs]
  			this_variant.del_nucs								= csv.first[:delnucs]
  			this_variant.filter_vcf							= csv.first[:"filter_(vcf)"]
  			
  			this_variant.ad											= csv.first[:"ad_(#{sample_id_lc})"]
  			this_variant.dp											= csv.first[:"dp_(#{sample_id_lc})"]
  			this_variant.gq											= csv.first[:"gq_(#{sample_id_lc})"]
  			this_variant.gt											= csv.first[:"gt_(#{sample_id_lc})"]
  			this_variant.pl											= csv.first[:"pl_(#{sample_id_lc})"]
  			
  			this_variant.parse_transcript()
  			variant_array.push(this_variant)
  			#puts this_variant.inspect
  		end
  	else
  		puts "ERROR :: variants file has no content :: SAMPLE ID : #{sample_id_lc}"
  		if File.exists?(file_name)
  			puts "File exists? #{File.exists?(file_name)}"
  			puts "File size? #{File.stat(file_name).size}"
  		else
  			puts "File exists? #{File.exists?(file_name)}"
  		end
  	end
  	return variant_array
  end
  
  def parse_sample_list(sample_file_path)
  	options = { :col_sep => "," }
  	sample_array = Array.new
  	
  	SmarterCSV.process( sample_file_path, options ) do |csv|
  		this_sample = Sample.new
  		this_sample.capture_number 	= csv.first[:capture_number]
  		this_sample.mody_number			= csv.first[:mody_number]
  		this_sample.ex_number				= csv.first[:ex_number]
  		this_sample.gender					= csv.first[:gender]
  		this_sample.profile					= csv.first[:profile]
  		this_sample.phenotype				= csv.first[:phenotype]
  		this_sample.sample_type			= csv.first[:sample_type]
  		this_sample.comment					= csv.first[:comment]
  		this_sample.parse_panel_version
  		this_sample.parse_sample_id
  		if this_sample.phenotype == nil
  			this_sample.parse_phenotype
  		end
  		
  		sample_array.push(this_sample)

  	end
  	
  	return sample_array
  end
  
  
  def parse_unwanted_variants(unwanted_file_path)
  	options = { :col_sep => "," }
  	unwanted_array = Array.new
  	
  	SmarterCSV.process( unwanted_file_path, options ) do |csv|
  		this_unwanted_variant = Variant.new

  		this_unwanted_variant.chromosome 						= csv.first[:chrom]
  		this_unwanted_variant.position							= csv.first[:start_pos]
  		this_unwanted_variant.gene 									= csv.first[:gene]		
  		this_unwanted_variant.wt_nuc								=	csv.first[:wtnuc]
  		this_unwanted_variant.var_nuc								=	csv.first[:varnuc]
  		this_unwanted_variant.ins_nucs							=	csv.first[:insnucs]
  		this_unwanted_variant.del_nucs							= csv.first[:delnucs]
  		this_unwanted_variant.reason_for_filtering 	= csv.first[:reason_for_filtering]
  		
  		unwanted_array.push(this_unwanted_variant)
  		#puts csv.first.inspect
  	end
  	
  	return unwanted_array
	end
	
	def parse_wanted_regions(wanted_file_path)
  	options = { :col_sep => "," }
  	wanted_array = Array.new
  	
  	SmarterCSV.process( wanted_file_path, options ) do |csv|
  		this_wanted_region = Region.new

  		this_wanted_region.chromosome 			= csv.first[:chromosome]
  		this_wanted_region.start_position		= csv.first[:start_position]
  		this_wanted_region.stop_position 		= csv.first[:stop_position]		

  		
  		wanted_array.push(this_wanted_region)
  		#puts csv.first.inspect
  	end
  	
  	return wanted_array
	end
	
	def parse_intervals(intervals_file_path)
  	options = { :col_sep => "\t" }
  	interval_array = Array.new
  	
  	SmarterCSV.process( intervals_file_path, options ) do |csv|
  		this_interval = Interval.new

  		this_interval.chromosome 			= csv.first[:chromosome]
  		this_interval.genomic_start		= csv.first[:genomic_start]
  		this_interval.genomic_end 		= csv.first[:genomic_end]		
  		this_interval.strand 					= csv.first[:strand]
  		this_interval.interval_name 	= csv.first[:interval_name]	
  		
  		interval_array.push(this_interval)
  		
  	end
  	
  	return interval_array
	end
	
	def load_phenotypes(transcript_file_path)
		phenotype_list = YAML.load_file("#{transcript_file_path}")
		return phenotype_list
	end
	
	def parse_metrics(input_file_path, metrics_line)
			counter = 0
			this_metric = Metric.new
			IO.foreach("#{input_file_path}") do |this_line| 
			
				if ( this_line.match /(^\r|^\n|^\r\n|^#)/ )
					counter = counter + 1
				else
					if counter == metrics_line.to_i
						
						metrics_array = this_line.split("\t")
						
						this_metric.variable_order.each_with_index do |var, index|
							file_var = metrics_array[index]
							if (file_var == "\n") || (file_var == "")
								file_var = nil 
							end
							this_metric.instance_variable_set("@#{var}", file_var)
						end
					end
					counter = counter + 1
				end
			end
			return this_metric
	end
	
	def populate_metrics(samples, this_batch, this_parser)
		processed_samples = Array.new
		samples.each do |this_sample|
			phenotype_metric_file_path = "#{this_batch.base_path}/#{this_batch.batch_id}/metrics/#{this_sample.panel_version}_#{this_sample.sample_id}_#{this_sample.gender.upcase!}_#{this_sample.phenotype}.phenotype.bait_capture_metrics"
			this_metric = this_parser.parse_metrics(phenotype_metric_file_path, this_batch.phenotype_metric_line)
			this_sample.add_metrics(Array.new.push(this_metric))
			processed_samples.push(this_sample)
		end
		return processed_samples
	end
  
  #Load the config YAML file and pass the settings to local variables
  this_batch = YAML.load_file('../configuration/config.yaml')
  
  #Init AlamutParser class
  parser = AlamutParser.new()
  
  parser.batch_id = this_batch.batch_id
  parser.base_path = this_batch.base_path
  parser.sample_list_path = this_batch.sample_list_path


  samples = parser.parse_sample_list(parser.sample_list_path)
  
  #populate phenotype specific metrics for each sample
  samples = parser.populate_metrics(samples, this_batch, parser)
  
  sample_store = SampleStore.new(samples)
  
  #build the workbook
  this_book = Spreadsheet::Workbook.new
  
  this_book = sample_store.process_samples(parser, this_batch, this_book)
 
  this_book = sample_store.samples_to_excel(this_book)
	

	#write out intervals that have been analysed by phenotype
	
	available_phenotypes = sample_store.available_phenotypes
	#puts available_phenotypes.inspect


	#available_phenotypes.each do |phenotype|
	#	intervals = parser.parse_intervals("#{parser.intervals_directory}/#{parser.panel_version}_#{phenotype}_diagnostic_ROI.tsv")
	#	this_sheet = this_book.create_worksheet :name => "#{phenotype} intervals"
	#	row_number = 0
	#	intervals.each do |this_interval|
	#		this_sheet.row(row_number).push "#{this_interval.chromosome}","#{this_interval.genomic_start}","#{this_interval.genomic_end}","#{this_interval.strand}","#{this_interval.interval_name}"
 	#		row_number = row_number + 1
	#	end
	#	
	#end
	
	this_book.write "#{parser.base_path}/#{parser.batch_id}/results/#{parser.batch_id}_variants.#{Time.now.strftime("%d-%m-%Y-%H%M%S")}.xls"
	
  
end#end AlamutParser class
