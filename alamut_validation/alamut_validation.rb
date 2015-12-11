class AlamutValidation
	require 'awesome_print'
  require 'yaml'
  require 'csv'
  require 'smarter_csv'
  require 'shellwords'
  require 'spreadsheet'
  require_relative '../shared/variant'
  require_relative '../shared/variant_store'
  require_relative '../shared/transcript_store'
  require_relative '../shared/batch'
  require_relative '../shared/panel'
  require_relative '../shared/gene'
  require_relative '../shared/phenotype'
  require_relative '../alamut_parser/alamut_parser'

  
  attr_accessor :batch_id, :base_path, :variants_directory, :intervals_directory, :transcript_file_path 
  attr_accessor :unwanted_file_path, :wanted_file_path, :sample_list_path, :panel_id, :panel_version
  
  def self.bash(command)
  	escaped_command = Shellwords.escape(command)
  	system "bash -c #{escaped_command}"
  end

  def self.parse_file(file_name, sample_id_lc)
  	options = { :col_sep => "\t" }
  	variant_array = Array.new 
  	puts file_name
  	if File.exists?(file_name) && ( File.stat(file_name).size > 0 )           
  		SmarterCSV.process( file_name, options ) do |csv|
  			this_variant = Variant.new
  			this_variant.chromosome 					= csv.first[:chrom]
  			this_variant.position							= csv.first[:pos]
  			this_variant.gene 								= csv.first[:gene]
  			this_variant.genotype							=	csv.first[:"gt_(#{sample_id_lc})"]
  			this_variant.transcript 					= csv.first[:transcript]
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
  			
  			variant_array.push(this_variant)
    	
  		end
  	else
  		puts "ERROR :: variants file has no content :: SAMPLE ID : #{sample_id_lc}"
  		if File.exists?
  			puts "File exists? #{File.exists?(file_name)}"
  			puts "File size? #{File.stat(file_name).size}"
  		else
  			puts "File exists? #{File.exists?(file_name)}"
  		end
  	end
  	return variant_array
  end
  
  
  #run control VCF through Alamut
  def self.run_alamut(run_type, alamut_path, panel_version, hgmd_user, hgmd_pass) 
  	bash("#{alamut_path} \
		--hgmdUser #{hgmd_user} --hgmdPasswd #{hgmd_pass} \
		--in controls/#{panel_version}_control.vcf  \
		--ann results/#{panel_version}_#{run_type}.alamut.alltrans.txt \
		--unann results/#{panel_version}_#{run_type}.alamut.alltrans.unannotated.txt \
		--alltrans \
		--ssIntronicRange 2 \
		--outputVCFInfo AC AF AN DP FS MQ MQ0 QD \
		--outputVCFGenotypeData AD DP GQ GT PL \
		--outputVCFQuality --outputVCFFilter")
		
		
  end
  
  def self.parse_transcripts(parser, transcript_file_path)
  		
  		  phenotype_list = parser.load_phenotypes(transcript_file_path)
			  phenotype_array = Array.new
			 
			  phenotype_list.each_pair do |phenotype, genes|
			  	gene_list = Array.new
			  	genes.each_pair do |gene_symbol, transcript|
			  		this_gene = Gene.new(gene_symbol, transcript)
			  		gene_list.push(this_gene)
			  	end
			  	this_phenotype = Phenotype.new(phenotype, gene_list)
			  	phenotype_array.push(this_phenotype) 	
			  end
			  
			  transcript_store = TranscriptStore.new(phenotype_array)
			  return transcript_store
  end
  
  def build_spreadsheet(selected_current_variants, sample_name)
  	#writes out selected variants data structure to an excel spreadsheet		
  	this_book = Spreadsheet::Workbook.new
  					
  	this_sheet = this_book.create_worksheet :name => "#{sample_name}"
  	row_number = 0
  	
  	#selected_current_variants["#{sample_name}"].first.first[1][0] equivalent to variants that have passed the selection criteris
  	#including transcript id
  	header_array = Array.new
 		this_sheet.row(row_number).push "SAMPLE ::", "#{sample_name}"
 		row_number = row_number + 1
 		
 		#write out annotated variants
 		selected_current_variants["#{sample_name}"].first.first[1][0].each do |this_selected_variant|
 			if header_array.empty?
 				header_array = this_selected_variant.variable_order
 				header_array.map!{ |element| element.to_s }
 				
 				header_array.each do |this_header|
 					  this_sheet.row(row_number).push this_header  
 				end
 			end
 			row_number = row_number + 1
 			
 			#output variables to spreadsheet in given order
 			this_selected_variant.variable_order.map {|var|  this_sheet.row(row_number).push  "#{this_selected_variant.send(var)}" }
 							
 		end#variants
 				
  			this_book.write "#{this_batch.base_path}/#{this_batch.batch_id}/results/#{this_batch.batch_id}_#{sample_name}_variants.xls"
  			return this_book
  end
  

  
  #Load the config YAML file and pass the settings to local variables
  this_batch = YAML.load_file('../configuration/config.yaml')
  alamut_path		=	this_batch.alamut_path
  panel_version	=	"v5"
  hgmd_user			=	this_batch.hgmd_user
  hgmd_pass			=	this_batch.hgmd_pass
  
  
  #run_alamut("control", alamut_path, panel_version, hgmd_user, hgmd_pass)
  
  #load control VCF Alamut output
  control_variants = parse_file("#{this_batch.base_path}/#{this_batch.batch_id}/scripts/alamut_validation/results/#{panel_version}_control.alamut.alltrans.txt", "control")
  
  #generate new alamut outputs from control VCF
  #run_alamut("current", alamut_path, panel_version, hgmd_user, hgmd_pass)
  
  #load current VCF Alamut output
  current_variants = parse_file("#{this_batch.base_path}/#{this_batch.batch_id}/scripts/alamut_validation/results/#{panel_version}_current.alamut.alltrans.txt", "current")
  
  			#select the panel
  			selected_panel = nil
  			this_batch.panels.each do |this_panel|
  			 if this_panel.panel_id == panel_version
  				selected_panel = this_panel
  			 end
  			end
  			
  			#Replicate the Alamut parser as much as possible
  			parser = AlamutParser.new
    		#Load the required transcripts for the panel
  			transcript_store = self.parse_transcripts(parser, selected_panel.transcript_file_path)
  			
  			all_transcripts = Hash.new
  			all_transcripts[selected_panel.panel_version] = transcript_store
  		 

  			#Collapse the transcript_store object to given a list of transcript IDs
  			transcripts = transcript_store.list_all_transcripts
  			
  			#Load blacklisted variants
  			unwanted_variants = parser.parse_unwanted_variants(selected_panel.unwanted_file_path)
  			#Load wanted regions
  			wanted_regions = parser.parse_wanted_regions(selected_panel.wanted_file_path)
  			special_gene_symbols = ["No alamut gene"]
  			
  			current_variants_store = VariantStore.new(current_variants)
  			selected_current_variants = Hash.new
  			tmp_sample_store = []
  			selected_current_variants["control"] = tmp_sample_store
  			selected_current_variants = current_variants_store.select_variants("current", "control", selected_current_variants, transcripts, unwanted_variants, wanted_regions, special_gene_symbols)
  			

  			puts selected_current_variants["control"].first.first[1][0].length
  				
  			puts selected_current_variants["control"].first.first[1][1].length
  			
  			
  			ap selected_current_variants["control"].first.first[1][0].first
  			#puts selected_current_variants["control"].first.first[1][1].first.inspect
  			
  			CSV.open('control_alamut_variants.csv', 'r') do |csv_obj|
  				header_array = []
  				selected_current_variants["control"].first.first[1][0].each do |variant|
  					if header_array.empty?
  						header_array = variant.variable_order
  						header_array.map!{ |element| element.to_s }
 				
  						csv_obj << header_array
  					end
  					variant_array = []
  					variant.variable_order.map {|var|  variant_array.push("#{variant.send(var)}") }
  					csv_obj << variant_array
  				end
  			end


 #compare
  puts "Control variants :: #{control_variants.length}"
  puts "Current variants :: #{current_variants.length}"
  
  unmatched_variants =[]
  matched_variants = []
  invalid_variants = []
  
  current_variants.each_with_index do |current_variant, current_index|
  	if current_variant.chromosome && current_variant.position
  		control_variants.each_with_index do |control_variant, control_index|
  			if current_variant.chromosome == control_variant.chromosome
  				if (current_variant.position == control_variant.position) && (current_variant.transcript == control_variant.transcript)
  					
  					matched_variants.push(current_variant)
  					current_variants.delete_at(current_index)
  				end
  			end
  		end
  	else
  		#No chromosome id and / or genomic position so invalid
  		invalid_variants.push(current_variant)
  		current_variants.delete_at(current_index)
  	end
  end

  puts matched_variants.length
  puts invalid_variants.length
  puts current_variants.length
  #note differences that would effect the ruleset
  
  CSV.open('matched_variants.csv', 'w') do |csv|
  	matched_variants.each do |variant|
  		puts variant.methods.sort
  	end
  end
  
end
