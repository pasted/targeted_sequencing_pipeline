class SampleStore
		require_relative '../shared/variant_store'
		require_relative '../shared/transcript_store'
		
		attr_accessor :samples_by_phenotype, :samples
		
		def initialize(samples)
			self.samples = samples
  	end
  	
  	def parse_transcripts(parser, transcript_file_path)
  		
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
  	
  	def select_panel(this_batch, this_sample)

  			selected_panel = nil
  		  this_batch.panels.each do |this_panel|
  		  	if this_panel.panel_id == this_sample.panel_version
  					selected_panel = this_panel
  				end
  			end
  			return selected_panel
    end

  	def process_samples(parser, this_batch, this_book)

  		#Special gene symbols
  		#catch any unnanotated gene symbols - can be expanded to include genes without transcripts
  		#or if you want to catch all transcripts for a gene
  		special_gene_symbols = ["No alamut gene"]
  		
  		phenotype_store = Hash.new
  		all_transcripts = Hash.new
  		
  		self.samples.each do |this_sample|
  		  
  			#select the relevant panel for each sample based on the panel version
  			this_panel = select_panel(this_batch, this_sample)


  			#Load the required transcripts for the panel
  			transcript_store = parse_transcripts(parser, this_panel.transcript_file_path)
  			
  			if !all_transcripts.keys.include?(this_panel.panel_version)
  				all_transcripts[this_panel.panel_version] = transcript_store
  		  end

  			#Collapse the transcript_store object to given a list of transcript IDs
  			transcripts = transcript_store.list_all_transcripts
  			
  			#Load blacklisted variants
  			unwanted_variants = parser.parse_unwanted_variants(this_panel.unwanted_file_path)

  			#Load wanted regions
  			wanted_regions = parser.parse_wanted_regions(this_panel.wanted_file_path)
  			
  			full_sample_id 	= this_sample.capture_number
  			phenotype		   	= this_sample.phenotype
  			gender					= this_sample.gender.upcase
  			sample_id 			= full_sample_id.split('_')[1]
  			#puts "#{full_sample_id}"
  			#puts "PHENOTYPE :: #{phenotype}"
  			#puts "SAMPLE ID :: #{sample_id}"
  			
  			
  			if !phenotype_store.has_key?(phenotype)
  				sample_store_array = Array.new
  				phenotype_store["#{phenotype}"] = sample_store_array
  			end
  			
  			#Load variants from alltrans alamut file
  			puts "#{this_batch.base_path}/#{parser.batch_id}/#{this_panel.variants_directory}/#{this_panel.panel_id}_#{sample_id}_#{gender}_#{phenotype}.alamut.alltrans.txt"
  			variants = parser.parse_file("#{this_batch.base_path}/#{parser.batch_id}/#{this_panel.variants_directory}/#{this_panel.panel_id}_#{sample_id}_#{gender}_#{phenotype}.alamut.alltrans.txt", "#{this_panel.panel_id}_#{sample_id}")
  			variant_store = VariantStore.new(variants)

  			#puts "variants number :: #{variants.length}"
  			
  			selected_variants = Array.new
  			not_selected_variants = Array.new
  			
  			phenotype_store = variant_store.select_variants(sample_id, phenotype, phenotype_store, transcripts, unwanted_variants, wanted_regions, special_gene_symbols)
  			#puts phenotype_store.inspect
  		end#samples loop
  		
  		all_transcripts.each_pair do |panel_version, transcript_store|
  			#puts "transcripts to excel :: #{panel_version}"
  			this_book = transcript_store.transcripts_to_excel(this_book, panel_version)
  		end
  		self.samples_by_phenotype = phenotype_store
  		return this_book
  	end
  	
    def available_phenotypes
    	available_phenotypes = Array.new
    	self.samples_by_phenotype.each_pair {|phenotype,samples| available_phenotypes.push(phenotype)}
    	return available_phenotypes
    end
    
    def samples_to_excel(workbook)
    	
    	  self.samples_by_phenotype.each_pair do |this_phenotype, this_sample_store|
  	
    	  	this_sheet = workbook.create_worksheet :name => "#{this_phenotype}"
    	  	row_number = 0
    	  	
    	  	puts "Number of samples with #{this_phenotype} phenotype :: #{this_sample_store.length}"
  	
 					this_sample_store.each do |samples|
 						
 						samples.each_pair do |this_sample_id, this_variant_array|
 							header_array = Array.new
 							ex_number_array = self.samples.select {|s| s.sample_id == this_sample_id }
 							puts ex_number_array.inspect
 							this_sheet.row(row_number).push "SAMPLE ::", this_sample_id, "EX NUMBER:", ex_number_array.first.ex_number
 							row_number = row_number + 1
 							
 							#puts "******"
 							#puts this_variant_array[0].inspect
 							#puts "Unwanted variants"
 							#this_variant_array[1].each do |this_unwanted_variant|
 							#	puts "Chromosome :: #{this_unwanted_variant.chromosome}"
 							#	puts "Position :: #{this_unwanted_variant.position}"
 							#	puts "Gene :: #{this_unwanted_variant.gene}"
 							#	
 							#	
 							#end
 							
 							this_variant_array[0].each do |this_selected_variant|
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
 							row_number = row_number + 2
 							
 						end#samples
 					end#sample_store
    	end#samples_by_phenotype
    	
    	return workbook
    end

end
