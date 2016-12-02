# @author Garan Jones
# @abstract SampleStore class: Object to help store and process Sample objects
class SampleStore
		require_relative 'variant_store'
		require_relative 'transcript_store'
		require_relative 'cnv_store'
		
		attr_accessor :samples_by_phenotype, :samples, :cnv_store
		
		def initialize(samples)
			self.samples = samples
			self.cnv_store = CnvStore.new
  	end
  	
    def capture_numbers
  		capture_number_array = self.samples.collect { |this_sample| this_sample.capture_number }
  		return capture_number_array
  	end
  	
  	
  		# @author Garan Jones
		  # Load in the list of transcript from a YAML file and construct a new TranscriptStore containing nested arrays of Phenotype and Gene objects
  		# @param parser [Object] A Parser instance
  		# @param transcript_file_path [String] The filepath to the yaml file containing the transcripts 
  		# @return [Object] TranscriptStore instance
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
  	
  		# @author Garan Jones
  		# Select the correct panel from available panels
  		# @param batch [Object] A Batch object containing details on the batch
  		# @param sample [Object] A Sample instance
  		# @return [Object] A Panel instance
  	def select_panel(batch, sample)

  		selected_panel = nil
  		batch.panels.each do |this_panel|
  		  	if this_panel.panel_id == sample.panel_version
  					selected_panel = this_panel
  			end
  		end
  		return selected_panel
    end

  	def process_samples(parser, this_batch, this_book, title_font)

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
  			unwanted_variants = parser.parse_unwanted_variants("#{this_batch.base_path}/#{this_batch.batch_id}/#{this_panel.unwanted_file_path}")

  			#Load wanted regions
  			wanted_regions = parser.parse_wanted_regions("#{this_batch.base_path}/#{this_batch.batch_id}/#{this_panel.wanted_file_path}")
  			
  			phenotype		= this_sample.phenotype
  			gender			= this_sample.gender.upcase

  			
  			if !phenotype_store.has_key?(phenotype)
  				sample_store_array = Array.new
  				phenotype_store["#{phenotype}"] = sample_store_array
  			end
  			
  			#Parse CNV files
  			if ( ["MALE","FEMALE"].include?(gender.upcase) ) && ( ["v501","v602","v603"].include?(this_panel.panel_id) )
  				cnv_array = parser.parse_cnvs("#{this_panel.panel_id}", "#{this_batch.base_path}/#{parser.batch_id}/cnv_analysis/#{this_sample.gender.downcase}/#{this_panel.panel_id.downcase}/results/Sample_#{this_panel.panel_id.downcase}_#{this_sample.ex_number}_#{gender}.realigned.bam.csv")
  				self.cnv_store.cnvs.store("#{this_sample.ex_number}", cnv_array)
  				self.cnv_store.process_cnvs(this_sample, this_batch, this_panel, parser)

  			else
  				puts "#{this_sample.capture_number} CNV file will not be parsed: check that gender(#{gender}) and panel (#{this_panel.panel_id}) meet criteria for running ExomeDepth"
  			end
  			

  			#Load variants from alltrans alamut file
  			variants = parser.parse_file("#{this_batch.base_path}/#{parser.batch_id}/#{this_panel.variants_directory}/#{this_panel.panel_id}_#{this_sample.ex_number}_#{gender}_#{phenotype}.alamut.alltrans.txt", "#{this_panel.panel_id.downcase}", "#{this_sample.ex_number.downcase}")
  			variant_store = VariantStore.new(variants)

  			selected_variants = Array.new
  			not_selected_variants = Array.new
  			
  			phenotype_store = variant_store.select_variants(this_sample.ex_number, phenotype, phenotype_store, transcripts, unwanted_variants, wanted_regions, special_gene_symbols)
  			
  		end#samples loop
  		
  		all_transcripts.each_pair do |panel_version, transcript_store|
  			this_book = transcript_store.transcripts_to_axlsx(this_book, panel_version)
  		end
  		self.samples_by_phenotype = phenotype_store
  		return self
  	end
  	
    def available_phenotypes
    	available_phenotypes = Array.new
    	self.samples_by_phenotype.each_pair {|phenotype,samples| available_phenotypes.push(phenotype)}
    	return available_phenotypes
    end
    
    def samples_to_axlsx(workbook, title_font)
    		spacer_font = workbook.styles.add_style :bg_color => "FF", :fg_color => "00", :sz => 12
  			not_selected_font = workbook.styles.add_style :bg_color => "CCCCCC", :fg_color => "00", :sz => 12
    	
    	  	self.samples_by_phenotype.each_pair do |this_phenotype, this_sample_store|
  	
    	  		#this_sheet = workbook.create_worksheet :name => "#{this_phenotype}"
    	  		workbook.add_worksheet(:name => "#{this_phenotype}") do |this_sheet|
    	  		#	row_number = 0
    	  	    	
    	  			puts "Number of samples with #{this_phenotype} phenotype :: #{this_sample_store.length}"
  	            	
 						this_sample_store.each do |samples|
 							
 							samples.each_pair do |this_ex_number, this_variant_array|
 								header_array = Array.new
 								cnv_header_array = Array.new
 								sample_array = self.samples.select {|s| s.ex_number == this_ex_number }
 								this_sample = sample_array.first
 								
 								this_sheet.add_row ["EX NUMBER:", "#{this_sample.ex_number}", "PANEL VERSION", "#{this_sample.panel_version}"], :style => title_font
 
 								this_sheet.add_row ["pct_target_bases_20x ::", "#{this_sample.pct_target_bases_20x}", "pct_target_bases_30x ::", "#{this_sample.pct_target_bases_30x}"], :style => title_font
                
 								if this_variant_array.length > 0
 									this_variant_array[0].each do |this_selected_variant|
 										if header_array.empty?
 											header_array = this_selected_variant.variable_order
 											header_array.map!{ |element| element.to_s }
 											this_sheet.add_row header_array
 										end
                	
 										#output variables to spreadsheet in given order
 										variant_array = Array.new
 										this_selected_variant.variable_order.map {|var|  variant_array.push( "#{this_selected_variant.send(var)}") }
 										this_sheet.add_row variant_array
 										
 									end#variants
 								else
 									this_sheet.add_row ["No variants passed selection criteria for this sample and profile."]
 								end
 								
 								this_sheet.add_row []
 								
 								this_sheet.add_row ["EX NUMBER:", "#{this_sample.ex_number}", "PANEL VERSION", "#{this_sample.panel_version}"], :style => title_font
 								
 								if self.cnv_store.cnvs.has_key?("#{this_sample.ex_number}")
 									this_cnv_array = self.cnv_store.cnvs.fetch("#{this_sample.ex_number}")
 									puts this_cnv_array.inspect
 									if this_cnv_array.length > 0
 										this_cnv_array.each do |this_selected_cnv|
 											if cnv_header_array.empty?
 												cnv_header_array = this_selected_cnv.variable_order
 												cnv_header_array.map!{ |element| element.to_s }
 												tmp_array = Array.new
 												#tmp_array.push("")
 												tmp_array.concat(cnv_header_array)
 												this_sheet.add_row tmp_array
 												#this_sheet.rows.last.cells[0].style = spacer_font

 											end
 											
 											#Output all cnvs for now
 												tmp_array = Array.new
 												#tmp_array.push("")
 												this_selected_cnv.variable_order.map {|var|  tmp_array.push  "#{this_selected_cnv.send(var)}" }
 												this_sheet.add_row tmp_array
 												#this_sheet.rows.last.cells[0].style = spacer_font
 											
 											#output variables to spreadsheet in given order
 											#if this_selected_cnv.wanted == true
 											#	tmp_array = Array.new
 											#	tmp_array.push("")
 											#	this_selected_cnv.variable_order.map {|var|  tmp_array.push  "#{this_selected_cnv.send(var)}" }
 											#	this_sheet.add_row tmp_array
 											#	this_sheet.rows.last.cells[0].style = spacer_font
 											#else
 											#	tmp_array = Array.new
 											#	tmp_array.push("")
 											#	this_selected_cnv.variable_order.map {|var|   tmp_array.push  "#{this_selected_cnv.send(var)}" }
 											#	this_sheet.add_row tmp_array, :style => not_selected_font
 											#	this_sheet.rows.last.cells[0].style = spacer_font
 											#end
 											
 										end
 									else
 										this_sheet.add_row ["","Sample had no CNVs within profile intervals"]
 									end
 								else
 									this_sheet.add_row ["","ExomeDepth not run for this sample"]
 								end
 								
 								this_sheet.add_row []
 								this_sheet.add_row []
 							end#samples
 						end#sample_store
 						
 					end#worksheet
 				end#samples_by_phenotype

    	return workbook
    end
    

end
