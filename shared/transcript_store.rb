class TranscriptStore
	attr_accessor :phenotypes
		
	def initialize(phenotypes)
		self.phenotypes = phenotypes
  	end
  	
  	def	list_all_transcripts
  		available_transcripts = Array.new
  		self.phenotypes.map{|this_phenotype| available_transcripts.concat(this_phenotype.all_transcripts)}
  		available_transcripts.uniq!
  		return available_transcripts
    end
  	
  	def transcripts_to_axlsx(workbook, panel_version)
  		
  		this_sheet = workbook.add_worksheet(:name => "#{panel_version} Transcripts")
  		#case panel_version
  		#when "V501"
  		#	title_font = workbook.styles.add_style :bg_color => "5CA56B", :fg_color => "FF", :sz => 12
  		#when "V5"
  		#	title_font = workbook.styles.add_style :bg_color => "4B2583", :fg_color => "FF", :sz => 12
  		#when "V602"
  		#	title_font = workbook.styles.add_style :bg_color => "7D8325", :fg_color => "FF", :sz => 12
  		#else
  		#	title_font = workbook.styles.add_style :bg_color => "727272", :fg_color => "FF", :sz => 12
  		#end
  		title_font = workbook.styles.add_style :bg_color => "5CA56B", :fg_color => "FF", :sz => 12
  		
    	self.phenotypes.each do |this_phenotype|
    		this_sheet.add_row ["#{this_phenotype.name}"], :style => title_font
    		  
    		this_phenotype.genes.each do |this_gene|
    			this_sheet.add_row ["#{this_gene.gene_symbol}", "#{this_gene.transcript}"]    			
    		end    		
    	end
    	return workbook
  	end

end
