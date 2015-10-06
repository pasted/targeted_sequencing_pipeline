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
    
    def transcripts_to_excel(workbook, panel_version)
    	this_sheet = workbook.create_worksheet :name => "#{panel_version} Transcripts"
    	row_number = 0
    	format = Spreadsheet::Format.new :color => :xls_color_42,
                                 :weight => :bold,
                                 :size => 12
                    
    	self.phenotypes.each do |this_phenotype|
    		this_sheet.row(row_number).push this_phenotype.name
    		this_sheet.row(row_number).default_format = format  
    		row_number = row_number + 1
    		this_phenotype.genes.each do |this_gene|
    			this_sheet.row(row_number).push this_gene.gene_symbol, this_gene.transcript
    			row_number = row_number + 1
    		end
    		row_number = row_number + 1
    	end
    	return workbook
  	end

end
