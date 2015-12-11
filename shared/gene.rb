class Gene
		attr_accessor :gene_symbol, :transcript
		
		def initialize(gene_symbol, transcript)
			self.gene_symbol = gene_symbol
			self.transcript = transcript
  	end
end
