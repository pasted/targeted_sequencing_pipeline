class Phenotype
	attr_accessor :name, :genes
		
	def initialize(name, genes)
		self.name = name
		self.genes = genes
  	end
  	
  	def all_transcripts
  		available_transcripts = self.genes.collect {|this_gene| this_gene.transcript}
  		return available_transcripts
  	end

end
