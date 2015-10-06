class Snp
	
	attr_accessor :chrom, :pos, :id, :ref, :alt, :qual, :filter, :ac, :af
	attr_accessor :alleles
	
	def initialize()
		self.alleles = []
	end
end

