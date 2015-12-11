class SnpStore
		require_relative 'snp'
		require_relative 'allele'
		
		attr_accessor :snps, :merged_snps
		
		def initialize(snps)
			self.snps = snps
  	end
  	
  	def process_merged_snps(snp_ids_to_merge)
  		
  		snps_to_merge = Array.new
 
  		self.snps.each do |this_snp| 			
  			if snp_ids_to_merge.include?("#{this_snp.id}") 
  				snps_to_merge.push(this_snp)
  			end
  		end
  		
  		merged_ids = []
  		alleles_to_merge = []
  		allele_array_lengths = []
  		snps_to_merge.each do |this_snp|
  			merged_ids.push(this_snp.id)
  			alleles_to_merge.push(this_snp.alleles)
  			allele_array_lengths.push(this_snp.alleles.length)
  		end
  	
  		number_of_samples = allele_array_lengths.first
  		count = 0
  		
  		merged_alleles = []
  		while count < number_of_samples  do
  			
   			sample_allele_left  = alleles_to_merge[0][count]
   			sample_allele_right = alleles_to_merge[1][count]
   			
   			merged_genotype = "#{sample_allele_left.genotype},#{sample_allele_right.genotype}"
   			if sample_allele_left.sample_id == sample_allele_right.sample_id
   			 merged_allele = Allele.new
   			 merged_allele.sample_id = "#{sample_allele_left.sample_id}"
   			 merged_allele.genotype = "#{merged_genotype}"
   			 merged_alleles.push(merged_allele)
   			else
   				puts "Sample mismatch in T1D merge alleles :: #{sample_allele_left.sample_id},#{sample_allele_left.sample_id}"
   			end
   			#puts this_sample.inspect
   			count +=1
   		end
  		
  		merged_snp = Snp.new
  		merged_snp.id = merged_ids.join(",")
  		merged_snp.alleles = merged_alleles
  		self.snps.push(merged_snp)
  	end

end
