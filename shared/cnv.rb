require 'active_support/core_ext/range'
class Cnv

	  attr_accessor :start_p, :end_p, :type, :nexons, :genome_start, :genome_end, :chromosome, :id
	  attr_accessor :bayes_factor, :reads_expected, :reads_observed, :reads_ratio, :exons, :wanted

	
	def variable_order
		variable_order = [:exons, :type, :bayes_factor, :reads_ratio, :reads_expected, :reads_observed]          
		variable_order = variable_order + [:nexons, :chromosome, :id, :genome_start, :genome_end, :start_p, :end_p]                                                                                                                                                                                   
		return variable_order
	end
	
	def is_wanted_region?(wanted_region_array)
		self.wanted = false
		wanted_region_array.each do |wanted_region|	
			wanted_range = Range.new(wanted_region.genomic_start.to_i, wanted_region.genomic_end.to_i)
			cnv_range = Range.new(self.genome_start.to_i, self.genome_end.to_i)
			
			#if ( ("#{self.chromosome}" == "#{wanted_region.chromosome}") && (self.genome_start.to_i.between?(wanted_region.genomic_start.to_i, wanted_region.genomic_end.to_i)) || (self.genome_end.to_i.between?(wanted_region.genomic_start.to_i, wanted_region.genomic_end.to_i)) )
			if ( ("#{self.chromosome}" == "#{wanted_region.chromosome}") && wanted_range.overlaps?(cnv_range) )
				self.wanted = true
			end

		end

	end
end
