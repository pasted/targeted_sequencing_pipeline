class Duplicate
	#Picard duplicates class

	attr_accessor :library, :unpaired_reads_examined, :read_pairs_examined, :unmapped_reads, :unpaired_read_duplicates, :read_pair_duplicates
	attr_accessor :read_pair_optical_duplicates, :percent_duplication, :estimated_library_size


	
	def variable_order
			variable_order = [:library, :unpaired_reads_examined, :read_pairs_examined, :unmapped_reads, :unpaired_read_duplicates, :read_pair_duplicates]          
			variable_order = variable_order + [:read_pair_optical_duplicates, :percent_duplication, :estimated_library_size]                                                                                                                                                                                   
			return variable_order
	end
		
	def print_attributes
			attr_array = Array.new
			self.instance_variables.map do |var|
				attr_array.push(self.instance_variable_get(var))
  		end
  			return attr_array
	end
		


end
