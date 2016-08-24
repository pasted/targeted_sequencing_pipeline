# @author Garan Jones
# Duplicate class: Duplicate class to store and manipulate information from Picard MarkDuplicates output
# @attr [String] library Name of library that MarkDuplicates was run against usually either the panel or phenotype intervals
# @attr [String] unpaired_reads_examined The number of mapped reads examined which did not have a mapped mate pair, either because the read is unpaired, or the read is paired to an unmapped mate.
# @attr [String] read_pairs_examined The number of mapped read pairs examined.
# @attr [String] unmapped_reads The total number of unmapped reads examined.
# @attr [String] unpaired_read_duplicates The number of fragments that were marked as duplicates.
# @attr [String] read_pair_duplicates The number of read pairs that were marked as duplicates.
# @attr [String] read_pair_optical_duplicates The number of read pairs duplicates that were caused by optical duplication. Value is always < READ_PAIR_DUPLICATES, which counts all duplicates regardless of source.
# @attr [String] percent_duplication 	The percentage of mapped sequence that is marked as duplicate.
# @attr [String] estimated_library_size The estimated number of unique molecules in the library based on PE duplication.
class Duplicate
	attr_accessor :library, :unpaired_reads_examined, :read_pairs_examined, :secondary_or_supplementary_rds, :unmapped_reads, :unpaired_read_duplicates, :read_pair_duplicates
	attr_accessor :read_pair_optical_duplicates, :percent_duplication, :estimated_library_size


		# @author Garan Jones
		# Specified order that the attributes will be printed out
  	# @return [Array<Symbol>] An array of attribute symbols
	def variable_order
			variable_order = [:library, :unpaired_reads_examined, :read_pairs_examined, :secondary_or_supplementary_rds, :unmapped_reads, :unpaired_read_duplicates, :read_pair_duplicates]          
			variable_order = variable_order + [:read_pair_optical_duplicates, :percent_duplication, :estimated_library_size]                                                                                                                                                                                   
			return variable_order
	end
		
		# @author Garan Jones
		# Print out the values of all existing instance variables
  	# @return [Array<String>] An array of attribute values
	def print_attributes
			attr_array = Array.new
			self.instance_variables.map do |var|
				attr_array.push(self.instance_variable_get(var))
  		end
  			return attr_array
	end
		


end
