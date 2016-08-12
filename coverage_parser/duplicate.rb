class Duplicate
	#GATK duplicates class
	
	attr_accessor :libarary, :unpaired_reads_examined, :read_pairs_examined, :secondary_or_supplementary_rds, :unmapped_reads
	attr_accessor :unpaired_read_duplicates, :read_pair_duplicates, :read_pair_optical_duplicates
	attr_accessor :percent_duplication, :estimated_library_size, :sample_id

end
