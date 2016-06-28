class Metric
	#Picard metrics class
	
	attr_accessor :bait_set, :genome_size, :bait_territory, :target_territory, :bait_design_efficiency, :total_reads
	attr_accessor :pf_reads, :pf_unique_reads, :pct_pf_reads, :pct_pf_uq_reads, :pf_uq_reads_aligned, :pct_pf_uq_reads_aligned
	attr_accessor :pf_uq_bases_aligned, :pf_uq_bases_aligned, :on_bait_bases, :near_bait_bases, :off_bait_bases, :on_target_bases
	attr_accessor :pct_selected_bases, :pct_off_bait, :on_bait_vs_selected, :mean_bait_coverage, :mean_target_coverage
	attr_accessor :pct_usable_bases_on_bait, :pct_usable_bases_on_target, :fold_enrichment, :zero_cvg_targets_pct, :fold_80_base_penalty
	attr_accessor :pct_target_bases_2x, :pct_target_bases_10x, :pct_target_bases_20x, :pct_target_bases_30x, :pct_target_bases_40x, :pct_target_bases_50x, :pct_target_bases_100x
	attr_accessor :hs_library_size, :hs_penalty_10x, :hs_penalty_20x, :hs_penalty_30x, :hs_penalty_40x, :hs_penalty_50x, :hs_penalty_100x
	attr_accessor :at_dropout, :gc_dropout, :ex_number, :library, :read_group
	
	
	def variable_order
			variable_order = [:bait_set, :genome_size, :bait_territory, :target_territory, :bait_design_efficiency, :total_reads, :pf_reads, :pf_unique_reads, :pct_pf_reads, :pct_pf_uq_reads, :pf_uq_reads_aligned, :pct_pf_uq_reads_aligned]          
			variable_order = variable_order + [:pf_uq_bases_aligned, :on_bait_bases, :near_bait_bases, :off_bait_bases, :on_target_bases, :pct_selected_bases, :pct_off_bait, :on_bait_vs_selected, :mean_bait_coverage, :mean_target_coverage]                                                                                                                                                                                   
			variable_order = variable_order + [:pct_usable_bases_on_bait, :pct_usable_bases_on_target, :fold_enrichment, :zero_cvg_targets_pct, :fold_80_base_penalty, :pct_target_bases_2x, :pct_target_bases_10x, :pct_target_bases_20x, :pct_target_bases_30x, :pct_target_bases_40x, :pct_target_bases_50x, :pct_target_bases_100x]                                                                                                                                                                 
			variable_order = variable_order + [:hs_library_size, :hs_penalty_10x, :hs_penalty_20x, :hs_penalty_30x, :hs_penalty_40x, :hs_penalty_50x, :hs_penalty_100x, :at_dropout, :gc_dropout, :ex_number, :library, :read_group]
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
