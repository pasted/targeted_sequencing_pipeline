class Metric
	#GATK metrics class
	
	attr_accessor :bait_set, :genome_size, :bait_territory, :target_territory, :bait_design_efficiency, :total_reads
	attr_accessor :pf_reads, :pf_unique_reads, :pct_pf_reads, :pct_pf_uq_reads, :pf_uq_reads_aligned, :pct_pf_uq_reads_aligned
	attr_accessor :pf_uq_bases_aligned, :pf_uq_bases_aligned, :on_bait_bases, :near_bait_bases, :off_bait_bases, :on_target_bases
	attr_accessor :pct_selected_bases, :pct_off_bait, :on_bait_vs_selected, :mean_bait_coverage, :mean_target_coverage
	attr_accessor :pct_usable_bases_on_bait, :pct_usable_bases_on_target, :fold_enrichment, :zero_cvg_targets_pct, :fold_80_base_penalty
	attr_accessor :pct_target_bases_2x, :pct_target_bases_10x, :pct_target_bases_20x, :pct_target_bases_30x, :pct_target_bases_40x, :pct_target_bases_50x, :pct_target_bases_100x
	attr_accessor :hs_library_size, :hs_penalty_10x, :hs_penalty_20x, :hs_penalty_30x, :hs_penalty_40x, :hs_penalty_50x, :hs_penalty_100x
	attr_accessor :at_dropout, :gc_dropout, :sample_id, :library, :read_group


end
