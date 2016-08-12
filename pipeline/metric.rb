# @author Garan Jones
# Metric class: Metric class to store and manipulate information from Picard CalculateHSMetrics output
# @attr [String] bait_set The name of the bait set used in the hybrid selection.
# @attr [String] genome_size The number of bases in the reference genome used for alignment.
# @attr [String] bait_territory The number of bases which have one or more baits on top of them.
# @attr [String] target_territory The unique number of target bases in the experiment where target is usually exons etc.
# @attr [String] bait_design_efficiency Target terrirtoy / bait territory.  1 == perfectly efficient, 0.5 = half of baited bases are not target. 
# @attr [String] total_reads The total number of reads in the SAM or BAM file examine.
# @attr [String] pf_reads The number of reads that pass the vendor's filter.
# @attr [String] pf_unique_reads The number of PF reads that are not marked as duplicates.
# @attr [String] pct_pf_reads PF reads / total reads.  The percent of reads passing filter.
# @attr [String] pct_pf_uq_reads PF Unique Reads / Total Reads.
# @attr [String] pf_uq_reads_aligned The number of PF unique reads that are aligned with mapping score > 0 to the reference genome.
# @attr [String] pct_pf_uq_reads_aligned PF Reads Aligned / PF Reads.
# @attr [String] pf_bases_aligned The number of PF unique bases that are aligned with mapping score > 0 to the reference genome.
# @attr [String] pf_uq_bases_aligned The number of bases in the PF aligned reads that are mapped to a reference base. Accounts for clipping and gaps.
# @attr [String] on_bait_bases The number of PF aligned bases that mapped to a baited region of the genome.
# @attr [String] near_bait_bases The number of PF aligned bases that mapped to within a fixed interval of a baited region, but not on a baited region.
# @attr [String] off_bait_bases The number of PF aligned bases that mapped to neither on or near a bait.
# @attr [String] on_target_bases The number of PF aligned bases that mapped to a targeted region of the genome.
# @attr [String] pct_selected_bases On+Near Bait Bases / PF Bases Aligned.
# @attr [String] pct_off_bait The percentage of aligned PF bases that mapped neither on or near a bait.
# @attr [String] on_bait_vs_selected The percentage of on+near bait bases that are on as opposed to near.
# @attr [String] mean_bait_coverage The mean coverage of all baits in the experiment.
# @attr [String] mean_target_coverage The mean coverage of targets.
# @attr [String] median_target_coverage The median coverage of targets.
# @attr [String] pct_usable_bases_on_bait The number of aligned, de-duped, on-bait bases out of the PF bases available.
# @attr [String] pct_usable_bases_on_target The number of aligned, de-duped, on-target bases out of the PF bases available.
# @attr [String] fold_enrichment The fold by which the baited region has been amplified above genomic background.
# @attr [String] zero_cvg_targets_pct The fraction of targets that did not reach coverage=1 over any base.
# @attr [String] pct_exc_dupe The fraction of aligned bases that were filtered out because they were in reads marked as duplicates.
# @attr [String] pct_exc_mapq The fraction of aligned bases that were filtered out because they were in reads with low mapping quality.
# @attr [String] pct_exc_baseq The fraction of aligned bases that were filtered out because they were of low base quality. 
# @attr [String] pct_exc_overlap The fraction of aligned bases that were filtered out because they were the second observation from an insert with overlapping reads.
# @attr [String] pct_exc_off_target The fraction of aligned bases that were filtered out because they did not align over a target base. 
# @attr [String] fold_80_base_penalty The fold over-coverage necessary to raise 80% of bases in "non-zero-cvg" targets to the mean coverage level in those targets.
# @attr [String] pct_target_bases_1x The percentage of all target bases achieving 1X or greater coverage.
# @attr [String] pct_target_bases_2x The percentage of all target bases achieving 2X or greater coverage.
# @attr [String] pct_target_bases_10x The percentage of all target bases achieving 10X or greater coverage. 
# @attr [String] pct_target_bases_20x The percentage of all target bases achieving 20X or greater coverage. 
# @attr [String] pct_target_bases_30x The percentage of all target bases achieving 30X or greater coverage.
# @attr [String] pct_target_bases_40x The percentage of all target bases achieving 40X or greater coverage.
# @attr [String] pct_target_bases_50x The percentage of all target bases achieving 50X or greater coverage.
# @attr [String] pct_target_bases_100x The percentage of all target bases achieving 100X or greater coverage.
# @attr [String] hs_library_size The estimated number of unique molecules in the selected part of the library.
# @attr [String] hs_penalty_10x The "hybrid selection penalty" incurred to get 80% of target bases to 10X. This metric should be interpreted as: if I have a design with 10 megabases of target, and want to get 0X coverage I need to sequence until PF_ALIGNED_BASES = 10^7 * 10 * HS_PENALTY_10X.
# @attr [String] hs_penalty_20x The "hybrid selection penalty" incurred to get 80% of target bases to 20X. This metric should be interpreted as: if I have a design with 10 megabases of target, and want to get 20X coverage I need to sequence until PF_ALIGNED_BASES = 10^7 * 20 * HS_PENALTY_20X.
# @attr [String] hs_penalty_30x The "hybrid selection penalty" incurred to get 80% of target bases to 30X. This metric should be interpreted as: if I have a design with 10 megabases of target, and want to get 30X coverage I need to sequence until PF_ALIGNED_BASES = 10^7 * 30 * HS_PENALTY_30X.
# @attr [String] hs_penalty_40x The "hybrid selection penalty" incurred to get 80% of target bases to 40X. This metric should be interpreted as: if I have a design with 10 megabases of target, and want to get 40X coverage I need to sequence until PF_ALIGNED_BASES = 10^7 * 40 * HS_PENALTY_40X.
# @attr [String] hs_penalty_50x The "hybrid selection penalty" incurred to get 80% of target bases to 50X. This metric should be interpreted as: if I have a design with 10 megabases of target, and want to get 50X coverage I need to sequence until PF_ALIGNED_BASES = 10^7 * 50 * HS_PENALTY_50X.
# @attr [String] hs_penalty_100x The "hybrid selection penalty" incurred to get 80% of target bases to 100X. This metric should be interpreted as: if I have a design with 10 megabases of target, and want to get 100X coverage I need to sequence until PF_ALIGNED_BASES = 10^7 * 100 * HS_PENALTY_100X.
# @attr [String] at_dropout A measure of how undercovered <= 50% GC regions are relative to the mean. For each GC bin [0..50] we calculate a = % of target territory, and b = % of aligned reads aligned to these targets. AT DROPOUT is then abs(sum(a-b when a-b < 0)). E.g. if the value is 5% this implies that 5% of total reads that should have mapped to GC<=50% regions mapped elsewhere.
# @attr [String] gc_dropout A measure of how undercovered >= 50% GC regions are relative to the mean. For each GC bin [50..100] we calculate a = % of target territory, and b = % of aligned reads aligned to these targets. GC DROPOUT is then abs(sum(a-b when a-b < 0)). E.g. if the value is 5% this implies that 5% of total reads that should have mapped to GC>=50% regions mapped elsewhere.
# @attr [String] het_snp_sensitivity The theoretical HET SNP sensitivity.
# @attr [String] het_snp_q The Phred Scaled Q Score of the theoretical HET SNP sensitivity.
# @attr [String] sample_id Sample id from the sample_list
# @attr [String] library Panel library
# @attr [String] read_group Read group
class Metric
	#Picard metrics class
	
	attr_accessor :bait_set, :genome_size, :bait_territory, :target_territory, :bait_design_efficiency, :total_reads
	attr_accessor :pf_reads, :pf_unique_reads, :pct_pf_reads, :pct_pf_uq_reads, :pf_uq_reads_aligned, :pct_pf_uq_reads_aligned, :pf_bases_aligned
	attr_accessor :pf_uq_bases_aligned, :pf_uq_bases_aligned, :on_bait_bases, :near_bait_bases, :off_bait_bases, :on_target_bases
	attr_accessor :pct_selected_bases, :pct_off_bait, :on_bait_vs_selected, :mean_bait_coverage, :mean_target_coverage, :median_target_coverage
	attr_accessor :pct_usable_bases_on_bait, :pct_usable_bases_on_target, :fold_enrichment, :zero_cvg_targets_pct	
	attr_accessor :pct_exc_dupe, :pct_exc_mapq, :pct_exc_baseq, :pct_exc_overlap, :pct_exc_off_target, :fold_80_base_penalty	
	attr_accessor :pct_target_bases_1x, :pct_target_bases_2x, :pct_target_bases_10x, :pct_target_bases_20x, :pct_target_bases_30x, :pct_target_bases_40x, :pct_target_bases_50x, :pct_target_bases_100x
	attr_accessor :hs_library_size, :hs_penalty_10x, :hs_penalty_20x, :hs_penalty_30x, :hs_penalty_40x, :hs_penalty_50x, :hs_penalty_100x
	attr_accessor :at_dropout, :gc_dropout, :het_snp_sensitivity, :het_snp_q, :sample_id, :library, :read_group
	

		# @author Garan Jones
		# Specified order that the attributes will be printed out
  	# @return [Array<Symbol>] An array of attribute symbols
	def variable_order
			variable_order = [:bait_set, :genome_size, :bait_territory, :target_territory, :bait_design_efficiency, :total_reads, :pf_reads, :pf_unique_reads, :pct_pf_reads, :pct_pf_uq_reads, :pf_uq_reads_aligned, :pct_pf_uq_reads_aligned, :pf_bases_aligned]          
			variable_order = variable_order + [:pf_uq_bases_aligned, :on_bait_bases, :near_bait_bases, :off_bait_bases, :on_target_bases, :pct_selected_bases, :pct_off_bait, :on_bait_vs_selected, :mean_bait_coverage, :mean_target_coverage, :median_target_coverage]                                                                                                                                                                                   
			variable_order = variable_order + [:pct_usable_bases_on_bait, :pct_usable_bases_on_target, :fold_enrichment, :zero_cvg_targets_pct, :pct_exc_dupe, :pct_exc_mapq, :pct_exc_baseq, :pct_exc_overlap, :pct_exc_off_target, :fold_80_base_penalty, :pct_target_bases_1x, :pct_target_bases_2x, :pct_target_bases_10x, :pct_target_bases_20x, :pct_target_bases_30x, :pct_target_bases_40x, :pct_target_bases_50x, :pct_target_bases_100x]                                                                                                                                                                 
			variable_order = variable_order + [:hs_library_size, :hs_penalty_10x, :hs_penalty_20x, :hs_penalty_30x, :hs_penalty_40x, :hs_penalty_50x, :hs_penalty_100x, :at_dropout, :gc_dropout, :het_snp_sensitivity, :het_snp_q, :sample_id, :library, :read_group]
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
