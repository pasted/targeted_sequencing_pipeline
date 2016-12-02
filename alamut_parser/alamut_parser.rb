class AlamutParser
  require 'yaml'
  require 'smarter_csv'
  require 'axlsx'
  require 'active_support'
  
  require_relative 'interval'
  require_relative 'region'
  require_relative '../shared/gene'
  require_relative '../shared/phenotype'
  require_relative '../shared/transcript_store'
  require_relative '../shared/sample_store'
  require_relative '../shared/variant'
  require_relative '../shared/sample'
  require_relative '../shared/batch'
  require_relative '../shared/panel'
  require_relative '../shared/metric'
  require_relative '../shared/cnv'

  
  attr_accessor :batch_id, :base_path, :variants_directory, :intervals_directory, :transcript_file_path 
  attr_accessor :unwanted_file_path, :wanted_file_path, :sample_list_path, :panel_id, :panel_version
  
  def parse_file(file_name, panel_version_lc, ex_number_lc)
  	options = { :col_sep => "\t" }
  	variant_array = Array.new 
  	
  	if File.exists?(file_name) && ( File.stat(file_name).size > 0 )
  	
  		SmarterCSV.process( file_name, options ) do |csv|
  			this_variant = Variant.new

  			this_variant.id										= csv.first[:id]
  			this_variant.assembly								= csv.first[:assembly]
  			this_variant.chromosome 							= csv.first[:chrom]
  			this_variant.position								= csv.first[:inputpos]
  			this_variant.genomic_dna_start						= csv.first[:gdnastart]
  			this_variant.genomic_dna_end						= csv.first[:gdnaend]
  			this_variant.complementary_dna_start				= csv.first[:cdnastart]
  			this_variant.complementary_dna_end					= csv.first[:cdnaend]
  			this_variant.gene_id								= csv.first[:geneid]
  			this_variant.gene 									= csv.first[:gene]
  			this_variant.gene_description						= csv.first[:genedesc]
  			this_variant.genotype								= csv.first[:"gt_(#{panel_version_lc}_#{ex_number_lc})"]
  			this_variant.full_transcript 						= csv.first[:transcript]
  			this_variant.transcript_length						= csv.first[:translen]
  			this_variant.var_type								= csv.first[:vartype]
  			this_variant.coding_effect							= csv.first[:codingeffect]
  			this_variant.var_location 							= csv.first[:varlocation]
  			this_variant.strand 								= csv.first[:strand]
  			this_variant.genomic_nomen 							= csv.first[:gnomen]
  			this_variant.cdna_nomen	 							= csv.first[:cnomen]
  			this_variant.protein_nomen 							= csv.first[:pnomen]
  			this_variant.protein								= csv.first[:protein]
  			this_variant.uniprot								= csv.first[:uniprot]
  			this_variant.alt_protein_nomen						= csv.first[:alt_pnomen]
  			this_variant.go_bio_process							= csv.first[:gobioprocess]
  			this_variant.go_cell_comp							= csv.first[:gocellcomp]
  			this_variant.go_mol_func							= csv.first[:gomolfunc]
  			this_variant.omim_id								= csv.first[:omimid]
  			this_variant.exon									= csv.first[:exon]
  			this_variant.intron									= csv.first[:intron]
  			this_variant.distance_nearest_splice_site 			= csv.first[:distnearestss]
  			this_variant.nearest_splice_site_type 				= csv.first[:nearestsstype]
  			this_variant.nearest_splice_site_change 			= csv.first[:nearestsschange]
  			this_variant.local_splice_effect					= csv.first[:localspliceeffect]
  			this_variant.local_ss_pos							= csv.first[:localss_pos]
  			this_variant.local_ss_wt_max_ent_score				= csv.first[:localss_wtmaxentscore]
  			this_variant.local_ss_wt_nns_score					= csv.first[:localss_wtnnsscore]
  			this_variant.local_ss_wt_hsf_score					= csv.first[:localss_wthsfscore]
  			this_variant.local_ss_var_max_ent_score				= csv.first[:localss_varmaxentscore]
  			this_variant.local_ss_var_nns_score					= csv.first[:localss_varnnsscore]
  			this_variant.local_ss_var_hsf_score					= csv.first[:localss_varhsfscore]
  			this_variant.wt_ssf_score							= csv.first[:wtssfscore]
  			this_variant.wt_max_ent_score						= csv.first[:wtmaxentscore]
  			this_variant.wt_nns_score							= csv.first[:wtnnsscore]
  			this_variant.wt_gs_score							= csv.first[:wtgsscore]
  			this_variant.wt_hsf_score							= csv.first[:wthsfscore]
  			this_variant.var_ssf_score							= csv.first[:varssfscore]
  			this_variant.var_max_ent_score						= csv.first[:varmaxentscore]
  			this_variant.var_nns_score							= csv.first[:varnnsscore]
  			this_variant.var_gs_score							= csv.first[:vargsscore]
  			this_variant.var_hsf_score							= csv.first[:varhsfscore]	
  			this_variant.branch_point_pos						= csv.first[:branchpointpos]
  			this_variant.branch_point_change					= csv.first[:branchpointchange]
  			this_variant.protein_domain_1						= csv.first[:proteindomain1]
  			this_variant.protein_domain_2						= csv.first[:proteindomain2]
  			this_variant.protein_domain_3						= csv.first[:proteindomain3]
  			this_variant.protein_domain_4						= csv.first[:proteindomain4]
  			this_variant.rs_id 									= csv.first[:rsid]
  			this_variant.rs_validated 							= csv.first[:rsvalidated]
  			this_variant.rs_suspect 							= csv.first[:rssuspect]
  			this_variant.rs_validations 						= csv.first[:rsvalidations]
  			this_variant.rs_validation_number					= csv.first[:rsvalidationnumber]
  			this_variant.rs_ancestral_allele					= csv.first[:rsancestralallele]
  			this_variant.rs_heterozygosity						= csv.first[:rsheterozygosity]
  			this_variant.rs_clinical_significance				= csv.first[:rsclinicalsignificance]
  			this_variant.rs_maf_allele							= csv.first[:rsmafallele]
  			this_variant.rs_maf									= csv.first[:rsmaf]
  			this_variant.rs_maf_count							= csv.first[:rsmafcount]
  			this_variant.genomes_1000_freq						= csv.first[:"1000g_af"]
  			this_variant.genomes_1000_afr_freq					= csv.first[:"1000g_afr_af"]
  			this_variant.genomes_1000_sas_freq  				= csv.first[:"1000g_sas_af"]
  			this_variant.genomes_1000_eas_freq  				= csv.first[:"1000g_eas_af"]
  			this_variant.genomes_1000_eur_freq  				= csv.first[:"1000g_eur_af"]
  			this_variant.genomes_1000_amr_freq  				= csv.first[:"1000g_amr_af"]
  			this_variant.exac_alt_freq_all       				= csv.first[:exacaltfreq_all]
  			this_variant.exac_alt_afr_freq						= csv.first[:exacaltfreq_afr]
  			this_variant.exac_alt_amr_freq						= csv.first[:exacaltfreq_amr]
  			this_variant.exac_alt_eas_freq						= csv.first[:exacaltfreq_eas]
  			this_variant.exac_alt_sas_freq						= csv.first[:exacaltfreq_sas]
  			this_variant.exac_alt_nfe_freq						= csv.first[:exacaltfreq_nfe]
  			this_variant.exac_alt_fin_freq						= csv.first[:exacaltfreq_fin]
  			this_variant.exac_alt_oth_freq						= csv.first[:exacaltfreq_oth]	
  			this_variant.exac_alt_count_all       				= csv.first[:exacaltcount_all]
  			this_variant.exac_alt_count_afr						= csv.first[:exacaltcount_afr]
  			this_variant.exac_alt_count_amr						= csv.first[:exacaltcount_amr]
  			this_variant.exac_alt_count_eas						= csv.first[:exacaltcount_eas]
  			this_variant.exac_alt_count_sas						= csv.first[:exacaltcount_sas]
  			this_variant.exac_alt_count_nfe						= csv.first[:exacaltcount_nfe]
  			this_variant.exac_alt_count_fin						= csv.first[:exacaltcount_fin]
  			this_variant.exac_alt_count_oth						= csv.first[:exacaltcount_oth]			
  			this_variant.exac_total_count_all       			= csv.first[:exactotalcount_all]
				this_variant.exac_total_count_afr					= csv.first[:exactotalcount_afr]
				this_variant.exac_total_count_amr					= csv.first[:exactotalcount_amr]
				this_variant.exac_total_count_eas					= csv.first[:exactotalcount_eas]
				this_variant.exac_total_count_sas					= csv.first[:exactotalcount_sas]
				this_variant.exac_total_count_nfe					= csv.first[:exactotalcount_nfe]
				this_variant.exac_total_count_fin					= csv.first[:exactotalcount_fin]
				this_variant.exac_total_count_oth					= csv.first[:exactotalcount_oth]			
				this_variant.exac_hom_freq_all       				= csv.first[:exachomfreq_all]
				this_variant.exac_hom_freq_afr						= csv.first[:exachomfreq_afr]
				this_variant.exac_hom_freq_amr						= csv.first[:exachomfreq_amr]
				this_variant.exac_hom_freq_eas						= csv.first[:exachomfreq_eas]
				this_variant.exac_hom_freq_sas						= csv.first[:exachomfreq_sas]
				this_variant.exac_hom_freq_nfe						= csv.first[:exachomfreq_nfe]
				this_variant.exac_hom_freq_fin						= csv.first[:exachomfreq_fin]
				this_variant.exac_hom_freq_oth						= csv.first[:exachomfreq_oth]			
				this_variant.exac_hom_count_all       				= csv.first[:exachomcount_all]
				this_variant.exac_hom_count_afr						= csv.first[:exachomcount_afr]
				this_variant.exac_hom_count_amr						= csv.first[:exachomcount_amr]
				this_variant.exac_hom_count_eas						= csv.first[:exachomcount_eas]
				this_variant.exac_hom_count_sas						= csv.first[:exachomcount_sas]
				this_variant.exac_hom_count_nfe						= csv.first[:exachomcount_nfe]
				this_variant.exac_hom_count_fin						= csv.first[:exachomcount_fin]
				this_variant.exac_hom_count_oth						= csv.first[:exachomcount_oth]
  			this_variant.exac_filter         	  				= csv.first[:exacfilter]
  			this_variant.exac_read_depth        				= csv.first[:exacreaddepth]
  			this_variant.clin_var_ids           				= csv.first[:clinvarids]
  			this_variant.clin_var_origins    	  				= csv.first[:clinvarorigins]
  			this_variant.clin_var_methods    	 	 			= csv.first[:clinvarmethods]
  			this_variant.clin_var_clin_signifs					= csv.first[:clinvarclinsignifs]
  			this_variant.clin_var_review_status					= csv.first[:clinvarreviewstatus]
  			this_variant.clin_var_phenotypes	  				= csv.first[:clinvarphenotypes]
  			this_variant.cosmic_ids         	  				= csv.first[:cosmicids]
  			this_variant.cosmic_tissues      	  				= csv.first[:cosmictissues]
  			this_variant.cosmic_freqs           				= csv.first[:cosmicfreqs]
  			this_variant.cosmic_sample_counts	  				= csv.first[:cosmicsamplecounts]
  			this_variant.esp_all_maf							= csv.first[:espallmaf]
				this_variant.esp_ref_ea_count  						= csv.first[:esprefeacount]
				this_variant.esp_ref_aa_count  						= csv.first[:esprefaacount]
				this_variant.esp_ref_all_count 						= csv.first[:esprefallcount]
				this_variant.esp_alt_ea_count  						= csv.first[:espalteacount]
				this_variant.esp_alt_aa_count  						= csv.first[:espaltaacount]
				this_variant.esp_alt_all_count 						= csv.first[:espaltallcount]
				this_variant.esp_ea_maf       						= csv.first[:espeamaf]
				this_variant.esp_aa_maf       						= csv.first[:espaamaf]
				this_variant.esp_ea_aaf       						= csv.first[:espeaaaf]
				this_variant.esp_aa_aaf       						= csv.first[:espaaaaf]
				this_variant.esp_all_aaf      						= csv.first[:espallaaf]
				this_variant.esp_avg_read_depth						= csv.first[:espavgreaddepth]
				this_variant.cosmic_ids								= csv.first[:cosmicids]
				this_variant.cosmic_tissues							= csv.first[:cosmictissues]
				this_variant.cosmic_freqs       					= csv.first[:cosmicfreqs]
				this_variant.cosmic_sample_counts					= csv.first[:cosmicsamplecounts]
				
				this_variant.subst_type								= csv.first[:substtype]
				this_variant.nuc_change   							= csv.first[:nucchange]
				this_variant.phast_cons   							= csv.first[:phastcons]
				this_variant.phylo_p      							= csv.first[:phylop]
				this_variant.wt_aa_1      							= csv.first[:wtaa_1]
				this_variant.wt_aa_3      							= csv.first[:wtaa_3]
				this_variant.wt_codon     							= csv.first[:wtcodon]
				this_variant.wt_codon_freq 							= csv.first[:wtcodonfreq]
				this_variant.var_aa_1     							= csv.first[:varaa_1]
				this_variant.var_aa_3     							= csv.first[:varaa_3]
				this_variant.var_codon    							= csv.first[:varcodon]
				this_variant.var_codon_freq							= csv.first[:varcodonfreq]
				this_variant.pos_aa       							= csv.first[:posaa]
				
				this_variant.blosum_45          					= csv.first[:blosum45]
				this_variant.blosum_62          					= csv.first[:blosum62]
				this_variant.blosum_80          					= csv.first[:blosum80]
				this_variant.wt_aa_composition   					= csv.first[:wtaacomposition]
				this_variant.var_aa_composition  					= csv.first[:varaacomposition]
				this_variant.wt_aa_polarity      					= csv.first[:wtaapolarity]
				this_variant.var_aa_polarity     					= csv.first[:varaapolarity]
				this_variant.wt_aa_volume        					= csv.first[:wtaavolume]
				this_variant.var_aa_volume       					= csv.first[:varaavolume]
			
			
  			this_variant.hgmd_phenotype 						= csv.first[:hgmdphenotype]
  			this_variant.hgmd_pub_med_id 						= csv.first[:hgmdpubmedid]
  			this_variant.hgmd_sub_category 						= csv.first[:hgmdsubcategory]
  			this_variant.n_orthos								= csv.first[:northos]
  			this_variant.conserved_orthos						= csv.first[:conservedorthos]
  			this_variant.conserved_dist_species					= csv.first[:conserveddistspecies]
  			this_variant.grantham_dist							= csv.first[:granthamdist]
  			this_variant.agv_gd_gv								= csv.first[:agvgdgv]
  			this_variant.agv_gd_gd								= csv.first[:agvgdgd]	
  			this_variant.agv_gd_class							= csv.first[:agvgdclass]
  			
  			this_variant.sift_prediction						= csv.first[:siftprediction]
  			this_variant.sift_weight							= csv.first[:siftweight]
  			this_variant.sift_median							= csv.first[:siftmedian]
  			this_variant.mapp_prediction  						= csv.first[:mappprediction]
  			this_variant.mapp_p_value      						= csv.first[:mapppvalue]
  			this_variant.mapp_p_value_median					= csv.first[:mapppvaluemedian]
  			
  			this_variant.wt_nuc									= csv.first[:wtnuc]
  			this_variant.var_nuc								= csv.first[:varnuc]
  			this_variant.ins_nucs								= csv.first[:insnucs]
  			this_variant.del_nucs								= csv.first[:delnucs]
  			
  			this_variant.quality_vcf							= csv.first[:"quality_(vcf)"]
  			this_variant.filter_vcf								= csv.first[:"filter_(vcf)"]
  			this_variant.ac 									= csv.first[:ac]
  			this_variant.af 									= csv.first[:af]
  			this_variant.an 									= csv.first[:an]
  			this_variant.dp 									= csv.first[:dp]
  			this_variant.fs 									= csv.first[:fs]
  			this_variant.mq 									= csv.first[:mq]
  			this_variant.mq_0									= csv.first[:mq0]
  			this_variant.qd 									= csv.first[:qd]
  			
  			this_variant.ad										= csv.first[:"ad_(#{panel_version_lc}_#{ex_number_lc})"]
  			this_variant.dp										= csv.first[:"dp_(#{panel_version_lc}_#{ex_number_lc})"]
  			this_variant.gq										= csv.first[:"gq_(#{panel_version_lc}_#{ex_number_lc})"]
  			this_variant.gt										= csv.first[:"gt_(#{panel_version_lc}_#{ex_number_lc})"]
  			this_variant.pl										= csv.first[:"pl_(#{panel_version_lc}_#{ex_number_lc})"]
  			
  			this_variant.parse_transcript()
  			variant_array.push(this_variant)
  			
  		end
  	else
  		puts "ERROR :: variants file has no content :: SAMPLE ID : #{ex_number_lc}"
  		if File.exists?(file_name)
  			puts "File exists? #{File.exists?(file_name)}"
  			puts "File size? #{File.stat(file_name).size}"
  		else
  			puts "File exists? #{File.exists?(file_name)}"
  		end
  	end
  	return variant_array
  end
  
  def parse_sample_list(sample_file_path)
  	options = { :col_sep => "," }
  	sample_array = Array.new
  	
  	SmarterCSV.process( sample_file_path, options ) do |csv|
  		this_sample = Sample.new
  		this_sample.capture_number 	= csv.first[:capture_number]
  		this_sample.mody_number			= csv.first[:mody_number]
  		this_sample.ex_number				= csv.first[:ex_number]
  		this_sample.gender					= csv.first[:gender]
  		this_sample.profile					= csv.first[:profile]
  		this_sample.phenotype				= csv.first[:phenotype]
  		this_sample.sample_type			= csv.first[:sample_type]
  		this_sample.comment					= csv.first[:comment]
  		this_sample.parse_panel_version
  		this_sample.parse_sample_id
  		if this_sample.phenotype == nil
  			this_sample.parse_phenotype
  		end
  		
  		sample_array.push(this_sample)

  	end
  	
  	return sample_array
  end
  
  
  def parse_unwanted_variants(unwanted_file_path)
  	options = { :col_sep => "," }
  	unwanted_array = Array.new
  	
  	SmarterCSV.process( unwanted_file_path, options ) do |csv|
  		this_unwanted_variant = Variant.new

  		this_unwanted_variant.chromosome 						= csv.first[:chrom]
  		this_unwanted_variant.position							= csv.first[:start_pos]
  		this_unwanted_variant.gene 									= csv.first[:gene]		
  		this_unwanted_variant.wt_nuc								=	csv.first[:wtnuc]
  		this_unwanted_variant.var_nuc								=	csv.first[:varnuc]
  		this_unwanted_variant.ins_nucs							=	csv.first[:insnucs]
  		this_unwanted_variant.del_nucs							= csv.first[:delnucs]
  		this_unwanted_variant.reason_for_filtering 	= csv.first[:reason_for_filtering]
  		
  		unwanted_array.push(this_unwanted_variant)

  	end
  	
  	return unwanted_array
	end
	
	def parse_wanted_regions(wanted_file_path)
  	options = { :col_sep => "," }
  	wanted_array = Array.new
  	
  	SmarterCSV.process( wanted_file_path, options ) do |csv|
  		this_wanted_region = Region.new

  		this_wanted_region.chromosome 			= csv.first[:chromosome]
  		this_wanted_region.start_position		= csv.first[:start_position]
  		this_wanted_region.stop_position 		= csv.first[:stop_position]		

  		
  		wanted_array.push(this_wanted_region)

  	end
  	
  	return wanted_array
	end
	
	def parse_bed_intervals(intervals_file_path)
  		
  		interval_array = Array.new
  		
  		File.open(intervals_file_path, "r") do |this_file|
  			this_file.each_line do |line|
  				
  				this_interval = Interval.new
  				interval_attributes = line.split("\t")
  				this_interval.chromosome 			= interval_attributes[0]
  				this_interval.genomic_start		= interval_attributes[1]
  				this_interval.genomic_end 		= interval_attributes[2]		
  				this_interval.strand 					= "+"
  				this_interval.interval_name 	= interval_attributes[3].strip!	
  			
  				interval_array.push(this_interval)
  			end
  			
  		end
  	
  		return interval_array
	end
	
	def parse_intervals(intervals_file_path)
  	options = { :col_sep => "\t" }
  	interval_array = Array.new
  	
  	SmarterCSV.process( intervals_file_path, options ) do |csv|
  		this_interval = Interval.new

  		this_interval.chromosome 			= csv.first[:chromosome]
  		this_interval.genomic_start		= csv.first[:genomic_start]
  		this_interval.genomic_end 		= csv.first[:genomic_end]		
  		this_interval.strand 					= csv.first[:strand]
  		this_interval.interval_name 	= csv.first[:interval_name]	
  		
  		interval_array.push(this_interval)
  		
  	end
  	
  	return interval_array
	end
	
	def parse_cnvs(panel, cnvs_file_path)
		options = { :col_sep => "," }
		cnv_array = Array.new
		puts cnvs_file_path
		if File.exists?(cnvs_file_path)
  			SmarterCSV.process( cnvs_file_path, options ) do |csv|
  				this_cnv = Cnv.new
  				puts csv.first.inspect
  				this_cnv.start_p 					= csv.first[:"start.p"]
  				this_cnv.end_p 						= csv.first[:"end.p"]
  				this_cnv.type 						= csv.first[:type]
  				this_cnv.nexons 					= csv.first[:nexons]
  				this_cnv.genome_start 		= csv.first[:start]
  				this_cnv.genome_end 			= csv.first[:end]
  				this_cnv.chromosome 			= csv.first[:chromosome]
  				this_cnv.id	 							= csv.first[:id]
  				this_cnv.bayes_factor			= csv.first[:bf]
  				this_cnv.reads_expected 	= csv.first[:"reads.expected"]
  				this_cnv.reads_observed 	= csv.first[:"reads.observed"]
  				this_cnv.reads_ratio 			= csv.first[:"reads.ratio"]
  				this_cnv.exons						= csv.first[:"#{panel}_exons"]
  				cnv_array.push(this_cnv)
  			end
  		end
  	
  		return cnv_array
	end
	
	def load_phenotypes(transcript_file_path)
		phenotype_list = YAML.load_file("#{transcript_file_path}")
		return phenotype_list
	end
	
	def parse_metrics(input_file_path, metrics_line)
			counter = 0
			this_metric = Metric.new
			IO.foreach("#{input_file_path}") do |this_line| 
			
				if ( this_line.match(/(^\r|^\n|^\r\n|^#)/) )
					counter = counter + 1
				else
					if counter == metrics_line.to_i
						
						metrics_array = this_line.split("\t")
						
						this_metric.variable_order.each_with_index do |var, index|
							file_var = metrics_array[index]
							if (file_var == "\n") || (file_var == "")
								file_var = nil 
							end
							this_metric.instance_variable_set("@#{var}", file_var)
						end
					end
					counter = counter + 1
				end
			end
			this_metric.sample_id = this_metric.sample
			return this_metric
	end
	
	def populate_metrics(samples, this_batch, this_parser)
		processed_samples = Array.new
		samples.each do |this_sample|
			phenotype_metric_file_path = "#{this_batch.base_path}/#{this_batch.batch_id}/metrics/#{this_sample.panel_version}_#{this_sample.ex_number}_#{this_sample.gender.upcase!}_#{this_sample.phenotype}.phenotype.bait_capture_metrics"

			this_metric = this_parser.parse_metrics(phenotype_metric_file_path, this_batch.phenotype_metric_line)

			this_sample.add_metrics(Array.new.push(this_metric))

			processed_samples.push(this_sample)
		end
		
		return processed_samples
	end
  
	def run_parser(parser)
		
		path = File.expand_path(__FILE__)
		base_path = path.split("/scripts/").first
		
  		#Load the config YAML file and pass the settings to local variables
  		this_batch = YAML.load_file("../configuration/config.yaml")
  		
  		#Init AlamutParser class
  		#parser = AlamutParser.new()
  		
  		parser.batch_id = this_batch.batch_id
  		parser.base_path = this_batch.base_path
  		parser.sample_list_path = this_batch.sample_list_path
  		
  		
  		samples = parser.parse_sample_list("../configuration/#{parser.sample_list_path}")
    	
  		#populate phenotype specific metrics for each sample
  		samples = parser.populate_metrics(samples, this_batch, parser)
  		
  		sample_store = SampleStore.new(samples)
  		
  		#build the workbook
  		#this_book = Spreadsheet::Workbook.new
  		axlsx_package = Axlsx::Package.new
  		this_book = axlsx_package.workbook
  		title_font = this_book.styles.add_style :b => true
   	
  		#this_book = sample_store.process_samples(parser, this_batch, this_book, title_font)
  		tmp_sample_store = sample_store.process_samples(parser, this_batch, this_book, title_font)
  		
  		
  		this_book = tmp_sample_store.samples_to_axlsx(this_book, title_font)
		
  		this_book.add_worksheet(:name => "Batch") do |this_sheet|
  			
  				this_sheet.add_row ["Batch Id", "#{this_batch.batch_id}"]
  				this_sheet.add_row ["FTP url", "#{this_batch.ftp_url}"]
  				this_sheet.add_row ["Base path", "#{this_batch.base_path}"]
  				this_sheet.add_row ["Sample list path", "#{this_batch.sample_list_path}"]
  				this_sheet.add_row ["Alamut path", "#{this_batch.alamut_path}"]
  				this_sheet.add_row ["Java path", "#{this_batch.java_path}"]
  				this_sheet.add_row ["Rscript path", "#{this_batch.rscript_path}"]
  				this_sheet.add_row ["BWA path", "#{this_batch.bwa_path}"]
  				this_sheet.add_row ["Picard path", "#{this_batch.picard_path}"]
  				this_sheet.add_row ["GATK version", "#{this_batch.gatk_version}"]
  				this_sheet.add_row ["GATK path", "#{this_batch.gatk_path}"]
  				this_sheet.add_row ["Reference path", "#{this_batch.reference_path}"]
  				this_sheet.add_row ["Common variant NKMI path", "#{this_batch.common_variants_nkmi_path}"]
  				this_sheet.add_row ["Common artefacts path", "#{this_batch.common_artefacts_path}"]
  				this_sheet.add_row ["dbSNP path", "#{this_batch.dbsnp_path}"]
  				this_sheet.add_row ["FTP url", "#{this_batch.ftp_url}"]
  				this_sheet.add_row ["Flowcell", "#{this_batch.flowcell}"]

  		end
  	
		#write out intervals that have been analysed by phenotype
		
		available_phenotypes = sample_store.available_phenotypes
  	
  	
		#available_phenotypes.each do |phenotype|
		#	intervals = parser.parse_intervals("#{parser.intervals_directory}/#{parser.panel_version}_#{phenotype}_diagnostic_ROI.tsv")
		#	this_sheet = this_book.create_worksheet :name => "#{phenotype} intervals"
		#	row_number = 0
		#	intervals.each do |this_interval|
		#		this_sheet.row(row_number).push "#{this_interval.chromosome}","#{this_interval.genomic_start}","#{this_interval.genomic_end}","#{this_interval.strand}","#{this_interval.interval_name}"
 		#		row_number = row_number + 1
		#	end
		#	
		#end

		axlsx_package.serialize "#{parser.base_path}/#{parser.batch_id}/results/#{parser.batch_id}_variants.#{Time.now.strftime("%d-%m-%Y-%H%M%S")}.xlsx"
		
		File.open("#{parser.base_path}/#{parser.batch_id}/results/#{parser.batch_id}_variants.yaml", 'w+') {|f| f.write(sample_store.to_yaml)}

		return sample_store

	end
	
	if __FILE__ == $0

		parser = AlamutParser.new()
		sample_store = parser.run_parser(parser)
		
	end
end#end AlamutParser class
