class CnvStore
	attr_accessor :cnvs
		
	def initialize()
			self.cnvs = Hash.new
	end
  
	def process_cnvs(this_sample, this_batch, this_panel, parser)
		phenotype_intervals_path="#{this_panel.intervals_directory}/#{this_sample.panel_version.downcase}_#{this_sample.phenotype}_variant_calling.bed"

		if File.exists?("#{phenotype_intervals_path}")
			intervals_array = parser.parse_bed_intervals("#{phenotype_intervals_path}")

			self.cnvs.each_pair do |this_ex_number, this_cnv_array|
				tmp_cnv_array = Array.new
				this_cnv_array.each do |this_cnv|
					this_cnv.is_wanted_region?(intervals_array)
					tmp_cnv_array.push(this_cnv)

				end
				self.cnvs.store(this_ex_number, tmp_cnv_array)
			end
			
		end

	end

end
