class Batch
	attr_accessor :batch_id, :ftp_url, :base_path, :sample_list_path, :panels, :alamut_path, :hgmd_user, :hgmd_pass 
	attr_accessor :sequencer_name, :java_path, :rscript_path, :bwa_path, :picard_path, :gatk_version, :gatk_path, :tmp_path
	attr_accessor :reference_path, :known_path, :common_variants_nkmi_path, :common_artefacts_path, :dbsnp_path, :ndm_snps_path, :type_one_snps_path
	attr_accessor :phenotype_metric_line, :overall_metric_line, :duplicate_metric_line, :flowcell
	

	def select_panel(panel_version)
		#should only ever be one panel for each version
		#however return the first panel from the selected_panels array
		panel_version.upcase!
		#puts panel_version
		selected_panels = self.panels.select{|p| p.panel_version == panel_version}
		return selected_panels.first
	end
end
