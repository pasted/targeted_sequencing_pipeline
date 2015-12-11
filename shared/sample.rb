class Sample

	  attr_accessor :capture_number, :mody_number, :ex_number, :gender, :profile, :phenotype, :sample_type, :comment
	  attr_accessor :mean_bait_coverage, :pct_target_bases_2x, :pct_target_bases_10x, :pct_target_bases_20x, :pct_target_bases_30x, :sample_id, :panel_version 
	  attr_accessor :sequencing_panel_version
	  
	  def parse_sample_id
	  	capture_number_array = self.capture_number.split('_')	
  		self.sample_id = "#{capture_number_array[1]}"
	  end
	  
	  def parse_panel_version
	  	#extract the part of the profile with either 1 or 3 integers
	  	self.panel_version = "#{self.profile.match(/v\d{3}|v\d{1}/)}"
	  	
	  	#check to see if the panel version has been parsed
	  	#if it has check to see if athe sequencing panel version has already been set
	  	#if not set it the same as the panel version
	  	
	  	if self.panel_version == nil
	  		puts "Error :: #{self.profile} :: not parsable to panel version"
	  	elsif self.sequencing_panel_version == nil
	  		self.sequencing_panel_version = self.panel_version
	  	end
	  	
	  end
	  
	  def parse_phenotype
	  	if self.panel_version != nil
	  		self.phenotype = self.profile.split(self.panel_version).last.strip
	  		#remove any non-phenotype string
	  		if self.phenotype =~ /\s/
	  			#puts ">>>>> #{self.phenotype}"
	  			self.phenotype = self.phenotype.split(" ").first
	  			#puts "#{self.phenotype}"
	  		end
	  	end
	  end
	  
	  def add_metrics(metrics_array)
	  	metrics_array.each do |this_metric|	
	  		if this_metric.sample_id == "#{self.panel_version}_#{self.sample_id}"
	  			self.mean_bait_coverage = this_metric.mean_bait_coverage
	  			self.pct_target_bases_2x = this_metric.pct_target_bases_2x
	  			self.pct_target_bases_10x = this_metric.pct_target_bases_10x
	  			self.pct_target_bases_20x = this_metric.pct_target_bases_20x
	  			self.pct_target_bases_30x = this_metric.pct_target_bases_30x
	  		end
	  	end
		end

end
