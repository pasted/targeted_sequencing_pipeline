class Sample

	  attr_accessor :capture_number, :mody_number, :ex_number, :gender, :analysis_type, :geneset, :sample_type, :comment
	  attr_accessor :mean_bait_coverage, :pct_target_bases_30x, :sample_id
	  
	  def parse_sample_id
	  	capture_number_array = self.capture_number.split('_')	
  		self.sample_id = "#{capture_number_array[0].upcase}_#{capture_number_array[2]}"
	  end
	  
	  def add_metrics(metrics_array)
	  	metrics_array.each do |this_metric|	
	  		if this_metric.sample_id == self.sample_id
	  			self.mean_bait_coverage = this_metric.mean_bait_coverage
	  			self.pct_target_bases_30x = this_metric.pct_target_bases_30x
	  		end
	  	end
		end

end
