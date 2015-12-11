class ParseMetrics
	
	require_relative 'wrapper'
	require_relative 'metric'
	require_relative 'duplicate'
	
	def parse_metrics(input_file_path, metrics_line)
			counter = 0
			this_metric = Metric.new
			IO.foreach("#{input_file_path}") do |this_line| 
			
				if ( this_line.match /(^\r|^\n|^\r\n|^#)/ )
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
			return this_metric
	end
	
	def parse_duplicates(input_file_path, metrics_line)
			counter = 0
			this_duplicate = Duplicate.new
			IO.foreach("#{input_file_path}") do |this_line| 
				
				if ( this_line.match /(^\r|^\n|^\r\n|^#)/ )
					counter = counter + 1
				else

					if counter == metrics_line.to_i

						duplicates_array = this_line.split("\t")
						
						this_duplicate.variable_order.each_with_index do |var, index|
							file_var = duplicates_array[index]
							if (file_var == "\n") || (file_var == "")
								file_var = nil 
							end

							this_duplicate.instance_variable_set("@#{var}", file_var.strip)
						end
					end
					counter = counter + 1
				end
			end
			return this_duplicate
	end
	
	def parse_batch_metrics(this_batch, samples)
		
			this_metric = ParseMetrics.new
			phenotype_metric_line = this_batch.phenotype_metric_line
			overall_metric_line = this_batch.overall_metric_line
			duplicate_metric_line = this_batch.duplicate_metric_line
	   
			phenotype_metric_hash = Hash.new
			overall_metric_array = Array.new
			duplicate_array = Array.new
			
			samples.each do |this_sample|
				phenotype_metric_file_path = "../../metrics/#{this_sample.panel_version.downcase}_#{this_sample.sample_id}_#{this_sample.gender.upcase}_#{this_sample.phenotype.upcase}.phenotype.bait_capture_metrics"
				overall_metric_file_path = "../../metrics/#{this_sample.panel_version.downcase}_#{this_sample.sample_id}_#{this_sample.gender.upcase}.overall.bait_capture_metrics"
				duplicates_file_path = "../../duplicates/#{this_sample.panel_version.downcase}_#{this_sample.sample_id}.duplicates"
				
				this_phenotype_metric = this_metric.parse_metrics(phenotype_metric_file_path, phenotype_metric_line)
				if phenotype_metric_hash.has_key?(this_sample.phenotype.upcase)
					this_phenotype_array = phenotype_metric_hash[this_sample.phenotype.upcase]
					this_phenotype_array.push(this_phenotype_metric)
					phenotype_metric_hash[this_sample.phenotype.upcase] = this_phenotype_array
				else
					this_phenotype_array = Array.new
					this_phenotype_array.push(this_phenotype_metric)
					phenotype_metric_hash[this_sample.phenotype.upcase] = this_phenotype_array
				end
				
				this_overall_metric = this_metric.parse_metrics(overall_metric_file_path, overall_metric_line)
				overall_metric_array.push(this_overall_metric)
				
				this_duplicate = this_metric.parse_duplicates(duplicates_file_path, duplicate_metric_line)
				duplicate_array.push(this_duplicate)
			end
			
			#write out phenotype metrics
			phenotype_metric_hash.each_pair do |phenotype, this_metric_array|
				csv_data = CSV.generate(col_sep: "\t") do |csv|
					csv << this_metric_array.first.variable_order
					this_metric_array.each do |this_metric|
						csv << this_metric.print_attributes
					end
				end
				# Writing the csv data back to the same file, (also specify UTF-8 format)
				File.open("../../metrics/#{this_batch.batch_id}_#{phenotype}.phenotype.metrics", 'w:UTF-8') { |file| file.write(csv_data)}
			end
			
			csv_data = CSV.generate(col_sep: "\t") do |csv|
					csv << overall_metric_array.first.variable_order
					overall_metric_array.each do |this_metric|
						csv << this_metric.print_attributes
					end
			end
			
			File.open("../../metrics/#{this_batch.batch_id}.overall.metrics", 'w:UTF-8') { |file| file.write(csv_data)}
			
			csv_data = CSV.generate(col_sep: "\t") do |csv|
					csv << duplicate_array.first.variable_order
					duplicate_array.each do |this_duplicate|
						csv << this_duplicate.print_attributes
					end
			end
				
			File.open("../../duplicates/#{this_batch.batch_id}.overall.duplicates", 'w:UTF-8') { |file| file.write(csv_data)}
			
			return true
	end

end
