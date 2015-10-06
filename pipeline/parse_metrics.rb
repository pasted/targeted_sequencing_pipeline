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
	


end
