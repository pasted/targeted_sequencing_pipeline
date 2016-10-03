require 'erb'
class RenderScript < ERB

	def self.template(batch, run_type)
		File.read("#{batch.base_path}/#{batch.batch_id}/scripts/templates/exome_depth_#{run_type}.r_script.erb")
	end

	def initialize(panel_version, gender, run_type, batch, options = {})
		@panel_version							= panel_version
		@gender 										= gender
		@batch 											= batch
		@run_type										= run_type
		@template										= options.fetch(:template, self.class.template(batch, run_type))
		
		super(@template)
	end
	
	def result
		super(binding)
	end
end
