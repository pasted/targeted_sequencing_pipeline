# @author Garan Jones
# @abstract Region class: Region class to hold the genomic coordinates with annotation that are of interest
# @attr [String] chromosome The chromosome name
# @attr [String] start_position The genomic start of the interval
# @attr [String] stop_position The genomic end of the interval
# @attr [String] annotation The annotation that relates to this region
class Region
		attr_accessor :chromosome, :start_position, :stop_position, :annotation
			
end
