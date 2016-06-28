# @author Garan Jones
# @abstract Interval class: Interval class to hold the genomic coordinates
# @attr [String] chromosome The chromosome name
# @attr [String] genomic_start The genomic start of the interval
# @attr [String] genomic_end The genomic end of the interval
# @attr [String] strand The strand that the interval is related too; not used
# @attr [String] interval_name The interval name, based on gene symbol_exon pattern
class Interval

	attr_accessor :chromosome, :genomic_start, :genomic_end, :strand, :interval_name


end
