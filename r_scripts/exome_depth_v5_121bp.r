library("GenomicRanges")
library("ExomeDepth")

#start clock
ptm <- proc.time()
 
data(Conrad.hg19)
data(exons.hg19)
data(exons.hg19.X)

exons.hg19.GRanges <- GRanges(seqnames = exons.hg19$chromosome,
                              IRanges(start=exons.hg19$start,end=exons.hg19$end),
                              names = exons.hg19$name)


exons.hg19.X.GRanges <- GRanges(seqnames = exons.hg19.X$chromosome,
                              IRanges(start=exons.hg19.X$start,end=exons.hg19.X$end),
                              names = exons.hg19.X$name)
                              
Alignment_folder <- '/mnt/Data4/targeted_sequencing/2015-05-18_1501931/assembly'
print(Alignment_folder)
setwd(Alignment_folder)
bed <- '/mnt/Data1/targeted_sequencing/intervals/v5/v5_CNV_121.bed'
reference.fasta <- '/mnt/Data1/resources/human_g1k_v37.fasta'

#get list of sample BAM files
bams.selected <- Sys.glob('v5_*_FEMALE_*.realigned.bam')
bais.selected <- Sys.glob('v5_*_FEMALE_*.realigned.bai')

print(bams.selected)

dir.base <- NULL
dir.sub <- NULL
dir.output <- NULL
dir.results <- NULL


	# Create counts dataframe for all BAMs
	my.counts <- getBamCounts(bed.file = bed, bam.files = bams.selected, index.files = bais.selected, include.chr = FALSE, referenceFasta = reference.fasta)
	ExomeCount.dafr <- as(my.counts[, colnames(my.counts)], 'data.frame')
	ExomeCount.dafr$chromosome <- gsub(as.character(ExomeCount.dafr$space),pattern = 'chr',replacement = '')
	
	# Create matrix
	ExomeCount.mat <- as.matrix(ExomeCount.dafr[, grep(names(ExomeCount.dafr), pattern = '*.bam')])
	
	nsamples <- ncol(ExomeCount.mat)
	samplenames <- colnames(ExomeCount.mat)
	
	dir.base <- "/mnt/Data4/targeted_sequencing/2015-05-18_1501931/cnv_analysis/CNVs/v5"
	save.image( sprintf("%s/exome_depth_count_v5_female.RData", dir.base) )
	
	for (i in 1:nsamples) {
  
            samplename <- samplenames[i]
   
            print(samplename)
  #### Create the aggregate reference set for this sample
            my.choice <- select.reference.set (test.counts = ExomeCount.mat[,i],
                                     reference.counts = ExomeCount.mat[,-i],
                                     bin.length = (ExomeCount.dafr$end - ExomeCount.dafr$start)/1000,
                                     n.bins.reduced = 10000)
            my.reference.selected <- apply(X = ExomeCount.mat[, my.choice$reference.choice, drop = FALSE],
                                 MAR = 1,
                                 FUN = sum)
            message('Now creating the ExomeDepth object')
            all.exons <- new('ExomeDepth',
                   test = ExomeCount.mat[,i],
                   reference = my.reference.selected,
                   formula = 'cbind(test, reference) ~ 1')
  ################ Now call the CNVs
            all.exons <- CallCNVs(x = all.exons,
                        transition.probability = 10^-4,
                        chromosome = ExomeCount.dafr$chromosome,
                        start = ExomeCount.dafr$start,
                        end = ExomeCount.dafr$end,
                        name = ExomeCount.dafr$names)
            if (nrow(all.exons@CNV.calls) > 0) {
      
                  all.exons <- AnnotateExtra(x = all.exons,
                              reference.annotation = Conrad.hg19.common.CNVs,
                              min.overlap = 0.5,
                              column.name = 'Conrad.hg19')
                  all.exons <- AnnotateExtra(x = all.exons,
                              reference.annotation = exons.hg19.GRanges,
                              min.overlap = 0.0001,
                              column.name = 'exons.hg19')
                 # all.exons <- AnnotateExtra(x = all.exons,
                 #             reference.annotation = exons.hg19.X.GRanges,
                 #             min.overlap = 0.0001,
                 #             column.name = 'exons.hg19.X')
                 
                  output.file <- paste(dir.base, '/Sample_', samplename, '.csv', sep = '')
                  print(output.file)
                  write.csv(file = output.file, x = all.exons@CNV.calls, row.names = FALSE)
   
            }
            
            
         }

#stop clock
finish <- proc.time() - ptm
print(finish)
summary(finish)
