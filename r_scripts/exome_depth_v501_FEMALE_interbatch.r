library("ExomeDepth")
library("GenomicRanges")
library("XML")
			
#start clock
ptm <- proc.time()

load("/mnt/Data1/resources/exome_depth_data/v501/v501_exons.RData")
load("/mnt/Data1/resources/exome_depth_data/v501/v501_controls.RData")

panel.version <- "v501"
sequencing.run <- "2016-07-04_2363"

#setup output XML DOM
xml.dom <- xmlOutputDOM()

dir.base <- sprintf("/mnt/Data4/targeted_sequencing/%s/cnv_analysis/female/%s", sequencing.run, panel.version)
dir.results <- sprintf("%s/results", dir.base)
dir.logs <- sprintf("%s/logs", dir.base)

v501_exons.GRanges <- GRanges(seqnames = v501_exons$chromosome,
                              IRanges(start=v501_exons$start,end=v501_exons$end),
                              names = v501_exons$name)
                              
alignment.folder <- sprintf("/mnt/Data4/targeted_sequencing/%s/assembly/", sequencing.run)

setwd(alignment.folder)

print(alignment.folder)

bed <- '/mnt/Data1/resources/exome_depth_data/v501/v501_CNV_121.bed'
reference.fasta <- '/mnt/Data1/resources/human_g1k_v37.fasta'



#get list of sample BAM files
bams.female <- Sys.glob('v501*_FEMALE.realigned.bam')
bais.female <- Sys.glob('v501*_FEMALE.realigned.bai')

bams.unknown <- Sys.glob('v501*_UNKNOWN.realigned.bam')
bais.unknown <- Sys.glob('v501*_UNKNOWN.realigned.bai')

bams.positive <- append(bams.female, bams.unknown)
bais.positive <- append(bais.female, bais.unknown)


FemaleExomeCounts.mat <- as.matrix(FemaleExomeCount.mat[, grep(names(FemaleExomeCount.mat), pattern = '*.bam')])

bams.selected <- c(bams.positive)
bais.selected <- c(bais.positive)
# Create counts dataframe for all BAMs
my.counts <- getBamCounts(bed.file = bed, bam.files = bams.selected, index.files = bais.selected, min.mapq = 20, include.chr = FALSE, referenceFasta = reference.fasta)
ExomeCount.dafr <- as(my.counts[, colnames(my.counts)], 'data.frame')
ExomeCount.dafr$chromosome <- gsub(as.character(ExomeCount.dafr$space),pattern = 'chr',replacement = '')
	
# Create matrix
ExomeCount.mat <- as.matrix(ExomeCount.dafr[, grep(names(ExomeCount.dafr), pattern = '*.bam')])
	
nsamples <- ncol(ExomeCount.mat)
samplenames <- colnames(ExomeCount.dafr)
	


	
dir.output <- file.path(dir.base)
	
dir.create(dir.output)
dir.create(file.path(dir.results))
save.image( sprintf("%s/female_exome_depth_count.RData", dir.output) )

	
for (i in 1:nsamples) {
  
          samplename <- samplenames[i + 6]
   
            print(samplename)
            #save sample name to XML
            xml.dom$addTag("sample", close=FALSE)
            xml.dom$addTag("sample_name", samplename)
            
 						#### Create the aggregate reference set for this sample
 						 my.choice <- select.reference.set (test.counts = ExomeCount.mat[,i],
 						                            reference.counts = FemaleExomeCounts.mat,
 						                            bin.length = (ExomeCount.dafr$end - ExomeCount.dafr$start)/1000,
 						                            n.bins.reduced = 10000)
 						 my.reference.selected <- apply(X = FemaleExomeCounts.mat[, my.choice$reference.choice, drop = FALSE],
 						                        MAR = 1,
 						                        FUN = sum)
 						 message('Now creating the ExomeDepth object')
 						 all.exons <- new('ExomeDepth',
 						          test = ExomeCount.mat[,i],
 						          reference = my.reference.selected,
 						          formula = 'cbind(test, reference) ~ 1')
 						 
 						print(my.choice$reference.choice)
 						#save reference choice to XML

 						xml.dom$addTag("reference_choice", my.choice$reference.choice)
 						xml.dom$addTag("summary_stats", my.choice$summary.stats)
 						xml.dom$closeTag()#close sample
 						 
 						################ Now call the CNVs
  
            all.exons <- CallCNVs(x = all.exons,
                        transition.probability = 10^-4,
                        chromosome = ExomeCount.dafr$chromosome,
                        start = ExomeCount.dafr$start,
                        end = ExomeCount.dafr$end,
                        name = ExomeCount.dafr$names)
            if (nrow(all.exons@CNV.calls) > 0) {
	    	          all.exons <- AnnotateExtra(x = all.exons,
                               reference.annotation = v501_exons.GRanges,
                               min.overlap = 0.0001,
                               column.name = 'v501_exons')

  
                  output.file <- paste(dir.results, '/Sample_', samplename, '.csv', sep = '')
                  write.csv(file = output.file, x = all.exons@CNV.calls, row.names = FALSE)

            }
            
            #out_plot_name <- paste(strsplit(samplename,'.bam'),'_ED_results.pdf',sep="")
            #
  					#if (nrow(all.exons@CNV.calls) > 0) {
  					#  for (row in row.names(all.exons@CNV.calls)){
  					#    data <- all.exons@CNV.calls[row,]
  					#    print(data)
  					#    
  					#    title_name <- paste(samplename, data$id,sep=" ")
  					#    pdf(paste(dir.results,'/plots/',gsub(":","-",data$id),"_",out_plot_name, sep=''))
  					#    plot(all.exons,sequence = data$chromosome, xlim = c(data$start-10000,data$end+10000), count.threshold = 20,main = title_name,with.gene=TRUE)
  					#    dev.off()
  					#  }
  					#}
		
} # Samples loop

dir.log <- sprintf("/mnt/Data4/targeted_sequencing/%s/cnv_analysis/female/v501/logs/%s_female_%s.log", sequencing.run, sequencing.run, panel.version)
saveXML(xml.dom$value(), file=dir.logs)

#stop clock
finish <- proc.time() - ptm
print(finish)
summary(finish)
