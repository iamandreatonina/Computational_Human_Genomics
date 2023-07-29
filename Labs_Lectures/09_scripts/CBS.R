library(DNAcopy)
folder = ""
cn <- read.table(file.path(folder,"SCNA.copynumber.called"),header=T) # loading the data

#pdf("~/Documents/HumanGenomics/09_SomaticCopyNumberCalling/Data/SegPlot.pdf")

plot(cn$raw_ratio,pch=".",ylim=c(-2.5,2.5)) # the data we are analyzing -> data with noise -> we can see by eye that there is a segmentation -> some values are from amplifications other from deletions
plot(cn$adjusted_log_ratio,pch=".",ylim=c(-2.5,2.5)) # same as above, just looking at the log2 p_adjusted values


CNA.object <-CNA(genomdat = cn$adjusted_log_ratio, # we generate a CNA object, as done in the other R script
                 chrom = cn$chrom,
                 maploc = cn$chr_start, data.type = 'logratio')

CNA.smoothed <- smooth.CNA(CNA.object) # perform a smoothing

# generate the segments, we are not specifiyng the alpha but we are using the default one.
segs <- segment(CNA.smoothed, min.width=2,
                undo.splits="sdundo", #undoes splits that are not at least this many SDs apart.
                undo.SD=3,verbose=1)
# the object segs is big

# lets plot the segmentation results
plot(segs,plot.type="w")

dev.off()

segs2 = segs$output
write.table(segs2, file=file.path(folder,"SCNA.copynumber.called.seg"), row.names=F, col.names=T, quote=F, sep="\t")
