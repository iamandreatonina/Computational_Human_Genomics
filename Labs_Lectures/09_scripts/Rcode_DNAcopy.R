##Load library: exec the commented lines to install

#source("http://bioconductor.org/biocLite.R")
#biocLite("DNAcopy")
library(DNAcopy)

##Generate a synthetic sample with some potential segments 

#number of snps for each segment
lengths <- c(10,10,10,10,10)
#mean value for each segment
values <- c(-0.1,-0.05,0.1,0.2,0.15)
#standard deviation of the noise
sdnoise <- 0.05

#genomic data = data + noise
genomdat <- rnorm(sum(lengths), sd=sdnoise) +  rep(values,lengths)
# with rep we repeat the values 10 times -> because lenghts is equal to 10.
# The genomic distribution is obtained using a normal distribution qith rnorm and adding some noise

#chromosomes 
chrs <- rep(1,sum(lengths))
maploc <- c(1:sum(lengths))

#generate a sample with seglength SNPs using normal distributions -> generate a toy sample. We generate an object that is used for further analysis
sample <- CNA(genomdat, #log2 values
              chrs, #chr
              maploc, #positions              
              data.type="logratio", #type of data
              sampleid="test") #name of the sample)
# calling sample we get the number of samples=1, number of probes=50 in datatype equal to log_ratio
# if we want to llok inside the object sample we can call View or str(sample). Sample is a dataframe

plot(maploc,genomdat,xlab="Chromosomal position",ylab="log2 ratio",main="Array data",
     pch=".",col="green",cex=4)
# we get the plt -> synthetically generated data -> random data that can be used to perform analysis on CNV and segmentation

#smooth the signal -> we can filter out some data
sample.smoothed <- smooth.CNA(sample) 

#segment alpha 0.1 -> we can perform the segmentation -> we generate a new object with the command segment -> in input it needs: the sample and alpga, a parameter needed for the segmentation taht rapresent the significative level
sample.segmented <- segment(sample.smoothed, #a CNA object
                            alpha=0.1, #significance levels for the test to accept change-points
                            nperm=100, #number of permutations used for p-value computation
                            verbose = 1) #print program progess
# the object sample.segmented is a test -> we have a table in which we have the start and end of the segments, the number of mark and the mean of teh segment 

# we can plot the segment -> as seen in the table, we have 4 segments
plot(sample.segmented,plot.type="w",ylab="log2 ratio",main="Array data segmented")

#segment alpha 1 -> changing the parameter alpha -> since here it is equal to 1 the significant level is less stringent -> more segmets! Alpha is the significative level to accept a change point
# alpha equal to 1 is too big!! We need a stringet alpha butnot too stringent though
sample.segmented.alpha1 <- segment(sample.smoothed, alpha=1,nperm=100)
plot(sample.segmented.alpha1,plot.type="w")

#segment alpha 0.01
sample.segmented.alpha0.01 <- segment(sample.smoothed, alpha=0.01,nperm=100)
plot(sample.segmented.alpha0.01,plot.type="w")


#undo -> we can undo the split / segments done
sdundo.sample.segmented <- segment(sample.smoothed,
                               undo.splits="sdundo", #undoes splits that are not at least this many SDs apart -> there is a threshold in the strandard deviation of the difference between segments that are/can be eliminated -> if high value for teh standard deviation (like 3) we obtain only one segment
                               undo.SD=3,verbose=1)

plot(sdundo.sample.segmented,plot.type="w")


#############################
##
## Load Real samples -> perform the analysis on real data from a .csv
log2data <- read.csv("AllLog2Sub.cn_medianCentered",sep="\t",header=F)
# log2data is a dataframe containing 9029 observation with info on the position on the genome (in the chromosome) and annotations -> total columns equal to 430

#select all samples  -> we left out the first three columns because they contains info on the position on teh chr not data of interest on the samples
samples_ind <- c(4:(dim(log2data)[2]))

#select N samples randomly
N<-200

# we sample a certain number of random values to selct random samples
samples_ind <- sample(samples_ind,length(samples_ind))[1:N]

genomedat <- cbind(apply(log2data[,samples_ind],2,unlist))

# run CNA on this data
samples <- CNA(genomedat, log2data[,2], log2data[,3], 
               data.type="logratio", 
               sampleid=paste("sample", c(1:(dim(genomedat)[2])),sep=""))

# we smooth out
samples.smoothed <- smooth.CNA(samples)

# we segment, indicating the threshold on alpha
samples.segmented <- segment(samples.smoothed, #a CNA object
                             alpha=0.05,#significance levels for the test to accept change-points
                             nperm=100, #number of permutations used for p-value computation
                             verbose=1) #print program progess

segFile <- samples.segmented$output # the output of the object segmented correspond to the segments -> we have 3580 segments

# save for IGV -> we can load this data on IGV
write.table(segFile,
            file="AllLog2Sub.seg",sep="\t",quote=F, 
            col.names = T,row.names=F)

# Distribution of log2_ratios of all segments detected in a sample set
hist(samples.segmented$output[,"seg.mean"],xlab="Sement Mean Log2 Ratio",breaks=100,
     main="CBS, alpha=0.05")

# we create function that filters the data for a specific position given in input
getSegments <- function(chr, initialPos, finalPos, segments){
  #select segments that intersect initialPos and finalPs
  goodSegs <- segments[which( (segments[,2] == chr) & (
    (segments[,3] >= initialPos & segments[,3] <= finalPos ) |
      (segments[,4] >= initialPos & segments[,4] <= finalPos ) |
      (segments[,3] <= initialPos & segments[,4] >= finalPos ) )),]
  
  return(goodSegs)  
}

#Distribution of one CNV across many individuals

# select gene in chromosome 15 in pos: 18748000,18758000
CHEK2P2segs<- getSegments(15,18748000,18758000,segFile)

# we can plot the data only for this gene
# from the histogram we can see if there are some delitions or insertion -> we see the distibution of teh reads therefore we can understand
hist(CHEK2P2segs[,"seg.mean"],xlab="Sement Mean Log2 Ratio",breaks=100,
     main="CBS, alpha=0.01")


OR4N4segs<- getSegments(15,19883836,19885179,segFile)
hist(OR4N4segs[,"seg.mean"],xlab="Sement Mean Log2 Ratio",breaks=400,
     main="CBS, alpha=0.01")

OR4N3Psegs<- getSegments(15,19914825,19915759,segFile)
hist(OR4N3Psegs[,"seg.mean"],xlab="Sement Mean Log2 Ratio",breaks=200,
     main="CBS, alpha=0.01",freq=T)


UGT2B17seg <- getSegments(4,69085498,69116840,segFile)
hist(UGT2B17seg[,"seg.mean"],xlab="Sement Mean Log2 Ratio",breaks=100,
     main="CNV across many individuals\nGene UGT2B17")



