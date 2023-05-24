library(data.table)
library(CLONETv2)
library(TPES)

setwd("~/Documents/HumanGenomics/10_PurityPloidyEstimation/")

normal = fread("Normal.csv",data.table=F) # WE ADD OUR DATA
normal$af = normal$altCount/normal$totalCount # WE CALCULATE MANUALLY TEH ALLELEIC FRACTION BASED ON TEH input data
tumor = fread("Tumor.csv",data.table=F)
tumor$af = tumor$altCount/tumor$totalCount

pileup.normal = normal[,c(1,2,4,5,14,8)] 
colnames(pileup.normal) = c("chr","pos","ref","alt","af","cov")

pileup.tumor = tumor[,c(1,2,4,5,14,8)]
colnames(pileup.tumor) = c("chr","pos","ref","alt","af","cov")

seg.tb <- fread("../09_SomaticCopyNumberCalling/Data/SCNA.copynumber.called.seg",data.table=F) ## we upload teh data of the segmentation we opertaed last lesson
# this info is organized like tgis: name of teh sample, position, mean of teh segments -> give us info on the segmentation and also info on the copy number for each region. The mean of teh sample teel us if we are losing or gaining an allele.
# we have info on teh segmentation and the coverage


# we put all the data togheter using clonet
bt <- compute_beta_table(seg.tb, pileup.tumor, pileup.normal) # create a beta table with teh info on teh segmentation and the two pileup
# in bt we have information for each segments -> clonet takes all teh samples and all the heterozygus snps for each segments. The AF is used to estimate beta.
# we have the beta, the number of snps and coverage. Plus a column named n_bita that correspod to teh beta of teh normal samples -> need to be close to 1.

## Compute ploidy table with default parameters
pl.table <- compute_ploidy(bt) # calculates the ploidy of our sample, starting from the bt object.
# this object tell us tehe ploidy of teh sampe, here it is 1.84, so basically 2 -> this is a diploid sample! It is rare to observe samples with 1.5 -> less dna than the expected

adm.table <- compute_dna_admixture(beta_table = bt, ploidy_table = pl.table)
# create an admixture table with an error associated to the admizture -> here we have an admixture of 0.8 -> means a sample with high quality, specifically high purity

# we can calculate allele specific copy number events. We give in input the beta table an the info on the ploidy, plus teh admixture table
allele_specific_cna_table <- compute_allele_specific_scna_table(beta_table = bt,
                                                                ploidy_table = pl.table, 
                                                                admixture_table = adm.table)
# give us a log2R corrected -> the software correct the observation of the covariets because now we now the purity and ploidy of teh samples, therefore now we can calculate exactly teh number of copies for each allele in each segment
# the estimation is not always the correct number -> we can understand for the allele if we have an amplification or a deletion -> example for sample1 we have cnA=6 -> gain of 5 copies and cnB = 0 -> deletion -> we have a complex event


# we can now plot the ploidy and check teh data
check.plot <- check_ploidy_and_admixture(beta_table = bt, ploidy_table = pl.table,
                                         admixture_table = adm.table)
print(check.plot)
# the plot is not so beautiful -> pretty much noisy -> we have someinfo on teh ploidy and the mizture, we have on the x axis the LogR info (=copy number info) while on the y axis we have the beta values (=info on the alleleic fraction). Each point is a different segment, th values are therefore calculated on teh whole segment.
# teh red dots are the point in teh space inwhich we expect specific events to fall, given that specfic purity and ploidy -> wild type events (=one copy of one allele, one copie of teh other one), abberating events -> we have an idea on where our abberation need to fall in 
# the better the clustering process, the more they are near teh red dots and more clustered.
# the red points take already into account teh admixture and purity
# basd on teh cloud of teh hemizigus cluster and specifically based on its position we can understand teh purity

# TPES -> do the same of what CLONET does but for SNVs data
# we get the snv data from the lesson 8 
snv.reads = fread("../08_SomaticVariantCalling/Data/somatic.pm",data.table=F)
snv.reads = snv.reads[which(snv.reads$somatic_status=="Somatic"),] # we filter teh events and take only the somatic one
snv.reads = snv.reads[,c("chrom","position","position","tumor_reads1","tumor_reads2")]
colnames(snv.reads) = c("chr","start","end","ref.count","alt.count")
snv.reads$sample = "Sample.1"
# snv_reads contains teh position , the ref allele and the alternative allele counts

TPES_purity(ID = "Sample.1", SEGfile = seg.tb,
            SNVsReadCountsFile = snv.reads,
            ploidy = pl.table,
            RMB = 0.47, maxAF = 0.6, minCov = 10, minAltReads = 10, minSNVs = 1)
# RMB is teh reference variance bias -> a technical bias we need to take in account when making this kind of estimation

# tpes report show us all the distribution for the SNVs we looked at -> based on teh distribution of the allelic fraction, it is possible to estimate teh purity
TPES_report(ID = "Sample.1", SEGfile = seg.tb,
            SNVsReadCountsFile = snv.reads,
            ploidy = pl.table,
            RMB = 0.47, maxAF = 0.6, minCov = 10, minAltReads = 10, minSNVs = 1)
