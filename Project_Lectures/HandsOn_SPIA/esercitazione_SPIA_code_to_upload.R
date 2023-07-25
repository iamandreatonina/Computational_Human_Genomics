

## HOW IS THE DATA?
# selected_snip.txt is a list of SNPs -> only the names
# input_genotype -> the data from the SPIA assay in vitro -> with SNPs number and cell lines
# somatic_genotype -> data from the SPIA -> SNPs and samples




rm(list=ls());verbose<-TRUE; options(max.print=1000);getOption("max.print")

# install.packages("SPIAssay")
## set working dir
# setwd("../Esercitazione_SPIA/")
# getwd()

library(SPIAssay)

## set param for high MAF SNPs -> the parameters of the probabilistic test -> see the SPIAassay.pdf and slides 2023_Human_genomics_slide_deck_5_8.pdf from slide 50
spia_parameters <- list(Pmm = 0.1, nsigma = 2, Pmm_nonM = 0.6, nsigma_nonM = 3, PercValidCall=0.7)
# Pmm and Pmm_notM set the gray area and where the median of the two distribution of paired and not paired sits.
# nsigma is the number 1, 2, 3, etc, not the sogma value but the multiplication value associated to the sigma ehrn calculating teh area of the gaussian for the probability of the distribution -> to define a specifci area

## load genotype data
genotype_data <- read.table("SPIA_input_genotype.tsv", header = TRUE, sep = "\t", as.is = TRUE ); 
dim(genotype_data); head(genotype_data); # print the dimention and the initial rows of the data

## load high MAF SNPs unique IDs -> we use the SPIA_selected_snp file -> only high quality snps
selected_snps <- read.delim("SPIA_selected_SNPs.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)[,1]
length(selected_snps); head(selected_snps); # total of snps = 145

##TASK 1
## select genotype data based on SNP IDs
genotype_data_selected <- genotype_data[selected_snps,] # select only the data fro the high quality snps
dim(genotype_data_selected);

## run test -> SPIATest using the parameters spia_parameters -> t computes SPIA distance and performs probabilistic test on a set of cell lines.
SPIA_selected <-  SPIATest(x = genotype_data_selected, row.names = FALSE, test.param = spia_parameters) 
## visualize distances -> SPIAPlots allows the user to rapidly visualize the result of the SPIA test.
SPIAPlot(SPIA_selected)
## check output matrix
head(SPIA_selected["SPIAresult"]);
##look at distances histogram
hist(as.numeric(SPIA_selected$SPIAresult[,3]),breaks = 100,xlim=c(0,1)) # to have a frequency -> distance goes from 0.3 to 0.8 -> none to 0 meaninig no identical cell lines. Mean around 0.65 -> cell lines are different from each other
## save output matrix
write.table(SPIA_selected["SPIAresult"],paste("output_","selected",".txt",sep=""),sep="\t",quote=FALSE,row.names=FALSE)

##subsampling of high MAF SNPs: 100 -> we take randomly and withut repetition (replace=FALSE) 100 snps 
c100<-sample(c(1:nrow(genotype_data_selected)),100,replace = FALSE);
# we do again the test 
SPIA_selected <-  SPIATest(x = genotype_data_selected[c100,], row.names = FALSE, test.param = spia_parameters) 
SPIAPlot(SPIA_selected) # the only thing that changes is the fact that they seems more distributed
## check output matrix
head(SPIA_selected["SPIAresult"]);
##look at distances histogram
hist(as.numeric(SPIA_selected$SPIAresult[,3]),breaks = 100,xlim=c(0,1)) # to have a frequency -> distance goes from 0.3 to 0.8 -> none to 0 meaninig no identical cell lines. Mean around 0.65 -> cell lines are different from each other
# the histogram changes -> changes the distribution


## TAST 1:  generate and observe the genetic distances obtained from n individuals’ pairs, by using 143, 100, 80, 60, 40, and 20 high MAF SNPs;

##subsampling of high MAF SNPs: 100 -> we take randomly and withut repetition (replace=FALSE) 100 snps 
c140<-sample(c(1:nrow(genotype_data_selected)),140,replace = FALSE);
# we do again the test 
SPIA_selected <-  SPIATest(x = genotype_data_selected[c140,], row.names = FALSE, test.param = spia_parameters) 
SPIAPlot(SPIA_selected) # the only thing that changes is the fact that they seems more distributed
## check output matrix
head(SPIA_selected["SPIAresult"]);
##look at distances histogram
hist(as.numeric(SPIA_selected$SPIAresult[,3]),breaks = 100,xlim=c(0,1)) # to have a frequency -> distance goes from 0.3 to 0.8 -> none to 0 meaninig no identical cell lines. Mean around 0.65 -> cell lines are different from each other
# the histogram changes -> changes the distribution 


##subsampling of high MAF SNPs: 100 -> we take randomly and withut repetition (replace=FALSE) 100 snps 
c80<-sample(c(1:nrow(genotype_data_selected)),80,replace = FALSE);
# we do again the test 
SPIA_selected <-  SPIATest(x = genotype_data_selected[c80,], row.names = FALSE, test.param = spia_parameters) 
SPIAPlot(SPIA_selected) # the only thing that changes is the fact that they seems more distributed
## check output matrix
head(SPIA_selected["SPIAresult"]);
##look at distances histogram
hist(as.numeric(SPIA_selected$SPIAresult[,3]),breaks = 100,xlim=c(0,1)) # to have a frequency -> distance goes from 0.3 to 0.8 -> none to 0 meaninig no identical cell lines. Mean around 0.65 -> cell lines are different from each other
# the histogram changes -> changes the distribution -> for the 

##subsampling of high MAF SNPs: 100 -> we take randomly and withut repetition (replace=FALSE) 100 snps 
c60<-sample(c(1:nrow(genotype_data_selected)),60,replace = FALSE);
# we do again the test 
SPIA_selected <-  SPIATest(x = genotype_data_selected[c60,], row.names = FALSE, test.param = spia_parameters) 
SPIAPlot(SPIA_selected) # the only thing that changes is the fact that they seems more distributed
## check output matrix
head(SPIA_selected["SPIAresult"]);
##look at distances histogram
hist(as.numeric(SPIA_selected$SPIAresult[,3]),breaks = 100,xlim=c(0,1)) # to have a frequency -> distance goes from 0.3 to 0.8 -> none to 0 meaninig no identical cell lines. Mean around 0.65 -> cell lines are different from each other


##subsampling of high MAF SNPs: 100 -> we take randomly and withut repetition (replace=FALSE) 100 snps 
c40<-sample(c(1:nrow(genotype_data_selected)),40,replace = FALSE);
# we do again the test 
SPIA_selected <-  SPIATest(x = genotype_data_selected[c40,], row.names = FALSE, test.param = spia_parameters) 
SPIAPlot(SPIA_selected) # the only thing that changes is the fact that they seems more distributed
## check output matrix
head(SPIA_selected["SPIAresult"]);
##look at distances histogram
hist(as.numeric(SPIA_selected$SPIAresult[,3]),breaks = 100,xlim=c(0,1)) # to have a frequency -> distance goes from 0.3 to 0.8 -> none to 0 meaninig no identical cell lines. Mean around 0.65 -> cell lines are different from each other
# the histogram changes -> changes the distribution -> for the 


##subsampling of high MAF SNPs: 20 -> NON FUNZIONAA -> NON PENSATO PER LAVORARE CON POCHI SNPs
c20<-sample(c(1:nrow(genotype_data_selected)),20,replace = FALSE);
# we do again the test 
param=list(Pmm = 0.01, nsigma = 1, Pmm_nonM = 0.01, nsigma_nonM = 2, PercValidCall=0.7)
SPIA_selected <-  SPIATest(x = genotype_data_selected[c20,], row.names = FALSE, test.param = param) 
SPIAPlot(SPIA_selected) # the only thing that changes is the fact that they seems more distributed
## check output matrix
head(SPIA_selected["SPIAresult"]);
##look at distances histogram
hist(as.numeric(SPIA_selected$SPIAresult[,3]),breaks = 100,xlim=c(0,1)) # to have a frequency -> distance goes from 0.3 to 0.8 -> none to 0 meaninig no identical cell lines. Mean around 0.65 -> cell lines are different from each other
# the histogram changes -> changes the distribution -> for the 



##TASK 2 generate and observe the genetic distances among m individuals using 100 unselected SNPs
SPIA_all <-  SPIATest(x = genotype_data, row.names = F, test.param = spia_parameters) 
SPIAPlot(SPIA_all)
c100<-sample(c(1:nrow(genotype_data)),100,replace = FALSE); # select 100 random SNPs among all the 2300 SNPs
genotype_data_random <- genotype_data[c100,]
SPIA_random <-  SPIATest(x = genotype_data_random, row.names = F, test.param = spia_parameters) 
SPIAPlot(SPIA_random)
hist(as.numeric(SPIA_random$SPIAresult[,3]),breaks = 100,xlim=c(0,1))
# avarage distance of 0.3 -> because the snp are not highly selected so the random selected snps could be snps that are common among all individualso we do not see high differencesnot good set of snps

##TASK 3 -> generate and observe the genetic distances among a set of samples that includes j normal/tumor pairs
## load somatic data
somatic_data <- read.table("SPIA_somatic_genotype.tsv", header = TRUE, sep = "\t", as.is = TRUE ); 
dim(somatic_data); head(somatic_data); # 334   7 # the data is divided in normal and tumor -> normal: 1_1, 1_2; tumor: 2_1, 2_2, etc till 2_5.
SPIA_somatic <-  SPIATest(x = somatic_data, row.names = F, test.param = spia_parameters)
SPIAPlot(SPIA_somatic) # we can see that the low distance points correspond to the comparison between normal samples or tumor samples, while the high distance is for the comparison among tumors and control samples
# look at SPIAresult
write.table(SPIA_somatic["SPIAresult"],paste("output_","somatic",".txt",sep=""),sep="\t",quote=FALSE,row.names=FALSE)
# in the 

##TASK 4
## chnage parameters
## set param for high MAF SNPs, exploring different values for Pmm, Pmm_nonM, nsigma, nsigma_nonM (one at the time)
spia_parameters <- list(Pmm = 0.2, nsigma = 2, Pmm_nonM = 0.6, nsigma_nonM = 3, PercValidCall=0.7)
# ...

