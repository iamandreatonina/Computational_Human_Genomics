folder = "~/Documents/HumanGenomics/07_AncestryAnalysis/Data/"
setwd(folder)

### 100 samples analysis -> upload teh results of the prevoius analysis in R using read.table -> we have the three files
ethseq.coord = read.table(file.path(folder,"ethseq/100s/Report.PCAcoord"),sep="\t",as.is=T,header=T) # each line is a sample, each coord is a dimention -> a colum
ethseq.report = read.table(file.path(folder,"ethseq/100s/Report.txt"),sep="\t",as.is=T,header=T) # report
smartpca.coord = read.table(file.path(folder,"smartpca/1000GP_Genotypes100s.pca.evec"),as.is=T) # each line is a sample, plus we have different variable dimentions -> coordinates in the space

## check samples order -> to be sure to compare sample 1 against sample 1 -> are the sample in the same position in all the files?
all(ethseq.coord[,1]==smartpca.coord[,1])
all(ethseq.report[,1]==smartpca.coord[,1])

## plot PCA space -> we can plot the PCA of teh two (SMARTPCA and EthSEQ) and compare them
par(mfrow=c(2,3),mar=c(4,4,4,1)) # to plot all toghether-> generate a window in which the plots will be included
plot(ethseq.coord$EV1,ethseq.coord$EV2,pch=19,xlab="PCA1",ylab="PCA2",main="EthSEQ") 
plot(ethseq.coord$EV2,ethseq.coord$EV3,pch=19,xlab="PCA2",ylab="PCA3",main="EthSEQ")
plot(ethseq.coord$EV1,ethseq.coord$EV3,pch=19,xlab="PCA1",ylab="PCA3",main="EthSEQ")
plot(smartpca.coord[,2],smartpca.coord[,3],pch=19,xlab="PCA1",ylab="PCA2",main="SMARTPCA")
plot(smartpca.coord[,3],smartpca.coord[,4],pch=19,xlab="PCA2",ylab="PCA3",main="SMARTPCA")
plot(smartpca.coord[,2],smartpca.coord[,4],pch=19,xlab="PCA1",ylab="PCA3",main="SMARTPCA")
# on the bottom we have the plots for the SMARTPCA, on the top we have the plots for the EthSEQ.
# For each we look at the first three dimentions -> explain the highes percentage of variance of teh data -> teh most informative dimentions -> the other dimentions explain less
# we can see in all the plots 4 clusters, both in EthSEQ and SMARTPCA -> cearer definition in EthSEQ though. The PCA component 1 is not exactly teh same of SMARTPCA component 1 but in teh end what is found is the same.

# comparing the PC component of the two tools
par(mfrow=c(1,3),mar=c(4,4,4,1))
plot(ethseq.coord$EV1,smartpca.coord[,2],pch=19,xlab="PCA1 (EthSEQ)",ylab="PCA1 (SMARTPCA)",main="EthSEQ/SMARTPCA")
plot(ethseq.coord$EV2,smartpca.coord[,3],pch=19,xlab="PCA2 (EthSEQ)",ylab="PCA2 (SMARTPCA)",main="EthSEQ/SMARTPCA")
plot(ethseq.coord$EV3,smartpca.coord[,4],pch=19,xlab="PCA3 (EthSEQ)",ylab="PCA3 (SMARTPCA)",main="EthSEQ/SMARTPCA")
# PC1 of the two tolls are the same, PC2 are slightly different, for PC3 we can see a correlation. Going up in the dimension (=less variance explained), we have less similarity between the two tools

## compare fastSTRUCTURE and EthSEQ report
faststructure = read.table(file.path(folder,"faststructure/1000GP_Genotypes100s.4.meanQ"),as.is=T,header=F)
ids =  read.table(file.path(folder,"faststructure/1000GP_Genotypes100s.fam"),as.is=T,header=F)
faststructure$id = ids[,2]

sel = ethseq.report$pop[which(ethseq.report$sample.id%in%faststructure[which(faststructure[,1]>=0.99),5])] # take the ids of the first population in EthSEQ and understand where these ids are in SMARTPCA -> are in the faststructure?
unique(sel)
sel = ethseq.report$pop[which(ethseq.report$sample.id%in%faststructure[which(faststructure[,1]<0.99),5])]
unique(sel)

table(sel) 

# to visualize lets use GGpolt
library(ggplot2)
library(reshape2)



# we can see that there are big chunks asscoiated to population 1 pr 2 and only in a small number there is 
