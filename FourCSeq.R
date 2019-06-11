###Prep
library (FourCSeq)
library(DESeq2)
library(ggplot2)
library(GenomicFeatures)
library(rtracklayer)
library(BSgenome)
library(AnnotationHub)
hub <- AnnotationHub()
query(hub, c("gtf","Homo sapiens"))
gtf <- hub[["AH7790"]]
tx <- makeTxDbFromGRanges(gtf)
gtf <- NULL


##Create metadata file
bamFilePath <- "C:/4C/kims_data/BWAalignment"
referenceGenomeFile <- "C:/4C/Final_Analysis/Data/human_g1k_v37.fasta"
primerFile <- "C:/4C/kims_data/4CPrimer_TTN_fw.fa"
ViewPointer<- "TTN_fw"
ConDition <- "Cardio"

#Create the actual metadata file   ####DATA BREACH####
metadata <- list(projectPath = "C:/4C/Final_Analysis/Data/",
                 fragmentDir = "Restriction_Fragments", #subfolder for saving restriction fragment data
                 referenceGenomeFile = referenceGenomeFile,
                 reSequence1 = "GATC", 
                 reSequence2 = "CATG", 
                 primerFile = primerFile,
                 bamFilePath = bamFilePath)

colData <- DataFrame(viewpoint = ViewPointer,
                     condition = factor(ConDition),
                     replicate = rep(c(1, 2, 3),
                                     1),
                     bamFile = c("reTTN_1.bam",
                                 "reTTN_2.bam",
                                 "reTTN_3.bam"),
                     sequencingPrimer="first")



fc <- FourC(colData, metadata)

fc <- addFragments(fc)

rowRanges(fc)

findViewpointFragments(fc)

fc <- addViewpointFrags(fc)

#### trim off the restriction site on the read. sometimes change this to 0 causes more reads to be aligned to the bins. 
#### Probably caused by data without restriction site. 
fc <- countFragmentOverlaps(fc, trim=4, minMapq=30)

fc <- combineFragEnds(fc)

#Function that alters the minCount variable if the function getZScores fails due to low quality data.
minnumber <- 45
boolFalse<-F
while(boolFalse==F)
{
  tryCatch({
    minnumber<- minnumber-5
    getZScores(fc, minCount = minnumber)
    boolFalse<-T
  }, error = function(e){
  }, finally={})
}



fcf <- getZScores(fc,minCount = minnumber)

#Here the minCount gets recorded so it is in the output file.
write.table(minnumber, "C:/4C/Enhancer_Output_data.txt", sep ="\t", dec=".", row.names = F, col.names = T, append=T)
zScore <- assay(fcf, "zScore")

#Quality assasment plots
hist(zScore[,1], breaks=100)
hist(zScore[,2], breaks=100)
hist(zScore[,3], breaks=100)


qqnorm(zScore[,1],
     main="Normal Q-Q Plot 1")
abline(a=0, b=1)

qqnorm(zScore[,2],
     main="Normal Q-Q Plot 2")
abline(a=0, b=1)

qqnorm(zScore[,3],
     main="Normal Q-Q Plot 3")
abline(a=0, b=1)

fcf <- addPeaks(fcf, zScoreThresh=3, fdrThresh=0.01) 

plotZScores(fcf[,1:3],
            txdb=tx,
            plotSingle = F)

#Storing the called peaks in an output file
object <- granges(fcf)
object <- as.data.frame(object)
dataframeZ <- as.data.frame(zScore)
total <-  cbind(dataframeZ, object) 
total.results <- total[total[,1] > 3 & total[,2] > 3 & total[,3] > 3,]
viewpoint.info <- read.table(file = "C:/4C/Final_Analysis/Data/Restriction_Fragments/primerFragments.txt", header=T, sep="\t", dec=".")

write.table(viewpoint.info, "C:/4C/Enhancer_Output_data.txt", sep ="\t", dec=".", row.names = F, col.names = T, append=T)

write.table(total.results, "C:/4C/Enhancer_Output_data.txt", sep ="\t", dec=".", row.names = F, col.names = T, append=T)
