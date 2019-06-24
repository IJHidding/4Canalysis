###Prep
## Loading all libraries used in the script
library (FourCSeq)
library(DESeq2)
library(ggplot2)
library(GenomicFeatures)
library(rtracklayer)
library(BSgenome)
library(AnnotationHub)

### Creating the reference genome file in TxDb and then emptying the variable to save memory
hub <- AnnotationHub()
query(hub, c("gtf","Homo sapiens"))
gtf <- hub[["AH7790"]]
tx <- makeTxDbFromGRanges(gtf)
gtf <- NULL


##Create metadata file
bamFilePath <- "Path/To/BamFiles"
referenceGenomeFile <- "Path/To/ReferenceGenome"
primerFile <- "Path/To/Primerfile.fa"
ViewPointer<- "ViewpointPrimername"
####Has to match the name in the file
ConDition <- "Condition if any"

#Create the actual metadata file   ####DATA BREACH####
metadata <- list(projectPath = "Path/To/All/Files",
                 fragmentDir = "Path/To/Subfolder/ForRestrictionfragments", #subfolder for saving restriction fragment data
                 referenceGenomeFile = referenceGenomeFile,
                 reSequence1 = "Restriction site used in the experiment 1", 
                 reSequence2 = "Restriction site used in the experiment 1", 
                 primerFile = primerFile,
                 bamFilePath = bamFilePath)

#Enter data based on number of replicates and conditions
colData <- DataFrame(viewpoint = ViewPointer,
                     condition = factor(ConDition),
                     replicate = rep(c(1, 2, 3),
                                     1),
                     bamFile = c("Name of the bamefile1",
                                 "Name of the bamefile2",
                                 "Name of the bamefile3"),
                     sequencingPrimer="first")


###Creating the objects 
fc <- FourC(colData, metadata)

fc <- addFragments(fc)

rowRanges(fc)

###Long step -- creating the fragment map.
findViewpointFragments(fc)

fc <- addViewpointFrags(fc)

#### trim off the restriction site on the read.
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
write.table(minnumber, "Path/To/Outputfile.txt", sep ="\t", dec=".", row.names = F, col.names = T, append=T)
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

##Adds info from the primerinfo
viewpoint.info <- read.table(file = "Path/To/PrimerFragments.txt", header=T, sep="\t", dec=".")

write.table(viewpoint.info, "Path/To/Outputfile.txt", sep ="\t", dec=".", row.names = F, col.names = T, append=T)

write.table(total.results, "Path/To/Outputfile.txt", sep ="\t", dec=".", row.names = F, col.names = T, append=T)
