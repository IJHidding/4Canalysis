#!/usr/bin/Rscript
library(peakC)
args <- commandArgs(trailingOnly = TRUE)
ViewpointGene <- args[1]
viewpoint <- args[2]
Outputfile <- args[3]
FileDirectory <- args[4]
Pattern <- args[5]
Outputtxt <- args[6]
Outputtxt 





#Pattern

viewpoint <- as.numeric(viewpoint)



files <- dir(path=FileDirectory, pattern=Pattern, full=TRUE)
data <- readMultiple(files, vp.pos=viewpoint)
res <- combined.analysis(data, vp.pos=viewpoint)
write.table(res$peak, Outputtxt, sep ="\t", dec=".", row.names = F, col.names = T, append=T)
setwd(Outputfile)
jpeg('Plot.jpg')
plot_C(res)

dev.off() 




