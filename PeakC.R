library(peakC)
data <- readMatrix("C:/4C/peakC/Gecombineerdefiles/BAG3_1.txt", vp.pos= 121416630, window = 700e3)
res <- single.analysis(data$data, vp.pos=121416630, qWd=2.5)
plot_C(res)
title(main="BAG3")

data <- readMatrix("C:/4C/peakC/Gecombineerdefiles/BAG3_2.txt", vp.pos= 121416630, window = 700e3)
res <- single.analysis(data$data, vp.pos=121416630, qWd=2.5)
plot_C(res)
title(main="BAG3")

data <- readMatrix("C:/4C/peakC/Gecombineerdefiles/BAG3_3.txt", vp.pos= 121416630, window = 700e3)
res <- single.analysis(data$data, vp.pos=121416630, qWd=2.5)
plot_C(res)
title(main="BAG3")


viewpoint <- 121416630
files <- dir(path="C:/4C/peakC/Gecombineerdefiles/", pattern="BAG3_[123].txt", full=TRUE)
data <- readMultiple(files, vp.pos=viewpoint)
res <- combined.analysis(data, vp.pos=viewpoint)
plot_C(res)
title(main="BAG3 with replicates")
```
```{r single experiment}
data <- readMatrix("C:/4C/peakC/Gecombineerdefiles/TTN_1.txt", vp.pos= 179673573, window = 700e3)
res <- single.analysis(data$data, vp.pos=179673573, qWd=2.5)
plot_C(res)
title(main="TTN")

data <- readMatrix("C:/4C/peakC/Gecombineerdefiles/TTN_2.txt", vp.pos= 179673573, window = 700e3)
res <- single.analysis(data$data, vp.pos=179673573, qWd=2.5)
plot_C(res)
title(main="TTN")

data <- readMatrix("C:/4C/peakC/Gecombineerdefiles/TTN_3.txt", vp.pos= 179673573, window = 700e3)
res <- single.analysis(data$data, vp.pos=179673573, qWd=2.5)
plot_C(res)
title(main="TTN")


viewpoint <- 179673573
files <- dir(path="C:/4C/peakC/Gecombineerdefiles/", pattern="TTN_[123].txt", full=TRUE)
data <- readMultiple(files, vp.pos=viewpoint)
res <- combined.analysis(data, vp.pos=viewpoint)
plot_C(res)
title(main="TTN with replicates")
```


```{r Combined data}
library(peakC)
viewpoint <- 121416630
files <- dir(path="C:/4C/peakC/Gecombineerdefiles/", pattern="BAG3_[123].txt", full=TRUE)
data <- readMultiple(files, vp.pos=viewpoint)
res <- combined.analysis(data, vp.pos=viewpoint)
head(res$peak)
plot_C(res)
title(main="BAG3")

viewpoint2 <- 179673573
files2 <- dir(path="C:/4C/peakC/Gecombineerdefiles/", pattern="TTN_[123].txt", full=TRUE)
data2 <- readMultiple(files2, vp.pos=viewpoint2)
res2 <- combined.analysis(data2, vp.pos=viewpoint2)
head(res2$peak)
plot_C(res2)
title(main="TTN")

viewpoint3 <- 201346777
files3 <- dir(path="C:/4C/peakC/Gecombineerdefiles/", pattern="TTN2_[123].txt", full=T)
data3 <- readMultiple(files3, vp.pos=viewpoint3)
res3 <- combined.analysis(data3, vp.pos=viewpoint3)
head(res3$peak)
plot_C(res3)
title(main="TNNT2")

viewpoint4<- 201341573
files4 <- dir(path="C:/4C/peakC/Gecombineerdefiles/", pattern="TTN2_2_[123].txt", full=TRUE)
data4 <- readMultiple(files4, vp.pos=viewpoint4)
res4 <- combined.analysis(data4, vp.pos=viewpoint4)
head(res4$peak)
plot_C(res4)
title(main="TNNT2_2")

viewpoint5 <- 23865593
files5 <- dir(path="C:/4C/peakC/Gecombineerdefiles/", pattern="MYH6_[123].txt", full=TRUE)
data5 <- readMultiple(files5, vp.pos=viewpoint5)
res5 <- combined.analysis(data5, vp.pos=viewpoint5)
head(res5$peak)
plot_C(res5)
title(main="MYH6")


viewpoint6 <- 23901458
files6 <- dir(path="C:/4C/peakC/Gecombineerdefiles/", pattern="MYH7_[123].txt", full=TRUE)
data6 <- readMultiple(files6, vp.pos=viewpoint6)
res6 <- combined.analysis(data6, vp.pos=viewpoint6)
head(res6$peak)
plot_C(res6, y.max = 5000)
title(main="MYH7")

viewpoint8 <- 156052727
files8 <- dir(path="C:/4C/peakC/Gecombineerdefiles/", pattern="LMNA1_[123].txt", full=TRUE)
data8 <- readMultiple(files8, vp.pos=viewpoint8)
res8 <- combined.analysis(data8, vp.pos=viewpoint8)
head(res8$peak)
plot_C(res8)
title(main="LMNA1")

viewpoint7 <- 156086255
files7 <- dir(path="C:/4C/peakC/Gecombineerdefiles/", pattern="LMNA2_[123].txt", full=TRUE)
data7 <- readMultiple(files7, vp.pos=viewpoint7)
res7 <- combined.analysis(data7, vp.pos=viewpoint7)
head(res7$peak)
plot_C(res7)
title(main="LMNA2")


#### 
#### 20 mei concept verslag
#### verslag begin juni
######Afstuderen 20 juni half 1 presentatie - kwart over 1
```


```{r }

library(peakC)
viewpoint <- 121416630
files <- dir(path="C:/4C/peakC/Gecombineerdefiles/", pattern="BAG3_[123].txt", full=TRUE)
data <- readMultiple(files, vp.pos=viewpoint)
res <- combined.analysis(data, vp.pos=viewpoint)
plot_C(res)
title(main="BAG3 peakC")

for (peaks in c(121357410,	121367756,	121403361,	121406094,	121407472,	121412863,	121413999, 121356399,	121367193,	121402787,	121405331,	121407302,121412102,	121413219)){
  FourCpeakinPeakC <- which(abs(res$dbR$frag_pos-peaks)==min(abs(res$dbR$frag_pos-peaks)))
  ActualPeakCpos <- c(ActualPeakCpos, res$dbR$frag_pos[FourCpeakinPeakC])
  
}

res$peak <- ActualPeakCpos
plot_C(res)
title(main="BAG3 FourCSeq")

######################################################################################################################
viewpoint2 <- 179673573
files2 <- dir(path="C:/4C/peakC/Gecombineerdefiles/", pattern="TTN_[123].txt", full=TRUE)
data2 <- readMultiple(files2, vp.pos=viewpoint2)
res2 <- combined.analysis(data2, vp.pos=viewpoint2)
plot_C(res2)
title(main="TTN PeakC")

#fourCpeaks <- c(179661092, 179685470, 179698789, 179702794, 179727086, 179729901, 179733185, 179741176, 179741528, 179745434, 179775810, 179785901)

for (peaks in c(179664586,	179664888,	179672190, 179673064, 179685470,	179698789,	179699999,	179702794,	179709915,	179720513,	179727086,	179729008,	179729901,	179730323,	179733185,	179741176,	179741528,	179745434,	179752851, 179664883,	179666902,	179673059,	179673587,	179685776,	179699367,	179700732,	179703358,	179710840,	179720777,	179727692,	179729887,	179730318,	179731515,	179734483,	179741523,	179741822,	179746531,	179752958)){
  FourCpeakinPeakC <- which(abs(res2$dbR$frag_pos-peaks)==min(abs(res2$dbR$frag_pos-peaks)))
  ActualPeakCpos <- c(ActualPeakCpos, res2$dbR$frag_pos[FourCpeakinPeakC])
  
}

res2$peak <- ActualPeakCpos
plot_C(res2)
title(main="TTN FourCSeq")
######################################################################################################################

viewpoint3 <- 201346777
files3 <- dir(path="C:/4C/peakC/Gecombineerdefiles/", pattern="TTN2_[123].txt", full=T)
data3 <- readMultiple(files3, vp.pos=viewpoint3)
res3 <- combined.analysis(data3, vp.pos=viewpoint3)
plot_C(res3)
title(main="TNNT2 peakC")

for (peaks in c()){
  FourCpeakinPeakC <- which(abs(res3$dbR$frag_pos-peaks)==min(abs(res3$dbR$frag_pos-peaks)))
  ActualPeakCpos <- c(ActualPeakCpos, res3$dbR$frag_pos[FourCpeakinPeakC])
  
}

res3$peak <- ActualPeakCpos
plot_C(res3)
title(main="TNNT2 FourCSeq")

######################################################################################################################

viewpoint4<- 201341573
files4 <- dir(path="C:/4C/peakC/Gecombineerdefiles/", pattern="TTN2_2_[123].txt", full=TRUE)
data4 <- readMultiple(files4, vp.pos=viewpoint4)
res4 <- combined.analysis(data4, vp.pos=viewpoint4)
plot_C(res4)
title(main="TNNT2_2 PeakC")

for (peaks in c(201336079,	201338149,	201339288,	201345697,	201337780,	201338717,	201339740,	201346772)){
  FourCpeakinPeakC <- which(abs(res4$dbR$frag_pos-peaks)==min(abs(res4$dbR$frag_pos-peaks)))
  ActualPeakCpos <- c(ActualPeakCpos, res4$dbR$frag_pos[FourCpeakinPeakC])
  
}
res4$peak <- ActualPeakCpos
plot_C(res4, y.max = 3000)
title(main="TNNT2_2 FourCSeq")

######################################################################################################################

viewpoint5 <- 23865593
files5 <- dir(path="C:/4C/peakC/Gecombineerdefiles/", pattern="MYH6_[123].txt", full=TRUE)
data5 <- readMultiple(files5, vp.pos=viewpoint5)
res5 <- combined.analysis(data5, vp.pos=viewpoint5)
plot_C(res5)
title(main="MYH6 PeakC")

for (peaks in c(23864623,	23867467,	23870458,	23873236,	23901694,	23904226,	23904489,	23907602,	23911007,	23914591,	23925331,	23864798,	23868026,	23871112, 23873576,	23902077,	23904484,	23907056,	23908184,	23913075,	23915127,	23925645)){
  FourCpeakinPeakC <- which(abs(res5$dbR$frag_pos-peaks)==min(abs(res5$dbR$frag_pos-peaks)))
  ActualPeakCpos <- c(ActualPeakCpos, res5$dbR$frag_pos[FourCpeakinPeakC])
  
}
res5$peak <- ActualPeakCpos
plot_C(res5, y.max = 3000)
title(main="MYH6 FourCSeq")
######################################################################################################################

viewpoint6 <- 23901458
files6 <- dir(path="C:/4C/peakC/Gecombineerdefiles/", pattern="MYH7_[123].txt", full=TRUE)
data6 <- readMultiple(files6, vp.pos=viewpoint6)
res6 <- combined.analysis(data6, vp.pos=viewpoint6)
plot_C(res6, y.max = 3000)
title(main="MYH7 PeakC")


for (peaks in c(23876080,	23902082,	23926758,	23877160,	23902260,	23927321)){
  FourCpeakinPeakC <- which(abs(res6$dbR$frag_pos-peaks)==min(abs(res6$dbR$frag_pos-peaks)))
  ActualPeakCpos <- c(ActualPeakCpos, res6$dbR$frag_pos[FourCpeakinPeakC])
  
}
res6$peak <- ActualPeakCpos
plot_C(res6, y.max = 3000)
title(main="MYH7 FourCSeq")

######################################################################################################################
viewpoint7 <- 156052727
files7 <- dir(path="C:/4C/peakC/Gecombineerdefiles/", pattern="LMNA1_[123].txt", full=TRUE)
data7 <- readMultiple(files7, vp.pos=viewpoint7)
res7 <- combined.analysis(data7, vp.pos=viewpoint7)
plot_C(res7)
title(main="LMNA1 PeakC")


for (peaks in c(156034859,	156038825,	156041749,	156055198,	156092798,	156101547,	156103541,	156103875,	156190661,	156193112,	156035708,	156039035,	156042805,	156056384,	156093829,	156101736,	156103870,	156104847,	156191555,	156193644)){
  FourCpeakinPeakC <- which(abs(res7$dbR$frag_pos-peaks)==min(abs(res7$dbR$frag_pos-peaks)))
  ActualPeakCpos <- c(ActualPeakCpos, res7$dbR$frag_pos[FourCpeakinPeakC])
  
}
res7$peak <- ActualPeakCpos
plot_C(res7, y.max = 3000)
title(main="LMNA1 FourCSeq")


######################################################################################################################
viewpoint8 <- 156086255
files8 <- dir(path="C:/4C/peakC/Gecombineerdefiles/", pattern="LMNA2_[123].txt", full=TRUE)
data8 <- readMultiple(files8, vp.pos=viewpoint8)
res8 <- combined.analysis(data8, vp.pos=viewpoint8)
plot_C(res8)
title(main="LMNA2 PeakC")

for (peaks in c(156055198,	156058836,	156061590,	156063165,	156063339,	156063783,	156065763,	156066068,	156066920,	156067167,	156069801,	156070767,	156077936,	156078239,	156079746,	156080547,	156081986,	156089149,	156089985,	156090948,	156091726,	156092798,	156093834,	156098438,	156101547,	156103541,	156166266,	156190661,	156056384,	156059005,	156062118,	156063310,	156063530,	156063927,	156066063,	156066915,	156067162,	156067759,	156070599,	156070986,	156078234,	156078340,	156080379, 156081333,	156082127,	156089843,	156090056,	156091418,	156091846,	156093829,	156094142,	156099384,	156101736,	156103870,	156166939,	156191555)){
  FourCpeakinPeakC <- which(abs(res8$dbR$frag_pos-peaks)==min(abs(res8$dbR$frag_pos-peaks)))
  ActualPeakCpos <- c(ActualPeakCpos, res8$dbR$frag_pos[FourCpeakinPeakC])
  
}
res8$peak <- ActualPeakCpos
plot_C(res8, y.max = 3000)
title(main="LMNA2 FourCSeq")
######################################################################################################################
