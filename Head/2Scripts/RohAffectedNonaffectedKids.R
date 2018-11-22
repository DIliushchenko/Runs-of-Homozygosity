rm(list=ls(all=TRUE))  # remove everything from R memory (old variables, datasets...) 

RohFull <- read.table("../../Body/2Derived/Roh.txt")

summary(RohFull$KB) 

pdf("../../Body/4Figures/RohAffectedNonaffectedKids.R.01.pdf")

VecOfThresholds = c('all','FirstQ','SecondQ','ThirdQ','FourthQ', 'longerthan90%')
for (i in 1:length(VecOfThresholds))
{ # i = 2
  if (VecOfThresholds[i] == 'all')     {ROH = RohFull}
  if (VecOfThresholds[i] == 'FirstQ')  {ROH = RohFull[RohFull$KB < quantile(RohFull$KB, 0.25),]}
  if (VecOfThresholds[i] == 'SecondQ') {ROH = RohFull[RohFull$KB < quantile(RohFull$KB, 0.5) & RohFull$KB >= quantile(RohFull$KB, 0.25),]}
  if (VecOfThresholds[i] == 'ThirdQ')  {ROH = RohFull[RohFull$KB < quantile(RohFull$KB, 0.75) & RohFull$KB >= quantile(RohFull$KB, 0.5),]}
  if (VecOfThresholds[i] == 'FourthQ') {ROH = RohFull[RohFull$KB >= quantile(RohFull$KB, 0.75),]}
  quantile(RohFull$KB, prob = seq(0, 1, length = 11), type = 5)
  if (VecOfThresholds[i] == 'longerthan90%') {ROH = RohFull[RohFull$KB >= quantile(RohFull$KB, 0.90),]}
  
  ##### get sum of ROH per each individuum, filter out non complete familie
  
  Agg = aggregate(ROH$KB, by = list(ROH$FidIid,ROH$FID,ROH$FamilyMember,ROH$Family,ROH$PHE), FUN = sum)
  names(Agg) = c('FidIid','FID','FamilyMember','Family','PHE','RohTotal')
  Agg = Agg[Agg$Family == 'CompleteFamily',]
  Agg = Agg[order(Agg$FID),]
  
  Agg = aggregate(ROH$KB, by = list(ROH$FidIid,ROH$FID,ROH$FamilyMember,ROH$Family,ROH$PHE,ROH$KidsInFamily,ROH$BothParentsAreHealthy), FUN = sum)
  names(Agg) = c('FidIid','FID','FamilyMember','Family','PHE','KidsInFamily','BothParentsAreHealthy','RohTotal')
  Agg = Agg[Agg$Family == 'CompleteFamily' & Agg$KidsInFamily == 'KidsAffectedAndNonAffected' & Agg$BothParentsAreHealthy == 'Yes',]
  
  AggHealthyKids = Agg[Agg$FamilyMember == 'Kid' & Agg$PHE == 1,]
  AggAffectedKids = Agg[Agg$FamilyMember == 'Kid' & Agg$PHE == 2,]
  
  AggHealthyKids = aggregate(AggHealthyKids$RohTotal, by = list(AggHealthyKids$FID), FUN = mean); names(AggHealthyKids)=c('FID','MeanHealthyKids')
  AggAffectedKids = aggregate(AggAffectedKids$RohTotal, by = list(AggAffectedKids$FID), FUN = mean); names(AggAffectedKids)=c('FID','MeanAffectedKids')
  
  Kids = merge(AggHealthyKids,AggAffectedKids)
  Kids$Contrast = Kids$MeanAffectedKids - Kids$MeanHealthyKids
  summary(Kids$Contrast)
  
  NumberOfAnalyzedFamilies = nrow(Kids)
  Title2 = paste('Number Of Analyzed Kids = ', NumberOfAnalyzedFamilies, sep = '')
  
  if (VecOfThresholds[i] == 'all')  {MinX = min(Kids$Contrast); MaxX = max(Kids$Contrast)}   
  
  a = wilcox.test(Kids$Contrast, mu = 0)
  Pvalue = as.numeric(a[3])
  ObservedMedian = round(median(Kids$Contrast),0)
  Title = paste(VecOfThresholds[i],': ObservedMedian = ',ObservedMedian, ', Pvalue = ',Pvalue, sep = ' ')
  Title = paste(Title,Title2, sep = '\n')
  
  hist(Kids$Contrast, breaks = 50, main = Title, xlab = 'Contrasts(MeanAffectedKidsMinusMeanHealthyKids, Kb)', 
       xlim = c(MinX,MaxX), col = 'grey')
  abline(v=0, lwd = 2, col = 'red')
}

dev.off()