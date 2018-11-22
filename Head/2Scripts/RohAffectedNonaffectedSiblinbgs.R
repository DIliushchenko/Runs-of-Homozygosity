rm(list=ls(all=TRUE))  # remove everything from R memory (old variables, datasets...) 

RohFull <- read.table("../../Body/2Derived/Roh.txt")

summary(RohFull$KB) 

## 51.89   559.34   909.33  3075.98  1945.79 91289.92 

pdf("../../Body/4Figures/RohAffectedNonaffectedSiblinbgs.R.01.pdf")

VecOfThresholds = c('all','FirstQ','SecondQ','ThirdQ','FourthQ', 'longer90%')
for (i in 1:length(VecOfThresholds))
{ # i = 2
if (VecOfThresholds[i] == 'all')     {ROH = RohFull}
if (VecOfThresholds[i] == 'FirstQ')  {ROH = RohFull[RohFull$KB < quantile(RohFull$KB, 0.25),]}
if (VecOfThresholds[i] == 'SecondQ') {ROH = RohFull[RohFull$KB < quantile(RohFull$KB, 0.5) & RohFull$KB >= quantile(RohFull$KB, 0.25),]}
if (VecOfThresholds[i] == 'ThirdQ')  {ROH = RohFull[RohFull$KB < quantile(RohFull$KB, 0.75) & RohFull$KB >= quantile(RohFull$KB, 0.5),]}
if (VecOfThresholds[i] == 'FourthQ') {ROH = RohFull[RohFull$KB >= quantile(RohFull$KB, 0.75),]}
if (VecOfThresholds[i] == 'longer90%') {ROH = RohFull[RohFull$KB >= quantile(RohFull$KB, 0.90),]}
    
##### get sum of ROH per each individuum, filter out non complete familie

Agg = aggregate(ROH$KB, by = list(ROH$FidIid,ROH$FID,ROH$FamilyMember,ROH$Family,ROH$PHE), FUN = sum)
names(Agg) = c('FidIid','FID','FamilyMember','Family','PHE','RohTotal')
Agg = Agg[Agg$Family == 'CompleteFamily',]
Agg = Agg[order(Agg$FID),]

###### compare parents with kids from each family
AggParents = Agg[Agg$FamilyMember == 'Father' | Agg$FamilyMember == 'Mother',]
AggParents = aggregate(AggParents$RohTotal, by = list(AggParents$FID), FUN = mean) 
names(AggParents)=c('FID','MeanParents')

AggKids = Agg[Agg$FamilyMember == 'Kid',]
AggKids = AggKids[AggKids$PHE == '1' | AggKids$PHE == '2',]
AggKids = aggregate(AggKids$RohTotal, by = list(AggKids$FID, AggKids$PHE), FUN = mean) 
names(AggKids)=c('FID','PHE','MeanKids')

AggKids2 = Agg[Agg$FamilyMember == 'Kid',]
AggKids2 = aggregate(AggKids2$RohTotal, by = list(AggKids2$FID), FUN = mean)
names(AggKids2) = c('FID', 'MeanKids')

ParentsKids = merge(AggParents,AggKids2)
ParentsKids$Contrasts = ParentsKids$MeanKids - ParentsKids$MeanParents
summary(ParentsKids$Contrasts)
NumberOfAnalyzedFamilies = nrow(ParentsKids)
Title2 = paste('Number Of Analyzed Kids = ', NumberOfAnalyzedFamilies, sep = '')

if (VecOfThresholds[i] == 'all')  {MinX = min(ParentsKids$Contrasts); MaxX = max(ParentsKids$Contrasts)}   

a = wilcox.test(ParentsKids$Contrasts, mu = 0)
Pvalue = as.numeric(a[3])
ObservedMedian = round(median(ParentsKids$Contrasts),0)
Title = paste(VecOfThresholds[i],':ObservedMedian =',ObservedMedian, 
',Pvalue =',Pvalue, sep = ' ')
Title = paste(Title,Title2, sep = '\n')


hist(ParentsKids$Contrasts, breaks = 50, main = Title, 
    xlab = 'Contrasts(MeanKidsMinusMeanParents, Kb)', xlim = c(MinX,MaxX), col = 'grey')
abline(v=0, lwd = 2, col = 'red')

}

dev.off()
