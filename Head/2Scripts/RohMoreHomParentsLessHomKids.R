
rm(list=ls(all=TRUE))  # remove everything from R memory (old variables, datasets...) 

ROH <- read.table("../../Body/2Derived/Roh.txt")

pdf("../../Body/4Figures/RohMoreHomParentsLessHomKids.R.01.pdf")

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
par(mfrow=c(3,1))
hist(ParentsKids$MeanParents, breaks = 30, col = 'grey', main = 'Mean Homozygosity of Parents, kb', xlim = c(0,600000))
hist(ParentsKids$MeanKids, breaks = 30, col = 'grey', main = 'Mean Homozygosity of Kids, kb', xlim = c(0,600000))
hist(ParentsKids$Contrasts, breaks = 30, col = 'grey', main = 'Contrasts (MeanOfKidsMinusMeanOfParents), kb')
abline(v=0, col = 'red', lwd = 2)
par(mfrow=c(1,1))
plot(ParentsKids$MeanParents,ParentsKids$Contrasts, xlab = 'Mean Homozygosity of Parents, kb', ylab = 'Contrasts (MeanOfKidsMinusMeanOfParents), kb', xlim = c(0,400000), main = 'all ROHs')
abline(h=0, col = 'red', lwd = 2)
cor.test(ParentsKids$MeanParents,ParentsKids$Contrasts, method = 'spearman')

######## repeat the same with very long ROHs
summary(ROH$KB)
ROH = ROH[ROH$KB >= quantile(ROH$KB,0.75),]

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

plot(ParentsKids$MeanParents,ParentsKids$Contrasts, xlab = 'Mean Homozygosity of Parents, kb', ylab = 'Contrasts (MeanOfKidsMinusMeanOfParents), kb', xlim = c(0,400000), main = 'ROHs longer than 70%')
abline(h=0, col = 'red', lwd = 2)
cor.test(ParentsKids$MeanParents,ParentsKids$Contrasts, method = 'spearman')

######## repeat the same with very long ROHs
summary(ROH$KB)
ROH = ROH[ROH$KB >= quantile(ROH$KB,0.9),]

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

plot(ParentsKids$MeanParents,ParentsKids$Contrasts, xlab = 'Mean Homozygosity of Parents, kb', ylab = 'Contrasts (MeanOfKidsMinusMeanOfParents), kb', xlim = c(0,400000), main = 'ROHs longer than 90%')
abline(h=0, col = 'red', lwd = 2)
cor.test(ParentsKids$MeanParents,ParentsKids$Contrasts, method = 'spearman')

dev.off()
