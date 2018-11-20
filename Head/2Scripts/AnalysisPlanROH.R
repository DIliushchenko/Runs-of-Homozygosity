rm(list=ls(all=TRUE))

user = 'DIMA'
if (user == 'DIMA') {dir = '/home/diliushchenko/Рабочий стол/ROH_ansar/'}
  
##### read tfam from each folder, add column 'FileName', bind all Tfam together, derive FamilyMember column and FidIid column
  VecOfTfam = c('Intellectual Disability/528.tfam',
                'New folder/VI2_pak_04052016_geno_t.tfam',
                'New folder1/MA_tFams.tfam',
                'New folder2/140317_VI_3rd_plate_geno_t.tfam',
                'New folder3/VI4_pak_06062017_geno_t.tfam')
  for (i in 1:length(VecOfTfam))
  { # i = 2
  file = paste(dir,VecOfTfam[i],sep=''); 
  Tfam = read.table(file); 
  names(Tfam) = c('FID','IID','FatherId','MotherId','Gender','Phenotype')
  Tfam$FID = as.character(Tfam$FID)
  Tfam$IID = as.character(Tfam$IID)
  Tfam$FileName = VecOfTfam[i]
  if (i == 1) {TFAM = Tfam}
  if (i > 1) {TFAM = rbind(TFAM,Tfam)}
  }
  
  TFAM$FamilyMember = 'Unknown'
  for (i in 1:nrow(TFAM))
  { 
    if (TFAM$FatherId[i] == 0 & TFAM$MotherId[i] == 0 & TFAM$Gender[i] == 1) {TFAM$FamilyMember[i] = 'Father'}
    if (TFAM$FatherId[i] == 0 & TFAM$MotherId[i] == 0 & TFAM$Gender[i] == 2) {TFAM$FamilyMember[i] = 'Mother'}
    if (TFAM$FatherId[i] != 0 & TFAM$MotherId[i] != 0) {TFAM$FamilyMember[i] = 'Kid'}
  }
  table(TFAM$FamilyMember)
  # TFAM[TFAM$FamilyMember == 'Unknown',] # check them by eye
  
  TFAM$FidIid = paste(TFAM$FID,TFAM$IID, sep = '|')
  
  TFAM$Family = 'Unknown'
  TFAM$KidsInFamily = 'Unknown'
  TFAM$BothParentsAreHealthy = 'Unknown'
  VecOfFam = unique(TFAM$FID); length(VecOfFam) # 166
  for (i in 1:length(VecOfFam))
  { # i = 1
    TEMP = TFAM[TFAM$FID == VecOfFam[i],]
    if (nrow(TEMP[TEMP$FamilyMember == 'Father',]) == 1 & nrow(TEMP[TEMP$FamilyMember == 'Mother',]) == 1 & nrow(TEMP[TEMP$FamilyMember == 'Kid',]) > 1)
      {
      TFAM[TFAM$FID == VecOfFam[i],]$Family = 'CompleteFamily';
      }
    if (nrow(TEMP[TEMP$FamilyMember == 'Kid' & TEMP$Phenotype == 1,]) > 0 & nrow(TEMP[TEMP$FamilyMember == 'Kid' & TEMP$Phenotype == 2,]) > 0)
     {
      TFAM[TFAM$FID == VecOfFam[i],]$KidsInFamily = 'KidsAffectedAndNonAffected';
     }
    if (nrow(TEMP[TEMP$FamilyMember == 'Father' & TEMP$Phenotype == 1,]) == 1 & nrow(TEMP[TEMP$FamilyMember == 'Mother' & TEMP$Phenotype == 1,]) == 1)
    {
      TFAM[TFAM$FID == VecOfFam[i],]$BothParentsAreHealthy = 'Yes';
    }
  }
  
  TFAM$KidsInFamilies = 'Unknown'
  VecOfFam = unique(TFAM$FID); length(VecOfFam) # 166
  for (i in 1:length(VecOfFam))
  { # i = 1
    TEMP = TFAM[TFAM$FID == VecOfFam[i],]
    if (nrow(TEMP[TEMP$KidFamilyMember == 'Father',]) == 1 & nrow(TEMP[TEMP$FamilyMember == 'Mother',]) == 1 & nrow(TEMP[TEMP$FamilyMember == 'Kid',]) > 1)
    {
      TFAM[TFAM$FID == VecOfFam[i],]$Family = 'CompleteFamily';
    }
  }
  
  
TFAM = TFAM[,grepl("FidIid|FamilyMember|Gender|FatherId|MotherId|Family|KidsInFamily|BothParentsAreHealthy", names(TFAM))]
    

##### read Hom from each folder, add column 'FileName', bind all Hom together, derive FidIid column
VecOfHom = c('Intellectual Disability/528_ROH.hom','New folder/VI_Pak_040516_corrected_RO.hom','New folder1/MA_tFams_new_RO.hom','New folder2/140317_VI_3rd_plate_geno_t_RO.hom','New folder3/VI4_pak_06062017_RO.hom')
for (i in 1:length(VecOfHom))
{ # i = 1
  file = paste(dir,VecOfHom[i],sep=''); 
  Hom = read.table(file, header = TRUE); 
  #names(Tfam) = c('FID','IID','FatherId','MotherId','Gender','Phenotype')
  Hom$FID = as.character(Hom$FID)
  Hom$IID = as.character(Hom$IID)
  Hom$FileName = VecOfHom[i]
  if (i == 1) {HOM = Hom}
  if (i > 1) {HOM = rbind(HOM,Hom)}
}
HOM$FidIid = paste(HOM$FID,HOM$IID, sep = '|')

##### merge HOM and TFAM by FidIid
nrow(HOM)
ROH = merge(HOM,TFAM, by = 'FidIid')
nrow(ROH)
ROH = ROH[order(ROH$Family,ROH$FID),]

##### save file
File = paste(dir,'Roh.txt', sep = '')
write.table(ROH,file = File)





################## SEPARATE SCRIPT - ANALYSES


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
hist(ParentsKids$Contrasts, breaks = 50)
abline(v=0, lwd = 2, col = 'red')
wilcox.test(ParentsKids$Contrasts, mu = 0)


ParentsKidsPHE = merge(AggParents,AggKids)
ParentsKidsPHE$Contrasts = ParentsKidsPHE$MeanKids - ParentsKidsPHE$MeanParents


g <- ggplot(ParentsKidsPHE, aes(MeanParents, Contrasts, color = ))+
  geom_point()
g + geom_vline(xintercept = median(ParentsKidsPHE$MeanParents), color ='red')
  
    
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
hist(Kids$Contrast, breaks= 50)
abline(v=0, lwd = 2, col = 'red')
wilcox.test(Kids$Contrast, mu = 0)

AggKids2 <- aggregate(Agg$RohTotal, by = list(Agg))


library(ggplot2)
ggplot()

