rm(list = ls(all = TRUE))  # remove everything from R memory (old variables, datasets...)

RohFull <- read.table("../../Body/2Derived/Roh.txt")

summary(RohFull$KB)

Agg = aggregate(
    RohFull$KB,
    by = list(
        RohFull$FidIid,
        RohFull$FID,
        RohFull$FamilyMember,
        RohFull$Family,
        RohFull$PHE,
        RohFull$KidsInFamily,
        RohFull$BothParentsAreHealthy
    ),
    FUN = sum
)
names(Agg) = c(
    'FidIid',
    'FID',
    'FamilyMember',
    'Family',
    'PHE',
    'KidsInFamily',
    'BothParentsAreHealthy',
    'RohTotal'
)
Agg = Agg[Agg$Family == 'CompleteFamily' &
              Agg$KidsInFamily == 'KidsAffectedAndNonAffected' &
              Agg$BothParentsAreHealthy == 'Yes'
          & Agg$FamilyMember == 'Kid' & Agg$PHE != -9, ]

Agg = Agg[order(Agg$FID),]

### all pairwise combinations between affected and non-affected siblings:

VecOfFamilies = unique(Agg$FID); length(VecOfFamilies) # 79
for (fam in (1:length(VecOfFamilies))) # go one by one through all families
{ # fam = 1
  NonAff = Agg[Agg$FID == VecOfFamilies[fam] & Agg$PHE == 1,]; NonAff = NonAff[,c(1,2,8)]; names(NonAff) = c('FidIid.NonAffKid','FID','RohTotal.NonAffKid') # keep only nonaffected from a given family
  Aff = Agg[Agg$FID == VecOfFamilies[fam] & Agg$PHE == 2,]; Aff = Aff[,c(1,8)]; names(Aff)=c('FidIid.AffKid','RohTotal.AffKid')              # keep only affected from a given family
  Fam = merge(NonAff,Aff); # merge all nonaffected with all nonaffected (all versus all = product)    
  Fam = Fam[,c(2,1,4,3,5)] # change the order of columns
  if (fam == 1) {Final = Fam}
  if (fam >  1) {Final = rbind(Final,Fam)}
}
write.table(Final, "../../Body/3Results/SiblingsRoh.PairwiseComparisonOfSiblingsKP.txt", quote = FALSE, row.names = FALSE)
