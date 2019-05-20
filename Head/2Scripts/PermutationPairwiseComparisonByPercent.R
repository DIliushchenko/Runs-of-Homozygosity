rm(list = ls(all = TRUE))  # remove everything from R memory (old variables, datasets...)

RohFull <- read.table("../../Body/2Derived/Roh.txt")

summary(RohFull$KB)

jpeg("../../Body/4Figures/PermutationKids.R.01.jpg")

for (i in 1:10000){
  
modified <- new.env()
for (i in 1:length(RohFull$IID)) {
  if (RohFull[i,]$IID %in% ls(modified)) {
    iidd <- RohFull[i,]$IID
    RohFull[i,]$PHE = modified[[as.character(iidd)]]
  }  else {
    PHE = sample(0:1, 1)
    print(PHE)
    
    RohFull[i,]$PHE = PHE
    
    iidd <- RohFull[i,]$IID
    modified[[as.character(iidd)]] <- PHE
  }
}

p <- c()
n <- c()

for (k in 1:100) {
  print("k")
  print(k)
  # i = 2
  ROH = RohFull[RohFull$KB >= quantile(RohFull$KB, k / 100),]
  
  Agg = aggregate(
    ROH$KB,
    by = list(
      ROH$FidIid,
      ROH$FID,
      ROH$FamilyMember,
      ROH$Family,
      ROH$PHE,
      ROH$KidsInFamily,
      ROH$BothParentsAreHealthy
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
            & Agg$FamilyMember == 'Kid' & Agg$PHE != -9,]
  
  Agg = Agg[order(Agg$FID), ]
  
  ### all pairwise combinations between affected and non-affected siblings:
  
  VecOfFamilies = unique(Agg$FID)
  length(VecOfFamilies) # 79
  for (fam in (1:length(VecOfFamilies)))
    # go one by one through all families
  {
    # fam = 1
    NonAff = Agg[Agg$FID == VecOfFamilies[fam] &
                   Agg$PHE == 1, ]
    NonAff = NonAff[, c(1, 2, 8)]
    names(NonAff) = c('FidIid.NonAffKid', 'FID', 'RohTotal.NonAffKid') # keep only nonaffected from a given family
    Aff = Agg[Agg$FID == VecOfFamilies[fam] &
                Agg$PHE == 2, ]
    Aff = Aff[, c(1, 8)]
    names(Aff) = c('FidIid.AffKid', 'RohTotal.AffKid')              # keep only affected from a given family
    Fam = merge(NonAff, Aff)
    # merge all nonaffected with all nonaffected (all versus all = product)
    Fam = Fam[, c(2, 1, 4, 3, 5)] # change the order of columns
    if (fam == 1) {
      Final = Fam
    }
    if (fam >  1) {
      Final = rbind(Final, Fam)
    }
  }
  write.table(Final)
  
  
  Final$Contrast = Final$RohTotal.NonAffKid - Final$RohTotal.AffKid
  summary(Final$Contrast)
  
  NumberOfAnalyzedFamilies = nrow(Final)
  Title2 = paste('Number Of Analyzed Kids = ', NumberOfAnalyzedFamilies, sep = '')
  
  a = wilcox.test(Final$Contrast, mu = 0)
  Pvalue = as.numeric(a[3])
  ObservedMedian = round(median(Final$Contrast), 0)
  p[k] <- k
  n[k] <- ObservedMedian
}

plot(p, n, xlab = "%", ylab = "ObservedMedian", "l", col = gray)
}

dev.off()
