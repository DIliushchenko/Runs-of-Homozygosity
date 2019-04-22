rm(list = ls(all = TRUE))  # remove everything from R memory (old variables, datasets...)

RohFull <- read.table("../../Body/2Derived/Roh.txt")

summary(RohFull$KB)

p <- c()
n <- c()

pdf("../../Body/4Figures/PairwiseComparisonByRercentKidsNonAffandKidsAff.R.01.pdf")

for (i in 1:100)
{
  # i = 2
  ROH = RohFull[RohFull$KB >= quantile(RohFull$KB, i / 100), ]
  
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
            & Agg$FamilyMember == 'Kid' & Agg$PHE != -9, ]
  
  Agg = Agg[order(Agg$FID),]
  
  ### all pairwise combinations between affected and non-affected siblings:
  
  VecOfFamilies = unique(Agg$FID)
  length(VecOfFamilies) # 79
  for (fam in (1:length(VecOfFamilies)))
    # go one by one through all families
  {
    # fam = 1
    NonAff = Agg[Agg$FID == VecOfFamilies[fam] &
                   Agg$PHE == 1,]
    NonAff = NonAff[, c(1, 2, 8)]
    names(NonAff) = c('FidIid.NonAffKid', 'FID', 'RohTotal.NonAffKid') # keep only nonaffected from a given family
    Aff = Agg[Agg$FID == VecOfFamilies[fam] &
                Agg$PHE == 2,]
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
  p[i] <- i
  n[i] <- ObservedMedian
  
  
  Title = paste(">",
                i,
                "%",
                ': ObservedMedian = ',
                ObservedMedian,
                ', Pvalue = ',
                Pvalue,
                sep = ' ')
  Title = paste(Title, Title2, sep = '\n')
  
  hist(
    Final$Contrast,
    breaks = 50,
    main = Title,
    xlab = 'Contrasts(MeanAffectedKidsMinusMeanHealthyKids, Kb)',
    col = 'grey'
  )
  abline(v = 0, lwd = 2, col = 'red')
}

dev.off()

plot(p, n, xlab = "%", ylab = "ObservedMedian", "l")
