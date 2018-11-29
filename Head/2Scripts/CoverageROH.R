rm(list=ls(all = TRUE))

##read our table from our directory
Roh = read.table("../../Body/2Derived/Roh.txt")   ### KP!!!

### create data for all windows 

for (j in 1:22) # KP!!! 
{ # j = 1 
  ROH = Roh[Roh$CHR == j,]
  VecOfUniqueBreakPoints = sort(unique(c(ROH$POS1,ROH$POS2)))
  VecOfStart = VecOfUniqueBreakPoints[(-length(VecOfUniqueBreakPoints))] # delete last
  VecEnd = VecOfUniqueBreakPoints[-1]  # delete first
  DataWind = data.frame(VecOfStart,VecEnd)
  DataWind$Coverage = 0
  DataWind$CHR = j
  for (i in 1 : nrow(DataWind))
  {
    Start = DataWind$VecOfStart[i]
    End = DataWind$VecEnd[i]
    coverage = nrow(ROH[ROH$POS1 <= Start & ROH$POS2 >= End,]) 
    DataWind$Coverage[i] = coverage
  }
  if (j == 1) {DataDrivedWindows = DataWind}                             # KP!!!!
  if (j >  1) {  DataDrivedWindows = rbind(DataDrivedWindows,DataWind)}  # KP!!!!
}

### write the table
names(DataDrivedWindows) <- c('StartofWindow', 'EndofWindow', 'coverage', 'CHR')
write.table(DataDrivedWindows,"../../Body/3Results/CoverageROHForWindows.txt")
