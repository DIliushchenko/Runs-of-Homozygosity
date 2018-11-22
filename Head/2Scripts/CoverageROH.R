rm(list=ls(all = TRUE))

##read our table from our directory
wd = getwd()
wd = paste(wd,'/Runs-of-Homozygosity/Body/2Derived/', sep="")
setwd(wd)
ROH = read.table('Roh.txt')

### create data for all windows 
DataDrivedWindows = data.frame()

for (j in 1:25) 
{ 
  ROH = ROH[ROH$CHR == j,]
  VecOfUniqueBreakPoints = sort(unique(c(ROH$POS1,ROH$POS2)))
  VecOfStart = VecOfUniqueBreakPoints[(-length(VecOfUniqueBreakPoints))]
  VecEnd = VecOfUniqueBreakPoints[-1]
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
  DataDrivedWindows = rbind(DataDrivedWindows,DataWind)
}
### write the table in our dir
names(DataDrivedWindows) <- c('StartofWindow', 'EndofWindow', 'coverage', 'CHR')
File = paste(wd,'CoverageROH.txt', sep = '')
write.table(CoverageCHR1, file = File)




