rm(list=ls(all=TRUE))

wd = getwd()
wd = paste(wd,'/Runs-of-Homozygosity/Body/2Derived/', sep="")
setwd(wd)
ROH = read.table('Roh.txt')
ROH = ROH[ROH$CHR == 1,]
VecOfUniqueBreakPoints = sort(unique(c(ROH$POS1,ROH$POS2)))

VecStart = VecOfUniqueBreakPoints[-length(VecOfUniqueBreakPoints)]
VecEnd = VecOfUniqueBreakPoints[-1]
DataDriwenWind = data.frame(VecStart,VecEnd)
DataDriwenWind$Coverage = 0
DataDriwenWind$CHR = 1

for (i in 1:nrow(DataDriwenWind))
{### i=1
  Start = DataDriwenWind$VecStart[i]
  End = DataDriwenWind$VecEnd[i]
  coverage = nrow(ROH[ROH$POS1 <= Start & ROH$POS2 >= End, ])
  DataDriwenWind$Coverage[i] = coverage

}


### plot of  windows' starts and segments of windows  
plot(DataDriwenWind$VecEnd, DataDriwenWind$Coverage, type = 'l', col = 'gray', xlab = 'Location of windows along CHR1' , ylab = 'Coverage of each window by ROH')
segments(x0 = DataDriwenWind$VecStart, x1 = DataDriwenWind$VecEnd, y0 = DataDriwenWind$Coverage, y1 = DataDriwenWind$Coverage, col = 'red')




