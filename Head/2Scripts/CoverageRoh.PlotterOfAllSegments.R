rm(list=ls(all=TRUE))

WindowsAll = read.table('../../Body/3Results/CoverageRoh.CoverageROH.WindowsForRoh.classAll.txt')
WindowsA = read.table('../../Body/3Results/CoverageRoh.CoverageROH.WindowsForRoh.classA.txt')
WindowsB = read.table('../../Body/3Results/CoverageRoh.CoverageROH.WindowsForRoh.classB.txt')
WindowsC = read.table('../../Body/3Results/CoverageRoh.CoverageROH.WindowsForRoh.classC.txt')

### PLOT WindowsAll in one page

pdf("../../Body/4Figures/CoverageRoh.PlotterOfAllSegments.R01.pdf", height = 400, width = 600) # Dima - why this -? , onefile = FALSE)
par(mfrow=c(22,1))
for (chr in 1:22)
 {
   #### chr=1
   WindowsAllForCHR = WindowsAll[WindowsAll$CHR==chr,]
   a = max(WindowsAllForCHR$StartofWindow)
   b = min(WindowsAllForCHR$EndofWindow)
   d = min(WindowsAllForCHR$coverage)
   f = max(WindowsAllForCHR$coverage)
   plot(NA, xlim=c(a,b), ylim=c(d,f), xlab='Location of Windows', ylab='Coverage')
   title(paste(chr), line = -30, cex.main = 50, adj = 0.01)
   for (i in 1:nrow(WindowsAllForCHR))
   {
     ### i =1
     segments(WindowsAllForCHR$StartofWindow[i],WindowsAllForCHR$coverage[i],WindowsAllForCHR$EndofWindow[i],WindowsAllForCHR$coverage[i],col = 'dark grey', lwd = 60) # Dima it was "WindowsForCHR$coverage" - wrong!!!
   }
}

### PLOT WindowsAll, A, B, C in one page

par(mfrow=c(22,1))
for (chr in 1:22)
{
  #### chr=1
  WindowsAllForCHR = WindowsAll[WindowsAll$CHR==chr,]
  WindowsAForCHR = WindowsA[WindowsA$CHR==chr,]
  WindowsBForCHR = WindowsB[WindowsB$CHR==chr,]
  WindowsCForCHR = WindowsC[WindowsC$CHR==chr,]
  a = max(WindowsAllForCHR$StartofWindow)
  b = min(WindowsAllForCHR$EndofWindow)
  d = min(WindowsAllForCHR$coverage)
  f = max(WindowsAllForCHR$coverage)  
  plot(NA, xlim=c(a,b), ylim=c(d,f), xlab='Location of Windows', ylab='Coverage')
  title(paste(chr), line = -30, cex.main = 50, adj = 0.01)
  for (i in 1:nrow(WindowsAllForCHR))   {segments(WindowsAllForCHR$StartofWindow[i],WindowsAllForCHR$coverage[i],WindowsAllForCHR$EndofWindow[i],WindowsAllForCHR$coverage[i],col = 'dark grey', lwd = 30)}
  for (i in 1:nrow(WindowsAForCHR))     {segments(WindowsAForCHR$StartofWindow[i],WindowsAForCHR$coverage[i],WindowsAForCHR$EndofWindow[i],WindowsAForCHR$coverage[i],col = 'red', lwd = 30)}
  for (i in 1:nrow(WindowsBForCHR))     {segments(WindowsBForCHR$StartofWindow[i],WindowsBForCHR$coverage[i],WindowsBForCHR$EndofWindow[i],WindowsBForCHR$coverage[i],col = 'green', lwd = 30)}
  for (i in 1:nrow(WindowsCForCHR))     {segments(WindowsCForCHR$StartofWindow[i],WindowsCForCHR$coverage[i],WindowsCForCHR$EndofWindow[i],WindowsCForCHR$coverage[i],col = 'blue', lwd = 30)}
}

### PLOT WindowsAll, A, B, C in separate pages for each chromosome

par(mfrow=c(1,1))
for (chr in 1:22)
{
  #### chr=1
  WindowsAllForCHR = WindowsAll[WindowsAll$CHR==chr,]
  WindowsAForCHR = WindowsA[WindowsA$CHR==chr,]
  WindowsBForCHR = WindowsB[WindowsB$CHR==chr,]
  WindowsCForCHR = WindowsC[WindowsC$CHR==chr,]
  a = max(WindowsAllForCHR$StartofWindow)
  b = min(WindowsAllForCHR$EndofWindow)
  d = min(WindowsAllForCHR$coverage)
  f = max(WindowsAllForCHR$coverage)  
  plot(NA, xlim=c(a,b), ylim=c(d,f), xlab='Location of Windows', ylab='Coverage')
  title(paste(chr), line = -30, cex.main = 50, adj = 0.01)
  for (i in 1:nrow(WindowsAllForCHR))   {segments(WindowsAllForCHR$StartofWindow[i],WindowsAllForCHR$coverage[i],WindowsAllForCHR$EndofWindow[i],WindowsAllForCHR$coverage[i],col = 'dark grey', lwd = 100)}
  for (i in 1:nrow(WindowsAForCHR))     {segments(WindowsAForCHR$StartofWindow[i],WindowsAForCHR$coverage[i],WindowsAForCHR$EndofWindow[i],WindowsAForCHR$coverage[i],col = 'red', lwd = 100)}
  for (i in 1:nrow(WindowsBForCHR))     {segments(WindowsBForCHR$StartofWindow[i],WindowsBForCHR$coverage[i],WindowsBForCHR$EndofWindow[i],WindowsBForCHR$coverage[i],col = 'green', lwd = 100)}
  for (i in 1:nrow(WindowsCForCHR))     {segments(WindowsCForCHR$StartofWindow[i],WindowsCForCHR$coverage[i],WindowsCForCHR$EndofWindow[i],WindowsCForCHR$coverage[i],col = 'blue', lwd = 100)}
}

dev.off()