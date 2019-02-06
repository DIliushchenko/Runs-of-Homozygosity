rm(list=ls(all = TRUE))

##read our table from our directory
Roh.All = read.table("../../Body/2Derived/Roh.txt")   

# from Pemberton 2012:
# less than 0.516 Mb - class A (out of Africa);
# between 0.516 Mb and 1.606 Mb - class B (population specific);
# more than 1.606 Mb - class C (pedigree-specific)

VecOfClasses = c('All','A','B','C')
# VecOfClasses = c('All')

for (class in VecOfClasses)
{
if (class == 'All')  {Roh = Roh.All;}
if (class == 'A')  {Roh = Roh.All[Roh.All$KB <= 516,];}
if (class == 'B')  {Roh = Roh.All[Roh.All$KB > 516 & Roh.All$KB <= 1606,];}
if (class == 'C')  {Roh = Roh.All[Roh.All$KB > 1606,];}

### create data for all windows 
for (j in 1:22) #
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
  if (j == 1) {DataDrivedWindows = DataWind}                             
  if (j >  1) {  DataDrivedWindows = rbind(DataDrivedWindows,DataWind)}  
}

### write the table
names(DataDrivedWindows) <- c('StartofWindow', 'EndofWindow', 'coverage', 'CHR')
outfile = paste('../../Body/3Results/CoverageRoh.CoverageROH.WindowsForRoh.class',class,'.txt',sep='')
write.table(DataDrivedWindows,outfile)
}
