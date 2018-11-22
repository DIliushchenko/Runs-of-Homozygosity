rm(list=ls(all=TRUE))  # remove everything from R memory (old variables, datasets...) 

test <- read.table("../../Body/1Raw/TestWind.txt")

## 1 to get sorted vector of unique breakpoints

VecOfUniqueBreakPoints = sort(unique(c(test$StartOfWindow,test$EndOfWindow))); length(VecOfUniqueBreakPoints)

## 2 form data-driven windows
VecStart = VecOfUniqueBreakPoints[-length(VecOfUniqueBreakPoints)]
VecEnd = VecOfUniqueBreakPoints[-1]
DataDrivenWindows = data.frame(VecStart,VecEnd)

## 3 estimate coverage for each window

DataDrivenWindows$Coverage = 0

for (i in 1:nrow(DataDrivenWindows))
 {
  Start = DataDrivenWindows$VecStart[i]
  End = DataDrivenWindows$VecEnd[i]
  Cov = nrow(test[test$StartOfWindow <= Start & test$EndOfWindow >= End,])
  DataDrivenWindows$Coverage[i] = Cov
}

### 4 plot segment
