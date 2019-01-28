rm(list=ls(all = TRUE))

## read coverage for each window:
Roh = read.table("../../Body/3Results/CoverageROHForWindows.txt")

## read genes:
Genes = read.table("../../Body/1Raw/gencode.v25.annotation.gtf.Genes.Shet.pLI.FIS.RVIS.GHIS.KnKs.GC.BrainSpecificRanking.Branch", header = TRUE)

## unzip recomb map
unzip("../../Body/1Raw/RecombMap/plink.GRCh37.map.zip", exdir = "../../Body/1Raw/RecombMap/")

for (chr in 1:22)
{ # chr = 1
  infile = paste('../../Body/1Raw/RecombMap/plink.chr',chr,'.GRCh37.map',sep = ''); RecMap = read.table(infile)
  RecMapStart = RecMap[-nrow(RecMap),c(3,4)]; names(RecMapStart)=c('RecStart','WindowStart')
  RecMapEnd = RecMap[-1,c(3,4)]; names(RecMapEnd)=c('RecEnd','WindowEnd')
  RecMap = cbind(RecMapStart,RecMapEnd); RecMap$Rec = RecMap$RecEnd-RecMap$RecStart
  plot(RecMap$RecStart,RecMap$Rec, pch = '.', main = chr)
  
  Genes1 = Genes[Genes$GeneChr == chr,]; Genes1$Coverage = 0
  Roh1 = Roh[Roh$CHR == chr,]
  for (j in 1:nrow(Genes1))
  { # j = 5
    Start = Genes1$GeneStart[j]; End = Genes1$GeneEnd[j];             #  169849631 169894267
    Roh2 = Roh1[Roh1$EndofWindow > Start & Roh1$StartofWindow < End,] #  169794359 170112399
    if (nrow(Roh2) == 1) {Genes1$Coverage[j] = Roh2$coverage}
    if (nrow(Roh2) >  1) 
      {
      for (line in 1:nrow(Roh2))
        { # line = 1
        if (Roh2$StartofWindow[line] < Start) {Roh2$StartofWindow[line] = Start}
        if (Roh2$EndofWindow[line] < End) {Roh2$EndofWindow[line] = End}
        }
      Roh2$Length = Roh2$EndofWindow - Roh2$StartofWindow
      Genes1$Coverage[j] = sum(as.numeric(Roh2$Length)*as.numeric(Roh2$coverage)) / sum(Roh2$Length)
      }
    }
  if (chr == 1) {FINAL = Genes1}
  if (chr >  1) {FINAL = rbind(FINAL,Genes1)}
} # warnings()