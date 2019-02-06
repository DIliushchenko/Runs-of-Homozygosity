rm(list=ls(all = TRUE))

## read genes:
Genes = read.table("../../Body/1Raw/gencode.v25.annotation.gtf.Genes.Shet.pLI.FIS.RVIS.GHIS.KnKs.GC.BrainSpecificRanking.Branch", header = TRUE)
Genes = Genes[Genes$GeneChr %in% c(1:22),]; nrow(Genes)

## unzip recomb map
unzip("../../Body/1Raw/RecombMap/plink.GRCh37.map.zip", exdir = "../../Body/1Raw/RecombMap/")

## read coverage for each window:

VecOfClasses = c('All','A','B','C')

for (class in VecOfClasses)
{
  if (class == 'All')  {Roh = read.table("../../Body/3Results/CoverageRoh.CoverageROH.WindowsForRoh.classAll.txt"); outfile = "../../Body/3Results/GeneRohRecomb.classAll.txt"}
  if (class == 'A')    {Roh = read.table("../../Body/3Results/CoverageRoh.CoverageROH.WindowsForRoh.classA.txt");   outfile = "../../Body/3Results/GeneRohRecomb.classA.txt"}
  if (class == 'B')    {Roh = read.table("../../Body/3Results/CoverageRoh.CoverageROH.WindowsForRoh.classB.txt");   outfile = "../../Body/3Results/GeneRohRecomb.classB.txt"}
  if (class == 'C')    {Roh = read.table("../../Body/3Results/CoverageRoh.CoverageROH.WindowsForRoh.classC.txt");   outfile = "../../Body/3Results/GeneRohRecomb.classC.txt"}
  
  for (chr in 1:22)
  { # chr = 1
    Genes1 = Genes[Genes$GeneChr == chr,]; 
    Genes1$RohCoverage = 0 
    Genes1$RecombActivity = 0
    
    ## get average recombination for each gene
    infile = paste('../../Body/1Raw/RecombMap/plink.chr',chr,'.GRCh37.map',sep = ''); RecMap = read.table(infile)
    RecMapStart = RecMap[-nrow(RecMap),c(3,4)]; names(RecMapStart)=c('RecStart','WindowStart')
    RecMapEnd = RecMap[-1,c(3,4)]; names(RecMapEnd)=c('RecEnd','WindowEnd')
    RecMap = cbind(RecMapStart,RecMapEnd); RecMap$Rec = RecMap$RecEnd-RecMap$RecStart
    # plot(RecMap$RecStart,RecMap$Rec, pch = '.', main = chr)
    
    for (j in 1:nrow(Genes1))
    { # j = 5
      Start = Genes1$GeneStart[j]; End = Genes1$GeneEnd[j];              #  169849631 169894267
      Rec = RecMap[RecMap$WindowEnd > Start & RecMap$WindowStart < End,] #  169794359 170112399
      if (nrow(Rec) == 1) {Genes1$RecombActivity[j] = Rec$Rec}
      if (nrow(Rec) >  1) 
      {
        for (line in 1:nrow(Rec))
        { # line = 1
          if (Rec$WindowStart[line] < Start) {Rec$WindowStart[line] = Start}
          if (Rec$WindowEnd[line] > End) {Rec$WindowEnd[line] = End}
        }
        Rec$Length = Rec$WindowEnd - Rec$WindowStart
        Genes1$RecombActivity[j] = sum(as.numeric(Rec$Length)*as.numeric(Rec$Rec)) / sum(Rec$Length)
      }
    }
  
    ## get average ROH coverage for each gene
    Roh1 = Roh[Roh$CHR == chr,]
    for (j in 1:nrow(Genes1))
    { # j = 5
      Start = Genes1$GeneStart[j]; End = Genes1$GeneEnd[j];             #  169849631 169894267
      Roh2 = Roh1[Roh1$EndofWindow > Start & Roh1$StartofWindow < End,] #  169794359 170112399
      if (nrow(Roh2) == 1) {Genes1$RohCoverage[j] = Roh2$coverage}
      if (nrow(Roh2) >  1) 
        {
        for (line in 1:nrow(Roh2))
          { # line = 1
          if (Roh2$StartofWindow[line] < Start) {Roh2$StartofWindow[line] = Start}
          if (Roh2$EndofWindow[line] > End) {Roh2$EndofWindow[line] = End}
          }
        Roh2$Length = Roh2$EndofWindow - Roh2$StartofWindow
        Genes1$RohCoverage[j] = sum(as.numeric(Roh2$Length)*as.numeric(Roh2$coverage)) / sum(Roh2$Length)
        }
    }
    
    if (chr == 1) {FINAL = Genes1}
    if (chr >  1) {FINAL = rbind(FINAL,Genes1)}
  } # warnings()
  
 write.table(FINAL, outfile, row.names = FALSE, quote = FALSE)
}

## delete all unziped recomb map files
files <- list.files("../../Body/1Raw/RecombMap/")
files = setdiff(files,c('plink.GRCh37.map.zip'))
for (i in 1:length(files)) 
{ # i = 1
  file = paste('../../Body/1Raw/RecombMap/',files[i],sep='')
  if (file.exists(file)) file.remove(file)
}
