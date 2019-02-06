rm(list=ls(all = TRUE))

##### read and prepare dataset for analyses:
RohAll = read.table("../../Body/3Results/GeneRohRecomb.classAll.txt", header = TRUE); RohAll$RohAllCoverage = RohAll$RohCoverage
RohA = read.table("../../Body/3Results/GeneRohRecomb.classA.txt", header = TRUE); RohA = RohA[grepl("EnsemblId|RohCoverage",colnames(RohA))]; names(RohA)=c('EnsemblId','RohACoverage')
RohB = read.table("../../Body/3Results/GeneRohRecomb.classB.txt", header = TRUE); RohB = RohB[grepl("EnsemblId|RohCoverage",colnames(RohB))]; names(RohB)=c('EnsemblId','RohBCoverage')
RohC = read.table("../../Body/3Results/GeneRohRecomb.classC.txt", header = TRUE); RohC = RohC[grepl("EnsemblId|RohCoverage",colnames(RohC))]; names(RohC)=c('EnsemblId','RohCCoverage')

Roh = merge(RohAll,RohA); Roh = merge(Roh,RohB); Roh = merge(Roh,RohC); 
cor.test(Roh$RohAllCoverage,Roh$RohACoverage+Roh$RohBCoverage+Roh$RohCCoverage) # perfect correlation as expected. 
Roh$FractionCoveredByC = Roh$RohCCoverage / Roh$RohAllCoverage
Roh$FractionCoveredByB = Roh$RohBCoverage / Roh$RohAllCoverage
Roh$FractionCoveredByA = Roh$RohACoverage / Roh$RohAllCoverage

write.table(Roh, "../../Body/3Results/GeneRohRecomb.classesABCAlls.txt")
