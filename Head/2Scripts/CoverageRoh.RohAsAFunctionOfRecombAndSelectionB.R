rm(list=ls(all = TRUE))

Roh = read.table("../../Body/3Results/GeneRohRecomb.classesABCAlls.txt")

pdf("../../Body/4Figures/CoverageRoh.RohAsAFunctionOfRecombAndSelectionB.R.01.pdf")

##### 1: correlations with each other:

Roh = Roh[Roh$RecombActivity > 0,]
cor.test(Roh$RohACoverage,Roh$RohBCoverage) # 0.12
cor.test(Roh$RohACoverage,Roh$RohCCoverage) # 0.07
cor.test(Roh$RohBCoverage,Roh$RohCCoverage) # 0.24

par(mfrow=c(2,3))
plot(log2(Roh$RohACoverage),log2(Roh$RohBCoverage))
plot(log2(Roh$RohACoverage),log2(Roh$RohCCoverage))
plot(log2(Roh$RohBCoverage),log2(Roh$RohCCoverage))

nrow(Roh[Roh$RohAllCoverage>0,]) # 52763
nrow(Roh[Roh$RohACoverage>0,])   # 17827
nrow(Roh[Roh$RohBCoverage>0,])   # 49524
nrow(Roh[Roh$RohCCoverage>0,])   # 52763

##### 2: recombination: => ALWAYS TAKE INTO MODEL RECOMBINATION!!!!!
cor.test(Roh[Roh$RohAllCoverage > 0,]$RohAllCoverage,Roh[Roh$RohAllCoverage > 0,]$RecombActivity, method = 'spearman') # -0.25
cor.test(Roh[Roh$RohACoverage > 0,]$RohACoverage,Roh[Roh$RohACoverage > 0,]$RecombActivity, method = 'spearman') # -0.16
cor.test(Roh[Roh$RohBCoverage > 0,]$RohBCoverage,Roh[Roh$RohBCoverage > 0,]$RecombActivity, method = 'spearman') # -0.257
cor.test(Roh[Roh$RohCCoverage > 0,]$RohCCoverage,Roh[Roh$RohCCoverage > 0,]$RecombActivity, method = 'spearman') # -0.134

plot(log2(Roh[Roh$RohACoverage > 0,]$RecombActivity),log2(Roh[Roh$RohACoverage > 0,]$RohACoverage))
plot(log2(Roh[Roh$RohBCoverage > 0,]$RecombActivity),log2(Roh[Roh$RohBCoverage > 0,]$RohACoverage))
plot(log2(Roh[Roh$RohCCoverage > 0,]$RecombActivity),log2(Roh[Roh$RohCCoverage > 0,]$RohACoverage))

######## 3: multiple models with recombination and selection:

##### -recomb + KnKsMouse (only C positively correlate with Kn/Ks)
a<-lm(Roh[Roh$RohAllCoverage > 0,]$RohAllCoverage ~ scale(Roh[Roh$RohAllCoverage > 0,]$RecombActivity) + scale(Roh[Roh$RohAllCoverage > 0,]$KnKsMouse)); summary(a)  # -1.9 + 1.4
a<-lm(Roh[Roh$RohACoverage > 0,]$RohACoverage ~ scale(Roh[Roh$RohACoverage > 0,]$RecombActivity) + scale(Roh[Roh$RohACoverage > 0,]$KnKsMouse)); summary(a)          # -0.4 + 0
a<-lm(Roh[Roh$RohBCoverage > 0,]$RohBCoverage ~ scale(Roh[Roh$RohBCoverage > 0,]$RecombActivity) + scale(Roh[Roh$RohBCoverage > 0,]$KnKsMouse)); summary(a)          # -0.7 + 0
a<-lm(Roh[Roh$RohCCoverage > 0,]$RohCCoverage ~ scale(Roh[Roh$RohCCoverage > 0,]$RecombActivity) + scale(Roh[Roh$RohCCoverage > 0,]$KnKsMouse)); summary(a)          # -1.04 + 1.32

a<-lm(Roh[Roh$RohAllCoverage > 0,]$RohAllCoverage ~ scale(Roh[Roh$RohAllCoverage > 0,]$RecombActivity)*scale(Roh[Roh$RohAllCoverage > 0,]$KnKsMouse)); summary(a)    # -2.28 + 0.7 + 006

# C cover more important genes (with low Knks) as compared to B and A:
a<-lm(Roh$KnKsMouse ~  scale(Roh$RohACoverage) + scale(Roh$RohBCoverage) + scale(Roh$RohCCoverage)); summary(a) # KnKsMouse = 0.31 - 0.01*RohBCoverage - 0.03*RohBCoverage - 0.08*RohCCoverage
a<-lm(Roh$KnKsMouse ~  scale(Roh$RohBCoverage) + scale(Roh$RohCCoverage)); summary(a) # KnKsMouse = 0.31 - 0.03*RohBCoverage - 0.08*RohCCoverage
a<-lm(Roh$KnKsMouse ~  scale(Roh$RecombActivity) + scale(Roh$RohACoverage) +  scale(Roh$RohBCoverage) + scale(Roh$RohCCoverage)); summary(a) # KnKsMouse = 0.31 - 0.03*RohBCoverage - 0.08*RohCCoverage
a<-lm(Roh$KnKsMouse ~  scale(Roh$RecombActivity) + scale(Roh$RohAllCoverage)); summary(a)

##### -recomb -ProbOfBeingLofIntolerant
a<-lm(Roh[Roh$RohAllCoverage > 0,]$RohAllCoverage ~ scale(Roh[Roh$RohAllCoverage > 0,]$RecombActivity) + scale(Roh[Roh$RohAllCoverage > 0,]$ProbOfBeingLofIntolerant)); summary(a)  # -2.27 - 0.6
a<-lm(Roh[Roh$RohACoverage > 0,]$RohACoverage ~ scale(Roh[Roh$RohACoverage > 0,]$RecombActivity) + scale(Roh[Roh$RohACoverage > 0,]$ProbOfBeingLofIntolerant)); summary(a)          # -0.41 - 0.63
a<-lm(Roh[Roh$RohBCoverage > 0,]$RohBCoverage ~ scale(Roh[Roh$RohBCoverage > 0,]$RecombActivity) + scale(Roh[Roh$RohBCoverage > 0,]$ProbOfBeingLofIntolerant)); summary(a)          # -0.73 - 0
a<-lm(Roh[Roh$RohCCoverage > 0,]$RohCCoverage ~ scale(Roh[Roh$RohCCoverage > 0,]$RecombActivity) + scale(Roh[Roh$RohCCoverage > 0,]$ProbOfBeingLofIntolerant)); summary(a)          # -1.23 - 0.3

a<-lm(Roh[Roh$RohAllCoverage > 0,]$RohAllCoverage ~ scale(Roh[Roh$RohAllCoverage > 0,]$RecombActivity)*scale(Roh[Roh$RohAllCoverage > 0,]$ProbOfBeingLofIntolerant)); summary(a)  # NO

##### -recomb -SelCoeffHet
a<-lm(Roh[Roh$RohAllCoverage > 0,]$RohAllCoverage ~ scale(Roh[Roh$RohAllCoverage > 0,]$RecombActivity) + scale(Roh[Roh$RohAllCoverage > 0,]$SelCoeffHet)); summary(a)  # -3.03 - 0.44
a<-lm(Roh[Roh$RohACoverage > 0,]$RohACoverage ~ scale(Roh[Roh$RohACoverage > 0,]$RecombActivity) + scale(Roh[Roh$RohACoverage > 0,]$SelCoeffHet)); summary(a)          # -0.48 - 0.38
a<-lm(Roh[Roh$RohBCoverage > 0,]$RohBCoverage ~ scale(Roh[Roh$RohBCoverage > 0,]$RecombActivity) + scale(Roh[Roh$RohBCoverage > 0,]$SelCoeffHet)); summary(a)          # -0.95 - 0
a<-lm(Roh[Roh$RohCCoverage > 0,]$RohCCoverage ~ scale(Roh[Roh$RohCCoverage > 0,]$RecombActivity) + scale(Roh[Roh$RohCCoverage > 0,]$SelCoeffHet)); summary(a)          # -1.69 - 0.19

##### -recomb -GenomeWideHaploinsufficiencyScore
a<-lm(Roh[Roh$RohAllCoverage > 0,]$RohAllCoverage ~ scale(Roh[Roh$RohAllCoverage > 0,]$RecombActivity) + scale(Roh[Roh$RohAllCoverage > 0,]$GenomeWideHaploinsufficiencyScore)); summary(a)  # -1.14 - 0
a<-lm(Roh[Roh$RohACoverage > 0,]$RohACoverage ~ scale(Roh[Roh$RohACoverage > 0,]$RecombActivity) + scale(Roh[Roh$RohACoverage > 0,]$GenomeWideHaploinsufficiencyScore)); summary(a)          # -0.35 - 0.74
a<-lm(Roh[Roh$RohBCoverage > 0,]$RohBCoverage ~ scale(Roh[Roh$RohBCoverage > 0,]$RecombActivity) + scale(Roh[Roh$RohBCoverage > 0,]$GenomeWideHaploinsufficiencyScore)); summary(a)          # -0.67 + 0
a<-lm(Roh[Roh$RohCCoverage > 0,]$RohCCoverage ~ scale(Roh[Roh$RohCCoverage > 0,]$RecombActivity) + scale(Roh[Roh$RohCCoverage > 0,]$GenomeWideHaploinsufficiencyScore)); summary(a)          # -0.27 - 0

# interactions
a<-lm(Roh[Roh$RohAllCoverage > 0,]$RohAllCoverage ~ scale(Roh[Roh$RohAllCoverage > 0,]$RecombActivity)*scale(Roh[Roh$RohAllCoverage > 0,]$GenomeWideHaploinsufficiencyScore)); summary(a)  # YES negative interaction
a<-lm(Roh[Roh$RohACoverage > 0,]$RohACoverage ~ scale(Roh[Roh$RohACoverage > 0,]$RecombActivity)*scale(Roh[Roh$RohACoverage > 0,]$GenomeWideHaploinsufficiencyScore)); summary(a)          # NO
a<-lm(Roh[Roh$RohBCoverage > 0,]$RohBCoverage ~ scale(Roh[Roh$RohBCoverage > 0,]$RecombActivity)*scale(Roh[Roh$RohBCoverage > 0,]$GenomeWideHaploinsufficiencyScore)); summary(a)          # NO
a<-lm(Roh[Roh$RohCCoverage > 0,]$RohCCoverage ~ scale(Roh[Roh$RohCCoverage > 0,]$RecombActivity)*scale(Roh[Roh$RohCCoverage > 0,]$GenomeWideHaploinsufficiencyScore)); summary(a)          # YES 

##### -recomb -FunctionalIndispensabilityScore
a<-lm(Roh[Roh$RohAllCoverage > 0,]$RohAllCoverage ~ scale(Roh[Roh$RohAllCoverage > 0,]$RecombActivity) + scale(Roh[Roh$RohAllCoverage > 0,]$FunctionalIndispensabilityScore)); summary(a)  # -2.33 + 0
a<-lm(Roh[Roh$RohACoverage > 0,]$RohACoverage ~ scale(Roh[Roh$RohACoverage > 0,]$RecombActivity) + scale(Roh[Roh$RohACoverage > 0,]$FunctionalIndispensabilityScore)); summary(a)          # -0.46 + 0
a<-lm(Roh[Roh$RohBCoverage > 0,]$RohBCoverage ~ scale(Roh[Roh$RohBCoverage > 0,]$RecombActivity) + scale(Roh[Roh$RohBCoverage > 0,]$FunctionalIndispensabilityScore)); summary(a)          # -0.79 + 0
a<-lm(Roh[Roh$RohCCoverage > 0,]$RohCCoverage ~ scale(Roh[Roh$RohCCoverage > 0,]$RecombActivity) + scale(Roh[Roh$RohCCoverage > 0,]$FunctionalIndispensabilityScore)); summary(a)          # -1.19 - 0

#### GcContent
a<-lm(Roh[Roh$RohAllCoverage > 0,]$RohAllCoverage ~ scale(Roh[Roh$RohAllCoverage > 0,]$RecombActivity) + scale(Roh[Roh$RohAllCoverage > 0,]$GcContent)); summary(a)  # -1.23 - 0.4
a<-lm(Roh[Roh$RohACoverage > 0,]$RohACoverage ~ scale(Roh[Roh$RohACoverage > 0,]$RecombActivity) + scale(Roh[Roh$RohACoverage > 0,]$GcContent)); summary(a)          # -0.41 + 0.35
a<-lm(Roh[Roh$RohBCoverage > 0,]$RohBCoverage ~ scale(Roh[Roh$RohBCoverage > 0,]$RecombActivity) + scale(Roh[Roh$RohBCoverage > 0,]$GcContent)); summary(a)          # -0.57 + 0.2
a<-lm(Roh[Roh$RohCCoverage > 0,]$RohCCoverage ~ scale(Roh[Roh$RohCCoverage > 0,]$RecombActivity) + scale(Roh[Roh$RohCCoverage > 0,]$GcContent)); summary(a)          # -0.46 - 0.64 (!!!!!!!!! WHY????)

# interaction
a<-lm(Roh[Roh$RohCCoverage > 0,]$RohCCoverage ~ scale(Roh[Roh$RohCCoverage > 0,]$RecombActivity)*scale(Roh[Roh$RohCCoverage > 0,]$GcContent)); summary(a)          # -0.46 - 0.64 (!!!!!!!!! WHY????)
#  scale(Roh[Roh$RohCCoverage > 0, ]$RecombActivity)                                              -0.59039    0.07140  -8.269  < 2e-16 ***
#  scale(Roh[Roh$RohCCoverage > 0, ]$GcContent)                                                   -0.64460    0.06409 -10.058  < 2e-16 ***
#  scale(Roh[Roh$RohCCoverage > 0, ]$RecombActivity):scale(Roh[Roh$RohCCoverage > 0, ]$GcContent)  0.25049    0.06094   4.110 3.96e-05 ***

#### backward stepwise multiple linear models starting with everyghing significant before:
a<-lm(Roh[Roh$RohAllCoverage > 0,]$RohAllCoverage ~ scale(Roh[Roh$RohAllCoverage > 0,]$RecombActivity) + scale(Roh[Roh$RohAllCoverage > 0,]$GcContent) + scale(Roh[Roh$RohAllCoverage > 0,]$ProbOfBeingLofIntolerant) + scale(Roh[Roh$RohAllCoverage > 0,]$SelCoeffHet) + scale(Roh[Roh$RohAllCoverage > 0,]$KnKsMouse)); summary(a)
a<-lm(Roh[Roh$RohAllCoverage > 0,]$RohAllCoverage ~ scale(Roh[Roh$RohAllCoverage > 0,]$RecombActivity) + scale(Roh[Roh$RohAllCoverage > 0,]$ProbOfBeingLofIntolerant) + scale(Roh[Roh$RohAllCoverage > 0,]$SelCoeffHet) + scale(Roh[Roh$RohAllCoverage > 0,]$KnKsMouse)); summary(a)
a<-lm(Roh[Roh$RohAllCoverage > 0,]$RohAllCoverage ~ scale(Roh[Roh$RohAllCoverage > 0,]$RecombActivity) + scale(Roh[Roh$RohAllCoverage > 0,]$ProbOfBeingLofIntolerant) + scale(Roh[Roh$RohAllCoverage > 0,]$KnKsMouse)); summary(a)
a<-lm(Roh[Roh$RohAllCoverage > 0,]$RohAllCoverage ~ scale(Roh[Roh$RohAllCoverage > 0,]$RecombActivity) + scale(Roh[Roh$RohAllCoverage > 0,]$KnKsMouse)); summary(a)  # -1.93 + 1.4

a<-lm(Roh[Roh$RohACoverage > 0,]$RohACoverage ~ scale(Roh[Roh$RohACoverage > 0,]$RecombActivity) + scale(Roh[Roh$RohACoverage > 0,]$GcContent) + scale(Roh[Roh$RohACoverage > 0,]$ProbOfBeingLofIntolerant) + scale(Roh[Roh$RohACoverage > 0,]$SelCoeffHet) + scale(Roh[Roh$RohACoverage > 0,]$KnKsMouse)); summary(a)
a<-lm(Roh[Roh$RohACoverage > 0,]$RohACoverage ~ scale(Roh[Roh$RohACoverage > 0,]$RecombActivity) + scale(Roh[Roh$RohACoverage > 0,]$GcContent) + scale(Roh[Roh$RohACoverage > 0,]$ProbOfBeingLofIntolerant) + scale(Roh[Roh$RohACoverage > 0,]$KnKsMouse)); summary(a)
a<-lm(Roh[Roh$RohACoverage > 0,]$RohACoverage ~ scale(Roh[Roh$RohACoverage > 0,]$RecombActivity) + scale(Roh[Roh$RohACoverage > 0,]$GcContent) + scale(Roh[Roh$RohACoverage > 0,]$ProbOfBeingLofIntolerant)); summary(a) # -0.4 + 0.4 - 0.6  

a<-lm(Roh[Roh$RohBCoverage > 0,]$RohBCoverage ~ scale(Roh[Roh$RohBCoverage > 0,]$RecombActivity) + scale(Roh[Roh$RohBCoverage > 0,]$GcContent) + scale(Roh[Roh$RohBCoverage > 0,]$ProbOfBeingLofIntolerant) + scale(Roh[Roh$RohBCoverage > 0,]$SelCoeffHet) + scale(Roh[Roh$RohBCoverage > 0,]$KnKsMouse)); summary(a)
a<-lm(Roh[Roh$RohBCoverage > 0,]$RohBCoverage ~ scale(Roh[Roh$RohBCoverage > 0,]$RecombActivity) + scale(Roh[Roh$RohBCoverage > 0,]$GcContent) + scale(Roh[Roh$RohBCoverage > 0,]$ProbOfBeingLofIntolerant) + scale(Roh[Roh$RohBCoverage > 0,]$KnKsMouse)); summary(a)
a<-lm(Roh[Roh$RohBCoverage > 0,]$RohBCoverage ~ scale(Roh[Roh$RohBCoverage > 0,]$RecombActivity) + scale(Roh[Roh$RohBCoverage > 0,]$GcContent) + scale(Roh[Roh$RohBCoverage > 0,]$KnKsMouse)); summary(a)
a<-lm(Roh[Roh$RohBCoverage > 0,]$RohBCoverage ~ scale(Roh[Roh$RohBCoverage > 0,]$RecombActivity) + scale(Roh[Roh$RohBCoverage > 0,]$GcContent)); summary(a) # -0.57 + 0.211

a<-lm(Roh[Roh$RohCCoverage > 0,]$RohCCoverage ~ scale(Roh[Roh$RohCCoverage > 0,]$RecombActivity) + scale(Roh[Roh$RohCCoverage > 0,]$GcContent) + scale(Roh[Roh$RohCCoverage > 0,]$ProbOfBeingLofIntolerant) + scale(Roh[Roh$RohCCoverage > 0,]$SelCoeffHet) + scale(Roh[Roh$RohCCoverage > 0,]$KnKsMouse)); summary(a)  # recomb and KnKs
a<-lm(Roh[Roh$RohCCoverage > 0,]$RohCCoverage ~ scale(Roh[Roh$RohCCoverage > 0,]$RecombActivity) + scale(Roh[Roh$RohCCoverage > 0,]$GcContent) + scale(Roh[Roh$RohCCoverage > 0,]$ProbOfBeingLofIntolerant) + scale(Roh[Roh$RohCCoverage > 0,]$KnKsMouse)); summary(a)  # recomb and KnKs
a<-lm(Roh[Roh$RohCCoverage > 0,]$RohCCoverage ~ scale(Roh[Roh$RohCCoverage > 0,]$RecombActivity) + scale(Roh[Roh$RohCCoverage > 0,]$GcContent) + scale(Roh[Roh$RohCCoverage > 0,]$ProbOfBeingLofIntolerant) + scale(Roh[Roh$RohCCoverage > 0,]$KnKsMouse)); summary(a)  # -1.2-0.75-0.22+1.47

dev.off()
