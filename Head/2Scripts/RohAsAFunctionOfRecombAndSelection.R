rm(list=ls(all = TRUE))

## read 
Roh = read.table("../../Body/3Results/GeneRohRecomb.txt", header = TRUE)
Roh = Roh[Roh$RohCoverage > 0,] # zero was default
Roh = Roh[Roh$RecombActivity > 0,] # zero was default

### ROH vs recombination
cor.test(Roh$RohCoverage,Roh$RecombActivity, method = 'spearman') # strong negative
plot(log2(Roh$RohCoverage),log2(Roh$RecombActivity)) # strong negative

### ROH vs selection
table(Roh$GeneType)
summary(Roh[Roh$GeneType != 'protein_coding',]$RohCoverage)
summary(Roh[Roh$GeneType == 'protein_coding',]$RohCoverage)
cor.test(Roh$RohCoverage,Roh$RecombActivity, method = 'spearman') # strong negative
plot(log2(Roh$RohCoverage),log2(Roh$RecombActivity)) # strong negative

cor.test(Roh$RohCoverage,Roh$SelCoeffHet, method = 'spearman') # nothing

##### tests for directions:
cor.test(Roh$SelCoeffHet,Roh$ProbOfBeingLofIntolerant, method = 'spearman') # strong positive


a<-lm(Roh$RohCoverage ~ Roh$RecombActivity + Roh$KnKsMouse); summary(a)  # positive
a<-lm(Roh$RohCoverage ~ Roh$RecombActivity + Roh$Branch); summary(a)  # strong positive

a<-lm(Roh$RohCoverage ~ Roh$RecombActivity + Roh$GcContent); summary(a)  # negative strong
a<-lm(Roh$RohCoverage ~ Roh$RecombActivity + Roh$ProbOfBeingLofIntolerant); summary(a)  # weak negative
a<-lm(Roh$RohCoverage ~ Roh$RecombActivity + Roh$SelCoeffHet); summary(a)               # weak negative

a<-lm(Roh$RohCoverage ~ Roh$RecombActivity + Roh$ResidualVariationIntoleranceScore); summary(a)  # nothing
a<-lm(Roh$RohCoverage ~ Roh$RecombActivity + Roh$GenomeWideHaploinsufficiencyScore); summary(a)  # nothing
a<-lm(Roh$RohCoverage ~ Roh$RecombActivity + Roh$FunctionalIndispensabilityScore); summary(a)  # nothing

### multiple backward linear model:
a<-lm(Roh$RohCoverage ~ Roh$RecombActivity + Roh$GcContent + Roh$KnKsMouse + Roh$Branch + Roh$ProbOfBeingLofIntolerant + Roh$SelCoeffHet); summary(a)
a<-lm(Roh$RohCoverage ~ Roh$RecombActivity + Roh$GcContent + Roh$KnKsMouse + Roh$Branch + Roh$SelCoeffHet); summary(a)
a<-lm(Roh$RohCoverage ~ Roh$RecombActivity + Roh$GcContent + Roh$KnKsMouse + Roh$Branch); summary(a)
a<-lm(Roh$RohCoverage ~ Roh$RecombActivity + Roh$GcContent + Roh$Branch); summary(a)

a<-lm(scale(Roh$RohCoverage) ~ scale(Roh$RecombActivity) + scale(Roh$GcContent) + scale(Roh$Branch)); summary(a)
a<-lm(scale(Roh$RohCoverage) ~ 0 + scale(Roh$RecombActivity) + scale(Roh$GcContent) + scale(Roh$Branch)); summary(a)  
#  scale(Roh$RecombActivity) -0.099733   0.010493  -9.505  < 2e-16 ***
#  scale(Roh$GcContent)      -0.035420   0.007831  -4.523 6.14e-06 ***
#  scale(Roh$Branch)          0.069151   0.007785   8.883  < 2e-16 ***

############# these fractions will change when we take into account only short (recomb is more important), middle and long (selection is more important) ROHs? 



