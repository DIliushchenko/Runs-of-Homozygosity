rm(list = ls(all = TRUE))  # remove everything from R memory (old variables, datasets...)

RohFull <- read.table("../../Body/2Derived/Roh.txt")

summary(RohFull$KB)

#Agg = aggregate(
#    RohFull$KB,
#    by = list(
#        RohFull$FidIid,
#        RohFull$FID,
#        RohFull$FamilyMember,
#        RohFull$Family,
#        RohFull$PHE
#    ),
#    FUN = sum
#)
#names(Agg) = c('FidIid', 'FID', 'FamilyMember', 'Family', 'PHE', 'RohTotal')
#Agg = Agg[Agg$Family == 'CompleteFamily', ]
#Agg = Agg[order(Agg$FID), ]

Agg = aggregate(
    RohFull$KB,
    by = list(
        RohFull$FidIid,
        RohFull$FID,
        RohFull$FamilyMember,
        RohFull$Family,
        RohFull$PHE,
        RohFull$KidsInFamily,
        RohFull$BothParentsAreHealthy
    ),
    FUN = sum
)
names(Agg) = c(
    'FidIid',
    'FID',
    'FamilyMember',
    'Family',
    'PHE',
    'KidsInFamily',
    'BothParentsAreHealthy',
    'RohTotal'
)
Agg = Agg[Agg$Family == 'CompleteFamily' &
              Agg$KidsInFamily == 'KidsAffectedAndNonAffected' &
              Agg$BothParentsAreHealthy == 'Yes'
          & Agg$FamilyMember == 'Kid' & Agg$PHE != -9, ]

Agg = Agg[order(Agg$FID),]


options(stringsAsFactors = FALSE)
full_list <- list(list("Family", "Id1", "Id2", "Difference"))
k <- 2

for (i in 1:nrow(Agg)) {
    elem1 <- Agg[i, ]
    
    if (is.na(elem1$FID)) {
        next
    }
    
    for (j in i + 1:nrow(Agg)) {
        elem2 <- Agg[j, ]
        
        if (is.na(elem2$FID)) {
            next
        }
        
        if (elem1$FID == elem2$FID) {
            if (elem1$PHE > elem2$PHE) {
                result <- elem1$RohTotal - elem2$RohTotal
                
                l <-
                    list(elem1$FID, elem1$FidIid, elem2$FidIid, result)
                full_list[[k]] <- l
                
                k <- k + 1
            }
            
            if (elem1$PHE < elem2$PHE) {
                result <- elem2$RohTotal - elem1$RohTotal
                
                l <-
                    list(elem1$FID, elem1$FidIid, elem2$FidIid, result)
                full_list[[k]] <- l
                
                k <- k + 1
            }
        }
    }
}

fr <- data.frame(t(sapply(full_list, c)))
View(fr)
