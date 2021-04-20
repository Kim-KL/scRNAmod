#########################################
# Add the quantifications of mRNA and m6A
#########################################

library(ggplot2)
library(dplyr)
library(yarrr)

mRNAfile = "FOV_mRNA.csv"
m6Afile = "FOV_m6A.csv"

setwd("~/Nanowell")

data <- read.csv("Filtered_Reg_FOV.csv", row.names = 1)
scFISH <- read.csv("scFISH.csv", row.names = 1)
# match the FOV number to the counting number on the FOV Matrix
scFISH[,1] <- scFISH[,1] - 1


data <- merge(data, scFISH, by.x="numFOV", by.y="ID")
write.csv(data, file="Filtered_Reg_FOV_scFISH.csv", row.names = FALSE)

scData <- data[data[,3]==1,]


mRNA <- read.csv(mRNAfile)
mRNA <- mRNA[c(FALSE,TRUE),]
mRNA$X <- (mRNA$X/2)-1
mRNA <- mRNA[, c("X", "Count")]

scData <- merge(scData, mRNA, by.x = "numFOV", by.y = "X")
names(scData)[dim(scData)[2]] <- "mRNA_Count"
scData$mRNA <- scData$mRNA_Count
scData <- scData[, -c(dim(scData)[2])]



m6A <- read.csv(m6Afile)
m6A <- m6A[c(FALSE,TRUE),]
m6A$X <- (m6A$X/2)-1
m6A <- m6A[, c("X", "Count")]

scData <- merge(scData, m6A, by.x = "numFOV", by.y = "X")
names(scData)[dim(scData)[2]] <- "m6A_Count"
scData$m6A <- scData$m6A_Count
scData <- scData[, -c(dim(scData)[2])]

write.csv(scData, file = "scData.csv", row.names = FALSE)





