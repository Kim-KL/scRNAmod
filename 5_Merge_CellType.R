###############################
# Add Single Cell Type Profiles
###############################

library(ggplot2)
library(yarrr)
library(dplyr)

setwd("~/Nanowell")


smScan = matrix(0, nrow = 31, ncol = 51)

posNum = 1
for (i in 1:31) {
  for (j in 1:51) {
    if((i%%2) == 0) {
      smScan[i, 52-j] = posNum
      #smScan[i, j] = posNum
    } else {
      smScan[i, j] = posNum
      #smScan[i, 52-j] = posNum
    }
    posNum = posNum + 1
  }
}

data <- as.data.frame(as.vector(smScan))
names(data)[1] <- "smScan"

# smScan <- read.csv(file = "smDensity.csv")


######################################################

# Summary of Nanowell Array (31 by 51)

######################################################

BestCor <- read.csv("BestCor.csv")
posNano <- as.numeric(BestCor[which.max(BestCor[,3]),])

nanoScanG_Int <- read.csv("nanoscanG_Int.csv")
nanoScanR_Int <- read.csv("nanoscanR_Int.csv")

RegiterNanoScanG_Int<- nanoScanG_Int[posNano[2]:(posNano[2]+30), posNano[1]:(posNano[1]+50)]
RegiterNanoScanR_Int <- nanoScanR_Int[posNano[2]:(posNano[2]+30), posNano[1]:(posNano[1]+50)]

data <- cbind(data, as.vector(as.matrix(RegiterNanoScanG_Int)))
data <- cbind(data, as.vector(as.matrix(RegiterNanoScanR_Int)))
names(data)[2:3] <- c("IntG", "IntR")


######################################################

# Summary of FOV Array (31 by 51)

######################################################


scData <- read.csv("scData.csv")
scData <- merge(data, scData, by.x = "smScan", by.y = "smScan")

write.csv(scData, file = "scData_Int.csv", row.names = FALSE)





