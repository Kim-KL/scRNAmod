####################################
# Arrange NanoWell_Info and FOV_Info
####################################

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


# Read BestCor 
BestCor <- read.csv("BestCor.csv")
posNano <- as.numeric(BestCor[which.max(BestCor[,3]),])

# Read number of cells in a well
nanoScan <- read.csv("nanoscan.csv")
nanoFOV <- nanoScan[posNano[2]:(posNano[2]+30), posNano[1]:(posNano[1]+50)]

data <- cbind(data, as.vector(as.matrix(nanoFOV)))
names(data)[2] <- "numCell"

# Read number of 'Green' cells in a well
nanoScanG <- read.csv("nanoscanG.csv")
nanoFOVG <- nanoScanG[posNano[2]:(posNano[2]+30), posNano[1]:(posNano[1]+50)]

data <- cbind(data, as.vector(as.matrix(nanoFOVG)))
names(data)[3] <- "GrnCell"

# Read number of 'Red' cells in a well
nanoScanR <- read.csv("nanoscanR.csv")
nanoFOVR <- nanoScanR[posNano[2]:(posNano[2]+30), posNano[1]:(posNano[1]+50)]

data <- cbind(data, as.vector(as.matrix(nanoFOVR)))
names(data)[4] <- "RedCell"

# Read number of captured RNA molecules in a FOV
smRNA <- read.csv("smDensity.csv")

data <- cbind(data, as.vector(as.matrix(smRNA)))
names(data)[5] <- "mRNA"

write.csv(data, file="RegFOV.csv")


# Read the index of filtered FOVs
index <- read.csv("Filtered_FOV_Index.csv", header = FALSE)
index <- index[,1:2]

data <- merge(data, index, by.x="smScan", by.y="V1")
names(data)[6] <- "numFOV"

write.csv(data, file="Filtered_Reg_FOV.csv")



