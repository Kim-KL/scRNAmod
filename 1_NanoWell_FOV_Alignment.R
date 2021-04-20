##################################################
# Allocate single cell profile into Nanowell array
# Aligment between Nanowell and FOV array
##################################################

## array for scRNA (31 by 51)
## List of Input Variables
setwd("~/Nanowell")
Threshold = 1000
inputfile = "FOV_Scan.csv"
## depending on bining and objective magnification
wellSize = 83
## a reference nanowell position on the plate image
refX = 14429
refY = 5518
## Cell Cut-Off Intensity (Grn = A488, Red = A647)
# GrnCutOff = 1000
# RedCutOff = 900




## Position_Intensity and Density
temp = read.csv(inputfile)
## - Density
temp <- temp[3:3164, c("Count")]
temp <- as.data.frame(temp)
temp <- temp[c(FALSE,TRUE),]
temp <- as.data.frame(temp)

smScan = matrix(0, nrow = 31, ncol = 51)

posNum = 1
for (i in 1:31) {
  for (j in 1:51) {
    if((i%%2) == 0) {
      ##      smScan[i, 52-j] = temp[posNum, 1]
      smScan[i, j] = temp[posNum, 1]
    } else {
      ##      smScan[i, j] = temp[posNum,1]
      smScan[i, 52-j] = temp[posNum, 1]
    }
    posNum = posNum + 1
  }
}

## Vertical Flip
smScan <- smScan[,51:1]
## Horizontal Flip
##smScan <- smScan[31:1,]

write.csv(smScan, file = "smDensity.csv", row.names = FALSE)


## cell occupancy array for scRNA (31 by 51)

smScan2 = matrix(0, nrow = 31, ncol = 51)

for (i in 1:31) {
  for(j in 1:51) {
    ## Single Molecule Density on each Position
    if(smScan[i,j]>Threshold) {
      smScan2[i,j] = 1
    } else {
    }
  }
}


## array for Nanowell Plate (128 by 385)

temp = read.csv("A405-Position.csv", row.names = 1)  # x,y coordinate for cell position
posCell =  temp

## 2-color intensity
temp = read.csv("A488-Intensity.csv", row.names = 1)    # A488 intensity of cell positions
posCell <- cbind(posCell, subset(temp, select = c(Mean)))
names(posCell)[3] <- paste("Green")

temp = read.csv("A555-Intensity.csv", row.names = 1)    # A555 intensity of cell positions
posCell <- cbind(posCell, subset(temp, select = c(Mean)))
names(posCell)[4] <- paste("Red")

write.csv(posCell, file = "Cell_Color_Intensity.csv", row.names = FALSE)

nanoScan = matrix(0, nrow = 140, ncol = 400)
nanoScanG_Int = matrix(0, nrow = 140, ncol = 400)
nanoScanR_Int = matrix(0, nrow = 140, ncol = 400)

refX = refX %% wellSize
refY = refY %% wellSize

write.csv(c(wellSize, refX, refY), file = "nanoScan_offset.csv")

for (i in 1:nrow(posCell)) {
  if(posCell[i,1] < refX){
    cNum = 1;
  } else {
    cNum = as.integer((posCell[i,1] - refX)/wellSize)
  }
  
  if(posCell[i,2] < refY){
    rNum = 1;
  } else {
    rNum = as.integer((posCell[i,2] - refY)/wellSize)
  }
  
  nanoScan[rNum, cNum] = nanoScan[rNum, cNum] + 1
  nanoScanG_Int[rNum, cNum] = nanoScanG_Int[rNum, cNum] + posCell[i,3]
  nanoScanR_Int[rNum, cNum] = nanoScanR_Int[rNum, cNum] + posCell[i,4]
}

write.csv(nanoScan, file = "nanoscan.csv", row.names = FALSE)
write.csv(nanoScanG_Int, file = "nanoscanG_Int.csv", row.names = FALSE)
write.csv(nanoScanR_Int, file = "nanoscanR_Int.csv", row.names = FALSE)



## Green array for Nanowell Plate (128 by 385)
nanoScanG = matrix(0, nrow = 140, ncol = 400)

for (i in 1:nrow(posCell)) {
  #if(posCell[i,3] > GrnCutOff) {
  if(posCell[i,3] > posCell[i,4]) {
    if(posCell[i,1] < refX){
      cNum = 1;
    } else {
      cNum = as.integer((posCell[i,1] - refX)/wellSize)
    }
    
    if(posCell[i,2] < refY){
      rNum = 1;
    } else {
      rNum = as.integer((posCell[i,2] - refY)/wellSize)
    }
    
    nanoScanG[rNum, cNum] = nanoScanG[rNum, cNum] + 1
  }
}

write.csv(nanoScanG, file = "nanoscanG.csv", row.names = FALSE)

## Red array for Nanowell Plate (128 by 385)
nanoScanR = matrix(0, nrow = 140, ncol = 400)

for (i in 1:nrow(posCell)) {
  if(posCell[i,4] > posCell[i,3]) {
    if(posCell[i,1] < refX){
      cNum = 1;
    } else {
      cNum = as.integer((posCell[i,1] - refX)/wellSize)
    }
    
    if(posCell[i,2] < refY){
      rNum = 1;
    } else {
      rNum = as.integer((posCell[i,2] - refY)/wellSize)
    }
    
    nanoScanR[rNum, cNum] = nanoScanR[rNum, cNum] + 1
  }
}

write.csv(nanoScanR, file = "nanoscanR.csv", row.names = FALSE)



## Digital (0 or 1) cell occupancy array for Nanowell Plate (140 by 400)
nanoScan2 <- 1*(nanoScan[,]>0)

## Alignment of smScan on nanoScan
smScan3 <- smScan

BestCor = "NULL"
for (i in 21:80) {
  for (j in 85:250) {
    temp <- nanoScan2[i:(i+30), j:(j+50)]
    BestCor = rbind(BestCor, c(j, i, cor(c(temp), c(smScan3))))
  }
}

BestCor <- BestCor[-c(1),]
BestCor[,3] <- formatC(as.numeric(BestCor[,3]), digits = 5, format = "f", mode = "double")

write.csv(BestCor, file = "BestCor.csv", row.names = FALSE)

max(BestCor[,3])
XY <- which.max(BestCor[,3])
XY <- as.numeric(BestCor[XY,])
temp <- nanoScan[XY[2]:(XY[2]+30), XY[1]:(XY[1]+50)]
write.csv(temp, file = "nanoScan_Crop.csv", row.names = FALSE)
temp <- nanoScanG[XY[2]:(XY[2]+30), XY[1]:(XY[1]+50)]
write.csv(temp, file = "nanoScanG_Crop.csv", row.names = FALSE)
temp <- nanoScanR[XY[2]:(XY[2]+30), XY[1]:(XY[1]+50)]
write.csv(temp, file = "nanoScanR_Crop.csv", row.names = FALSE)









