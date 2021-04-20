###################
# Load SeqFISH data
###################

# Read SeqFISH Gene Name
setwd("~/Nanowell")

# List of SeqFISH Genes
smProbeFile = "List_smP.csv"
# Number of FOVs
NumFOV = 469


# Load the list of SeqFISH Genes
nameFISH <- read.csv(smProbeFile)
nameFISH <- as.character(nameFISH[,2])
nameFISH <- nameFISH[1:9]

setwd("~/SeqFISH_dxdy")

# Number of FISH Cycles
FISH <- length(nameFISH)
nFISH <- FISH + 3

# Initialize Data.Frame
data <- data.frame(matrix(ncol = (nFISH+FISH+1), nrow = 0))

for (i in 1:NumFOV) {
  fname <- paste(toString(i), "_SeqDot.csv", sep = "")
  
  if(file.exists(fname)){
    # Read data
    SeqDot <- read.table(fname, sep = ",")
    names(SeqDot)[1:3] <- c("xcent", "ycent", "m6A")
    names(SeqDot)[4:nFISH] <- nameFISH
    
    # Exclude Fiducial Markers
    SeqDot <- cbind(SeqDot, rowSums(SeqDot[,4:nFISH]))
    names(SeqDot)[(nFISH + 1)] <- "SumFISH"
    SeqDot <- SeqDot[SeqDot$SumFISH<2,] # No fiducial markers
    
    # Single Cell Gene Expression Table
    SumDot <- colSums(SeqDot)
    SeqDot2 <- SeqDot[SeqDot$m6A==1,]      # fiducial marker
    names(SeqDot2) <- paste(names(SeqDot2), "_m6A", sep="")    # fiducial marker
    SumDot2 <- colSums(SeqDot2)            # fiducial marker
    #data <- rbind(data, SumDot[2:(nFISH + 1)])
    data <- rbind(data, c(SumDot[2:(nFISH + 1)], SumDot2[4:(nFISH + 1)]))    # fiducial marker
    data[nrow(data),1] <- i
  }
}

names(data)[1:2] <- c("ID", "m6A")
names(data)[3:(nFISH-1)] <- nameFISH
names(data)[nFISH] <- "SumFISH"
names(data)[(nFISH+1):(nFISH+FISH)] <- paste(nameFISH, "_m6A", sep = "")  # fiducial marker
names(data)[(nFISH+FISH+1)] <- "SumFISH_m6A"   # fiducial marker


## add dxdy correlation coefficient

refPos = 2
# reference position

dname <- getwd()
wdname <- paste(dname, "/Position ", toString(refPos), sep = "")
setwd(wdname)
fname <- "dxdy.csv"
refdxdy <- read.table(fname, sep = ",")

data <- cbind(data, dxdy_coeff=0)

for (i in 1:NumFOV) {
  wdname <- paste(dname, "/Position ", toString(i), sep = "")
  setwd(wdname)
  
  if(file.exists(fname)){
    # Read dxdy & Corr Coef calculation
    dxdy <- read.table(fname, sep = ",")
    dxdy_cor <- cor(c(as.matrix(refdxdy[,1:2])), c(as.matrix(dxdy[,1:2])))
    data[i, ncol(data)] <- dxdy_cor
  }
}

## Save data
setwd("~/Nanowell")
write.csv(data, file = "scFISH.csv")
