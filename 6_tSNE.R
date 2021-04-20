############################################
# tSNE Mapping using SeqFISH gene expression
############################################

options(max.print = 500)
options(stringsAsFactors = FALSE)
options(scipen = 999)

library("Rtsne")
library("heatmap3")
library(ggplot2)
library(gplots)
library(dplyr)
library(tidyr)
library(scales)
library(reshape2)

scaleMinMax <- function(x, min=0, max=1, z=NULL, keepwithin=TRUE) {
  if(!is.null(z)) {   # use z-score normalization instead
    min <- mean(x) - z*sd(x)
    max <- mean(x) + z*sd(x)
  }
  x <- (x-min) / (max-min)
  if(keepwithin) {
    x[x < 0] <- 0
    x[x > 1] <- 1
  }
  x
}

colCustom <- function(x, z=NULL, colors=c("blue", "#FFFFFF", "red")) {   # use grey to red as default
  if(is.null(z)) {   # just scale to min and max value
    x <- scaleMinMax(x, min(x, na.rm=TRUE), max(x, na.rm=TRUE))
  } else if(length(z) == 1) {   # zscore
    x <- scaleMinMax(x, z = z, keepwithin = TRUE)
  } else {  # scale to min max as provided
    x <- scaleMinMax(x, min = z[1], max = z[2], keepwithin = TRUE)
  }
  
  m <- is.na(x)
  x[m] <- 0.5
  r <- colorRamp(colors)(x)
  y <- apply(r, 1, function(rr) rgb(rr[1], rr[2], rr[3], maxColorValue = 255))
  y[m] <- NA
  y
}


## load data
setwd("~/Nanowell")

d <- read.csv("scData_Int.csv")
str(d)
table(table(d$numCell))

d <- d[d$dxdy_coeff>0.99,]

# normalize/scale
i <- colnames(d)[10:18]
d[10:18] <- d[10:18]+1
x <- log2(t(as.matrix(d[, i[1:9]]/d$mRNA*100+1)))
y <- x - rowMeans(x)


y.cor1 <- cor(t(y))
y.cor1.hc <- as.dendrogram(hclust(as.dist(1-y.cor1)))

y.cor2 <- cor(y)
y.cor2.hc <- as.dendrogram(hclust(as.dist(1-y.cor2)))


# tsne
set.seed(7)
y.cor2.tsne <- Rtsne(1-y.cor2, perplexity = 10, pca = FALSE, is_distance = TRUE)$Y
plot(y.cor2.tsne, col=l, asp = 1)

