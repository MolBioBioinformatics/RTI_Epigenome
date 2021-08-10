##############################################################
# RTI calculation and Random Forrest Model of Predicting RTI #
##############################################################
library(rpart)
library(rpart.plot)
library(randomForest)
library(gbm)
library(caret)
library(MASS)
library(ISLR)
library(pROC)

### RTI calculation ####
cpm.qnorm.loess <- read.table("dat/BrdU.CPM.loess.csv") # BrdU were count over 50Kb bin, quantile normalized and LOESS smoothed as described in [Marchal et al., 2018]
bin.bed <- read.table("dat/50Kb.bin.bed")

bin.crd <- rowMeans(bin.bed[,2:3])/1e6
calc_RTI <- function(dat){
  dat2 <- dat/rowSums(dat)
  tmp <- ((1*dat2[,1]+2*dat2[,2]+3*dat2[,3]+4*dat2[,4])-1)/(4-1)
  return(tmp)
}
RTI.ct1 <- calc_RTI(cpm.qnorm.loess[,4:7])
RTI.ct2 <- calc_RTI(cpm.qnorm.loess[,11:14])
RTI.ct <- (RTI.ct1+RTI.ct2)/2
RTI.4A1 <- calc_RTI(cpm.qnorm.loess[,18:21])
RTI.4A2 <- calc_RTI(cpm.qnorm.loess[,25:28])
RTI.4A <- (RTI.4A1+RTI.4A2)/2

RTI.mtx <- cbind(RTI.ct1,RTI.ct2,RTI.4A1,RTI.4A2)
diffRT <- rowMeans(RTI.mtx[,3:4]) - rowMeans(RTI.mtx[,1:2])

### Random Forrest model predicting ∆RTI from broad histone marks ###
tab.Ctrl.Broad <- as.matrix(read.table("dat/bin.Ctrl.enrich.tab",header=T)[,-(1:3)]) # Input normalized log2 enrichment of all broad histone marks in each 50Kb bin
tab.deltaBroad <- as.matrix(read.table("dat/bin.diff.broad.tab",header=T)[,-(1:3)]) # difference log2 enrichment between Ctrl and KDM4A-OE of all broad histone marks in each 50Kb bin

na.idx <- which(is.na(diffRT)) # bins without BrdU count
inp <- data.frame(diffRT=diffRT, tab.Ctrl.Broad,tab.deltaBroad)[-na.idx,]

inp_idx <- sample(1:nrow(inp), floor(nrow(inp)/2))
inp_trn <- inp[inp_idx,]
inp_tst <- inp[-inp_idx,]

diffRT_forest <- randomForest(diffRT ~ ., data = inp_trn, mtry = 4, importance = TRUE, ntrees = 500)
diffRT_pred <- predict(diffRT_forest, newdata = inp_tst)

png(sprintf("Predict.diffRT.broad.png"),6,6,units="in",res=300)
par(mar=c(4,4,2,2))
par(mgp=c(2,0.5,0))
plot(diffRT_pred, inp_tst$diffRT,col=densCols(diffRT_pred, inp_tst$diffRT,nbin=512),
     xlab = "Predicted ∆RTI", ylab = "Observed ∆RTI",main = "",
     pch = 16,cex=0.3)
dev.off()



