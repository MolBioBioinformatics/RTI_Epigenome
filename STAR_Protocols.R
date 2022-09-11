#' Load required packages 
library(preprocessCore)
library(runner)
library(tidyverse)
library(randomForest)
library(bedtoolsr)
library(ComplexHeatmap)
library(dplyr)
library(circlize)

#' ### Replication Timing (RT) data QC 
#' Load Replication Timing (RT) genomic bin coverage.
#' RT_rawCount.txt counted 
tab <- read.table("dat/RT_rawCount.txt", header=TRUE)
colnames(tab) <- gsub(".bam","",colnames(tab))
rownames(tab) <- paste("bin",1:dim(tab)[1],sep="")

RT.bin.cnt <- as.matrix(tab[,-(1:3)])
bin.bed <- tab[,1:3]
bin_crd <- rowMeans(bin.bed[,2:3])
bin.bed$binID <- rownames(bin.bed)
RT.bin.cpm <- t(t(RT.bin.cnt)/(colSums(RT.bin.cnt)/1e6))

#' Genomic coverage scatter plot for visual inspection of consistency between replicates
plot(RT.bin.cpm[,"RPE.ct1.S1"],
     RT.bin.cpm[,"RPE.ct2.S1"],
     col=densCols(log(RT.bin.cpm[,"RPE.ct1.S1"]+0.1),log(RT.bin.cpm[,"RPE.ct2.S1"]+0.1)),
     xlab="RPE.ct1.S1",ylab="RPE.ct2.S1",
     cex=0.5,pch=16,
     log="xy") %>% suppressWarnings()

#' Calculate Pearson correlation between replicates
r <- cor(log10(RT.bin.cpm[,"RPE.ct1.S1"]+0.01), log10(RT.bin.cpm[,"RPE.ct2.S1"]+0.01))
message(sprintf("Pearson R=%.3f",r))

#' Principal Components Analysis over all Repli-seq samples
pca <- prcomp(t(log10(RT.bin.cpm+0.01)))
plot(pca$x[,1],pca$x[,2],type="n",xlab="PC1",ylab="PC2",xlim=c(-200,200))
text(pca$x[,1],pca$x[,2],rownames(pca$x))

#' ### Repli-seq data quantile normalization and LOESS smoothing
#' Genome wide quantile normalization and LOESS smoothing
qn <- function(inp){
  out <- normalize.quantiles(inp)
  rownames(out) <- rownames(inp)
  colnames(out) <- colnames(inp)
  return(out)
}
RT.bin.S1.norm <- qn(RT.bin.cpm[,c("RPE.ct1.S1","RPE.ct2.S1","RPE.4A1.S1","RPE.4A2.S1")])
RT.bin.S2.norm <- qn(RT.bin.cpm[,c("RPE.ct1.S2","RPE.ct2.S2","RPE.4A1.S2","RPE.4A2.S2")])
RT.bin.S3.norm <- qn(RT.bin.cpm[,c("RPE.ct1.S3","RPE.ct2.S3","RPE.4A1.S3","RPE.4A2.S3")])
RT.bin.S4.norm <- qn(RT.bin.cpm[,c("RPE.ct1.S4","RPE.ct2.S4","RPE.4A1.S4","RPE.4A2.S4")])

RT.bin.norm <- cbind(RT.bin.S1.norm, RT.bin.S2.norm, RT.bin.S3.norm, RT.bin.S4.norm)
RT.bin.loess <- RT.bin.norm #initialization 

for(i in 1:dim(RT.bin.norm)[2]){
  for(chr in unique(bin.bed[,1])){
    idx <- which(bin.bed[,1]==chr)
    x <- bin.bed[idx,2]
    y <- RT.bin.norm[idx,i]
    lspan <-500000/(max(bin_crd[idx])-min(bin_crd[idx]))
    Rpla <- loess(y ~ x, span=lspan) %>% suppressWarnings()
    RT.bin.loess[idx,i] <- Rpla$fitted
  }
}
RT.bin.loess[RT.bin.loess<0]=0

#' Example of RT coverage before and after pre-processing.
idx <- bin.bed %>% filter(chr=="chr1" & start>10e6 & end<15e6) %>% rownames()
plot(bin_crd[idx]/1e6,RT.bin.cpm[idx,"RPE.ct1.S1"],type="l",col="grey",xlab="Coordinate (Mb)",ylab="CPM",main="CT1.S1@chr1:10Mb-11Mb")
lines(bin_crd[idx]/1e6,RT.bin.loess[idx,"RPE.ct1.S1"],col="red")
legend(x="topright",col=c("grey","red"),lty=1,legend = c("raw CPM","LOESS CPM"))

#' QC and quantile normalization steps can be applied to ChIP-seq coverage over same genomic bins.

#' ### Sub-clustering genomic regions with the same chromatin state by their replication timing
#' Chromatin states were learn using ChromHMM from combination of histone marks.
#' To understand the general replication behavior of regions under specific
#' chromatin states, it is informative to define and analyze subgroups of
#' regions sharing similar replication timing within each chromatin state.

chromHMM.bed <- read.table("dat/broad.chromatin.state.bed",sep="\t")[,1:4]
colnames(chromHMM.bed) <- c("chr","start","end","state")

target.state <- "act"
target.state.bins <- bedtoolsr::bt.intersect(bin.bed,chromHMM.bed[chromHMM.bed$state==target.state,],wa=T)[,4]
state.RT <- log10(RT.bin.norm[target.state.bins,c("RPE.ct1.S1","RPE.ct1.S2","RPE.ct1.S3","RPE.ct1.S4")]+0.1)
km <- kmeans(state.RT,3)

ComplexHeatmap::Heatmap(state.RT,
                        use_raster=TRUE,
                        col=c("white","grey10"),
                        show_row_names = FALSE,
                        show_column_names = T,
                        show_row_dend = FALSE,
                        cluster_columns=FALSE,
                        cluster_rows = F,
                        split = km$cluster,
                        column_labels = c("S1","S2","S3","S4"),
                        cluster_row_slices=T,
                        name='log10(RPE)')
#' Using cell cycle H3K36me3 ChIP-seq coverage as an example to investigate
#' above histone modification temporal dynamics at each time point (G1, ES, LS,
#' G2/M) through cell cycle.
#' 
#' Load H3K36me3 and Input ChIP-seq genomic coverage were normalized as shown before
H3K36me3.bin.cpm <- read.table("dat/H3K36me3.qn.txt",header=T,check.names = F)
Input.bin.cpm <- read.table("dat/Input.qn.txt",header=T,check.names = F)

#' Calculate log2 ChIP over Input enrichment
pseudo=0.1
H3K36me3.bin.enrich <- log2(H3K36me3.bin.cpm+pseudo)-log2(Input.bin.cpm+pseudo)
rownames(H3K36me3.bin.enrich) <- rownames(bin.bed)
H3K36me3.ctrl.enrich <- (H3K36me3.bin.enrich[,c("ct1-G1-K36me3","ct1-ES-K36me3","ct1-LS-K36me3","ct1-G2-K36me3")]+
                             H3K36me3.bin.enrich[,c("ct2-G1-K36me3","ct2-ES-K36me3","ct2-LS-K36me3","ct2-G2-K36me3")])/2

#' Calculate enrichment dynamic through cell cycle
H3K36me3.cycle <- H3K36me3.ctrl.enrich-rowMeans(H3K36me3.ctrl.enrich)

#' Visualize H3K36me3 cell cycle dynamic in each chromatin subset defined by RT above
col <- colorRampPalette(c("purple","black","yellow"))
ComplexHeatmap::Heatmap(H3K36me3.cycle[target.state.bins,] %>% as.matrix(),
                        use_raster=TRUE,
                        show_row_names = FALSE,
                        show_column_names = FALSE,
                        show_row_dend = FALSE,
                        cluster_columns=FALSE,
                        cluster_rows=FALSE,
                        col=colorRamp2(seq(-0.4,0.4,length.out=100), col(100)),
                        split = km$cluster)


#' ### Calculate Replication Timing Index (RTI) and RTI change
#' Calculate RTI based on Repli-seq data at multiple time points across cell
#' cycle in all samples (2 replicates of control cells and 2 replicates of
#' over-expression cells)
calc_rt_index <- function(dat) {
  dat.norm <- dat/rowSums(dat)
  sum <- rep(0,dim(dat.norm)[1])
  for(i in 1:dim(dat.norm)[2]){
    sum = sum+i*dat.norm[,i]
  }
  RTI = (sum-1)/(dim(dat.norm)[2]-1) # scale from 0 to 1 instead or 1 to 4
  return(RTI)
}
RTI.ctrl.1 <- calc_rt_index(RT.bin.loess[,c("RPE.ct1.S1", "RPE.ct1.S2", "RPE.ct1.S3", "RPE.ct1.S4")])
RTI.ctrl.2 <- calc_rt_index(RT.bin.loess[,c("RPE.ct2.S1", "RPE.ct2.S2", "RPE.ct2.S3", "RPE.ct2.S4")])
RTI.oe.1 <- calc_rt_index(RT.bin.loess[,c("RPE.4A1.S1", "RPE.4A1.S2", "RPE.4A1.S3", "RPE.4A1.S4")])
RTI.oe.2 <- calc_rt_index(RT.bin.loess[,c("RPE.4A2.S1", "RPE.4A2.S2", "RPE.4A2.S3", "RPE.4A2.S4")])
RTI.ctrl <- rowMeans(cbind(RTI.ctrl.1, RTI.ctrl.2))
RTI.oe <- rowMeans(cbind(RTI.oe.1, RTI.oe.2))

#' Quantify the magnitude and statistical significance of RTI differences between biological conditions
RTI.diff <- rowMeans(cbind(RTI.oe.1, RTI.oe.2)) - rowMeans(cbind(RTI.ctrl.1, RTI.ctrl.2))
# RTI.diff[is.na(RTI.diff)] <- 0
RTI.rep.sd <- sd(c(RTI.oe.1- RTI.oe.2, RTI.ctrl.1- RTI.ctrl.2),na.rm=T)
RTI.diff <- data.frame(bin.bed,
                       diff=RTI.diff,
                       P.value=2*pnorm(-abs(RTI.diff/RTI.rep.sd)))

#' List top differential RT genomic bins
RTI.diff[order(RTI.diff$P.value),] %>% head()

#' Example of a differential RT region
idx <- bin.bed %>% filter(chr=="chr5" & start>3e6 & end<8e6) %>% rownames()
plot(bin_crd[idx]/1e6,1-RTI.ctrl[idx],type="l",col="black",xlab="Coordinate (Mb)",ylab="RTI",
     main="RTI@chr5: 3Mb-8Mb",ylim=c(0,1),yaxt="n")
axis(2,0:1,1:0,las=2)
lines(bin_crd[idx]/1e6,1-RTI.oe[idx],col="red")
legend(x="topright",col=c("black","red"),lty=1,legend = c("Control","Over-Expression"))

#' Output differential RT regions to a bed file.
RTI.diff.bins <- bin.bed[RTI.diff %>% filter(abs(diff)>0.05 & P.value<0.01) %>% select(binID) %>% as.matrix(),]
write.table(bedtoolsr::bt.merge(RTI.diff.bins[order(RTI.diff.bins$chr,RTI.diff.bins$start),]),
            "sample_output/diff.RT.bed",
            quote=F,row.names = F,col.names = F,sep="\t")


#' ### Train and evaluate computational models to predict RTI from the
#' combination of epigenomic and transcriptomic data
tab.Ctrl.Broad <- as.matrix(read.table("dat/bin.Ctrl.enrich.tab",header=T)[,-(1:3)]) # Input normalized log2 enrichment of all broad histone marks in each 50Kb bin
tab.sync.Ctrl.Broad <- as.matrix(read.table("dat/bin.sync.ctrl.enrich.tab",header=T)[,-(1:3)]) # Input normalized log2 enrichment of all broad histone marks across cell cycle

tab.deltaBroad <- as.matrix(read.table("dat/bin.diff.broad.tab",header=T)[,-(1:3)]) # difference log2 enrichment between Ctrl and KDM4A-OE of all broad histone marks in each 50Kb bin
tab.sync.deltaBroad <- as.matrix(read.table("dat/bin.diff.broad.cycle.tab",header=T)[,-(1:3)]) # difference log2 enrichment between Ctrl and KDM4A-OE of all broad histone marks in each 50Kb bin across cell cycle
tab.diffPeak <- as.matrix(read.table("dat/bin.diff.peak.cnt.tab",header=T)[,-(1:3)]) # number of differential peaks between Ctrl and KDM4A-OE in each 50Kb bin

#' Split training and testing sets
rand.chr.list <- sample(unique(bin.bed[,1]))
inp_idx <- c()
for(chr in rand.chr.list){
  message(chr)
  inp_idx <- c(inp_idx,which(bin.bed[,1]==chr))
  if(length(inp_idx) > floor(nrow(bin.bed)/4)){
    break
  }
}
dat <- data.frame(RT=RTI.ctrl, tab.sync.Ctrl.Broad) 
dat_trn = dat[-inp_idx,] %>% filter(!is.na(RT)) # Training 
dat_tst = dat[inp_idx,] %>% filter(!is.na(RT)) # Test

#' Train and predict using linear model
fit <- lm(RT ~ .,data = dat_trn)
pred <- predict(fit, newdata = dat_tst)
plot(dat_tst$RT, pred, ylab = "Predicted", xlab = "Observed",
     col=densCols(pred, dat_tst$RT),cex=0.5,pch=16,main="RTI")

#' Train and predict using random forrest model
fit <- randomForest(RT ~ ., data = dat_trn, mtry = round(sqrt(dim(dat_trn)[2]-1)), importance = TRUE, ntrees = 500)
pred <- predict(fit, newdata = dat_tst)
plot(dat_tst$RT, pred, ylab = "Predicted", xlab = "Observed",
     col=densCols(pred, dat_tst$RT),cex=0.5,pch=16,main="RTI")

#' ### Train and evaluate computational models to predict RTI change
#' from the combination of epigenomic and transcriptomic data
dat <- data.frame(diffRT=RTI.diff$diff, tab.sync.deltaBroad,tab.diffPeak)
dat_trn = dat[-inp_idx,] %>% filter(!is.na(diffRT))# Training 
dat_tst = dat[inp_idx,] %>% filter(!is.na(diffRT))# Test

#' Train and predict using linear model
fit <- lm(diffRT ~ .,data = dat_trn)
pred <- predict(fit, newdata = dat_tst)
plot(dat_tst$diffRT, pred, ylab = "Predicted", xlab = "Observed",
     col=densCols(pred, dat_tst$diffRT),
     cex=0.5,pch=16,main="OE vs Control RTI change")

#' Train and predict using random forrest model
fit <- randomForest(diffRT ~ ., data = dat_trn, mtry = round(sqrt(dim(dat_trn)[2]-1)), importance = TRUE, ntrees = 500)
pred <- predict(fit, newdata = dat_tst)
plot(dat_tst$diffRT, pred, ylab = "Predicted", xlab = "Observed",
     col=densCols(pred, dat_tst$diffRT),
     cex=0.5,pch=16,main="OE vs Control RTI change")

#' ### Define genomic regions with specific types of local replication patterns
#' Using genome-wide RTI values calculated at the previous step, define genomic
#' locations of initiation zones, termination sites, constant timing regions.
w <- 10
var.cutoff<-1e-3
bs <- head(bin.bed[,3]-bin.bed[,2],1)
dist <- floor(w/2)*bs
rti = RTI.ctrl
RT_local_pattern <- function(rti,output.name){
  win.stat <- data.frame(max=runner(rti,k=w,f=max),
                         min=runner(rti,k=w,f=min),
                         var=runner(rti,k=w,f=var),
                         cnt.noneNA=runner(rti,k=w,f=function(x){sum(!is.na(x))}),
                         cnt.chr=runner(bin.bed[,1],k=w,f=function(x){length(unique(x))}))
  win.stat <- win.stat[-(1:floor(w/2)),] # sliding window center at current position
  win.stat[win.stat$cnt.noneNA<w | win.stat$cnt.chr>1,c("max","min","var")] <- NA
  bin.min.bed <- bin.bed[which(rti[1:dim(win.stat)[1]]==win.stat$min & win.stat$var>var.cutoff),]
  bin.max.bed <- bin.bed[which(rti[1:dim(win.stat)[1]]==win.stat$max & win.stat$var>var.cutoff),]
  bin.const.bed <- bin.bed[which(win.stat$var<var.cutoff),]
  
  bin.min.bed %>% bedtoolsr::bt.merge(d=dist) %>%
    bedtoolsr::bt.window(b=bin.const.bed,w=dist,v=T) %>%
    write.table(file=sprintf("%s.IZ.bed",output.name),row.names = F,col.names = F,sep="\t",quote=F)
  bin.max.bed %>% bedtoolsr::bt.merge(d=dist) %>%
    bedtoolsr::bt.window(b=bin.const.bed,w=dist,v=T) %>%
    write.table(file=sprintf("%s.TS.bed",output.name),row.names = F,col.names = F,sep="\t",quote=F)
  bedtoolsr::bt.merge(bin.const.bed) %>%
    write.table(file=sprintf("%s.CTR.bed",output.name),row.names = F,col.names = F,sep="\t",quote=F)
}
RT_local_pattern(RTI.ctrl,"Ctrl.RT")


