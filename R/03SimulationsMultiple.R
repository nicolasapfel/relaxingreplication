rm(list=ls())
source("01Functions_repo.R")

m <- 1000
n <- 1000
L <- 20

set.seed(42)
Z <- matrix(rnorm(n*L,0,1), ncol=L)
Zc <- cbind(1,Z)

# I. P=1 ---- 

gammam <- rep(1, 20)

epsvar <- 0.5
a_all <- 1:20

# 1. ONE ENDOGENOUS REGRESSOR

comb <- combinations(L,1, repeats.allowed=F)
Len <- choose(L,1)

nr.cores <- detectCores() - 1

# Initiate cluster
cl <- makeCluster(nr.cores)

registerDoRNG(152, once=FALSE)
GD <- foreach(s=1:18, .combine=rbind) %dopar%
  once(s,P=1)
stopCluster(cl)
save.image('../temp/Onem100n10000.RData')
load('../temp/Onem100n10000.RData')
GD <- data.frame(matrix(unlist(GD), ncol=8))
colnames(GD) <- c("biasn1", "biasor1", "biasada1", "sdn1", "sdor1", "sdada1", "nInvada", "freqAllada")
biasn <- abs(GD$biasn1)
biasor <- abs(GD$biasor1)
biasada <- abs(GD$biasada1)
nInvada <- GD$nInvada
freqAllada <- GD$freqAllada
freqAllhc <- GD$freqAllhc
GD <- data.frame(cbind(biasn,biasor,biasada,nInvada,freqAllada))

par(mfrow=c(1,3))
par(mar=c(5,5,1,1), xpd=F)
plot(1:18, biasn, type="p",pch=1, lwd=2, col="grey28", ylim=c(0,10), xlab="# of invalid IVs", ylab="Sum of MAE")#, + , "|",hat(beta[2]), "|"
points(1:18, biasor, pch=2, cex=0.85, col="grey28", lty=5, lwd=2)
points(1:18, biasada, pch=4, cex=0.85, col="black", lwd=1.8)
abline(v=10, lty=3)
par(mar=c(5,5,1,1), xpd=F)
plot(1:18, rep(0,18), type="p",pch=1, lwd=2, main ="One regressor", col="grey28", ylim=c(0,20), xlab="# of invalid IVs", ylab="# of IVs selected as invalid")#, + , "|",hat(beta[2]), "|"
points(1:18, 1:s, pch=2, cex=0.85, col="grey28", lty=5, lwd=2, xlab="# of invalid IVs", ylab="# of IVs selected as invalid")
points(1:18, nInvada, pch=4, cex=0.85, col="black", lty=5, lwd=1.8, xlab="# of invalid IVs", ylab="# of IVs selected as invalid")
abline(v=10, lty=3)
plot(1:18, rep(0,18), type="p",pch=1, lwd=2, col="grey28", ylim=c(0,1), xlab="# of invalid IVs", ylab="Freq. with which all invalid selected as invalid")#, + , "|",hat(beta[2]), "|"
points(1:18, rep(1,18), pch=2, cex=0.85, col="grey28", lty=5, lwd=2, xlab="# of invalid IVs", ylab="# of regressors selected as invalid")
points(1:18, freqAllada, pch=4, cex=0.85, col="black", lty=5, lwd=1.8, xlab="# of invalid IVs", ylab="# of regressors selected as invalid")
abline(v=10, lty=3)

plot(0,0)
legend("topright", c("naive", "oracle", "ada"), col=c("grey28", "grey28", "black"), pch=c(1,2,4), lwd=2, lty=0)

# II.P=2 ----

rm(list=setdiff(ls(),c("once", "Z", "Zc", "MCmad", "MCsd", "nr.cores", "epsvar", "a_all")))

m <- 1000
n <- 1000
L <- 20

gamma1 <- seq(0.05, 1, 0.05)
gamma2 <- seq(1, 0.05, by=-0.05)
gammam <- cbind(gamma1, gamma2)
a_all <- 1:20

# 1. TWO ENDOGENOUS REGRESSORS

comb <- combinations(L,2, repeats.allowed=F)
Len <- choose(L,2)

cl <- makeCluster(nr.cores)

registerDoRNG(152, once=FALSE)
GD <- foreach(s=1:18, .combine=rbind) %dopar%
  once(s,P=2)
stopCluster(cl)

GD <- data.frame(matrix(unlist(GD), ncol=14))
colnames(GD) <- c("biasn1", "biasn2", "biasor1","biasor2", "biasada1", "biasada2", "sdn1", "sdn2", "sdor1", "sdor2", "sdada1", "sdada2", "nInvada", "freqAllada")
biasn <- abs(GD$biasn1) + abs(GD$biasn2)
biasor <- abs(GD$biasor1) + abs(GD$biasor2)
biasada <- abs(GD$biasada1) + abs(GD$biasada2)
nInvada <- GD$nInvada
freqAllada <- GD$freqAllada

#Create .pdf graph
par(mfrow=c(1,3))
par(mar=c(5,5,1,1), xpd=F)
plot(1:18, biasn, type="p",pch=1, lwd=2, col="grey28", ylim=c(0,20), xlab="# of invalid IVs", ylab="Sum of MAE")#, + , "|",hat(beta[2]), "|"
points(1:18, biasor, pch=2, cex=0.85, col="grey28", lty=5, lwd=2)
points(1:18, biasada, pch=4, cex=0.85, col="black", lwd=1.8)
abline(v=5, lty=3)
par(mar=c(5,5,1,1), xpd=F)
plot(1:18, rep(0,18), type="p",pch=1, lwd=2, main ="Two regressors", col="grey28", ylim=c(0,20), xlab="# of invalid IVs", ylab="# of IVs selected as invalid")#, + , "|",hat(beta[2]), "|"
points(1:18, 1:s, pch=2, cex=0.85, col="grey28", lty=5, lwd=2, xlab="# of invalid IVs", ylab="# of regressors selected as invalid")
points(1:18, nInvada, pch=4, cex=0.85, col="black", lty=5, lwd=1.8, xlab="# of invalid IVs", ylab="# of regressors selected as invalid")
abline(v=5, lty=3)
plot(1:18, rep(0,18), type="p",pch=1, lwd=2, col="grey28", ylim=c(0,1), xlab="# of invalid IVs", ylab="Freq. with which all invalid selected as invalid")#, + , "|",hat(beta[2]), "|"
points(1:18, rep(1,18), pch=2, cex=0.85, col="grey28", lty=5, lwd=2, xlab="# of invalid IVs", ylab="# of regressors selected as invalid")
points(1:18, freqAllada, pch=4, cex=0.85, col="black", lty=5, lwd=1.8, xlab="# of invalid IVs", ylab="# of regressors selected as invalid")
abline(v=5, lty=3)

# III.P=3 ----

rm(list=setdiff(ls(),c("once", "Z", "Zc", "MCmad", "MCsd", "nr.cores", "epsvar", "a_all")))

m <- 1000
n <- 1000
L <- 20

gamma1 <- seq(0.05,1,0.05)
gamma2 <- seq(1,0.05,-0.05)
gamma3 <- rep(c(0.05,0.1, 0.15, 2), L/4)
gammam <- cbind(gamma1, gamma2, gamma3)
a_all <- 1:20

comb <- combinations(L,3, repeats.allowed=F)
Len <- choose(L,3)

cl <- makeCluster(nr.cores)

registerDoRNG(152, once=FALSE)
GD <- foreach(s=1:17, .combine=rbind) %dopar%
  once(s,P=3)
stopCluster(cl)

GD <- data.frame(matrix(unlist(GD), ncol=20))
colnames(GD) <- c("biasn1", "biasn2","biasn3", "biasor1","biasor2","biasor3", "biasada1", "biasada2","biasada3", "sdn1", "sdn2","sdn3", "sdor1", "sdor2", "sdor3", "sdada1", "sdada2","sdada3", "nInvada", "freqAllada")
biasn <- abs(GD$biasn1) + abs(GD$biasn2) + abs(GD$biasn3)
biasor <- abs(GD$biasor1) + abs(GD$biasor2) + abs(GD$biasor3)
biasada <- abs(GD$biasada1) + abs(GD$biasada2) + abs(GD$biasada3)
nInvada <- GD$nInvada
freqAllada <- GD$freqAllada

par(mfrow=c(1,3))
par(mar=c(5,5,1,1), xpd=F)
plot(1:17, biasn, type="p",pch=1, lwd=2, col="grey28", ylim=c(0,20), xlab="# of invalid IVs", ylab="Sum of MAE")#, + , "|",hat(beta[2]), "|"
points(1:17, biasor, pch=2, cex=0.85, col="grey28", lty=5, lwd=2)
points(1:17, biasada, pch=4, cex=0.85, col="black", lwd=1.8)
abline(v=3, lty=3)
par(mar=c(5,5,1,1), xpd=F)
plot(1:17, rep(0,17), type="p",pch=1, lwd=2, main ="Three regressors",  col="grey28",ylim=c(0,20), xlab="# of invalid IVs", ylab="# of IVs selected as invalid")#, + , "|",hat(beta[2]), "|"
points(1:17, 1:s, pch=2, cex=0.85, col="grey28", lty=5, lwd=2, xlab="# of invalid IVs", ylab="# of IVs selected as invalid")
points(1:17, nInvada, pch=4, cex=0.85, col="black", lty=5, lwd=1.8, xlab="IVs", ylab="# of IVs selected as invalid")
abline(v=3, lty=3)
plot(1:17, rep(0,17), type="p",pch=1, lwd=2, col="grey28", ylim=c(0,1), xlab="# of invalid IVs", ylab="Freq. with which all invalid selected as invalid")#, + , "|",hat(beta[2]), "|"
points(1:17, rep(1,17), pch=2, cex=0.85, col="grey28", lty=5, lwd=2, xlab="# of invalid IVs", ylab="# of IVs selected as invalid")
points(1:17, freqAllada, pch=4, cex=0.85, col="black", lty=5, lwd=1.8, xlab="# of invalid IVs", ylab="# of IVs selected as invalid")
abline(v=3, lty=3)

#End----