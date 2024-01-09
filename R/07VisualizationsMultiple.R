# Upper bound ---- 
#Calculates the upper bounds mentioned for the multiple regressor case 
# for several parameter combinations

rm(list=ls())

K <- 2
G <- 30
L <- 50
S <- L - G
R <- 40

phyper(1,G,S,K, lower.tail = F)^R

dhyper(0,G,S,K)
dhyper(1,G,S,K)
dhyper(2,G,S,K)
dhyper(4,G,S,K)

#q=1 at least two drawn from largest group
#m cardinality of largest group. Vary this. Should be larger than 
#A - m = k
#k number of centers we end up with

# Limit ----

rm(list=ls())

limitG <- function(CandidateIVn){
  gLim <- (1 + sqrt(1 + 4*0.5*(CandidateIVn-1)*CandidateIVn))/(2*CandidateIVn)
}

Cands <- matrix(1:10000, nrow=10000)
minG <- apply(Cands,MARGIN=1,FUN=limitG)
plot(Cands, minG, type="l", ylim=c(0.7,0.71))
abline(h=0.7071, col="red")

# Fraction valid ----

rm(list=ls())

L <- 100

validFrac <- function(g){choose(g, k)/choose(L,k)}

par(mfrow=c(1,3))
par(mar=c(5, 4, 2, 1))
k <- 1
x <- seq(-100,0,by=1)
vf <- validFrac(x)
plot(x=x, y=-vf, xlab="", main="P=1", ylab="Fraction of models using valid IVs", xaxt="n",yaxt="n", type="l")
axis(side = 2, at = c(0,0.5,0.9))
axis(side = 1, at = c(-100, -75, -50, 0), labels=c(0,25,50,100))
abline(h=0.5)
abline(v=-50, lty=2)

k <- 2
x <- seq(0,100,by=0.1)
vf <- validFrac(x)
plot(x=-x, y = vf, main="P=2", xlab="Number of invalid IVs", ylab="", yaxt  = "n", xaxt="n", type="l")
axis(side = 2, at = c(0,0.5,0.9))
axis(side = 1, at = c(-100, -71, -50, 0), labels=c(0,29,50,100))
abline(h=0.5)
abline(v=-71, lty=2)

k <- 3
vf <- validFrac(x)
plot(x=-x, y = vf, main="P=3", ylab="", xlab="", yaxt  = "n", xaxt="n", type="l")
axis(side = 2, at = c(0,0.5,0.9))
axis(side = 1, at = c(-100, -80, -50, 0), labels=c(0,20,50,100))
abline(h=0.5)
abline(v=-80, lty=2)

L <- 20
k <- 2
g <- 14

choose(g,k)/choose(L,k)
