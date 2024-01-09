# Toy illustrations

# Illustrations for CIM in paper
par(mar=c(2,4,1,1))
m <- matrix(c(1:4, 1,1.12,1.15,1.4,1.1,1.2,1.28,1.5,1.2,1.3,1.45,1.6), nrow=4)
plot(m[,1], m[,3], xlim=c(1,4), ylim=c(1,1.65), pch=19, ylab = "CIs", xaxt="n",xlab = "")
segments(x0=1:4, x1=1:4, y0=m[,4], y1=m[,2], lwd=3)

par(mar=c(2,4,1,1))
m <- matrix(c(1:4, 1,1.12,1.15,1.4,1.1,1.2,1.28,1.5,1.2,1.3,1.45,1.6), nrow=4)
plot(m[,1], m[,3], xlim=c(1,4), ylim=c(1,1.65), pch=19, ylab = "CIs", xaxt="n",xlab = "")
segments(x0=1:4, x1=1:4, y0=m[,4], y1=m[,2], lwd=3)
segments(x0=2,x1=0,y0=1.12, y1=1.12, col="red", lwd=3)
text("2", x=1.5,y=1.2)

m <- matrix(c(1:4, 1,1.12,1.15,1.4,1.1,1.2,1.28,1.5,1.2,1.3,1.45,1.6), nrow=4)
plot(m[,1], m[,3], xlim=c(1,4), ylim=c(1,1.65), pch=19, ylab = "CIs", xaxt="n",xlab = "")
segments(x0=1:4, x1=1:4, y0=m[,4], y1=m[,2], lwd=3)
segments(x0=3,x1=0,y0=1.15, y1=1.15, col="red", lwd=3)
text("3", x=2.5,y=1.2)

m <- matrix(c(1:4, 1,1.12,1.15,1.4,1.1,1.2,1.28,1.5,1.2,1.3,1.45,1.6), nrow=4)
plot(m[,1], m[,3], xlim=c(1,4), ylim=c(1,1.65), pch=19, ylab = "CIs", xaxt="n",xlab = "")
segments(x0=1:4, x1=1:4, y0=m[,4], y1=m[,2], lwd=3)
segments(x0=4,x1=0,y0=1.4, y1=1.4, col="red", lwd=3)
text("2", x=3.5,y=1.45)

par(mar=c(2,4,1,1))
m <- cbind(1:5,c(-0.12,-0.05,0.1,0.4,0.6), c(-0.3,-0.15,-0.05,0.35,0.55), c(0.06,0.05,0.25,0.45,0.65))
plot(m[,1], m[,2], xlim=c(1,5), ylim=c(-0.28,0.8), pch=19, ylab = "CIs", xaxt="n",xlab = "", frame.plot = F)
axis(1,at=1:5, labels=c("A", "B","C", "D", "E"))
abline(h=0)
segments(x0=1:5, x1=1:5, y0=m[,3], y1=m[,4], lwd=3)
segments(x0=2,x1=0,y0=-0.15, y1=-0.15, col="grey46", lwd=2, lty=2)
text("2", x=1.5,y=-0.2)

plot(m[,1], m[,2], xlim=c(1,5), ylim=c(-0.3,0.8), pch=19, ylab="", yaxt="n", xaxt="n",xlab = "", frame.plot = F)
axis(1,at=1:5, labels=c("A", "B","C", "D", "E"))
abline(h=0)
segments(x0=1:5, x1=1:5, y0=m[,3], y1=m[,4], lwd=3)
segments(x0=3,x1=0,y0=-0.05, y1=-0.05, col="grey46", lwd=2, lty=2)
text("3", x=2.5,y=-0.1)

plot(m[,1], m[,2], xlim=c(1,5), ylim=c(-0.3,0.8), pch=19, ylab="", yaxt="n", xaxt="n",xlab = "", frame.plot = F)
axis(1,at=1:5, labels=c("A", "B","C", "D", "E"))
abline(h=0)
segments(x0=1:5, x1=1:5, y0=m[,3], y1=m[,4], lwd=3)
segments(x0=4,x1=0,y0=0.35, y1=0.35, col="grey46", lwd=2, lty=2)
text("1", x=3.5,y=0.3)

plot(m[,1], m[,2], xlim=c(1,5), ylim=c(-0.3,0.8), pch=19, ylab="", yaxt="n", xaxt="n",xlab = "", frame.plot = F)
axis(1,at=1:5, labels=c("A", "B","C", "D", "E"))
abline(h=0)
segments(x0=1:5, x1=1:5, y0=m[,3], y1=m[,4], lwd=3)
segments(x0=5,x1=0,y0=0.55, y1=0.55, col="grey46", lwd=2, lty=2)
text("1", x=4.5,y=0.5)

# All in one plot
par(mfrow=c(1,4))
par(mar=c(2,2.2,1,1))
m <- cbind(1:5,c(-0.12,-0.05,0.1,0.4,0.6), c(-0.3,-0.15,-0.05,0.35,0.55), c(0.06,0.05,0.25,0.45,0.65))
plot(m[,1], m[,2], xlim=c(1,5), ylim=c(-0.28,0.8), pch=19, ylab = "CIs", xaxt="n",xlab = "", frame.plot = F)
axis(1,at=1:5, labels=c("A", "B","C", "D", "E"))
abline(h=0)
segments(x0=1:5, x1=1:5, y0=m[,3], y1=m[,4], lwd=3)
segments(x0=2,x1=0,y0=-0.15, y1=-0.15, col="grey46", lwd=2, lty=2)
text("2", x=1.5,y=-0.2)
text("a)", x=1.4, y=0.8, cex=2)

plot(m[,1], m[,2], xlim=c(1,5), ylim=c(-0.28,0.8), pch=19, ylab="", yaxt="n", xaxt="n",xlab = "", frame.plot = F)
axis(1,at=1:5, labels=c("A", "B","C", "D", "E"))
abline(h=0)
segments(x0=1:5, x1=1:5, y0=m[,3], y1=m[,4], lwd=3)
segments(x0=3,x1=0,y0=-0.05, y1=-0.05, col="grey46", lwd=2, lty=2)
text("3", x=2.5,y=-0.1)
text("b)", x=1.4, y=0.8, cex=2)

plot(m[,1], m[,2], xlim=c(1,5), ylim=c(-0.28,0.8), pch=19, ylab="", yaxt="n", xaxt="n",xlab = "", frame.plot = F)
axis(1,at=1:5, labels=c("A", "B","C", "D", "E"))
abline(h=0)
segments(x0=1:5, x1=1:5, y0=m[,3], y1=m[,4], lwd=3)
segments(x0=4,x1=0,y0=0.35, y1=0.35, col="grey46", lwd=2, lty=2)
text("1", x=3.5,y=0.3)
text("c)", x=1.4, y=0.8, cex=2)

plot(m[,1], m[,2], xlim=c(1,5), ylim=c(-0.28,0.8), pch=19, ylab="", yaxt="n", xaxt="n",xlab = "", frame.plot = F)
axis(1,at=1:5, labels=c("A", "B","C", "D", "E"))
abline(h=0)
segments(x0=1:5, x1=1:5, y0=m[,3], y1=m[,4], lwd=3)
segments(x0=5,x1=0,y0=0.55, y1=0.55, col="grey46", lwd=2, lty=2)
text("1", x=4.5,y=0.5)
text("d)", x=1.4, y=0.8, cex=2)

# ALasso Illustration
m_true <- cbind(c(0,0,0,0.5,0.5), c(0.1,0.11,0.12,0.1,0.11))
m_est <- cbind(c(-0.12,-0.05,0.05,0.4,0.6), c(0.1, 0.1, 0.1, 0.1, 0.1))

par(mfrow=c(1,1))
par(mar=c(4,0,0,0))
plot(m_est, ylim=c(0.1,0.11), xlim = c(-0.3,0.8),col="darkgrey", yaxt="n", xaxt="n", pch=16, cex=2, type="p", xlab=expression(beta), ylab="", frame.plot = F)
points(m_est, col="black", pch=1, cex=2, lwd=2)
lines(cbind(c(0,0),c(0,0.104)), lty=2)

text(expression(hat(beta)[A]), x=-0.09,y=0.1001)
text(expression(hat(beta)[B]), x=-0.02,y=0.1001)
text(expression(hat(beta)[C]), x=0.08,y=0.1001)

text(expression(hat(beta)[D]), x=0.43,y=0.1001)
text(expression(hat(beta)[E]), x=0.63,y=0.1001)

text(bquote({beta[A]==beta[B]}==beta[C]), x=0,y=0.105)

axis(1,at=c(-0.2, 0,0.05,0.5, 0.8), labels=c(-0.2, bquote(beta[0]==0), expression(~~~~~~~~~~~~~hat(beta)[m] == 0.05), expression(beta[0] + c == 0.5), 0.8))
