# Simulations in paper, for single regressor

rm(list=ls())
# Number of cores:
nr.cores <- detectCores() - 1

# Initiate cluster
cl <- makeCluster(nr.cores)

registerDoParallel(cl)

source("01Functions_repo.R")

M=500; s <- 3; L <- 10; alpha <- c(rep(0.2, s), rep(0, L-s))
res <- simulone(M=M, s=s, L=L, alpha=alpha, gammaval=0.6,title="Majority valid", ylimup=0.42, breaks=c(0, 0.1, 0.2, 0.3, 0.4), yticklab=c("      0","0.1", "0.2", "0.3", "0.4"), xlab="Strong IV, Weak violation", ylimF=23, breaksF=c(0,5,10,15,20), yticklabF=c("      0","5", "10", "15", "20"))
save.image("../temp/SE_maj_a02_g06.RData")
load("../temp/SE_maj_a02_g06.RData")
mnoac <- res$mnoac

mnoac$Method <- factor(mnoac$Method, levels = c("Oracle", "Standard", "AL", "CIM"))
class(mnoac$Method)

# Plotting
# Black-white
m1_main <- ggplot(mnoac, aes(x = ev, y = mad, color = Method, size = Method, linetype = Method)) +
  geom_path() +  # Use geom_path instead of geom_line
  labs(x = "ev", y = "mad")
m1_leg <- m1_main +  
  scale_color_manual(values = c("Standard" = "black", "Oracle" = "darkgray", "AL" = "black", "CIM" = "black")) +
  scale_size_manual(values = c("Standard" = 1, "Oracle" = 3, "AL" = 2, "CIM" = 1)) + 
  scale_linetype_manual(values = c("Standard" = "solid", "Oracle" = "solid", "AL" = "dotted", "CIM" = "dotdash")) +
  theme_classic() + 
  ggtitle("Majority valid") + 
  scale_y_continuous(limits=c(0,0.22), name ="Median Absolute Deviation", breaks=c(0,0.05,0.1,0.15,0.2), labels=c("0","0.05", "0.1", "0.15", "0.2")) +
  theme(plot.margin = unit(c(0,0,-8,0), "pt")) 
m1_leg
m1_main <- m1_main +  scale_color_manual(values = c("Standard" = "black", "Oracle" = "darkgray", "AL" = "black", "CIM" = "black")) +
  scale_size_manual(values = c("Standard" = 1, "Oracle" = 3, "AL" = 2, "CIM" = 1)) + 
  scale_linetype_manual(values = c("Standard" = "solid", "Oracle" = "solid", "AL" = "dotted", "CIM" = "dotdash")) +
  scale_x_continuous(name ="", breaks=c(400, 1200, 2000, 2800, 3600, 4400, 5200, 6000), labels=c("400", "1200", "2000", "2800", "3600", "4400", "5200", "6000")) +
  theme_classic() + ggtitle("Majority valid") + scale_y_continuous(limits=c(0,0.22), name ="Median Absolute Deviation", breaks=c(0,0.05,0.1,0.15,0.2), labels=c("      0","0.05", "0.1", "0.15", "0.2")) +
  theme(plot.margin = unit(c(0,0,-8,0), "pt"), legend.position = "none") 
m1_main

# Colour
m1_col_main <- ggplot(mnoac, aes(x = ev, y = mad, color = Method, size = Method, linetype = Method)) +
  geom_line() +
  labs(x = "ev", y = "mad")
m1_col_main <- m1_main +  scale_color_manual(values=c("orangered3","grey","gold2","dodgerblue3")) +
  scale_size_manual(values = c("Standard" = 1, "Oracle" = 3, "AL" = 2, "CIM" = 1)) + 
  scale_linetype_manual(values = c("Standard" = "solid", "Oracle" = "solid", "AL" = "dotted", "CIM" = "dotdash")) +
  scale_x_continuous(name ="", breaks=c(400, 1200, 2000, 2800, 3600, 4400, 5200, 6000), labels=c("400", "1200", "2000", "2800", "3600", "4400", "5200", "6000")) +
  theme_classic() + ggtitle("Majority valid") + scale_y_continuous(limits=c(0,0.22), name ="Median Absolute Deviation", breaks=c(0,0.05,0.1,0.15,0.2), labels=c("      0","0.05", "0.1", "0.15", "0.2")) +
  theme(plot.margin = unit(c(0,0,-8,0), "pt"), legend.position = "none") 
m1_col_main

# nInv
# Black-white
m1_nInv <- ggplot(mnoac, aes(x = ev, y = nInv, color = Method, size = Method, linetype = Method)) +
  geom_line() +
  labs(x = "ev", y = "nInv")
m1_nInv <- m1_nInv +  scale_color_manual(values = c("Standard" = "black", "Oracle" = "darkgray", "AL" = "black", "CIM" = "black")) +
  scale_size_manual(values = c("Standard" = 1, "Oracle" = 3, "AL" = 2, "CIM" = 1)) + 
  scale_linetype_manual(values = c("Standard" = "solid", "Oracle" = "solid", "AL" = "dotted", "CIM" = "dotdash")) +
  scale_x_continuous(name ="", breaks=c(400, 1200, 2000, 2800, 3600, 4400, 5200, 6000), labels=c("400", "1200", "2000", "2800", "3600", "4400", "5200", "6000")) +
  theme_classic() + scale_y_continuous(name="# invalid", limits=c(0,8.2), breaks=c(0, 2, 4, 6, 8), labels=c("      0","2","4", "6", "8")) +
  theme(plot.margin = unit(c(0,0,-8,0), "pt"), legend.position = "none") 
m1_nInv

# Color
m1_col_nInv <- ggplot(mnoac, aes(x = ev, y = nInv, color = Method, size = Method, linetype = Method)) +
  geom_line() +
  labs(x = "ev", y = "nInv")
m1_col_nInv <- m1_col_nInv +  scale_color_manual(values=c("orangered3","grey","gold2","dodgerblue3")) +
  scale_size_manual(values = c("Standard" = 1, "Oracle" = 3, "AL" = 2, "CIM" = 1)) + 
  scale_linetype_manual(values = c("Standard" = "solid", "Oracle" = "solid", "AL" = "dotted", "CIM" = "dotdash")) +
  scale_x_continuous(name ="", breaks=c(400, 1200, 2000, 2800, 3600, 4400, 5200, 6000), labels=c("400", "1200", "2000", "2800", "3600", "4400", "5200", "6000")) +
  theme_classic() + ggtitle("Majority valid") + scale_y_continuous(name="# invalid", limits=c(0,8.2), breaks=c(0, 2, 4, 6, 8), labels=c("      0","2","4", "6", "8")) +
  theme(plot.margin = unit(c(0,0,-8,0), "pt"), legend.position = "none") 
m1_col_nInv

# Frequency
m1_freq <- ggplot(mnoac, aes(x = ev, y = freq, color = Method, size = Method, linetype = Method)) +
  geom_line() +
  labs(x = "ev", y = "nInv")
m1_freq <- m1_freq +  scale_color_manual(values = c("Standard" = "black", "Oracle" = "darkgray", "AL" = "black", "CIM" = "black")) +
  scale_size_manual(values = c("Standard" = 1, "Oracle" = 3, "AL" = 2, "CIM" = 1)) + 
  scale_linetype_manual(values = c("Standard" = "solid", "Oracle" = "solid", "AL" = "dotted", "CIM" = "dotdash")) +
  scale_x_continuous(name ="", breaks=c(400, 1200, 2000, 2800, 3600, 4400, 5200, 6000), labels=c("400", "1200", "2000", "2800", "3600", "4400", "5200", "6000")) +
  theme_classic() + ylab("Frequency") +
  theme(plot.margin = unit(c(0,0,-8,0), "pt"), legend.position = "none") 
m1_freq

# Color
m1_col_freq <- ggplot(mnoac, aes(x = ev, y = freq, color = Method, size = Method, linetype = Method)) +
  geom_line() +
  labs(x = "ev", y = "nInv")
m1_col_freq <- m1_col_freq +  scale_color_manual(values=c("orangered3","grey","gold2","dodgerblue3")) +
  scale_size_manual(values = c("Standard" = 1, "Oracle" = 3, "AL" = 2, "CIM" = 1)) + 
  scale_linetype_manual(values = c("Standard" = "solid", "Oracle" = "solid", "AL" = "dotted", "CIM" = "dotdash")) +
  scale_x_continuous(name ="", breaks=c(400, 1200, 2000, 2800, 3600, 4400, 5200, 6000), labels=c("400", "1200", "2000", "2800", "3600", "4400", "5200", "6000")) +
  theme_classic() + ylab("Frequency") +
  theme(plot.margin = unit(c(0,0,-8,0), "pt"), legend.position = "none") 
m1_col_freq

# Coverage
# Black-white
m1_cover <- ggplot(mnoac, aes(x = ev, y = cover, color = Method, size = Method, linetype = Method)) +
  geom_line() +
  labs(x = "ev", y = "cover")
m1_cover <- m1_cover +  scale_color_manual(values = c("Standard" = "black", "Oracle" = "darkgray", "AL" = "black", "CIM" = "black")) +
  scale_size_manual(values = c("Standard" = 1, "Oracle" = 3, "AL" = 2, "CIM" = 1)) + 
  scale_linetype_manual(values = c("Standard" = "solid", "Oracle" = "solid", "AL" = "dotted", "CIM" = "dotdash")) +
  scale_x_continuous(name ="", breaks=c(400, 1200, 2000, 2800, 3600, 4400, 5200, 6000), labels=c("400", "1200", "2000", "2800", "3600", "4400", "5200", "6000")) +
  theme_classic() +scale_y_continuous(limits=c(0,1.1), "Coverage", breaks=c(0,0.2,0.4,0.6,0.8,0.95), labels=c("      0","0.2", "0.4", "0.6", "0.8", "0.95")) +
  theme(plot.margin = unit(c(0,0,-8,0), "pt"), legend.position = "none") 
m1_cover

# Colour
m1_col_cover <- ggplot(mnoac, aes(x = ev, y = cover, color = Method, size = Method, linetype = Method)) +
  geom_line() +
  labs(x = "ev", y = "cover")
m1_col_cover <- m1_col_cover +  scale_color_manual(values=c("orangered3","grey","gold2","dodgerblue3")) +
  scale_size_manual(values = c("Standard" = 1, "Oracle" = 3, "AL" = 2, "CIM" = 1)) + 
  scale_linetype_manual(values = c("Standard" = "solid", "Oracle" = "solid", "AL" = "dotted", "CIM" = "dotdash")) +
  scale_x_continuous(name ="", breaks=c(400, 1200, 2000, 2800, 3600, 4400, 5200, 6000), labels=c("400", "1200", "2000", "2800", "3600", "4400", "5200", "6000")) +
  theme_classic() + ggtitle("Majority valid") + scale_y_continuous(limits=c(0,1.1), "Coverage", breaks=c(0,0.2,0.4,0.6,0.8,0.95), labels=c("      0","0.2", "0.4", "0.6", "0.8", "0.95")) +
  theme(plot.margin = unit(c(0,0,-8,0), "pt"), legend.position = "none") 
m1_col_cover

mcol_main <- plot_grid( m1_main, m1_cover, m1_nInv, m1_freq, 
                        align = 'h', 
                        hjust = 0, 
                        nrow = 4)
mrow_main <- plot_grid( m1_main, m1_cover, m1_nInv, m1_freq,
                        align = 'h', 
                        hjust = 0, 
                        nrow = 1)
mrow_col <- plot_grid( m1_col_main, m1_col_cover, m1_col_nInv, m1_col_nInv,
                       align = 'h', 
                       hjust = 0, 
                       nrow = 1)

m1_tF_main <- data.frame(res$mnoac_F, seq(400,6000, by=400))
colnames(m1_tF_main) <- c("tF", "ev")
m1_tF_main <- ggplot(m1_tF_main, aes(x=ev, y=tF))
m1_tF_main <- m1_tF_main + geom_line(colour="black", size=1, linetype="solid") + theme_classic() + theme(plot.margin = unit(c(0,0,0,0), "pt"), legend.position = "none")
m1_tF_main <- m1_tF_main + scale_x_continuous(name ="", breaks=c(400, 1200, 2000, 2800, 3600, 4400, 5200, 6000), labels=c("400", "1200", "2000", "2800", "3600", "4400", "5200", "6000"))
m1_tF_main <- m1_tF_main + scale_y_continuous(limits=c(0,23), name ="F", breaks=c(0,5,10,15,20), labels=c("      0","5", "10", "15", "20"))
m1_tF_main
# Saved this manually via Export in RStudio

#####################
# Plurality setting #
#####################
# Number of cores:
nr.cores <- detectCores() - 1

# Initiate cluster
cl <- makeCluster(nr.cores)

registerDoParallel(cl)

M <- 500; s <- 6; L <- 10; alpha1 <- c(0.1, 0.1, 0.2, 0.2, 0.3, 0.3); alpha2 <- rep(0, L-s)
alpha <- c(alpha1, alpha2)
res_plu <- simulone(M=M, s=s, L=L, alpha=alpha, gammaval=0.6,title="Plurality valid", ylimup=1.33, breaks=c(0, 0.3, 0.6, 0.9, 1.2), yticklab=c("      0","0.3", "0.6", "0.9", "1.2"), xlab="Strong IV, Weak violation", ylimF=23, breaksF=c(0,5,10,15,20), yticklabF=c("      0","5", "10", "15","20"))
save.image("../temp/SE_plu_a02_g06.RData")
load("../temp/SE_plu_a02_g06.RData")
pnoac <- res_plu$mnoac

pnoac$Method <- factor(pnoac$Method, levels = c("Standard", "Oracle", "AL", "CIM"))
class(pnoac$Method)

# Black-white
p1_main <- ggplot(pnoac, aes(x = ev, y = mad, color = Method, size = Method, linetype = Method)) +
  geom_line() +
  labs(x = "ev", y = "mad")
p1_main <- p1_main +  scale_color_manual(values = c("Standard" = "black", "Oracle" = "darkgray", "AL" = "black", "CIM" = "black")) +
  scale_size_manual(values = c("Standard" = 1, "Oracle" = 3, "AL" = 2, "CIM" = 1)) + 
  scale_linetype_manual(values = c("Standard" = "solid", "Oracle" = "solid", "AL" = "dotted", "CIM" = "dotdash")) +
  scale_x_continuous(name ="", breaks=c(400, 1200, 2000, 2800, 3600, 4400, 5200, 6000), labels=c("400", "1200", "2000", "2800", "3600", "4400", "5200", "6000")) +
  theme_classic() + ggtitle("Plurality valid") + scale_y_continuous(limits=c(0,0.22), name ="Median Absolute Deviation", breaks=c(0,0.05,0.1,0.15,0.2), labels=c("      0","0.05", "0.1", "0.15", "0.2")) +
  theme(plot.margin = unit(c(0,0,-8,0), "pt"), legend.position = "none") 
p1_main

# Colour
p1_col_main <- ggplot(pnoac, aes(x = ev, y = mad, color = Method, size = Method, linetype = Method)) +
  geom_line() +
  labs(x = "ev", y = "mad")
p1_col_main <- p1_main +  scale_color_manual(values=c("orangered3","grey","gold2","dodgerblue3")) +
  scale_size_manual(values = c("Standard" = 1, "Oracle" = 3, "AL" = 2, "CIM" = 1)) + 
  scale_linetype_manual(values = c("Standard" = "solid", "Oracle" = "solid", "AL" = "dotted", "CIM" = "dotdash")) +
  scale_x_continuous(name ="", breaks=c(400, 1200, 2000, 2800, 3600, 4400, 5200, 6000), labels=c("400", "1200", "2000", "2800", "3600", "4400", "5200", "6000")) +
  theme_classic() + ggtitle("Majority valid") + scale_y_continuous(limits=c(0,0.22), name ="Median Absolute Deviation", breaks=c(0,0.05,0.1,0.15,0.2), labels=c("      0","0.05", "0.1", "0.15", "0.2")) +
  theme(plot.margin = unit(c(0,0,-8,0), "pt"), legend.position = "none") 
p1_col_main

# nInv
# Black-white
p1_nInv <- ggplot(pnoac, aes(x = ev, y = nInv, color = Method, size = Method, linetype = Method)) +
  geom_line() +
  labs(x = "ev", y = "nInv")
p1_nInv <- p1_nInv +  scale_color_manual(values = c("Standard" = "black", "Oracle" = "darkgray", "AL" = "black", "CIM" = "black")) +
  scale_size_manual(values = c("Standard" = 1, "Oracle" = 3, "AL" = 2, "CIM" = 1)) + 
  scale_linetype_manual(values = c("Standard" = "solid", "Oracle" = "solid", "AL" = "dotted", "CIM" = "dotdash")) +
  scale_x_continuous(name ="", breaks=c(400, 1200, 2000, 2800, 3600, 4400, 5200, 6000), labels=c("400", "1200", "2000", "2800", "3600", "4400", "5200", "6000")) +
  theme_classic() + scale_y_continuous(name="# invalid", limits=c(0,8.2), breaks=c(0, 2, 4, 6, 8), labels=c("      0","2","4", "6", "8")) +
  theme(plot.margin = unit(c(0,0,-8,0), "pt"), legend.position = "none") 
p1_nInv

# Color
p1_col_nInv <- ggplot(pnoac, aes(x = ev, y = nInv, color = Method, size = Method, linetype = Method)) +
  geom_line() +
  labs(x = "ev", y = "nInv")
p1_col_nInv <- p1_col_nInv +  scale_color_manual(values=c("orangered3","grey","gold2","dodgerblue3")) +
  scale_size_manual(values = c("Standard" = 1, "Oracle" = 3, "AL" = 2, "CIM" = 1)) + 
  scale_linetype_manual(values = c("Standard" = "solid", "Oracle" = "solid", "AL" = "dotted", "CIM" = "dotdash")) +
  scale_x_continuous(name ="", breaks=c(400, 1200, 2000, 2800, 3600, 4400, 5200, 6000), labels=c("400", "1200", "2000", "2800", "3600", "4400", "5200", "6000")) +
  theme_classic() + ggtitle("Plurality valid") + scale_y_continuous(name="# invalid", limits=c(0,8.2), breaks=c(0, 2, 4, 6, 8), labels=c("      0","2","4", "6", "8")) +
  theme(plot.margin = unit(c(0,0,-8,0), "pt"), legend.position = "none") 
p1_col_nInv

# Frequency
p1_freq <- ggplot(pnoac, aes(x = ev, y = freq, color = Method, size = Method, linetype = Method)) +
  geom_line() +
  labs(x = "ev", y = "nInv")
p1_freq <- p1_freq +  scale_color_manual(values = c("Standard" = "black", "Oracle" = "darkgray", "AL" = "black", "CIM" = "black")) +
  scale_size_manual(values = c("Standard" = 1, "Oracle" = 3, "AL" = 2, "CIM" = 1)) + 
  scale_linetype_manual(values = c("Standard" = "solid", "Oracle" = "solid", "AL" = "dotted", "CIM" = "dotdash")) +
  scale_x_continuous(name ="", breaks=c(400, 1200, 2000, 2800, 3600, 4400, 5200, 6000), labels=c("400", "1200", "2000", "2800", "3600", "4400", "5200", "6000")) +
  theme_classic() + ylab("Frequency") +
  theme(plot.margin = unit(c(0,0,-8,0), "pt"), legend.position = "none") 
p1_freq

# Color
p1_col_freq <- ggplot(pnoac, aes(x = ev, y = freq, color = Method, size = Method, linetype = Method)) +
  geom_line() +
  labs(x = "ev", y = "nInv")
p1_col_freq <- p1_col_freq +  scale_color_manual(values=c("orangered3","grey","gold2","dodgerblue3")) +
  scale_size_manual(values = c("Standard" = 1, "Oracle" = 3, "AL" = 2, "CIM" = 1)) + 
  scale_linetype_manual(values = c("Standard" = "solid", "Oracle" = "solid", "AL" = "dotted", "CIM" = "dotdash")) +
  scale_x_continuous(name ="", breaks=c(400, 1200, 2000, 2800, 3600, 4400, 5200, 6000), labels=c("400", "1200", "2000", "2800", "3600", "4400", "5200", "6000")) +
  theme_classic() + ylab("Frequency") +
  theme(plot.margin = unit(c(0,0,-8,0), "pt"), legend.position = "none") 
p1_col_freq

# Coverage
# Black-white
p1_cover <- ggplot(pnoac, aes(x = ev, y = cover, color = Method, size = Method, linetype = Method)) +
  geom_line() +
  labs(x = "ev", y = "cover")
p1_cover <- p1_cover +  scale_color_manual(values = c("Standard" = "black", "Oracle" = "darkgray", "AL" = "black", "CIM" = "black")) +
  scale_size_manual(values = c("Standard" = 1, "Oracle" = 3, "AL" = 2, "CIM" = 1)) + 
  scale_linetype_manual(values = c("Standard" = "solid", "Oracle" = "solid", "AL" = "dotted", "CIM" = "dotdash")) +
  scale_x_continuous(name ="", breaks=c(400, 1200, 2000, 2800, 3600, 4400, 5200, 6000), labels=c("400", "1200", "2000", "2800", "3600", "4400", "5200", "6000")) +
  theme_classic() +scale_y_continuous(limits=c(0,1.1), "Coverage", breaks=c(0,0.2,0.4,0.6,0.8,0.95), labels=c("      0","0.2", "0.4", "0.6", "0.8", "0.95")) +
  theme(plot.margin = unit(c(0,0,-8,0), "pt"), legend.position = "none") 
p1_cover

# Colour
p1_col_cover <- ggplot(pnoac, aes(x = ev, y = cover, color = Method, size = Method, linetype = Method)) +
  geom_line() +
  labs(x = "ev", y = "cover")
p1_col_cover <- p1_col_cover +  scale_color_manual(values=c("orangered3","grey","gold2","dodgerblue3")) +
  scale_size_manual(values = c("Standard" = 1, "Oracle" = 3, "AL" = 2, "CIM" = 1)) + 
  scale_linetype_manual(values = c("Standard" = "solid", "Oracle" = "solid", "AL" = "dotted", "CIM" = "dotdash")) +
  scale_x_continuous(name ="", breaks=c(400, 1200, 2000, 2800, 3600, 4400, 5200, 6000), labels=c("400", "1200", "2000", "2800", "3600", "4400", "5200", "6000")) +
  theme_classic() + ggtitle("Plurality valid") + scale_y_continuous(limits=c(0,1.1), "Coverage", breaks=c(0,0.2,0.4,0.6,0.8,0.95), labels=c("      0","0.2", "0.4", "0.6", "0.8", "0.95")) +
  theme(plot.margin = unit(c(0,0,-8,0), "pt"), legend.position = "none") 
p1_col_cover

pcol_main <- plot_grid( p1_main, p1_cover, p1_nInv, p1_freq, 
                        align = 'h', 
                        hjust = 0, 
                        nrow = 4)
prow_main <- plot_grid( p1_main, p1_cover, p1_nInv, p1_freq,
                        align = 'h', 
                        hjust = 0, 
                        nrow = 1)
prow_col <- plot_grid( p1_col_main, p1_col_cover, p1_col_nInv, p1_col_freq,
                       align = 'h', 
                       hjust = 0, 
                       nrow = 1)

# Pdfs

# Main graph- column
legend_m <- get_legend(m1_leg + theme(legend.margin=margin(t=-4, r=0, b=0, l=0, unit="cm"), legend.position="bottom"))
p <- plot_grid(mcol_main, pcol_main, nrow=1, ncol = 2, rel_heights = c(1, .2))
p <- plot_grid(p, legend_m, nrow=2, ncol = 1, rel_heights = c(1, .2))
pdf("../gfx/SS_100_main.pdf", height=12, width=8)
p
dev.off()

# Main graph- row
p <- plot_grid(mrow_main, prow_main, nrow=2, ncol = 1, rel_heights = c(1, 1))
p <- plot_grid(p, legend_m, nrow=2, ncol = 1, rel_heights = c(1, .2))
pdf("../gfx/SS_100_row.pdf", height=8, width=14)
p
dev.off()

