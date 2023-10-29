library(plotrix);
library(tidyverse)
library(scales) # to access break formatting functions
library(latex2exp)

muc <- 1e-10
mus <- 3.64e-9

divB <- 10^seq(-10,0,0.01)
divV <- (divB[-1]*divB[-length(divB)])^0.5

#tau  <- 1.93e7; pair <-  "Enterobacter_asburiae_vs_Enterobacter_roggenkampii"; nfiles <- 1
#tau  <- 1.49e8; pair <- "Citrobacter_vs_Enterobacter_roggenkampii"; nfiles <- 4
#tau  <- 4.01e7; pair <- "Enterobacter_hormaechei_vs_Enterobacter_roggenkampii"; nfiles <- 7
#tau  <- 1.03e8/3; pair <- "Enterobacter_asburiae_vs_Enterobacter_hormaechei"; nfiles <- 7
tau  <- 1.00e8; pair <- "Escherichia_coli_vs_Salmonella_enterica"; nfiles <- 10


counts <- 0
for (i in 1:nfiles)
{
  counts <- counts + read.table(file = paste0("/home/misha/Downloads/segmut_",i,".L") ,header=FALSE, row.names=NULL,sep="\t",skip = 0)$V1
}

divBnew <- 10^seq(-10,0,0.1)
divVnew <- (divBnew[-1]*divBnew[-length(divBnew)])^0.5
counts <- weighted.hist(divV,counts,breaks=divBnew,plot=FALSE)$counts/diff(divBnew)/sum(counts)

# plot(log10(divVnew[counts>0]/tau),log10(counts[counts>0]))
# lines(log10(divVnew[counts>0]/tau),17+log10(divVnew[counts>0]/tau))


mindens <- 3e-3

mu <- divVnew[counts>=mindens]/tau
dens <- counts[counts>=mindens]

densLin <- dens[mu>=muc & mu<=mus]
muLin <- mu[mu>=muc & mu<=mus]

fit <- lm(log10(densLin) - log10(muLin) ~ 1)

densmax <- 10^fit$coefficients*mus
densmin <- 10^fit$coefficients*muc


# dev.new(width=4, height=4)
pdf(paste0("/home/misha/Documents/Development/HAT/plots/mu_dist/",pair,".pdf"),width=4, height=4)
ggplot(
  data=data.frame(divVnew=mu,dens=dens), aes(x=divVnew, y=dens))  + 
  geom_point(color="black",shape=1) + 
  scale_x_log10(TeX("$\\mu_i$"), breaks = 10^seq(-12,0,1),
                labels = trans_format("log10", math_format(10^.x)), 
                sec.axis = sec_axis(~ . *tau, name = TeX("$\\theta_i$"),breaks = 10^seq(-12,0,1),
                                    labels = trans_format("log10", math_format(10^.x)))) +
  scale_y_log10("probability density",breaks = 10^seq(-6,9,1),
                labels = trans_format("log10", math_format(10^.x))) +
  theme(panel.border = element_rect(color = "black", linewidth = 0.5, fill=NA),aspect.ratio=1,axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),text = element_text(size = 20)) +
  geom_line(data=data.frame(x = rep(mus,2), y = c(mindens, densmax)), aes(x=x, y=y), color="blue", linetype="solid", size=0.5) +
  geom_line(data=data.frame(x = c(muc,muc), y = c(mindens, densmin)), aes(x=x, y=y), color="blue", linetype="solid", size=0.5) +
  geom_line(data=data.frame(x = c(muc,muc,mus,mus), y = c(densmin,densmin,densmax, densmax)), aes(x=x, y=y), color="blue", linetype="solid", size=0.5)
dev.off()











