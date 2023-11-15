# install.packages("/home/msheinman/Development/HAT/src/SimulateHAT", repos = NULL, type = "source",force=TRUE)
rm(list = ls())
# load libraries----
library(SimulateHAT)
library(RcppArmadillo)
library(Rcpp)
library(data.table)
library(igraph)
library(ggpubr)
require(MASS)
library(matrixStats)
library(stringi)
library(tidyverse)
library(plotrix)
library(DistributionUtils)
library(scales)
library(latex2exp)
library(ggExtra)
library(gridExtra)
library(RColorBrewer) #need colors to make heatmaps
library(DescTools)
library(nloptr)
library(parallel)
library(foreach)
library(iterators)
library(doParallel)
library(EnvStats)
library(metaheuristicOpt)
library(expint)
library(phytools)
library(castor)
library(ComplexHeatmap)
library(circlize)
library(expint)
library(lme4)


rB <- c(-0.5,2.5,seq(5.5,35.5,3),10^seq(log10(35.5)+0.1,5.5,0.1))
rV <- 0.5*(rB[-length(rB)]+rB[-1])
species <- c("Escherichia_coli",
             "Klebsiella_pneumoniae",
             "Salmonella_enterica",
             "Enterobacter_asburiae",
             "Escherichia_fergusonii",
             "Escherichia_albertii",
             "Raoultella",
             "Citrobacter",
             "Cronobacter",
             "Serratia",
             "Enterobacter_roggenkampii",
             "Enterobacter_hormaechei",
             "Vibrio")

# species <- c(
#   "Simulation1e7",
#   "Simulation2e7",
#   "Simulation3e7",
#   "Simulation4e7",
#   "Simulation5e7",
#   "Simulation6e7",
#   "Simulation7e7",
#   "Simulation8e7",
#   "Simulation9e7",
#   "Simulation1e8",
#   "Simulation2e8",
#   "Simulation3e8",
#   "Simulation4e8",
#   "Simulation5e8",
#   "Simulation6e8",
#   "Simulation7e8",
#   "Simulation8e8",
#   "Simulation9e8"
#   )
# species <- c(species,paste0(species,"_init"))


# species <- c("Saccharomyces_cerevisiae","Saccharomyces_pastorianus")

ncores <- 20
registerDoParallel(20)
pairV <- c()
tauV <- c()
K1V <- c()
rho1V <- c()
rho2V <- c()
s1 <- 1; s2 <- 2

dat <- foreach (s1 = 1:(length(species)-1), .inorder=FALSE, .combine=rbind) %dopar%
{
  foreach (s2 = (s1+1):length(species), .inorder=FALSE, .combine=rbind) %dopar%
  {
  dat <- data.frame()
  species1 <- species[s1]
  species2 <- species[s2]
  files <- sort(Sys.glob(file.path(paste0("/scratch/ws1/msheinman-msheinman/nucmer/",species1,"_vs_",species2,"/data_*.L"))))
  filePairs <- sort(Sys.glob(file.path(paste0("/scratch/ws1/msheinman-msheinman/nucmer/",species1,"_vs_",species2,"/listOfPairs.csv"))))
  if (length(files)==0 & s1!=s2)
  {
    species2 <- species[s1]
    species1 <- species[s2]
    files <- sort(Sys.glob(file.path(paste0("/scratch/ws1/msheinman-msheinman/nucmer/",species1,"_vs_",species2,"/data_*.L"))))
    filePairs <- sort(Sys.glob(file.path(paste0("/scratch/ws1/msheinman-msheinman/nucmer/",species1,"_vs_",species2,"/listOfPairs.csv"))))
    if (length(files)==0)
    {
      print(paste0(species1,"_vs_",species2," NOT FOUND!"))
    }
  }
  if (length(files)>0)
  {
    print(paste0(species1,"_vs_",species2))

    Lall <- foreach (f = files, .combine=rbind, .inorder=TRUE) %do% {read.table(file = f ,header=TRUE, row.names=NULL,sep="\t",skip = 0)}
    if (nrow(Lall)>0)
    {
      dataPairs <- read.table(file = filePairs ,header=FALSE, row.names=NULL,sep="\t",skip = 0)
      Lall <- Lall[!is.na(Lall[,7]),]
      nE <- colSums(Lall[,7:ncol(Lall)])
      nE <- nE/nrow(dataPairs)
      mE <- nE/diff(rB)
      dat <- data.frame(species1=species1,species2=species2,name=paste0(species1,"_vs_",species2),
                        # div=sum(Lall$tau*Lall$L)/sum(Lall$L),
                        # div=mean(Lall$tau),
                        div=sum(Lall[,7:ncol(Lall)])/sum(Lall$L),
                        divSD=sd(Lall$tau),
                        L=mean(Lall$L),L1=mean(dataPairs[,5]),L2=mean(dataPairs[,6]),
                        Lfrac=mean(Lall$L)/(mean(dataPairs[,5])+mean(dataPairs[,5]))*2, 
                        LfracSD=sd(Lall$L/rowMins(cbind(dataPairs[,6],dataPairs[,5]))),
                        n1=length(unique(Lall[,3])),n2=length(unique(Lall[,4])))
      dat <- cbind(dat,t(mE))
      dat$L0 <- mean(rowMins(as.matrix(dataPairs[,c(5,6)])))
    }
    else print(paste0(species1,"_vs_",species2," EMPTY!"))

  }
  dat
  }
}
nPairs <- nrow(dat)

dat <- dat[order(dat$div),]
rTh <- rep(1e6,nPairs)

# dat <- dat[!dat$name=="Enterobacter_asburiae_vs_Cronobacter",]; nPairs <- nrow(dat)
# dat <- dat[!dat$name=="Cronobacter_vs_Enterobacter_roggenkampii",]; nPairs <- nrow(dat)
# dat <- dat[!dat$name=="Cronobacter_vs_Enterobacter_hormaechei",]; nPairs <- nrow(dat)


# add times from TimeTree
tree <- read.newick("/home/msheinman/Development/HAT/data/external/4052_Marin_2016.nwk")
for (i in 1:nrow(dat))
{
  species1 <- dat$species1[i]
  species2 <- dat$species2[i]
  Ind1 <- which(str_detect(tree$tip.label,species1))
  if (length(Ind1)==0)
  {
    print(species1)
  }
  Ind2 <- which(str_detect(tree$tip.label,species2))
  if (length(Ind2)==0)
  {
    print(species2)
  }
  if (length(Ind1)>0 & length(Ind2)>0)
    dat$TimeTree[i] <- get_pairwise_distances(tree, tree$tip.label[Ind1[1]], tree$tip.label[Ind2[1]] , as_edge_counts=FALSE, check_input=TRUE)
  else
    dat$TimeTree[i] <- 0
}

r <- rV
dr <- 0.1
rm <- r - dr
rp <- r + dr


a <- 1
d <- 0.25
mus <- 8.9e-11*200*941000/4.6e6 #  Mutation Rate Inferred From Synonymous Substitutions in a Long-Term Evolution Experiment With Escherichia coli
muc <- 1.0e-2/100e6 #Calibrating bacterial evolution
# muc <- mus/180 #Calibrating bacterial evolution

mus/muc
par <- log10(c(muc))

registerDoParallel(nPairs)

# global_fit <- function(par)
{
  muc <- 10^par[1]
  print(10^par)

  results <- foreach (i = 1:nPairs) %dopar%
  {
    mE <- as.numeric(dat[i,13+(1:length(rV))-1])
    L0 <- dat$L0[i]
    LLlocal <- function(par1)
    {
      tau <- as.numeric(10^par1[1]);
      rho <- as.numeric(10^par1[2])
      mua <- min(d/tau,mus)

      # mc <- (exp(-(mua/mus) - r*(mua + muc)*tau)*(-(exp(muc*(1/mus + r*tau))*(mua + mus + r*mua*mus*tau)) + exp(mua*(1/mus + r*tau))*(muc + mus + r*muc*mus*tau)))/((muc + mus)*(1 + r*mus*tau)^2)
      # mcm <- (exp(-(mua/mus) - rm*(mua + muc)*tau)*(-(exp(muc*(1/mus + rm*tau))*(mua + mus + rm*mua*mus*tau)) + exp(mua*(1/mus + rm*tau))*(muc + mus + rm*muc*mus*tau)))/((muc + mus)*(1 + rm*mus*tau)^2)
      # mcp <- (exp(-(mua/mus) - rp*(mua + muc)*tau)*(-(exp(muc*(1/mus + rp*tau))*(mua + mus + rp*mua*mus*tau)) + exp(mua*(1/mus + rp*tau))*(muc + mus + rp*muc*mus*tau)))/((muc + mus)*(1 + rp*mus*tau)^2)

      mc <- 2*((1 + r*mua*tau)/exp(r*mua*tau) - (1 + r*muc*tau)/exp(r*muc*tau))/(r^2*(muc^2 - mus^2)*tau^2)
      mcm <- 2*((1 + rm*mua*tau)/exp(rm*mua*tau) - (1 + rm*muc*tau)/exp(rm*muc*tau))/(rm^2*(muc^2 - mus^2)*tau^2)
      mcp <- 2*((1 + rp*mua*tau)/exp(rp*mua*tau) - (1 + rp*muc*tau)/exp(rp*muc*tau))/(rp^2*(muc^2 - mus^2)*tau^2)


      # mc <- -(((1 + a)*(r*tau)^(-1 - a)*((-exp(-(r*mua*tau)) + exp(-(r*muc*tau)))*muc^a*(r*tau)^a + gammainc(1 + a,r*mua*tau) - gammainc(1 + a,r*muc*tau)))/(a*muc^a*(muc - mus) + mus*(-muc^a + mus^a)))
      # mcm <- -(((1 + a)*(rm*tau)^(-1 - a)*((-exp(-(rm*mua*tau)) + exp(-(rm*muc*tau)))*muc^a*(rm*tau)^a + gammainc(1 + a,rm*mua*tau) - gammainc(1 + a,rm*muc*tau)))/(a*muc^a*(muc - mus) + mus*(-muc^a + mus^a)))
      # mcp <- -(((1 + a)*(rp*tau)^(-1 - a)*((-exp(-(rp*mua*tau)) + exp(-(rp*muc*tau)))*muc^a*(rp*tau)^a + gammainc(1 + a,rp*mua*tau) - gammainc(1 + a,rp*muc*tau)))/(a*muc^a*(muc - mus) + mus*(-muc^a + mus^a)))

      mc <- L0*(mcm + mcp - 2*mc)/dr^2
      mc[is.na(mc)] <- 0

      if (tau<d/mus)
      {
        mh <- (2*(-exp(-(r*muc*tau)) + exp(-(r*mus*tau)) + r*(-muc + mus)*tau))/(r^2*(-muc^2 + mus^2)*tau)
        mhm <- (2*(-exp(-(rm*muc*tau)) + exp(-(rm*mus*tau)) + rm*(-muc + mus)*tau))/(rm^2*(-muc^2 + mus^2)*tau)
        mhp <- (2*(-exp(-(rp*muc*tau)) + exp(-(rp*mus*tau)) + rp*(-muc + mus)*tau))/(rp^2*(-muc^2 + mus^2)*tau)
      }      else
      {
        mh <- (-2*(-exp(-(r*muc*tau)) + r*(-muc + mus)*tau + (1 + r*(d - mus*tau))/exp(r*d)))/(r^2*(muc^2 - mus^2)*tau)
        mhm <- (-2*(-exp(-(rm*muc*tau)) + rm*(-muc + mus)*tau + (1 + rm*(d - mus*tau))/exp(rm*d)))/(rm^2*(muc^2 - mus^2)*tau)
        mhp <- (-2*(-exp(-(rp*muc*tau)) + rp*(-muc + mus)*tau + (1 + rp*(d - mus*tau))/exp(rp*d)))/(rp^2*(muc^2 - mus^2)*tau)
      }
      mh <- L0*rho*(mhm + mhp - 2*mh)/dr^2

      mT <- mc + mh
      # mT <- mT/sum(mT*rV*diff(rB))*sum(mE*rV*diff(rB))
      # mT <- mT/max(mT)*max(mE)
      return(mean( ((mT-mE)/(mE + mT))^2))
    }

    res <- metaOpt(LLlocal, optimType = "MIN", algorithm = c("HS"), 2,
                   rangeVar=t(as.matrix(data.frame(lower=log10(c(5e6,1e-14)),
                                                   upper = log10(c(8e9,1.01e-8))))),
                   control = list(numPopulation=10000,maxIter=10), seed = 3)
    res$par <- res$result

    res <- Nelder_Mead(par=res$par,
                lower=log10(c(1e6,1e-16)),
                upper=log10(c(1e10,1.01e-7)),
                fn=LLlocal); res$value <- res$fval

    goodnessFits <- res$value

    tauFit <- 10^res$par[1]
    tau <- tauFit # !!!!!!!!
    rho <- 10^res$par[2]
    div <- dat[i,"div"]
    species1 <- dat[i,"species1"]
    species2 <- dat[i,"species2"]

    

    mua <- min(d/tau,mus)


    # mc <- (exp(-(mua/mus) - r*(mua + muc)*tau)*(-(exp(muc*(1/mus + r*tau))*(mua + mus + r*mua*mus*tau)) + exp(mua*(1/mus + r*tau))*(muc + mus + r*muc*mus*tau)))/((muc + mus)*(1 + r*mus*tau)^2)
    # mcm <- (exp(-(mua/mus) - rm*(mua + muc)*tau)*(-(exp(muc*(1/mus + rm*tau))*(mua + mus + rm*mua*mus*tau)) + exp(mua*(1/mus + rm*tau))*(muc + mus + rm*muc*mus*tau)))/((muc + mus)*(1 + rm*mus*tau)^2)
    # mcp <- (exp(-(mua/mus) - rp*(mua + muc)*tau)*(-(exp(muc*(1/mus + rp*tau))*(mua + mus + rp*mua*mus*tau)) + exp(mua*(1/mus + rp*tau))*(muc + mus + rp*muc*mus*tau)))/((muc + mus)*(1 + rp*mus*tau)^2)

    mc <- 2*((1 + r*mua*tau)/exp(r*mua*tau) - (1 + r*muc*tau)/exp(r*muc*tau))/(r^2*(muc^2 - mus^2)*tau^2)
    mcm <- 2*((1 + rm*mua*tau)/exp(rm*mua*tau) - (1 + rm*muc*tau)/exp(rm*muc*tau))/(rm^2*(muc^2 - mus^2)*tau^2)
    mcp <- 2*((1 + rp*mua*tau)/exp(rp*mua*tau) - (1 + rp*muc*tau)/exp(rp*muc*tau))/(rp^2*(muc^2 - mus^2)*tau^2)

    # mc <- -(((1 + a)*(r*tau)^(-1 - a)*((-exp(-(r*mua*tau)) + exp(-(r*muc*tau)))*muc^a*(r*tau)^a + gammainc(1 + a,r*mua*tau) - gammainc(1 + a,r*muc*tau)))/(a*muc^a*(muc - mus) + mus*(-muc^a + mus^a)))
    # mcm <- -(((1 + a)*(rm*tau)^(-1 - a)*((-exp(-(rm*mua*tau)) + exp(-(rm*muc*tau)))*muc^a*(rm*tau)^a + gammainc(1 + a,rm*mua*tau) - gammainc(1 + a,rm*muc*tau)))/(a*muc^a*(muc - mus) + mus*(-muc^a + mus^a)))
    # mcp <- -(((1 + a)*(rp*tau)^(-1 - a)*((-exp(-(rp*mua*tau)) + exp(-(rp*muc*tau)))*muc^a*(rp*tau)^a + gammainc(1 + a,rp*mua*tau) - gammainc(1 + a,rp*muc*tau)))/(a*muc^a*(muc - mus) + mus*(-muc^a + mus^a)))

    mc <- L0*(mcm + mcp - 2*mc)/dr^2
    mc[is.na(mc)] <- 0

    if (tau<d/mus)
    {
      mh <- (2*(-exp(-(r*muc*tau)) + exp(-(r*mus*tau)) + r*(-muc + mus)*tau))/(r^2*(-muc^2 + mus^2)*tau)
      mhm <- (2*(-exp(-(rm*muc*tau)) + exp(-(rm*mus*tau)) + rm*(-muc + mus)*tau))/(rm^2*(-muc^2 + mus^2)*tau)
      mhp <- (2*(-exp(-(rp*muc*tau)) + exp(-(rp*mus*tau)) + rp*(-muc + mus)*tau))/(rp^2*(-muc^2 + mus^2)*tau)
    }    else
    {
      mh <- (-2*(-exp(-(r*muc*tau)) + r*(-muc + mus)*tau + (1 + r*(d - mus*tau))/exp(r*d)))/(r^2*(muc^2 - mus^2)*tau)
      mhm <- (-2*(-exp(-(rm*muc*tau)) + rm*(-muc + mus)*tau + (1 + rm*(d - mus*tau))/exp(rm*d)))/(rm^2*(muc^2 - mus^2)*tau)
      mhp <- (-2*(-exp(-(rp*muc*tau)) + rp*(-muc + mus)*tau + (1 + rp*(d - mus*tau))/exp(rp*d)))/(rp^2*(muc^2 - mus^2)*tau)
    }
    mh <- L0*rho*(mhm + mhp - 2*mh)/dr^2

    mT <- mc + mh
    
    # mT <- mT/sum(mT*rV*diff(rB))*dat$L[i]

    mexp <- L0*div^2*exp(-div*rV)
    p <- ggplot(data=data.frame(rV=rV[mE>0],mE=mE[mE>0],mT=mT[mE>0],mc=mc[mE>0],mh=mh[mE>0],mexp=mexp[mE>0]),aes(rV[mE>0], mE[mE>0])) +
      # geom_vline(xintercept = 1/muc/tau,size=2, linetype="solid", color = "grey") +
      geom_point(aes(rV[mE>0], mE[mE>0]),col="black",size=0.03) +
      geom_line(aes(rV[mE>0], mc[mE>0]),col="blue",size=0.03) +
      geom_line(aes(rV[mE>0], mh[mE>0]),col="red",size=0.03) +
      geom_line(aes(rV[mE>0], mT[mE>0]),col="black", size=0.25) +
      # geom_line(aes(rV[mE>0], mexp[mE>0]),col="grey",size=0.03) +
      # geom_line(data = data.frame(rVdense=rV,mTdense=mT), aes(rVdense, mTdense),col="black", size=0.3, inherit.aes = TRUE) +
      scale_x_continuous(limits = c(1e0, 1e5),trans='log10',breaks = 10^(-10:10), labels = trans_format("log10", math_format(10^.x))) + expand_limits(x = c(1e-1, 1e5)) +
      scale_y_continuous(limits = c(min(mE[mE>0])/1.1, 5.3*max(c(mE,mT))),trans='log10',breaks = 10^seq(-10,10,2), labels = trans_format("log10", math_format(10^.x))) +
      labs(x="r", y = "m(r)", width = 10) +
      annotate(geom="text", x=rV[max(which(mE>0))]/10.16, y=max(c(mE,mT)), label=TeX(paste("$\\tau = ",formatC(tau, format = "e", digits = 2),"$")),size=2)+
      annotate(geom="text", x=rV[max(which(mE>0))]/10.16, y=max(c(mE,mT))/10, label=TeX(paste("$\\rho = ",formatC(rho, format = "e", digits = 2),"$")),size=2)+
      annotate(geom="text", x=rV[max(which(mE>0))]/10.16, y=max(c(mE,mT))/100, label=TeX(paste("$\\theta = ",round(div,4),"$")),size=2)+
      annotate(geom="text", x=0.4e2, y=min(mE[mE>0])*100, label=paste0(paste(strsplit(species1,"_")[[1]],collapse=' ')," (",dat$n1[i],")"),size=1.8,fontface = "italic")+
      annotate(geom="text", x=0.4e2, y=min(mE[mE>0])*10, label="vs.",size=1.8,fontface = "italic")+
      annotate(geom="text", x=0.4e2, y=min(mE[mE>0]), label=paste0(paste(strsplit(species2,"_")[[1]],collapse=' ')," (",dat$n2[i],")"),size=1.8,fontface = "italic")+
      ggtitle(paste0(paste(strsplit(species1,"_")[[1]],collapse=' ')," (",dat$n1[i],") vs. ",paste(strsplit(species2,"_")[[1]],collapse=' ')," (",dat$n2[i],")"))+
      theme(plot.title = element_text(size=3.7, face="italic"),
            aspect.ratio=1,
            panel.border = element_rect(colour = "black", fill=NA),
            panel.background = element_blank(),
            axis.text = element_text(size=6,face="italic"),
            axis.text.x=element_text(size=7,face="italic"),
            axis.text.y=element_text(size=7,face="italic"))

    list(plot=p,tau=tau,rho=rho,L0=L0,goodnessFits=goodnessFits)
  }

  goodnessFits <- 0
  for (i in (1:nPairs))
  {
    goodnessFits <- goodnessFits + results[[i]][["goodnessFits"]]/nPairs
    dat$tauFit[i] <- results[[i]][["tau"]]
    dat$rho[i] <- results[[i]][["rho"]]
    dat$L0[i] <- results[[i]][["L0"]]
  }

  dat$tau <- dat$tauFit

  dat <- dat[order(dat$tau),]

  tauV <- 10^seq(7,9.5,0.01)

  deltacV <- ifelse (tauV<d/mus,
                     (1+a)/(2+a)*(mus^(2+a)-muc^(2+a))/(mus^(1+a)-muc^(1+a))*tauV,
                     (1+a)/(2+a)*((d/tauV)^(2+a)-muc^(2+a))/((d/tauV)^(1+a)-muc^(1+a))*tauV)

  dat$deltac <- ifelse (dat$tau<d/mus,
                     (1+a)/(2+a)*(mus^(2+a)-muc^(2+a))/(mus^(1+a)-muc^(1+a))*dat$tau,
                     (1+a)/(2+a)*((d/dat$tau)^(2+a)-muc^(2+a))/((d/dat$tau)^(1+a)-muc^(1+a))*dat$tau)

  deltahV <- (a+1)/(a+2)*ifelse (tauV<d/mus,
                     (mus^(a+2)-muc^(a+2))/ (mus^(a+1)-muc^(a+1))*tauV/2,
                     d/2*(2*d^(a+1)+a/d*(muc*tauV)^(a+2)-(a+2)*d*(mus*tauV)^a)/(d^(a+1)+a*(muc*tauV)^(a+1)-(a+1)*d*(mus*tauV)^a)
                     )

  dat$deltah <- (a+1)/(a+2)*ifelse (dat$tau<d/mus,
                                 (mus^(a+2)-muc^(a+2))/ (mus^(a+1)-muc^(a+1))*dat$tau/2,
                                 d/2*(2*d^(a+1)+a/d*(muc*dat$tau)^(a+2)-(a+2)*d*(mus*dat$tau)^a)/(d^(a+1)+a*(muc*dat$tau)^(a+1)-(a+1)*d*(mus*dat$tau)^a)
  )

  deltacV <- 3/4*(1-exp(-deltacV*4/3))
  dat$deltac <- 3/4*(1-exp(-dat$deltac*4/3))
  deltahV <- 3/4*(1-exp(-deltahV*4/3))
  dat$deltah <- 3/4*(1-exp(-dat$deltah*4/3))

  LcV <- ifelse (tauV<d/mus, 1, ((d/tauV)^(a+1)-muc^(a+1))/(mus^(a+1)-muc^(a+1)))
  dat$Lc <- ifelse (dat$tau<d/mus, 1, ((d/dat$tau)^(a+1)-muc^(a+1))/(mus^(a+1)-muc^(a+1)))
  dat$Lh <- dat$rho*ifelse (dat$tau<d/mus, dat$tau, d/mus+1/a/mus/(mus^(a+1)-muc^(a+1))*(a*dat$tau*muc^(a+1)*(d/dat$tau-mus)-d*mus*(d^a/dat$tau^a-mus^a)))

  dat$delta <- (dat$Lh*dat$deltah+dat$Lc*dat$deltac)/(dat$Lc+dat$Lh)

  dat$tauKnown <- sapply(1:nrow(dat),function(i){str <- strsplit(dat$species1[i],"Simulation")[[1]][2];  as.numeric(strsplit(str,"_")[[1]][1])})


  pdf(paste0("/scratch/ws1/msheinman-msheinman/nucmer/plots/pars.pdf"),width=2.5,height=2.5)
  {
    p <- ggplot(dat,aes(x=tau,y=L/L0)) +
      # geom_errorbar(aes(ymin=Lfrac-LfracSD, ymax=Lfrac+LfracSD), width=.2, position=position_dodge(0.05)) +
      scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x), limit=c(6e-6, 1.5),labels = scales::trans_format("log10", scales::math_format(10^.x))) +
      scale_x_log10(labels = scales::trans_format("log10", scales::math_format(10^.x)), limit=c(0.9e7, 1.5e9), breaks=10^(7:9)) +
      geom_line(aes(x=tau,y=rho*tau),col="grey", size=1,linetype = "dotted") +
      geom_line(aes(x=tau,y=Lh),col="red", size=0.5) +
      # geom_line(aes(x=tau,y=Lc+Lh),col="black") +
      geom_line(data=data.frame(tauV=tauV,LcV=LcV),aes(x=tauV,y=LcV),col="blue") +
      geom_point(shape=1, fill=NA, color="black", size=2) +
      # xlab(TeX(paste("$\\tau"))) + ylab(TeX(paste("$\\frac{L}{L_0}"))) +
      theme(panel.border = element_rect(colour = "black", fill=NA),aspect.ratio=1,axis.title.y = element_text(angle = 0, vjust = 0.5))
    print(p)

    p <- ggplot(dat,aes(x=tau,y=div),) +
      # geom_errorbar(aes(ymin=div-divSD, ymax=div+divSD), width=.2, position=position_dodge(0.05)) +
      scale_x_log10(labels = scales::trans_format("log10", scales::math_format(10^.x)), limit=c(0.9e7, 1.5e9), breaks=10^(7:9)) +
      # geom_line(aes(x=tau,y=delta))
      geom_line(data=data.frame(tauV=tauV,deltacV=deltacV),aes(x=tauV,y=deltacV),col="blue") +
      geom_line(data=data.frame(tauV=tauV,deltahV=deltahV),aes(x=tauV,y=deltahV),col="red") +
      geom_point(shape=0, fill=NA, color="black", size=2) +
      # geom_line(aes(x=tau,y=deltac),col="black") +
      ylim(0, 0.2) + 
      # xlab(TeX(paste("$\\tau"))) + ylab(TeX(paste("$\\theta"))) +
      theme(panel.border = element_rect(colour = "black", fill=NA),aspect.ratio=1,axis.title.y = element_text(angle = 0, vjust = 0.5))
    print(p)

    p <- ggplot(dat,aes(x=tau,y=rho)) +
      geom_point(shape=4, fill=NA, color="black", size=2) +
      theme(panel.border = element_rect(colour = "black", fill=NA),aspect.ratio=1,axis.title.y = element_text(angle = 0, vjust = 0.5)) +
      # xlab(TeX(paste("$\\tau"))) + ylab(TeX(paste("$\\rho"))) +
      # geom_line(aes(x=tau,y=1e6/tau^2),col="blue") +
      scale_x_log10(labels = scales::trans_format("log10", scales::math_format(10^.x)), limit=c(0.9e7, 1.5e9), breaks=10^(7:9)) +
      scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),labels = scales::trans_format("log10", scales::math_format(10^.x)))
      
      print(p)

      p <- ggplot(dat,aes(x=tau,y=L0)) +
        geom_point(shape=1, fill=NA, color="black", size=1) +
        scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),labels = scales::trans_format("log10", scales::math_format(10^.x))) +
        scale_x_log10(labels = scales::trans_format("log10", scales::math_format(10^.x)), limit=c(0.9e7, 1e9), breaks=10^(7:9)) +
        xlab(TeX(paste("$\\tau"))) + ylab(TeX(paste("$\\frac{L}{L_0}"))) +
        theme(panel.border = element_rect(colour = "black", fill=NA),aspect.ratio=1,axis.title.y = element_text(angle = 0, vjust = 0.5))
      print(p)
      p <- ggplot(dat,aes(x=tau,y=tauKnown)) +
        geom_point(shape=1, fill=NA, color="black", size=1) +
        scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),labels = scales::trans_format("log10", scales::math_format(10^.x))) +
        scale_x_log10(labels = scales::trans_format("log10", scales::math_format(10^.x)), limit=c(0.9e7, 1e9), breaks=10^(7:9)) +
        xlab(TeX(paste("$\\tau"))) + ylab(TeX(paste("$\\frac{L}{L_0}"))) +
        geom_line(aes(x=tau,y=tau),col="blue") +
        theme(panel.border = element_rect(colour = "black", fill=NA),aspect.ratio=1,axis.title.y = element_text(angle = 0, vjust = 0.5))
      print(p)
  }
  dev.off()

  distance_matrix <- matrix(0,nrow=length(species),ncol=length(species))
  colnames(distance_matrix) <- species
  rownames(distance_matrix) <- species
  div_matrix <- distance_matrix
  for (i in 1:nrow(dat))
  {
    species1 <- dat$species1[i]
    species2 <- dat$species2[i]
    distance_matrix[species1,species2] <- as.numeric(dat$tau[i])
    distance_matrix[species2,species1] <- as.numeric(dat$tau[i])
    div_matrix[species2,species1] <- as.numeric(dat$div[i])
    div_matrix[species1,species2] <- as.numeric(dat$div[i])
  }

  # distance_matrix <- distance_matrix[rownames(distance_matrix)!="Vibrio",]
  # distance_matrix <- distance_matrix[,colnames(distance_matrix)!="Vibrio"]
  # distance_matrix <- distance_matrix[rownames(distance_matrix)!="Cronobacter",]
  # distance_matrix <- distance_matrix[,colnames(distance_matrix)!="Cronobacter"]

  # distance_matrix["Enterobacter_asburiae","Cronobacter"] <- 1
  # distance_matrix["Cronobacter","Enterobacter_asburiae"] <- 1
  # distance_matrix["Enterobacter_roggenkampii","Cronobacter"] <- 1
  # distance_matrix["Cronobacter","Enterobacter_roggenkampii"] <- 1
  # distance_matrix["Enterobacter_hormaechei","Cronobacter"] <- 1
  # distance_matrix["Cronobacter","Enterobacter_hormaechei"] <- 1

  pdf(paste0("/scratch/ws1/msheinman-msheinman/nucmer/plots/tree.pdf"),width=8,height=4)
  {
    # tree <- phangorn::upgma(distance_matrix)
    tree <- as.phylo(hclust(as.dist(distance_matrix), method = "average"))
    plot(tree,font=3,cex=0.5,direction = "leftwards")
    axisPhylo(side = 1, root.time = NULL, backward = TRUE)
    x <- as.vector(as.dist(distance_matrix))
    y <- as.vector(as.dist(cophenetic(tree)))

    p <- ggplot(data=data.frame(x=(y)/2, y=(x)/2),aes(x=x,y=y))
    p <- p + geom_point(data=data.frame(x=(y)/2, y=(x)/2),aes(x=x,y=y)) +
      scale_y_log10(labels = scales::trans_format("log10", scales::math_format(10^.x)), breaks=c(1e7,1e8,1e9), limit=c(1e7, 1e9)) +
      scale_x_log10(labels = scales::trans_format("log10", scales::math_format(10^.x)), breaks=c(1e7,1e8,1e9), limit=c(1e7, 1e9)) +
      geom_line(data=data.frame(x=(y)/2, y=(x)/2),aes(x=x,y=x),col="blue") +
      theme(panel.border = element_rect(colour = "black", fill=NA),aspect.ratio=1,axis.title.y = element_text(angle = 0, vjust = 0.5),text = element_text(size = 20))
    print(p)


    plot((y)/2, (y)/2,lty="solid", log="xy",xaxt = "n", yaxt = "n", type="l")
    for (liney in unique(y/2)) {lines(c(liney,liney),c(1e5,1e10),lty="dashed",col="grey")}
    points((y)/2*(1+0*runif(length(y))), (x)/2, pch=20, cex=1)
    axis(1, at = c(1e7*seq(1,9,1),1e8*seq(1,9,1)))
    axis(2, at = c(1e7*seq(1,9,1),1e8*seq(1,9,1)), labels=format(c(1e7*seq(1,9,1),1e8*seq(1,9,1)), scientific=TRUE), hadj=0.9, cex.axis=0.8, las=2)


    # tree <- phangorn::upgma(distance_matrix)
    tree <- as.phylo(hclust(as.dist(div_matrix), method = "average"))
    plot(tree,font=3,cex=0.5,direction = "leftwards")
    axisPhylo(side = 1, root.time = NULL, backward = TRUE)
    x <- as.vector(as.dist(div_matrix))
    y <- as.vector(as.dist(cophenetic(tree)))

    p <- ggplot(data=data.frame(x=(y)/2, y=(x)/2),aes(x=x,y=y))
    p <- p + geom_point(data=data.frame(x=(y)/2, y=(x)/2),aes(x=x,y=y)) +
      scale_y_log10(labels = scales::trans_format("log10", scales::math_format(10^.x)), breaks=c(1e7,1e8), limit=c(1e7, 5e8)) +
      scale_x_log10(labels = scales::trans_format("log10", scales::math_format(10^.x)), breaks=c(1e7,1e8), limit=c(1e7, 5e8)) +
      geom_line(data=data.frame(x=(y)/2, y=(x)/2),aes(x=x,y=x),col="blue") +
      theme(panel.border = element_rect(colour = "black", fill=NA),aspect.ratio=1,axis.title.y = element_text(angle = 0, vjust = 0.5),text = element_text(size = 20))
    print(p)


    plot((y)/2, (y)/2,lty="solid", log="xy",xaxt = "n", yaxt = "n", type="l")
    for (liney in unique(y/2)) {lines(c(liney,liney),c(1e5,1e10),lty="dashed",col="grey")}
    points((y)/2*(1+0*runif(length(y))), (x)/2, pch=20, cex=1)
    axis(1, at = c(1e7*seq(1,9,1),1e8*seq(1,9,1)))
    axis(2, at = c(1e7*seq(1,9,1),1e8*seq(1,9,1)), labels=format(c(1e7*seq(1,9,1),1e8*seq(1,9,1)), scientific=TRUE), hadj=0.9, cex.axis=0.8, las=2)



    # lines(log10(x), log10(x), col="red")
    # abline(lm(y~x), col="red")
    tree <- ape::nj(distance_matrix)
    x <- as.vector(as.dist(distance_matrix))
    y <- as.vector(as.dist(cophenetic(tree)))
    tree <- root(tree, out = "Vibrio")
    tree <- midpoint.root(tree)

    plot(tree,font=3,cex=0.5,direction = "leftwards")
    axisPhylo(side = 1, root.time = NULL, backward = TRUE)
    plot(log10(y/2), log10(x/2), pch=20, cex=1)
    lines(log10(y/2), log10(y/2))

    p <- ggplot(data=data.frame(x=y/2, y=x/2),aes(x=x,y=y))
    p <- p + geom_point(data=data.frame(x=y/2, y=x/2),aes(x=x,y=y)) +
      scale_y_log10(labels = scales::trans_format("log10", scales::math_format(10^.x)), breaks=c(1e7,1e8), limit=c(1e7, 5e8)) +
      scale_x_log10(labels = scales::trans_format("log10", scales::math_format(10^.x)), breaks=c(1e7,1e8), limit=c(1e7, 5e8)) +
      geom_line(data=data.frame(x=(y)/2, y=(x)/2),aes(x=x,y=x),col="blue") +
      theme(panel.border = element_rect(colour = "black", fill=NA),aspect.ratio=1,axis.title.y = element_text(angle = 0, vjust = 0.5),text = element_text(size = 20))
    print(p)


    p <- Heatmap(distance_matrix,
            col=colorRamp2(c(min(as.dist(distance_matrix)), max(as.dist(distance_matrix))), c("white", "red")),
            cluster_columns = hclust(as.dist(distance_matrix),method="average"),
            cluster_rows = hclust(as.dist(distance_matrix),method="average"))
    print(p)
  }
  dev.off()

  # dL <- mean((dat$Lc-dat$L)^2/(dat$L)^2)^0.5
  ddelta <- mean((dat$deltac-dat$div)^2/(dat$div)^2)^0.5
  print(c(ddelta,goodnessFits))
  goodnessFits <- 0

  pdf(paste0("/scratch/ws1/msheinman-msheinman/nucmer/plots/m.pdf"),width=2,height=2)
  for (i in (1:nPairs))
  {
    print(results[[i]][["plot"]])
  }
  dev.off()

  return(goodnessFits)
}

# par <- log10(c(d,muc))
# res <- Nelder_Mead(par=par,
#              lower=log10(c(0.19,muc/10)),
#              upper=log10(c(0.25,muc*10)),
#              fn=global_fit
#              )
#
# par <- log10(muc)
# optim(par=par, fn=global_fit,
#       method = c("Brent"),
#       lower = log10(muc/10), upper = log10(muc*10),
#       control = list(), hessian = FALSE)
#-------------------------------------------------------------------------------------------------------------------------


datdiv <- foreach (s1 = 1:(length(species)-1), .inorder=FALSE, .combine=rbind) %dopar%
{
  foreach (s2 = (s1+1):length(species), .inorder=FALSE, .combine=rbind) %dopar%
  {
    dat <- data.frame()
    species1 <- species[s1]
    species2 <- species[s2]
    files <- sort(Sys.glob(file.path(paste0("/scratch/ws1/msheinman-msheinman/nucmer/",species1,"_vs_",species2,"/data_*.L"))))
    filePairs <- sort(Sys.glob(file.path(paste0("/scratch/ws1/msheinman-msheinman/nucmer/",species1,"_vs_",species2,"/listOfPairs.csv"))))

    if (length(files)==0 & s1!=s2)
    {
      species2 <- species[s1]
      species1 <- species[s2]
      files <- sort(Sys.glob(file.path(paste0("/scratch/ws1/msheinman-msheinman/nucmer/",species1,"_vs_",species2,"/data_*.L"))))
      filePairs <- sort(Sys.glob(file.path(paste0("/scratch/ws1/msheinman-msheinman/nucmer/",species1,"_vs_",species2,"/listOfPairs.csv"))))

      if (length(files)==0)
      {
        print(paste0(species1,"_vs_",species2," NOT FOUND!"))
      }
    }
    if (length(files)>0)
    {
      print(paste0(species1,"_vs_",species2))
      Lall <- foreach (f = files, .combine=rbind, .inorder=TRUE) %do% {read.table(file = f ,header=TRUE, row.names=NULL,sep="\t",skip = 0)}
      if (nrow(Lall)>0)
      {
        dataPairs <- read.table(file = filePairs ,header=FALSE, row.names=NULL,sep="\t",skip = 0)
        dat <- data.frame(species1=species1,species2=species2,name=paste0(species1," vs. ",species2),
                          div=Lall$tau,divmedian=median(Lall$tau),
                          Lfrac=Lall$L/min(c(mean(dataPairs[,5]),mean(dataPairs[,6]))),
                          Lfracmedian=median(Lall$L/min(c(mean(dataPairs[,5]),mean(dataPairs[,6])))))
      }
      else print(paste0(species1,"_vs_",species2," EMPTY!"))
    }
    dat
  }
}

datdiv$divmedian = factor(datdiv$divmedian , level = sort(unique(datdiv$divmedian)) )
datdiv$Lfracmedian = factor(datdiv$Lfracmedian , level = sort(unique(datdiv$Lfracmedian)) )


pdf(paste0("/scratch/ws1/msheinman-msheinman/nucmer/plots/div.pdf"),width=12,height=12)
p <- ggplot(datdiv, aes(x=reorder(name,div,na.rm = TRUE), y=div,fill=divmedian)) +
  geom_violin(adjust = 1,linewidth=0.05) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=10),panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5),legend.position="none")
print(p)
p <- ggplot(datdiv, aes(x=reorder(name,-Lfrac,na.rm = TRUE), y=Lfrac,fill=Lfracmedian)) +
  geom_violin(adjust = 2,linewidth=0.05) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=10),panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5),legend.position="none") +
  scale_y_continuous(breaks = c(0.001,0.01,0.1,1), trans = "log10")
print(p)
dev.off()

