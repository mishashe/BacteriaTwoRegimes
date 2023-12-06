#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
ncores <- as.numeric(args[1])
number <- as.numeric(args[2])
tau  <- args[3]


 
# tau <- "1e8"
# number <- 100
# ncores <- 101

print(tau)
library(expint)
library(tidyverse)
library(scales)
library(ggplot2)
library(gridExtra)
library(doParallel)
library(foreach)
library(seqinr)

system(paste0("mkdir /scratch/ws1/msheinman-sheinman/Simulation",tau,"/"))
system(paste0("mkdir /scratch/ws1/msheinman-sheinman/Simulation",tau,"_init/"))
system(paste0("mkdir /scratch/ws1/msheinman-sheinman/lastz/Simulation",tau,"_init_vs_Simulation",tau,"/"))

rB <- c(seq(-0.5,35.5,3),10^seq(log10(35.5)+0.1,5.5,0.1))
rV <- 0.5*(rB[-length(rB)]+rB[-1])
L <- 1e7
setwd("/home/msheinman/Development/HAT/src/")
Rcpp::sourceCpp(paste0("/home/msheinman/Development/HAT/src/SimulateDivergence.cpp"))
# Genome0 <- sample(c('a','c','t','g'),L,replace = TRUE)
Genome0 <- read.fasta(file = "/scratch/ws1/msheinman-sheinman/Escherichia_coli/GCF_022453585.1_ASM2245358v1_chr.fna")$NZ_CP092647.1
L <- length(Genome0)
registerDoParallel(ncores)

print("start simulating genomes")

write.fasta(Genome0, names=paste0("Simulation_",0), file.out=paste0("/scratch/ws1/msheinman-sheinman/Simulation",tau,"_init/Simulation_",0,"_chr.fna"), 
            open = "w", nbchar = 60, as.string = FALSE)
mEr <- foreach(a=1:number, .inorder=FALSE, .combine='+') %dopar%
{
  Genome <- SimulateDivergence(Genome0[1:L],a, as.numeric(tau))
  # Ind <- sample(1:L,round(L/500),replace=FALSE)
  # Ind <- c(Ind+1,Ind+2,Ind+3,Ind+4,Ind+5,Ind+6,Ind+7,Ind+8,Ind+9,Ind+10,Ind+11,Ind+12,Ind+13,Ind+14,Ind+15,Ind+16,Ind+17,Ind+18,Ind+19,Ind+20)
  # Genome[Ind] <- Genome0[Ind]
  write.fasta(Genome, names=paste0("Simulation_",a), file.out=paste0("/scratch/ws1/msheinman-sheinman/Simulation",tau,"/Simulation_",a,"_chr.fna"), 
                                                open = "w", nbchar = 60, as.string = FALSE)
  # r <- diff(sort(which(Genome0!=Genome)))-1
  # hist(r,breaks=rB,plot=FALSE)$counts
}/number

