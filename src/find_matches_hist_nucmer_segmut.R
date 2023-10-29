#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
species1 <- args[1]
species2 <- args[2]
ncores <- as.numeric(args[3])
batch <- as.numeric(args[4])
NumberBatches <- as.numeric(args[5])
Sys.sleep(1+3*batch)
print(paste0("Doing ", species1," vs. ",species2," with ",ncores," cores. Batch number ",batch," from ",NumberBatches," batches." ))
outdir <- paste0("/scratch/ws1/msheinman-msheinman/",species1,"_vs_",species2)
system(paste0("mkdir ",outdir,"/"))
setwd(outdir)
system(paste0("cd ",outdir,"/"))
print(0)
divB <- 10^seq(-10,0,0.01)
howoften <- 1

# library(seqinr);
# ncores <- 101
# batch <- 1
# NumberBatches <- 1
# species1 <- "Escherichia_coli"
# species2 <- "Salmonella_enterica"
# outdir <- paste0("/scratch/ws1/msheinman-msheinman/",species1,"_vs_",species2)
# system(paste0("mkdir ",outdir,"/"))

# run this
# install.packages("/home/msheinman/Software/", repos = NULL, type = "source",force=TRUE)
# species1 <- "Escherichia_coli"
# species1 <- "Enterobacter_asburiae"
# species1 <- "Escherichia_fergusonii"
# species2 <- "Salmonella_enterica"
# species1 <- "Escherichia_albertii"
# species1 <- "Klebsiella_pneumoniae"
# species2 <- "Raoultella"
# ncores <- 128

# import libraries----
library(stringr);print("loaded all libraries")
library(parallel);print("loaded all libraries")
library(foreach);print("loaded all libraries")
library(iterators);print("loaded all libraries")
library(doParallel);print("loaded all libraries")
library(data.table);print("loaded all libraries")
library(seqinr);print("loaded all libraries")
library(plotrix);print("loaded all libraries")
library(stringi);print("loaded all libraries")
library(optimx);print("loaded all libraries")
library(splitstackshape);print("loaded all libraries")
library(segmut);print("loaded all libraries")
registerDoParallel(ncores);print("loaded all libraries")


# function to remove plasmids----
RemovePlasmids <- function(dir)
{
  # setwd(dir)
  files <- Sys.glob(file.path(dir, "*_genomic.fna"))

  foo <- foreach (file = files) %dopar%
  {
    # system(paste0("awk '/^>/{if(N==2)exit;++N;} {print;}' ",file," > ",str_replace(file,"_genomic.fna","_chr.fna")))
    # system(paste0("awk '/^>/ {ok=index($0,'chromosome');} {if(ok) print;}' ",file," > ",str_replace(file,"_genomic.fna","_chr.fna")))
    system(paste0("awk '/^>/{x = /plasmid/==0;}(x)' ",file," > ",str_replace(file,"_genomic.fna","_still_genomic.fna")), wait=TRUE)
    system(paste0("awk 'NR <= 1 || !/^>/' ", str_replace(file,"_genomic.fna","_still_genomic.fna"), " > ", str_replace(file,"_genomic.fna","_still2_genomic.fna")))
    system(paste0("tr A-Z a-z < ",str_replace(file,"_genomic.fna","_still2_genomic.fna") ," > ",str_replace(file,"_genomic.fna","_chr.fna")))
  }
  system(paste0("rm ",dir,"/*_genomic.fna"))
  return()
}

system(paste0("mkdir ",outdir,"/temp/"))
system(paste0("mkdir ",outdir,"/plots/"))
# system(paste0("rm -rf ",outdir))

# make list of pairs----
if (batch==0) #prepare
{
  print("make list of pairs")
  dir1 <- paste0("/scratch/ws1/msheinman-msheinman/",species1); RemovePlasmids(dir1)
  files1 <- sort(Sys.glob(file.path(dir1, "*_chr.fna")))
  files1 <- sample(files1,replace=FALSE)

  dir2 <- paste0("/scratch/ws1/msheinman-msheinman/",species2); RemovePlasmids(dir2)
  files2 <- sort(Sys.glob(file.path(dir2, "*_chr.fna")))
  files2 <- sample(files2,replace=FALSE)


  Sys.sleep(1)
  registerDoParallel(ncores)
  print("record lengths")
  print(paste0("/scratch/ws1/msheinman-msheinman/",species1))
  FastaLengths1 <- foreach(i=1:length(files1), .inorder=TRUE, .combine=c) %dopar%
  {
    getLength(read.fasta(file = files1[i]))
  }

  Sys.sleep(1)
  FastaLengths2 <- foreach(i=1:length(files2), .inorder=TRUE, .combine=c) %dopar%
  {
    getLength(read.fasta(file = files2[i]))
  }
  Sys.sleep(1)
  print("recorded lengths")
  FastaHeaders1 <- foreach(i=1:length(files1), .inorder=TRUE, .combine=c) %dopar%
  {
    headers <- readLines(con = files1[i], n = -1L, ok = TRUE, warn = TRUE, encoding = "unknown", skipNul = FALSE)
    headers <- headers[substr(headers,1,1)==">"]
    headers <- sapply(headers,function(h){strsplit(h," ")[[1]][1]})
    headers <- str_replace(headers,">","")
    headers
  }
  FastaHeaders2 <- foreach(i=1:length(files2), .inorder=TRUE, .combine=c) %dopar%
  {
    headers <- readLines(con = files2[i], n = -1L, ok = TRUE, warn = TRUE, encoding = "unknown", skipNul = FALSE)
    headers <- headers[substr(headers,1,1)==">"]
    headers <- sapply(headers,function(h){strsplit(h," ")[[1]][1]})
    headers <- str_replace(headers,">","")
    headers
  }
  Sys.sleep(1)
  print("recorded headers")
  listOfFiles1 <- data.frame(ind=1:length(files1),file=files1,header=FastaHeaders1,Length=FastaLengths1)
  write.table(listOfFiles1,file = paste0(outdir,"/listOfFiles1.csv"),append=FALSE,row.names=FALSE,col.names=FALSE,sep = "\t",quote=FALSE)
  listOfFiles2 <- data.frame(ind=1:length(files2),file=files2,header=FastaHeaders2,Length=FastaLengths2)
  write.table(listOfFiles2,file = paste0(outdir,"/listOfFiles2.csv"),append=FALSE,row.names=FALSE,col.names=FALSE,sep = "\t",quote=FALSE)

  listOfPairs <- as.data.frame(expand.grid(1:length(files1),1:length(files2)))
  colnames(listOfPairs) <- c("ind1","ind2")
  listOfPairs$header1 <- listOfFiles1$header[listOfPairs$ind1]
  listOfPairs$header2 <- listOfFiles2$header[listOfPairs$ind2]
  listOfPairs$Length1 <- listOfFiles1$Length[listOfPairs$ind1]
  listOfPairs$Length2 <- listOfFiles2$Length[listOfPairs$ind2]
  listOfPairs$file1 <- listOfFiles1$file[listOfPairs$ind1]
  listOfPairs$file2 <- listOfFiles2$file[listOfPairs$ind2]
  listOfPairs <- listOfPairs[listOfPairs$file1!=listOfPairs$file2,]
  write.table(listOfPairs,file = paste0(outdir,"/listOfPairs.csv"),append=FALSE,row.names=FALSE,col.names=FALSE,sep = "\t",quote=FALSE)
} else # make comparisons----
{
  rbindlist2 <- function(list1,list2) {return(list(rbind(list1[[1]],list2[[1]]),rbind(list1[[2]],list2[[2]])))}

  listOfPairs <- NA
  listOfPairs <- tryCatch(
  {
    read.table(paste0(outdir,"/listOfPairs.csv"), header = FALSE, sep = "\t",row.names=NULL)
  },
  error=function(cond){print(paste0("failed (batch ",batch,")"));print(cond);return(NA)},
  warning=function(cond){print(paste0("failed (batch ",batch,")"));print(cond);return(NA)}
    )

  if(any(is.na(listOfPairs)))
  {
    print(paste0("failed (batch ",batch,")"))
    fawefasfasfas <- asedfasedfas
  }

  colnames(listOfPairs) <- c("ind1","ind2","header1","header2","Length1","Length2","file1","file2")

  breaksBatches <- round(seq(1,nrow(listOfPairs)+1,length.out=NumberBatches+1))
  listOfPairs <- listOfPairs[breaksBatches[batch]:(breaksBatches[batch+1]-1),]
  gc()
  registerDoParallel(ncores)
  iV <- c(seq(1,nrow(listOfPairs),by=ncores*howoften),nrow(listOfPairs)+1)
  rB <- c(seq(-0.5,35.5,3),10^seq(log10(35.5)+0.1,5.5,0.1))
  datsum <- rep(0,length(divB)-1)
  for (iVt in 1:(length(iV)-1))
  {
    print(paste0(iVt," from ",(length(iV)-1)," batch ",batch," ",Sys.time()))
    dat <- foreach (i1 = iV[iVt]:(iV[iVt+1]-1), .combine='+', .inorder=TRUE) %dopar%
    {
      Length1 <- listOfPairs$Length1[i1]
      header1 <- listOfPairs$header1[i1]
      Length2 <- listOfPairs$Length2[i1]
      header2 <- listOfPairs$header2[i1]
      ind1 <- listOfPairs$ind1[i1]
      ind2 <- listOfPairs$ind2[i1]
      #system(paste0("/scratch/ws1/msheinman-msheinman/mummer4.0/bin/nucmer --threads 1 --mum --prefix=",outdir,"/temp/Batch_",ind1,"_",ind2," ",listOfPairs$file1[i1]," ",listOfPairs$file2[i1]), intern = TRUE, wait=TRUE)
      #alignment <- system(paste0("/scratch/ws1/msheinman-msheinman/mummer4.0/bin/show-aligns -w ",1+max(c(Length1,Length2))," ",outdir,"/temp/Batch_",ind1,"_",ind2,".delta ", header1," ",header2), intern = TRUE, wait=TRUE)
      #alignment <- substr(alignment,12,nchar(alignment)) # removing empty space in the front
      alignment <- system(paste0("/home/msheinman/Software/lastz-distrib/bin/lastz ",listOfPairs$file1[i1],"[unmask] ",listOfPairs$file2[i1],"[unmask] --ambiguous=iupac --strand=both --gfextend --chain --gapped --format=maf"), intern = TRUE, wait=TRUE)
      
      
      if (length(alignment)==0) return(data.frame())
      # replace long indels by 1bp indel to make sure that (1 indel = 1 mutation)
      Ls <- c()
      number_muts <- c()
      r <- c()
      mutsS <- list()
      for (i in seq(1,length(alignment)-1,1))
      {
        if (substr(alignment[[i]],1,1) %in% c("s") & substr(alignment[[i+1]],1,1) %in% c("s"))
        {
          str1 <- str_split(alignment[[i]]," ")[[1]]; str1 <- str1[length(str1)]; str1 <- str_split(str1,"")[[1]]
          str2 <- str_split(alignment[[i+1]]," ")[[1]]; str2 <- str2[length(str2)]; str2 <- str_split(str2,"")[[1]]
          Ind <- which(str1[1:(length(str1)-1)]=="-" & str1[2:length(str1)]=="-")
          while (length(Ind)>0)
          {
            str1 <- str1[-Ind]
            str2 <- str2[-Ind]
            Ind <- which(str1[1:(length(str1)-1)]=="-" & str1[2:length(str1)]=="-")
          }
          Ind <- which(str2[1:(length(str2)-1)]=="-" & str2[2:length(str2)]=="-")
          while (length(Ind)>0)
          {
            str1 <- str1[-Ind]
            str2 <- str2[-Ind]
            Ind <- which(str2[1:(length(str2)-1)]=="-" & str2[2:length(str2)]=="-")
          }
          Ls <- c(Ls,length(str1))
          number_muts <- c(number_muts,sum(str1!=str2))
          r <- c(r,diff(c(0,which(str1!=str2),length(str1)+1))-1)
          mutsS[[length(mutsS)+1]] <- which(str1!=str2)
        }
      }
      divergence <- sum(number_muts)/sum(Ls)

      Ks <- c()
      taus <- c()
      nmutsS <- c()
      start_time <- Sys.time()
      k <- 0
      for (i in order(-Ls))
      {
        k <- k+1
        muts <- mutsS[[i]]
        L <- Ls[[i]]
        if (length(muts)>0)
        {
          breaks0L <- improve(muts, L,c(0,L))
        }
        else
        {
          breaks0L <- c(0,L)
        }
        for (j in 1:(length(breaks0L)-1))
        {
          Ks <- c(Ks, breaks0L[j+1]-breaks0L[j]+1)
          nmutsS <- c(nmutsS, sum(muts<breaks0L[j+1] & muts>breaks0L[j]))
        }
      }
      taus <- (2+nmutsS)/Ks
      p <- weighted.hist(taus,w=Ks, breaks=divB,plot=FALSE)
      return(p$counts)
    }
    datsum <- datsum + dat*1.000000001
    write.table(datsum,file = paste0(outdir,"/segmut_",batch,".L"),append=FALSE,row.names=FALSE,col.names=FALSE,sep = "\t",quote=FALSE)
  }
  gc()
  print(paste0("Batch ",batch," completed on ",Sys.time()))
}



















