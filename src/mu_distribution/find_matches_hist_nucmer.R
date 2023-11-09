#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
species1 <- args[1]
species2 <- args[2]
ncores <- as.numeric(args[3])
batch <- as.numeric(args[4])
NumberBatches <- as.numeric(args[5])
Sys.sleep(1+3*batch)
print(paste0("Doing ", species1," vs. ",species2," with ",ncores," cores. Batch number ",batch," from ",NumberBatches," batches." ))
outdir <- paste0("/scratch/ws1/msheinman-sheinman/nucmer/",species1,"_vs_",species2)
system(paste0("mkdir ",outdir,"/"))
setwd(outdir)
system(paste0("cd ",outdir,"/"))
print(0)
 
# library(seqinr);
# ncores <- 101
# batch <- 1
# NumberBatches <- 1
# species1 <- "Simulation1e8_init"
# species2 <- "Simulation1e8"
# outdir <- paste0("/scratch/ws1/msheinman-sheinman/",species1,"_vs_",species2)
# system(paste0("mkdir ",outdir,"/"))

# run this
# install.packages("/home/msheinman/Software/segmut/", repos = NULL, type = "source",force=TRUE)
# species1 <- "Escherichia_coli"
# species1 <- "Enterobacter_asburiae"
# species1 <- "Escherichia_fergusonii"
# species2 <- "Salmonella_enterica"
# species1 <- "Escherichia_albertii"
# species1 <- "Klebsiella_pneumoniae"
# species2 <- "Raoultella"
# ncores <- 128

# import libraries----
library(stringr);
library(parallel);
library(foreach);
library(iterators);
library(doParallel);
library(data.table);
library(seqinr);
library(plotrix);
library(stringi);
library(optimx);
library(splitstackshape);
library(segmut)
registerDoParallel(ncores);


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
  dir1 <- paste0("/scratch/ws1/msheinman-sheinman/",species1); RemovePlasmids(dir1)
  files1 <- sort(Sys.glob(file.path(dir1, "*_chr.fna")))
  files1 <- sample(files1,replace=FALSE)

  dir2 <- paste0("/scratch/ws1/msheinman-sheinman/",species2); RemovePlasmids(dir2)
  files2 <- sort(Sys.glob(file.path(dir2, "*_chr.fna")))
  files2 <- sample(files2,replace=FALSE)


  Sys.sleep(1)
  registerDoParallel(ncores)
  print("record lengths")
  print(paste0("/scratch/ws1/msheinman-sheinman/",species1))
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
  iV <- c(seq(1,nrow(listOfPairs),by=ncores*10),nrow(listOfPairs)+1)
  rB <- c(seq(-0.5,35.5,3),10^seq(log10(35.5)+0.1,5.5,0.1))
  write.table(matrix(c("Ind1","Ind2","AccNum1","AccNum2","L","tau",0.5*(rB[-length(rB)]+rB[-1])),nrow=1),
              file = paste0(outdir,"/data_",batch,".L"),append=FALSE,row.names=FALSE,col.names=FALSE,sep = "\t",quote=FALSE)
  # write.table(matrix(c("Length"),nrow=1),file = paste0(outdir,"/coords.L"),append=FALSE,row.names=FALSE,col.names=FALSE,sep = "\t",quote=FALSE)
  for (iVt in 1:(length(iV)-1))
  {
    print(paste0(iVt," from ",(length(iV)-1)," batch ",batch," ",Sys.time()))
    dat <- foreach (i1 = iV[iVt]:(iV[iVt+1]-1), .combine=rbind, .inorder=TRUE) %dopar%
    { 
      Length1 <- listOfPairs$Length1[i1]
      header1 <- listOfPairs$header1[i1] 
      Length2 <- listOfPairs$Length2[i1]
      header2 <- listOfPairs$header2[i1]
      ind1 <- listOfPairs$ind1[i1]
      ind2 <- listOfPairs$ind2[i1]
      system(paste0("/scratch/ws1/msheinman-sheinman/mummer4.0/bin/nucmer --threads 1 --mum --prefix=",outdir,"/temp/Batch_",ind1,"_",ind2," ",listOfPairs$file1[i1]," ",listOfPairs$file2[i1]), intern = TRUE, wait=TRUE)
      # coords <- system(paste0("/home/msheinman/Software/mummer4.0/bin/show-coords -c -d -r -l -o -T  ",outdir,"/temp/Batch_",ind1,"_",ind2,".delta"), intern = TRUE, wait=TRUE)
      # coords <- as.data.frame(coords[5:length(coords)])
      # coords <- as.data.frame(cSplit(coords,splitCols=colnames(coords), "\t"))
      # coords <- data.frame(round(coords[,5]+coords[,6])/2)
      alignment <- system(paste0("/scratch/ws1/msheinman-sheinman/mummer4.0/bin/show-aligns -w ",1+max(c(Length1,Length2))," ",outdir,"/temp/Batch_",ind1,"_",ind2,".delta ", header1," ",header2), intern = TRUE, wait=TRUE)
      alignment <- substr(alignment,12,nchar(alignment)) # removing empty space in the front
      if (length(alignment)==0) return(data.frame())
      # replace long indels by 1bp indel to make sure that (1 indel = 1 mutation)
      Ls <- c()
      number_muts <- c()
      r <- c()
      mutsS <- list()
      for (i in seq(1,length(alignment)-1,1))
      {
        if (substr(alignment[[i]],1,1) %in% c("a","c","t","g","A","C","T","G") & substr(alignment[[i+1]],1,1) %in% c("a","c","t","g","A","C","T","G"))
        {
          str1 <- str_split(alignment[[i]],"")[[1]]
          str2 <- str_split(alignment[[i+1]],"")[[1]]
          strmuts <- str_split(alignment[[i+2]],"")[[1]]
          Ind <- which(str1[1:(length(str1)-1)]=="." & str1[2:length(str1)]==".")
          while (length(Ind)>0)
          {
            str1 <- str1[-Ind]
            str2 <- str2[-Ind]
            strmuts <- strmuts[-Ind]
            Ind <- which(str1[1:(length(str1)-1)]=="." & str1[2:length(str1)]==".")
          }
          Ind <- which(str2[1:(length(str2)-1)]=="." & str2[2:length(str2)]==".")
          while (length(Ind)>0)
          {
            str1 <- str1[-Ind]
            str2 <- str2[-Ind]
            strmuts <- strmuts[-Ind]
            Ind <- which(str2[1:(length(str2)-1)]=="." & str2[2:length(str2)]==".")
          }
          Ls <- c(Ls,length(strmuts))
          number_muts <- c(number_muts,sum(strmuts=="^"))
          r <- c(r,diff(c(0,which(strmuts=="^"),length(strmuts)+1))-1)
          mutsS[[length(mutsS)+1]] <- which(strmuts=="^")
        }
      }
      divergence <- sum(number_muts)/sum(Ls)
   
      counts <- hist(r,breaks=rB,plot=FALSE)$counts
      dat <- data.frame(Ind1=ind1,Ind2=ind2,AccNum1=header1,AccNum2=header2,L=sum(Ls),div=divergence)
      dat <- cbind(dat,matrix(counts,nrow=1))
      system(paste0("rm ",outdir,"/temp/Batch_",ind1,"_",ind2,".delta "),wait=TRUE,intern = TRUE)
      return(dat)
    }
    write.table(dat,file = paste0(outdir,"/data_",batch,".L"),append=TRUE,row.names=FALSE,col.names=FALSE,sep = "\t",quote=FALSE)
  }
  gc()
  print(paste0("Batch ",batch," completed on ",Sys.time()))
}



















