library(stringr)
library(dplyr);
library(plyr);
library(phylotools)
geneid <- read.csv("geneid.csv",header = F)
geneid <- geneid[order(geneid[,2],decreasing = F),]
chr_wd <- "gene_sam"
fa_wd <-"gene_fa"
sam_filelist <- list.files(path = chr_wd,pattern = ".*.txt")
sam_filelist <- str_replace_all(sam_filelist, "\r", "")
fasta_filelist <- list.files(path = fa_wd,pattern = ".*.fa")
fasta_filelist <- str_replace_all(fasta_filelist, "\r", "")
gene_bed_4 <- read.csv("gene_bed_4.csv",header = T)
gene_file <- matrix(0,nrow = length(fasta_filelist),ncol = 2)
for(i in seq(1,length(sam_filelist)))
{
  loc <- which(geneid[,2]==sam_filelist[[i]])
  gene_file[i,1] <- geneid[loc,1]
  gene_file[i,2] <- geneid[loc,2]
}
gene_file <- as.data.frame(gene_file)
gene_bed_4 <- merge(gene_file,gene_bed_4,by.x="V1",by.y = "SYMBOL")
gene_bed_4 <- gene_bed_4[match(gene_file$V2, gene_bed_4$V2),]
total_transfererror <- matrix(0,nrow = length(gene_file[,2]),ncol = 3)
rownames(total_transfererror) <- gene_file[,1]
colnames(total_transfererror) <- c("Bases sequenced","Errors detected","Error rate")
total_transfererror <- as.data.frame(total_transfererror)
depth <- read.table("bedtools_fildepth.txt")
#sam_CIGAR
sam_CIGAR <- function(x){
  x$V6 <- gsub("S","-S-",x$V6)
  x$V6 <- gsub("M","-M-",x$V6)
  x$V6 <- gsub("D","-D-",x$V6)
  x$V6 <- gsub("I","-I-",x$V6)
  x$V6 <- gsub("N","-N-",x$V6)
  x$V6 <- gsub("H","-H-",x$V6)
  x$V6 <- gsub("P","-P-",x$V6)
  x$V6 <- gsub("=","-=-",x$V6)
  x$V6 <- gsub("X","-X-",x$V6)
  for(i in seq(1,nrow(x)))
  {
    if(x[i,6]!="*")
    {
      CIGAR <-  unlist(strsplit(x[i,6],split = '-'))
      len <- length(CIGAR)
      if(len==2)
      {x[i,6] <- paste(CIGAR[[1]],CIGAR[[2]],sep = "-")}
      else
      {
        x[i,6] <- CIGAR[[1]]
        for(j in seq(1,len-1))
        {
          x[i,6] <- paste(x[i,6],CIGAR[[j+1]],sep = "-")
        }
      }
    }
  }
  return(x)
}
#aftererror
aftererror <- function(x,exon_start){
  after_sam1 <- x[x$V4>=exon_start,]
  if(nrow(after_sam1)!=0)
  {
    after_sam2 <- after_sam1[after_sam1$V6!="*",]
    colnames(after_sam2)[1] <- c("readname");colnames(after_sam2)[4] <- "matchloc";colnames(after_sam2)[6] <- c("CIGAR");colnames(after_sam2)[18] <- c("NM")
    after_error <- matrix(0,nrow = nrow(after_sam2),ncol = 4)
    rownames(after_error) <- after_sam2[,1]
    colnames(after_error) <- c("M","XM","match","dismatch")
    for(i in seq(1,nrow(after_sam2)))
    {
      CIGAR <-  unlist(strsplit(after_sam2[i,6],split = '-'))
      M <- which(CIGAR=="M")
      if(length(M)==1)
      {after_error[i,1] <- as.numeric(CIGAR[[M-1]])}
      else
      {
        for(j in seq(1,length(M)))
        {after_error[i,1] <- after_error[i,1]+as.numeric(CIGAR[[M[[j]]-1]])}
      }
      V14 <- unlist(strsplit(after_sam2[i,14],split = ':'))
      if("XM" %in% V14)
      {after_error[i,2] <- as.numeric(V14[[3]])}
      V15 <- unlist(strsplit(after_sam2[i,15],split = ':'))
      if("XM" %in% V15)
      {after_error[i,2] <- as.numeric(V15[[3]])}
      after_error[i,4] <- after_error[i,2]
      after_error[i,3] <- after_error[i,1]-after_error[i,4]
    }
  }
  else{after_error <- matrix(0,nrow = 1,ncol = 4)}
  return(after_error)
}
#beforeerror
beforeerror <- function(x,exon_start){
  #before
  before_sam1 <- x[x$V4<exon_start,]
  if(nrow(before_sam1)!=0)
  {
    before_sam2 <- before_sam1[before_sam1$V6!="*",]
    colnames(before_sam2)[1] <- c("readname");colnames(before_sam2)[4] <- "matchloc";colnames(before_sam2)[6] <- c("CIGAR");colnames(before_sam2)[18] <- c("NM")
    before_error <- matrix(0,nrow = nrow(before_sam2),ncol = 4)
    rownames(before_error) <- before_sam2[,1]
    colnames(before_error) <- c("M","XM","match","dismatch")
    for(i in seq(1,nrow(before_sam2)))
    {
      CIGAR <-  unlist(strsplit(before_sam2[i,6],split = '-'))
      M <- which(CIGAR=="M")
      if(length(M)==1)
      {before_error[i,1] <- as.numeric(CIGAR[[M-1]])}
      else
      {
        for(j in seq(1,length(M)))
        {before_error[i,1] <- before_error[i,1]+as.numeric(CIGAR[[M[[j]]-1]])}
      }
      V14 <- unlist(strsplit(before_sam2[i,14],split = ':'))
      if("XM" %in% V14)
      {before_error[i,2] <- as.numeric(V14[[3]])}
      V15 <- unlist(strsplit(before_sam2[i,15],split = ':'))
      if("XM" %in% V15)
      {before_error[i,2] <- as.numeric(V15[[3]])}
      before_error[i,4] <- before_error[i,2]
      before_error[i,3] <- before_error[i,1]-before_error[i,4]
    }
  }
  else{before_error <- matrix(0,nrow = 1,ncol = 4)}
  return(before_error)
}
#depth
depthfile <- function(x){
  depth <- x[,-1]
  a <- which(depth[,1]==exon_start)
  b <- which(depth[,1]==exon_end)
  depth <- depth[a:b,]
  colnames(depth) <- c("location","depth")
  return(depth)
}
#MDtag_error
find_AGCT<-function(index,x){
  a<-c()
  k <- 1
  for (i in index)
  {
    if(i %in% x)
    { 
      b <- grep(i,x,value=F)
      if(length(b) <2)
      {a[k] <- b;k <- k+1}
      else
      {
        for(j in seq(1,length(b)))
        {a[k] <- b[[j]];k <- k+1}
      }
    }
  }
  a <- sort(a)
  return(a)
}
MD_tag <- function(x,after_error){
  after_sam1 <- x[x$V4>=exon_start,]
  after_sam2 <- after_sam1[after_sam1$V6!="*",]
  after_mismatch_sam <- after_sam2[after_error[,2]!=0,]
  if(nrow(after_mismatch_sam)!=0)
  {
    after_mismatch_sam <- after_mismatch_sam[,c(1,4,6,10,18,19)]
    after_mismatch_sam <- cbind(after_mismatch_sam,0)
    for(i in seq(1,nrow(after_mismatch_sam)))
    {
      a <- unlist(strsplit(after_mismatch_sam[i,5],split = ':'))
      b <- unlist(strsplit(after_mismatch_sam[i,6],split = ':'))
      if("MD" %in% a)
      {after_mismatch_sam[i,7] <- after_mismatch_sam[i,5]}
      else
      {after_mismatch_sam[i,7] <- after_mismatch_sam[i,6]}
    }
    after_mismatch_sam <- after_mismatch_sam[,-c(5,6)]
    after_mismatch_sam <- cbind(after_mismatch_sam,0)
    colnames(after_mismatch_sam)[[5]] <- "V19"
    after_mismatch_sam$V19 <- gsub("MD:Z:","",after_mismatch_sam$V19)
    colnames(after_mismatch_sam) <- c("readname","start","CIGAR","reads","MD:Z:","MD_num")
    find_AGCT<-function(index,x){
      a<-c()
      k <- 1
      for (i in index)
      {
        if(i %in% x)
        { 
          b <- grep(i,x,value=F)
          if(length(b) <2)
          {a[k] <- b;k <- k+1}
          else
          {
            for(j in seq(1,length(b)))
            {a[k] <- b[[j]];k <- k+1}
          }
        }
      }
      a <- sort(a)
      return(a)
    }
    #x<-a
    #index<- c("A","D","H")
    #f(index,x)
    index <- c("A","G","C","T","N")
    #after_mismatch_sam[i,4] <- "23G^AGT23G^AGT32G^CG44G"
    for(i in seq(1,nrow(after_mismatch_sam)))
    {
      if(grepl("\\^",after_mismatch_sam[i,5]))
      {
        temporary <-  unlist(strsplit(after_mismatch_sam[i,5],split = ''))
        loc <- grep("\\^",temporary)
        agct <- find_AGCT(index,temporary)
        d <- length(temporary)
        if(length(loc)==1)
        {
          e <- loc+1
          while(e %in% agct)
          {
            e <- e+1
          }
          a <- str_sub(after_mismatch_sam[i,5],1,loc-1)
          b <- str_sub(after_mismatch_sam[i,5],e,d)
          after_mismatch_sam[i,5] <- paste(a,b,sep = "")
        }
        else
        {
          a <- str_sub(after_mismatch_sam[i,5],1,loc[[1]]-1)
          for(k in seq(1,length(loc)))
          {
            e <- loc[[k]]+1
            while(e %in% agct)
            {
              e <- e+1
            }
            if(k != length(loc))
            {b <- str_sub(after_mismatch_sam[i,5],e,loc[[k+1]]-1)}
            else {b <- str_sub(after_mismatch_sam[i,5],e,d)}
            a <- paste(a,b,sep = "")
          }
          after_mismatch_sam[i,5] <- a
        }
      }
    }
    
    for(i in seq(1,nrow(after_mismatch_sam)))
    {
      MD <- after_mismatch_sam[i,5]
      MD_tag <- unlist(strsplit(MD,split = ''))
      agct <- find_AGCT(index,MD_tag)
      after_mismatch_sam[i,6] <- MD_tag[[1]]
      if(length(agct)==1)
      {
        if(agct[[1]]==2)
        {after_mismatch_sam[i,6] <- as.numeric(MD_tag[[1]])}
        else
        {
          for(j in seq(2,agct-1))
          {
            after_mismatch_sam[i,6] <- as.numeric(paste(after_mismatch_sam[i,6],MD_tag[[j]],sep = ""))
          }
        }
      }
      else
      {
        tep <- c()
        n <- 1
        for(j in seq(1,length(agct)))
        {
          if(j ==1)
          {
            if(agct[[j]]!=2)
            {
              num1 <- after_mismatch_sam[i,6]
              for(k in seq(2,agct[[j]]-1))
              {
                num1 <- as.numeric(paste(num1,MD_tag[[k]],sep = ""))
              }
              tep[n] <- as.numeric(num1)
              n <- n+1
            }
            else
            {tep[n] <- MD_tag[[1]];n <- n+1}
          }
          else
          {
            m <- agct[[j-1]]+1
            if(m+1!=agct[[j]])
            {
              num2 <- MD_tag[[m]]
              for(k in seq(m+1,agct[[j]]-1))
              {
                num2 <- as.numeric(paste(num2,MD_tag[[k]],sep = ""))
              }
              tep[n] <- as.numeric(num2)
              n <- n+1
            }
            else
            {tep[n] <- as.numeric(MD_tag[[m]]);n < n+1}
          }
        }
        dd <- tep[1]
        for(d in seq(2,length(tep)))
        {dd <- paste(dd,tep[d],sep = "-")}
        after_mismatch_sam[i,6] <- dd
      }
    }
    for(i in seq(1,nrow(after_mismatch_sam)))
    {
      CIGAR <- unlist(strsplit(after_mismatch_sam[i,3],split = '-'))
      MD_tag <- as.numeric(unlist(strsplit(after_mismatch_sam[i,6],split = '-')))
      if("N" %in% CIGAR)
      {
        a <- which(CIGAR=="N")
        e <- 0
        for(k in seq(1,length(a)))
        {
          ee <- a[[k]]
          e <- e+as.numeric(CIGAR[[ee-3]])
        }
        if(length(MD_tag)==1)
        {
          if(MD_tag+1 < (e+1))
          {
            b <- which(rownames(after_error)==after_mismatch_sam[i,1])
            after_error[b,2] <- after_error[b,2]-1
          }
        }
        else
        {
          if(MD_tag[[1]]+1 <  (e+1))
          {
            b <- which(rownames(after_error)==after_mismatch_sam[i,1])
            after_error[b,2] <- after_error[b,2]-1
          }
          c <- MD_tag[[1]]+1
          for(j in seq(2,length(MD_tag)))
          {
            d <- c+MD_tag[[j]]+1
            if(d <  (e+1))
            {
              b <- which(rownames(after_error)==after_mismatch_sam[i,1])
              after_error[b,2] <- after_error[b,2]-1
            }
          }
        }
      }
      else
      {
        e <- exon_start-after_mismatch_sam[i,2]
        if(length(MD_tag)==1)
        {
          if(MD_tag+1 < (e+1))
          {
            b <- which(rownames(after_error)==after_mismatch_sam[i,1])
            after_error[b,2] <- after_error[b[[1]],2]-1
          }
        }
        else
        {
          if(MD_tag[[1]]+1 <  (e+1))
          {
            b <- which(rownames(after_error)==after_mismatch_sam[i,1])
            after_error[b,2] <- after_error[b,2]-1
          }
          c <- MD_tag[[1]]+1
          for(j in seq(2,length(MD_tag)))
          {
            d <- c+MD_tag[[j]]+1
            if(d <  (e+1))
            {
              b <- which(rownames(after_error)==after_mismatch_sam[i,1])
              after_error[b,2] <- after_error[b,2]-1
            }
          }
        }
      }
      
    }
    after_mismatch_sam <- cbind(after_mismatch_sam,0)
    for(i in seq(1,nrow(after_mismatch_sam)))
    {
      CIGAR <- unlist(strsplit(after_mismatch_sam[i,3],split = '-'))
      span <- exon_end - after_mismatch_sam[i,2]
      for(j in seq(2,length(CIGAR),2))
      {
        if(CIGAR[[j]]!="S")
        {
          if(CIGAR[[j]]=="M")
          {
            if(as.numeric(CIGAR[[j-1]]) < span)
            {scribe <- as.numeric(CIGAR[[j-1]])}else{scribe <- span}
            MD_error <-  unlist(strsplit(after_mismatch_sam[i,6],split = '-'))
            MD_distance <- 0
            for(t in seq(1,length(MD_error)))
            {
              MD_distance <- MD_distance +as.numeric(MD_error[[t]])+1 
              if(MD_distance >scribe)
              {after_mismatch_sam[i,7] <- as.numeric(after_mismatch_sam[i,7])+1}
            }
          }
          else{
            a <- 0
            if(j == 2)
            {a <- a+as.numeric(CIGAR[[j-1]])}
            if(j >2)
            {
              for(k in seq(1,j-1,2))
              {a <- a+as.numeric(CIGAR[[j-1]]) }
            }
            if(a > span)
            {break}
          }
        }
      }
    }
  }
  else(after_mismatch_sam <- matrix(0,nrow = 1,ncol = 7))
  return(after_mismatch_sam)
}

mistakebase <- function(x,exon_start,exon_end,after_error,fasta_seq,fildepth){
  after_sam1 <- x[x$V4>=exon_start,]
  after_sam2 <- after_sam1[after_sam1$V6!="*",]
  mismatch_sam <- after_sam2[after_error[,2]!=0,]
  error_minus_base <- 0
  mutbase <- 0
  if(nrow(mismatch_sam)==0){gene_TR <- 0}
  if(nrow(mismatch_sam)!=0)
  {  mismatch_sam <- mismatch_sam[,c(1,4,6,10,18,19)]
  colnames(mismatch_sam) <- c("readname","start","CIGAR","sequence")
  error_base <- matrix(0,nrow = nrow(mismatch_sam),ncol = ncol(fasta_seq))
  rownames(error_base) <- mismatch_sam[,1]
  error_base <- rbind(fasta_seq,error_base)
  rownames(error_base)[[1]] <- "Refseq"
  colnames(error_base) <- seq(exon_start,exon_end)
  if(nrow(mismatch_sam)>1)
  {
    for(i in seq(2,nrow(mismatch_sam)))
    {
      reads <- unlist(strsplit(mismatch_sam[i-1,4],split = ''))
      read_start <- mismatch_sam[i-1,2]
      a <- 1 
      b <- which(colnames(error_base)==read_start)
      CIGAR <- unlist(strsplit(mismatch_sam[i-1,3],split = '-'))
      for(j in seq(2,length(CIGAR),2))
      {
        if(CIGAR[[j]]=="S")
        {a <- a+as.numeric(CIGAR[[j-1]])}
        if(CIGAR[[j]]=="N")
        {
          c <- exon_start +as.numeric(CIGAR[[j-1]])
          if(c > exon_end)
          {break}
          else
          {b <- b+as.numeric(CIGAR[[j-1]])}
        }
        if(CIGAR[[j]]=="M")
        {
          for(k in seq(1,as.numeric(CIGAR[[j-1]])))
          {
            if(a ==1)
            {
              if(k < length(reads)+1)
              {read <- reads[[k]]}
              if(b+k-1 < ncol(error_base)+1)
              {
                Ref <- error_base[1,b+k-1]
                if(read != Ref)
                {
                  error_base[i,b+k-1] <- read
                  if(k<3||k>(length(reads)-3))
                  {error_minus_base <- error_minus_base+1;error_base[i,b+k-1] <- 0}
                }
              }
            }
            else
            {
              if(a+k-1 < length(reads)+1)
              {read <- reads[[a+k-1]]}
              if(b+k-1 < ncol(error_base)+1)
              {
                Ref <- error_base[1,b+k-1]
                if(read != Ref)
                {
                  error_base[i,b+k-1] <- read
                  if(a+k<3||a+k>(length(reads)-3))
                  {error_minus_base <- error_minus_base+1;error_base[i,b+k-1] <- 0}
                }
              }
            }
          }
          a <- a+as.numeric(CIGAR[[j-1]])
          b <- b+as.numeric(CIGAR[[j-1]])
        }
        if(CIGAR[[j]]=="I")
        {a <- a+as.numeric(CIGAR[[j-1]])}
        if(CIGAR[[j]]=="D")
        {b <- b+as.numeric(CIGAR[[j-1]])}
      }
    }
  }
  else
  {
    i <- 1
    reads <- unlist(strsplit(mismatch_sam[i,4],split = ''))
    read_start <- mismatch_sam[i,2]
    a <- 1 
    b <- which(colnames(error_base)==read_start)
    CIGAR <- unlist(strsplit(mismatch_sam[i,3],split = '-'))
    for(j in seq(2,length(CIGAR),2))
    {
      if(CIGAR[[j]]=="S")
      {a <- a+as.numeric(CIGAR[[j-1]])}
      if(CIGAR[[j]]=="M")
      {
        for(k in seq(1,as.numeric(CIGAR[[j-1]])))
        {
          if(a ==1)
          {
            if(k < length(reads)+1)
            {read <- reads[[k]]}
            if(b+k-1 < ncol(error_base)+1)
            {
              Ref <- error_base[1,b+k-1]
              if(read != Ref)
              {
                error_base[i,b+k-1] <- read
                if(k<3||k>(length(reads)-3))
                {error_minus_base <- error_minus_base+1;error_base[i,b+k-1] <- 0}
              }
            }
          }
          else
          {
            if(a+k-1 < length(reads)+1)
            {read <- reads[[a+k-1]]}
            if(b+k-1 < ncol(error_base)+1)
            {
              Ref <- error_base[1,b+k-1]
              if(read != Ref)
              {
                error_base[i,b+k-1] <- read
                if(a+k<3||a+k>(length(reads)-3))
                {error_minus_base <- error_minus_base+1;error_base[i,b+k-1] <- 0}
              }
            }
          }
        }
      }
      if(CIGAR[[j]]=="I")
      {a <- a+as.numeric(CIGAR[[j-1]])}
      if(CIGAR[[j]]=="D")
      {b <- b+as.numeric(CIGAR[[j-1]])}
    }
  }
  #mutbase
  sum_errorbase <- matrix(0,nrow = 4,ncol = ncol(fasta_seq))
  rownames(sum_errorbase) <- c("A","G","C","T")
  sum_errorbase <- rbind(fasta_seq,sum_errorbase)
  rownames(sum_errorbase)[[1]] <- "Refseq"
  colnames(sum_errorbase) <- seq(exon_start,exon_end)
  for(i in seq(1,ncol(error_base)))
  {
    a <- which("A" == error_base[2:nrow(error_base),i])
    sum_errorbase[2,i] <- length(a)
    g <- which("G" == error_base[2:nrow(error_base),i])
    sum_errorbase[3,i] <- length(g)
    c <- which("C" == error_base[2:nrow(error_base),i])
    sum_errorbase[4,i] <- length(c)
    t <- which("T" == error_base[2:nrow(error_base),i])
    sum_errorbase[5,i] <- length(t)
  }
  mut <- t(sum_errorbase)
  mut[,1] <- rownames(mut)
  mut_depth <- as.data.frame(merge(fildepth,mut,by.x = "location",by.y = "Refseq"))
  mut_depth$fil <- mut_depth$depth*0.5
  for(i in seq(1,nrow(mut_depth)))
  {
    for(j in seq(3,6))
    {if(as.numeric(mut_depth[i,j])>mut_depth[i,7])
    {mutbase <-mutbase +as.numeric(mut_depth[i,j]) }
    }
  }
  #pergene
  pergene <- matrix(1,nrow = ncol(sum_errorbase),ncol = 1)
  rownames(pergene) <- colnames(sum_errorbase)
  for(i in seq(1,ncol(sum_errorbase)))
  {
    if(sum(as.numeric(sum_errorbase[2:5,i]))!= 0 )
    {
      pergene[i,1] <- 1-sum(as.numeric(sum_errorbase[2:5,i]))/nrow(mismatch_sam)
    }
  }
  gene_TR <- 1-prod(as.numeric(pergene[,1]))
  }
  return(list(error_minus_base,mutbase,sum_errorbase,gene_TR))
}
totalminbase <- function(after_error,before_error){
  length(after_error[,1])*4+length(before_error[,1])*4
}
#error_trand
codon_matrix <- function(sum_errorbase,changecodon,c,gene_bed_4,exon_start){
  codon <- expand.grid(
    a = c("A","G","C","T"),
    b = c("A","G","C","T"),
    c = c("A","G","C","T")) %>% mutate(codon = paste(a,b,c,sep = ""))
  codon_codon <- matrix(0,nrow = 64,ncol = 64)
  rownames(codon_codon) <- codon$codon
  colnames(codon_codon) <- codon$codon
  codon_codon <- as.data.frame(codon_codon)
  ref_start <- as.numeric(gene_bed_4[c,5])
  loc_start <- exon_start- ref_start +1
  loc_start_codon <- loc_start%%3
  if(loc_start_codon==0){t <- 2}
  if(loc_start_codon==1){t <- 1}
  if(loc_start_codon==2){t <- 3}
  minuscodon <- (ncol(sum_errorbase))%%3
  for(i in seq(t,ncol(sum_errorbase)-minuscodon,3))
  {
    if(i != ncol(sum_errorbase) && i != ncol(sum_errorbase)-1)
    {
      if(sum(as.numeric(sum_errorbase[2:5,i])) != 0)
      {
        error_loc <- which(as.numeric(sum_errorbase[2:5,i])!= 0) +1
        for(n in seq(1,length(error_loc)))
        {
          tt <- error_loc[[n]]
          codon_orignal <- paste(sum_errorbase[1,i],sum_errorbase[1,i+1],sum_errorbase[1,i+2],sep = "")
          codon_after <- paste(rownames(sum_errorbase)[[tt]],sum_errorbase[1,i+1],sum_errorbase[1,i+2],sep = "")
          a <- which(rownames(codon_codon)==codon_orignal)
          b <- which(colnames(codon_codon)==codon_after)
          codon_codon[a,b] <- codon_codon[a,b] +1
        }
      }
      if(sum(as.numeric(sum_errorbase[2:5,i+1])) != 0)
      {
        error_loc <- which(as.numeric(sum_errorbase[2:5,i+1])!= 0) +1
        for(n in seq(1,length(error_loc)))
        {
          tt <- error_loc[[n]]
          codon_orignal <- paste(sum_errorbase[1,i],sum_errorbase[1,i+1],sum_errorbase[1,i+2],sep = "")
          codon_after <- paste(sum_errorbase[1,i],rownames(sum_errorbase)[[tt]],sum_errorbase[1,i+2],sep = "")
          a <- which(rownames(codon_codon)==codon_orignal)#before
          b <- which(colnames(codon_codon)==codon_after)#after
          codon_codon[a,b] <- codon_codon[a,b] +1
        }
      }
      if(sum(as.numeric(sum_errorbase[2:5,i+2])) != 0)
      {
        error_loc <- which(as.numeric(sum_errorbase[2:5,i+2])!= 0) +1
        for(n in seq(1,length(error_loc)))
        {
          tt <- error_loc[[n]]
          codon_orignal <- paste(sum_errorbase[1,i],sum_errorbase[1,i+1],sum_errorbase[1,i+2],sep = "")
          codon_after <- paste(sum_errorbase[1,i],sum_errorbase[1,i+1],rownames(sum_errorbase)[[tt]],sep = "")
          a <- which(rownames(codon_codon)==codon_orignal)
          b <- which(colnames(codon_codon)==codon_after)
          codon_codon[a,b] <- codon_codon[a,b] +1
        }
      }
    }
  }
  for(i in seq(1,nrow(codon_codon)))
  {
    if(sum(codon_codon[i,]) != 0)
    {
      loc_after_codon <- which(codon_codon[i,] !=0)
      for(j in seq(1,length(loc_after_codon)))
      {
        m <- loc_after_codon[[j]]
        codon_ch <- paste(rownames(codon_codon)[[i]],colnames(codon_codon)[[m]],sep = "-")
        loc_ch_codon <- which(rownames(changecodon)==codon_ch)
        changecodon[loc_ch_codon,c] <- changecodon[loc_ch_codon,c] +as.numeric(codon_codon[i,m])
      }
    }
  }
  changecodon <- as.numeric(changecodon[,c])
  return(changecodon)
}
errorloc <- function(sum_errorbase,geneid,exonloc){
    
  a <- sum_errorbase[2:5, ]
  a <- apply(a, 2, as.numeric)  
  b <- sum(colSums(a) != 0) 
  if(b != 0)
  {
    gene_errorloc <- matrix(0,nrow = b,ncol = 5)
    colnames(gene_errorloc) <- c("CHROM","POS","REF","ALT","reads")
    gene_errorloc[,1] <- exonloc[[1]]
    for(i in seq(1,b))
    {
      loc_1 <- which(colSums(a) != 0)[[i]]
      gene_errorloc[i,2] <- colnames(a)[[loc_1]]
      gene_errorloc[i,3] <- sum_errorbase[1,loc_1]
      loc_2 <- which.max(a[,loc_1])
      gene_errorloc[i,4] <- rownames(sum_errorbase)[[loc_2+1]]
      gene_errorloc[i,5] <- sum_errorbase[loc_2+1,loc_1]
    }
    gene_errorloc <- as.data.frame(gene_errorloc)
  }
  else{gene_errorloc <- 0}
  return(gene_errorloc)
}


transfer <- c("T-G","T-C","T-A","G-T","G-C","G-A","C-T","C-G","C-A","A-T","A-G","A-C")
codons <- c("AAA","GAA","CAA","TAA","AGA","GGA","CGA","TGA","ACA","GCA","CCA","TCA","ATA","GTA","CTA","TTA","AAG","GAG","CAG","TAG","AGG","GGG","CGG","TGG","ACG","GCG","CCG",
            "TCG","ATG","GTG","CTG","TTG","AAC","GAC","CAC","TAC","AGC","GGC","CGC","TGC","ACC","GCC","CCC","TCC","ATC","GTC","CTC","TTC","AAT","GAT","CAT","TAT","AGT","GGT",
            "CGT","TGT","ACT","GCT","CCT","TCT","ATT","GTT","CTT","TTT")
codon_change <- expand.grid(
  a = codons,
  b = codons)%>% mutate(change = paste(a,b,sep = "-"))
changecodon <- matrix(0,nrow = nrow(codon_change),ncol = length(gene_file[,1]))
rownames(changecodon) <- codon_change$change
colnames(changecodon) <- gene_file[,1]
changecodon <- as.data.frame(changecodon)
TR_pergene <- data.frame(err1=rep(0,length(gene_file[,1])))
rownames(TR_pergene) <- gene_file[,1]
gene_errorloc_list <- list()
for(x in seq(1,length(sam_filelist)))
{
  sum_after_match <-0
  sum_before_match <-0
  sum_after_error <-0
  sum_error_minus_base <- 0
  sum_mutbase <- 0
  sum_total_minus_base <- 0
  sum_MD_outscripe_error <- 0
  samloc <- paste(chr_wd,sam_filelist[[x]],sep = "/")
  samloc <- paste(samloc,"\r",sep = "")
  sam <- read.table(samloc,header = F,fill  = T,col.names = paste("V", 1:23, sep = ""))
  fastaloc <- paste(fa_wd,fasta_filelist[[x]],sep = "/")
  fastaloc <- paste(fastaloc,"\r",sep = "")
  fasta <- read.fasta(fastaloc)
  fasta_seq <- toupper(t(as.matrix(unlist(strsplit(fasta[1,2],split = '')))))
  exonloc <- unlist(strsplit(sam_filelist[[x]],split = '\\.t'))
  exonloc <- unlist(strsplit(exonloc[[1]],split = '_'))
  exon_start <- as.numeric(exonloc[[2]])
  exon_end <- as.numeric(exonloc[[3]])
  fildepth <- depthfile(depth)
  sam <- sam_CIGAR(sam)
  after_error <- aftererror(sam,exon_start)
  sum_after_match <- sum_after_match+sum(after_error[,1])
  before_error <- beforeerror(sam,exon_start)
  sum_before_match <- sum_before_match+sum(before_error[,1])
  sum_after_error <- sum_after_error+sum(after_error[,4])
  after_mismatch_sam <- MD_tag(sam,after_error)
  MD_outscripe_error <- sum(after_mismatch_sam[,7])
  sum_MD_outscripe_error <- sum_MD_outscripe_error+MD_outscripe_error
  mistake_Base <- mistakebase(sam,exon_start,exon_end,after_error,fasta_seq,fildepth)
  sum_error_minus_base <- sum_error_minus_base+mistake_Base[[1]]
  sum_mutbase <- sum_mutbase+mistake_Base[[2]]
  total_minus_base <- totalminbase(after_error,before_error)
  sum_total_minus_base <- sum_total_minus_base +total_minus_base 
  sum_errorbase <- mistake_Base[[3]]
  TR_pergene[x,1] <- mistake_Base[[4]]
  changecodon[,x] <- codon_matrix(sum_errorbase,changecodon,x,gene_bed_4,exon_start)
  gene_errorloc <- errorloc(sum_errorbase,geneid,exonloc)
  gene_errorloc_list[[x]] <- gene_errorloc
  write.csv(changecodon,"changecodon.csv",quote = F)
  write.csv(TR_pergene,"TR_pergene.csv",quote = F)
  names(gene_errorloc_list)[[x]] <- geneid[x,1]
  save(gene_errorloc_list,file="gene_errorloc_list.RData")
  total_transfererror[x,1] <- (sum_after_match+sum_before_match-sum_total_minus_base)
  total_transfererror[x,2] <- (sum_after_error-sum_mutbase-sum_error_minus_base-sum_MD_outscripe_error)
  total_transfererror[x,3] <- total_transfererror[x,2]/total_transfererror[x,1]
  history_wd <- paste(chr_wd,"error_history.RData",sep = "/")
  save(sum_after_match,sum_before_match,sum_total_minus_base,sum_after_error,sum_mutbase,sum_error_minus_base,file = history_wd)
}
write.csv(total_transfererror,"gene_transfererror.csv",quote = F)
write.csv(changecodon,"drosophila_changecodon.csv",quote = F)
write.csv(TR_pergene,"drosophila_TR_pergene.csv",quote = F)