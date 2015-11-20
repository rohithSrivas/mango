#' aligns reads using bowtie
#'
#' This function aligns reads using bowtie
#'  
#' @param fastq full path to fastq file
#' @param bowtiepath full path to bowtie
#' @param bowtieref full path to bowtie reference
#' @param samtools path full path to samtools
#' @param arguments arguments to give to bowtie
#' @param verbose boolean whether or not to print command
#' @param nlines the number of lines to look at to determine the scoring method
#' @param num.threads The number of parallel search threads to use for alignment
#' @export
#' 
alignBowtie <- function(fastq,output,bowtiepath,bowtieref,samtoolspath,
                        shortreads,num.threads=1,verbose=TRUE)
{
  
  # choose alignment parameters
  # note- "-m 1" ensures that only uniquely mapped reads are reported.
  bowtievar=paste("-S -v 0 -k 1 --chunkmbs 500 -S --mapq 40 -m 1 -p ",num.threads,sep="")
  if (shortreads == FALSE)
  {
    bowtievar=paste("-S -n 2 -l 50 -k 1 --chunkmbs 500 -S --mapq 40 -m 1 --best -p ",num.threads,sep="")
  }
  
  # determine the illumina score encoding
  illuminascore = findScore(fastq,nlines =10000) 
  
  # create the output to samtools pipe
  samtoolsPipe <- " | samtools view -b - > "
  
  # form full command
  bowtiecommand = paste (bowtiepath, bowtieref, fastq, illuminascore, bowtievar, samtoolsPipe)
  
  # print command if desirec
  if (verbose ==TRUE)
  {
    print (bowtiecommand)
  }
  
  # execute command
  system(bowtiecommand)
}

