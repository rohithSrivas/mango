\name{runAlign2RawSignal}
\alias{runAlign2RawSignal}
\title{Produce signal track data using align2rawsignal}
\usage{
runAlign2RawSignal(input.bam,temp.filtered.bam,path.to.bamutils,output.mat.file,temp.output.bedgraph.file,output.bw.file,path.to.mcr,path.to.align2rawsignal,bedtoolsgenome,fragLength,chrDir,mapDir,verbose)
}
\arguments{
  \item{input.bam}{full path to input bam file}
  
  \item{temp.filtered.bam}{full path to temporary bam file filtered for very short or very long reads}
  
  \item{path.to.bamutils}{full path to the NSGUtils bamutils program}

  \item{output.mat.file}{output file containing signal track data in MATLAB mat format}

  \item{temp.output.bedgraph.file}{temporary output file which contains the bedgraph file [this file will be deleted]}

  \item{output.bw.file}{output file containing signal track data in BigWig file}

  \item{path.to.mcr}{full path to the matlab compile runtime [MCR v17]}

  \item{path.to.align2rawsignal}{full path to the align2rawsingal binary file}
  
  \item{bedtoolsgenome}{path to the bedtools genome file (name of contigs and the size of each contig)}
  
  \item{fragLength}{the average fragment length estimated from phantompeakqual toolkit}
  
  \item{chrDir}{path to directory containing individual chromsome fasta files}
  
  \item{mapDir}{path to directory containing mappability files}
  
  \item{verbose}{Flag to determine whether to output the executed commands to file}
}
\description{
This function produces signal track data in both MAT and BW formats using align2rawsignal
}