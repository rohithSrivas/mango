\name{alignBowtie}
\alias{alignBowtie}
\title{aligns reads using bowtie}
\usage{
alignBowtie(fastq, output, bowtiepath, bowtieref, shortreads, num.threads=1, verbose = TRUE)
}
\arguments{
  \item{fastq}{full path to fastq file}

  \item{bowtiepath}{full path to bowtie}

  \item{bowtieref}{full path to bowtie reference}

  \item{arguments}{arguments to give to bowtie}

  \item{verbose}{boolean whether or not to print command}

  \item{nlines}{the number of lines to look at to determine
  the scoring method}
  
  \item{num.threads}{The number of parallel search threads to use for alignment}
}
\description{
This function aligns reads using bowtie
}

