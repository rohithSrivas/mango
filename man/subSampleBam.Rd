\name{subSampleBam}
\alias{subSampleBam}
\title{Subsample BAM files produced from the Mango pipeline}
\usage{
subSampleBam(bamFile1, bamFile2, outputBamFile1, outputBamFile2, sampleFrac) {
}
\arguments{
  \item{bamFile1}{full path to first paired end BAM file}
  
  \item{bamFile2}{full path to second paired end BAM file}
  
  \item{outputBamFile1}{full path to first subsampled paired end output BAM file}

  \item{outputBamFile2}{full path to second subsampled paired end output BAM file}

  \item{sampleFrac}{fraction of reads to keep [must be decimal between 0 and 1]}
}
\description{
This function subsamples both paired end BAM files to the desired percentage.
}