\name{runAlign2RawSignal}
\alias{runAlign2RawSignal}
\title{Assess ChIP quality using phantompeakquals toolkit}
\usage{
runPhantomPeakQual(input.bam,output.results.file,output.plot.file,path.to.phantom.script,num.threads=1)
}
\arguments{
  \item{input.bam}{full path to input bam file}

  \item{output.results.file}{output file containing ENCODE chip quality metrics and est. fragment length}

  \item{output.results.plot.file}{output plot containing results of strand cross-correlation}

  \item{path.to.phantom.script}{path to run_spp.R [script for running phantom peak quals]}

  \item{num.threads}{Number of parallel threads to utilize for computing strand cross-correlation}
}
\description{
This function assesses ChIP quality using the phantompeakquals toolkit; it produces estimates of RSC, NSC, and estimated fragment length.
}