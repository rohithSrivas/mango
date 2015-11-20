#' aligns reads using bowtie
#'
#' This function assesses ChIP quality using the phantompeakquals toolkit
#'  
#' @param input.bam full path to input bam file
#' @param output.results.file output file containing ENCODE chip quality metrics and est. fragment length
#' @param output.results.plot.file output plot containing results of strand cross-correlation
#' @param path.to.phantom.script path to run_spp.R [script for running phantom peak quals]
#' @param num.threads Number of parallel threads to utilize for computing strand cross-correlation
#' @export
#' 
runPhantomPeakQual <- function(	input.bam,
								output.results.file,
								output.plot.file,
								path.to.phantom.script,
                        		num.threads=1)
{
  
  # setup command
  qual.command <- paste(	path.to.phantom.script,
	  						" -c=",input.bam,
							" -savp=",output.plot.file,
							" -out=",output.results.file,
							" -p=",num.threads,sep="")
  
  
  # print command invoked
  print (qual.command)
  
  # execute command
  system(paste("Rscript ",qual.command,sep=""))
}

