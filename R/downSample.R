#' downsample BAM file using picard tools
#'
#' This function downsamples a BAM file using the picard tools
#'  
#' @param bamInputFile full path to the input BAM file
#' @param bamOutputFile full path to the output BAM file
#' @param picardToolsPath full path to picard tools
#' @param prob fraction of reads to retain
#' @param verbose if true, the command invoked is printed to the console
#' @export
#'
downSampleBam <- function(bamInputFile, bamOutputFile, picardToolsPath,prob,verbose)
{
	#Setup command
	command <- paste("java -jar -Xmx2g ",picardToolsPath,"DownsampleSam.jar I=",bamInputFile," O=",bamOutputFile," P=",prob," R=123",sep="")
	
	# print command if desired
	if (verbose ==TRUE) {
		print (command)
	}

	# execute command
	system(command)
}