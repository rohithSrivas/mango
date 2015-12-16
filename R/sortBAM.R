#' sort BAM file using samtools
#'
#' This function sorts reads in a BAM file by either name or coordinates using samtools
#'  
#' @param bamInputFile full path to the input BAM file
#' @param bamOutputFile full path to the output BAM file
#' @param samtools path full path to samtools
#' @param by.name if true, sorting will be done by read name. otherwise by coordinate.
#' @param num.threads The number of parallel search threads to use for alignment
#' @param verbose if true, the command invoked is printed to the console
#' @export
#' 
sortBAM <- function(bamInputFile,bamOutputFile,samtools,by.name=TRUE,num.threads=1,verbose=TRUE)
{
	# setup temporary file names
	temp.extension <- ifelse(dirname(bamInputFile)==".","temp_sorted",paste(dirname(bamInputFile),"/temp_sorted",sep=""))

	# setup samtools sort parameters
	samtoolsoptions=paste("sort -O 'bam' -T ",temp.extension,"-@ ",num.threads,sep="")
	if (by.name == TRUE)
	{
	samtoolsoptions=paste("sort -n -O 'bam' -T ",temp.extension," -@ ",num.threads,sep="")
	}

	# form full command
	samtoolscommand = paste (samtools, samtoolsoptions, bamInputFile,">",bamOutputFile)

	# print command if desirec
	if (verbose ==TRUE)
	{
	print (samtoolscommand)
	}

	# execute command
	system(samtoolscommand)
}

