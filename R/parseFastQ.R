parseFASTQs <- function(	fastq1,
							fastq2,
							basename,
							minlength,
							maxlength,
							keepempty,
							linker1,
							linker2,)
{
	#Step 1: Check if file is zipped or not
    is.fastq1.zipped <- any(str_detect(fastq1,c("\\.gz","\\.gzip")))
    is.fastq2.zipped <- any(str_detect(fastq2,c("\\.gz","\\.gzip")))
  
	if(is.fastq1.zipped!=is.fastq2.zipped) {
  	  print ("Both FASTQ files need to be either .GZIP or not zipped!")
  	  stop()
    }
    print ("finding linkers")
	
	#Step 2: Split and find linkers
    #If zipped, no need to unzip, just process it
	parsingresults = NULL
    if(is.fastq1.zipped) {
    	parsingresults = parseFastq_gzip	(	fastq1=fastq1,
                								fastq2=fastq2,
	                					  		basename = basename,
	                					  		minlength = minlength,
	                					  		maxlength = maxlength, 
	                					  		keepempty = keepempty,
	  											verbose=TRUE,
	                					  		linker1=linker1,
	                					  		linker2=linker2)
	} else {
	  	parsingresults = parseFastq			(	fastq1=fastq1,
			              						fastq2=fastq2,
			              					  	basename = basename,
			              					  	minlength = minlength,
			              					  	maxlength = maxlength, 
			              					  	keepempty = keepempty,
			              					  	linker1=linker1,
			              					  	linker2=linker2)
	}
	
	return(parsingresults)
}