processAndFilterPETS <- function(	outname,
									downSample)
{
	#Step 1: Set various filenames
    bedpefile          = paste(outname ,".bedpe",sep="")
    bedpefilesort      = paste(outname ,".sort.bedpe",sep="")
    bedpefilesortrmdup = paste(outname ,".sort.rmdup.bedpe",sep="")
    bam1 = ifelse(downSample<1.0,paste(outname ,"_1.same_downSampled_",downSample,".sorted.bam",sep=""),paste(outname ,"_1.same.sorted.bam",sep=""))
    bam2 = ifelse(downSample<1.0,paste(outname ,"_2.same_downSampled_",downSample,".sorted.bam",sep=""),paste(outname ,"_2.same.sorted.bam",sep=""))
	
	#Step 2: Constrct BEDPE file from the two BAM files
    print ("building bedpe")
    if (file.exists(bedpefile)) {
		file.remove(bedpefile)
	}
    buildBedpefromBam(bam1 =bam1, bam2 = bam2, bedpefile = bedpefile)
	
	#Step 3: Remove duplicate PETs.
	#To speed things up, PETs are split into regions (chromsomes/bins) and then PETs with similar start/end coordinate are 
	#discarded.
    print ("removing duplicate PETs")
    distancesplit = 10000000
    rmdupresults = removeDups(bedpefile,outname,distancesplit)
	
	#Step 4: Cleanup unnecessary files
    for (f in (5:length(rmdupresults)))
    {
      if (file.exists(rmdupresults[f])==TRUE){file.remove(rmdupresults[f])} 
    }
	
	return(rmdupresults)
}