alignReads <- function(	outname,
						bowtieref,
						shortreads,
						numThreads,
						downSample=1.0)
{
	#Step 0: Set various filenames
    fastq1 = 			paste(outname ,"_1.same.fastq",sep="")
    fastq2 = 			paste(outname ,"_2.same.fastq",sep="")
    bam1   = 			paste(outname ,"_1.same.bam",sep="")
    bam2   = 			paste(outname ,"_2.same.bam",sep="")
    bam1.sorted = 		paste(outname ,"_1.same.sorted.bam",sep="")
    bam2.sorted = 		paste(outname ,"_2.same.sorted.bam",sep="")
	
	#Step 1: Align each FASTQ file separately
	print ("aligning reads")
    count1 <- alignBowtie(fastq=fastq1,output=bam1,bowtiepath=bowtiepath,bowtieref=bowtieref,samtoolspath=samtoolspath,shortreads,num.threads=numThreads)
    count2 <- alignBowtie(fastq=fastq2,output=bam2,bowtiepath=bowtiepath,bowtieref=bowtieref,samtoolspath=samtoolspath,shortreads,num.threads=numThreads)
	
	#Step 2: Sort by name to ensure consistent pairing of PETs
    print("sorting reads")
    sortBAM(bamInputFile=bam1,bamOutputFile=bam1.sorted,samtools=samtoolspath,by.name=TRUE,num.threads=numThreads,verbose=TRUE)
    sortBAM(bamInputFile=bam2,bamOutputFile=bam2.sorted,samtools=samtoolspath,by.name=TRUE,num.threads=numThreads,verbose=TRUE)
	
	#Step 3: If requested perform downsampling
    if(downSample < 1.0) 
    {
		checkRequired(opt,c("picardtoolspath"))

		print ("performing downsampling")
		bam1.ds = 		paste(outname ,"_1.same_downSampled_",downSample,".bam",sep="")
		bam2.ds = 		paste(outname ,"_2.same_downSampled_",downSample,"bam",sep="")

		bam1.ds.sorted = 	paste(outname ,"_1.same.downSampled_",downSample,".sorted.bam",sep="")
		bam2.ds.sorted = 	paste(outname ,"_2.same.downSampled_",downSample,".sorted.bam",sep="")

		#Down sample
		downSampleBam(bam1,bam1.ds,downSample,verbose=TRUE)
		downSampleBam(bam2,bam2.ds,downSample,verbose=TRUE)

		#Sort by name
		sortBAM(bamInputFile=bam1.ds,bamOutputFile=bam1.ds.sorted,samtools=samtoolspath,by.name=TRUE,num.threads=numThreads,verbose=TRUE)
		sortBAM(bamInputFile=bam1.ds,bamOutputFile=bam1.ds.sorted,samtools=samtoolspath,by.name=TRUE,num.threads=numThreads,verbose=TRUE)
		
		#Remove unsorted BAM files
		file.remove(bam1.ds)
		file.remove(bam2.ds)
	}
	
	#Step 4: Perform clean-up
	if(file.exists(fastq1)) {file.remove(fastq1)}
	if(file.exists(fastq2)) {file.remove(fastq2)}
	if(file.exists(bam1)) {file.remove(bam1)}
	if(file.exists(bam2)) {file.remove(bam2)}
	
	return(c(count1,count2))
}