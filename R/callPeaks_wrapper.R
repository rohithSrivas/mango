callPeaks.wrapper <- function(	bedtoolsgenome,
								outname,
								MACS_qvalue=0.05,
								peakslop=500,
								peakinput,
								MACS_shiftsize=NULL,
								blacklist)
{
	#Step 1: Setup and establish filenames
    bedpefilesortrmdup = paste(outname ,".rmdup.bedpe",sep="")
    tagAlignfile       = paste(outname,".tagAlign",sep="")
    peaksfile          = paste(outname,"_peaks.narrowPeak",sep="")
    peaksfileslop      = paste(outname,"_peaks.slopPeak",sep="")
	
	#Step 2: Check if we even need to call peaks using MACS2, i.e., has user supplied a list of peaks to use?
    if (peakinput != "NULL") {
      peaksfile = peakinput
    } else {
    	#Step 2.1: Build a TAG align file for peak calling (reverse strands)
	    print ("building tagAlign file")
	    if (file.exists(tagAlignfile)){file.remove(tagAlignfile)}
	    buildTagAlign(bedpefilesortrmdup,tagAlignfile)
		
		#Step 2.2: Call peaks
	    print ("calling peaks")
	    callpeaks(	macs2path=macs2path,
					tagAlignfile=tagAlignfile,
					peaksfile=outname,
					qvalue=MACS_qvalue,
					MACS_shiftsize=MACS_shiftsize)
    }
	
	#Step 3: Perform the following operations:
	#	(1) Extend peaks by 'peakslop' amount in either direction.
	#	(2) Merge overlapping peaks into a single interval
	#	(3) Remove peaks overlapping with regions in the blacklist file
	print ("extending peaks")
    peakcounts = extendpeaks(	peaksfile=peaksfile,
  	  						peaksfileslop=peaksfileslop,
  							bedtoolspath=bedtoolspath,
               			   	bedtoolsgenome=bedtoolsgenome,
  							peakslop=peakslop,
  							blacklist=blacklist)
	
	return(peakcounts)
}