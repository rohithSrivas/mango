#' produce signal track data using align2rawsignal
#'
#' This function produces signal track data using the align2rawsignal program
#'  
#' @param input.bam full path to input bam file
#' @param output.mat.file output file containing signal track data in MATLAB mat format
#' @param temp.output.bedgraph.file temporary output file which contains the bedgraph file [this file will be deleted]
#' @param output.bw.file output file containing signal track data in BigWig file
#' @param path.to.mcr full path to the matlab compile runtime [MCR v17]
#' @param path.to.align2rawsignal full path to the align2rawsingal binary file
#' @param bedtoolsgenome path to the bedtools genome file (name of contigs and the size of each contig)
#' @param fragLength the average fragment length estimated from phantompeakqual toolkit
#' @param chrDir path to directory containing individual chromsome fasta files
#' @param mapDir path to directory containing mappability files
#' @export
#' 
runAlign2RawSignal <- function(	input.bam,
								output.mat.file,
								temp.output.bedgraph.file,
								output.bw.file,
								path.to.mcr,
								path.to.align2rawsignal,
								bedtoolsgenome,
								fragLength,
								chrDir,
								mapDir)
{
	
	# setup environmental path variables command
  	env.command <- paste(	"MCRROOT=",path.to.mcr,"\n",
  							"LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRROOT}/runtime/glnxa64\n",
  							"LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRROOT}/bin/glnxa64\n"
							"LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRROOT}/sys/os/glnxa64\n"
							"MCRJRE=${MCRROOT}/sys/java/jre/glnxa64/jre/lib/amd64\n"
							"LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRJRE}/native_threads\n"
							"LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRJRE}/server\n"
							"LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRJRE}\n"
							"XAPPLRESDIR=${MCRROOT}/X11/app-defaults\n"
							"export LD_LIBRARY_PATH\n"
							"export XAPPLRESDIR\n"
							"export MCR_CACHE_ROOT=",tempdir(),"\n",sep="")
	system(env.command)
	
	# setup actual run command
	command.mat <- paste(	path.to.align2rawsignal,
							" -i=\"",input.bam,"\" ",
							"-s=\""+chrDir+"\" ",
							"-u=\""+mapDir+"\" ",
							"-o=\""+output.mat.file+"\" ",
							"-of=\"mat\" ",
							"-n=5 ",
							"-l=",fragLength," ",
							"-mm=35",sep="")
  
 	 # print command invoked
	 print (command.mat)
	 
	 # execute command
	 system(command.mat)
	 
	 # setup actual run command
	 command.bg <- paste(	path.to.align2rawsignal,
 							" -i=\"",input.bam,"\" ",
 							"-s=\""+chrDir+"\" ",
 							"-u=\""+mapDir+"\" ",
 							"-o=\""+temp.output.bedgraph.file+"\" ",
 							"-of=\"bg\" ",
 							"-n=5 ",
 							"-l=",fragLength," ",
 							"-mm=35",sep="")
  
  	 # print command invoked
 	 print (command.bg)
	 
 	 # execute command
 	 system(command.bg)
	 
	 #convert bedgraph to bigWig file
	 command.bw <- paste("bedGraphToBigWig ",temp.output.bedgraph.file," ",bedtoolsgenome," ",output.bw.file,sep="")
	 print (command.bw)
	 system(command.bw)
	 
	 #clean up by removing temporary bedgraph file
	 file.remove(temp.output.bedgraph.file)
}
