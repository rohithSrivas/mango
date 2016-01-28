callInteractions_wrapper <- function(	outname,
										bedtoolsgenome,
										bedtoolspath=NULL,
										distcutrangemin=1000,
										distcutrangemax=100000,
										biascut=0.05,
										maxinteractingdist=1000000,
										numofbins=50,
										FDR=0.05,
										minPETS=2,
										chrominclude="NULL",
										chromexclude="chrM,chrY",
										reportallpairs=TRUE,
										corrMethod="BH",
										MHT="all",
										extendreads=120,
										verboseoutput=FALSE,
										saveRoutput=TRUE,
										set.manual.cutoff=NA)
{
	#Step 1: Establish various filenames required
    tagAlignfile       	= paste(outname,".tagAlign",sep="")
    tagAlignfileExt    	= paste(outname ,".tagAlign.extended.bed",sep="")
    temppeakoverlap    	= paste(outname ,".temppeakoverlap.bed",sep="")
    peaksfile          	= paste(outname ,"_peaks.narrowPeak",sep="")
    peaksfileslop      	= paste(outname ,"_peaks.slopPeak",sep="")
    peaksfileslopdepth 	= paste(outname ,"_peaks.slopPeak.depth",sep="")
    bedpefilesortrmdup 	= paste(outname ,".rmdup.bedpe",sep="") 
    distancefile       	= paste(outname ,".distance",sep="")
    distancecutpdf     	= paste(outname ,".distance.pdf",sep="")
    modelspdf          	= paste(outname ,".models.pdf",sep="")
    allpairsfile       	= paste(outname ,".interactions.all.mango",sep="")
    fdrpairsfile       	= paste(outname ,".interactions.fdr.mango",sep="")
	routputfile			= paste(outname,".allInfo.Rdata",sep="")
	
	#Step 1.5: If bedtoolspath is not provided auotmagically find from path
	if(is.null(bedtoolspath)){
		bedtoolspath  = DefinePaths(c("bedtools"))[1]
	}
	
	#Step 2: Count number of reads per peak; this performs the following operations
	#	(1) Extend reads 'extendreads' length on 3' end.
	#	(2) Intersect the extended tagAlign file with the extended peakset file (one line per overlap)
	#	(3) Convert this intersect file to a peak depth file, i.e., peak location<tab>peak name<tab># of reads
    if (file.exists(tagAlignfileExt) ==TRUE){file.remove(tagAlignfileExt)}
    if (file.exists(temppeakoverlap) ==TRUE){file.remove(temppeakoverlap)}
    DeterminePeakDepths(	bedtools=bedtoolspath,bedtoolsgenome=bedtoolsgenome,extendreads=extendreads,tagAlignfile=tagAlignfile,
                    		tagAlignfileExt=tagAlignfileExt,peaksfileslop=peaksfileslop,temppeakoverlap=temppeakoverlap,
							peaksfileslopdepth=peaksfileslopdepth)
    if (file.exists(tagAlignfileExt) ==TRUE){file.remove(tagAlignfileExt)}
    if (file.exists(temppeakoverlap) ==TRUE){file.remove(temppeakoverlap)}
		
	
	#Step 3: Determine self-ligation distance cutoff. This is done via a two-step process
	# Sub-Step 1: 	Generate a distance file wherein each line represents a single PET with the distance between PETs and whether 
	#			or not they are on the same strand or not. This only outputs PETs within a certain distance range
	#			defined 
	distancecutoff <- NULL
	
	if(is.na(set.manual.cutoff))
	{
		print ("determining self-ligation distance de-novo")
	    makeDistanceFile(bedpefilesortrmdup,distancefile,distcutrangemin,distcutrangemax)
	
		# Sub-Step 2:	(1) Use the distance file bin PETs into 50 bins between distcutmin and distcutmax
		#				(2) Calculate the percent of PETs on opposing strands (self-ligation products)[P]
		#				(3) Self-ligation distance cutoff is determined as the largest bin for which P >= biastCut
	    distancecutoff = calcDistBias(	distancefile,
										distancecutpdf=distancecutpdf,
	                                  	range=c(distcutrangemin,distcutrangemax),
	                                  	biascut= biascut)
		print (paste("self-ligation cutoff =",distancecutoff))
	} else {
		print ("setting user-specified self-ligation distance cutoff")
		distancecutoff = set.manual.cutoff
		print (paste("self-ligation cutoff =",distancecutoff))
	}
	
	
	#Step 4: Begin to assemble information about putatitve peak pairs
	
	#Step 4.1: 	Split PETs by chromosome (for computational ease), extend by "extendread" amount in 3' direction and intersect
	#			with list of peaks, and then generate putatitve peak pairs on a chromsome by chromosome bases, this is a file
	#			with the following information
	#			Column 1: Chromosome # for peak 1
	#			Column 2: Start coordinate for peak 1
	#			Column 3: End coordinate for peak 1
	#			Column 4: Chromosome # for peak 2
	#			Column 5: Start coordinate for peak 2
	#			Column 6: End Coordinate for peak 2
	#			Column 7: Arbitrary/unique string name for peak pair (peak1,":",peak2)
	#			Column 8: Arbitrary/unique string name for peak1
	#			Column 9: Arbitrary/unique string name for peak2
	#			Column 10: Read depth at peak 1
	#			Column 11: Read depth at peak 2
	#			Column 12: # of PETs connecting peak 1 and peak 2
	#			Column 13: Distance between peak 1 and peak 2 [defined as start2-end1]
    print ("grouping PETs into interactions")
    chromosomes = groupPairs(bedpefilesortrmdup=bedpefilesortrmdup,
                             outname=outname,
                             peaksfile=peaksfileslop,
                             bedtoolspath = bedtoolspath,
                             bedtoolsgenome = bedtoolsgenome,
                             extendreads=extendreads,peaksfileslopdepth=peaksfileslopdepth,
                             verbose=FALSE)
     
	 # filter out unwanted chromosomes
     originalchroms = chromosomes
					 
	#Step 4.2: Keep peak pairs located on "chrominclude" and remove peak pairs located on "chromexclude"
    bedtoolsgenomeinfo = read.table(bedtoolsgenome,header=FALSE,sep="\t")
    chromosomes = bedtoolsgenomeinfo[,1]
    chromosomes = chromosomes[grep("_",chromosomes,invert=TRUE)]
    if(chrominclude[1] != "NULL")
    {
      chromosomes = unlist(strsplit(chrominclude,split=","))
    }
  
    if (chromexclude[1] !=  "NULL")
    {
      chromosomestpremove = unlist(strsplit(chromexclude,split=","))
      chromosomes = chromosomes[-which(chromosomes %in% chromosomestpremove)] 
    }
	
	#Step 4.3: After inclusion/exclusion of certain chromosomes merge peak-pair list across all chromsomes
	print("merging peak pairs across all chromosomes")
	putpairs = combineputativepairs(chromosomes,outname)
	
	#Step 4.4: Calculate a new distance metric which is the difference in midpoints between the two peaks
	putpairs$distances = abs( (putpairs[,2] + putpairs[,3] ) / 2 - (putpairs[,5] + putpairs[,6] ) / 2  )
	
	#Step 4.5: Filter out peak pairs that fall below self-ligation cutoff ('distancecutoff) or are above maxinteractindist
	putpairs = putpairs[which(putpairs$distances < maxinteractingdist & putpairs$distances > distancecutoff),]
	
	#Step 4.6: Calculate a combined peak depth metric (either the sum or product of the two individual peak depths)
	putpairs$depths = calcDepths(putpairs[,10:11],type="product")
	
	#Step 4.7: Estimate p-value of observing the # of PETs seen between each peak-pair. This is done as follows:
	#	(1) Estimate numerous individual probabilities:
	#		(1.1) Probability of observing a PET at various peak pair distances
	#		(1.2) Probability of observing a PET at various peak pair depths
	#		(1.3) Probability of observing any peak pair at various distances
	#		(1.4) Probability of observing any peak pair at various peak pair depths
	#	(2) Utilize Bayes formula to combine these into probability of observing a PET between each peak pair.
	#	(3) Utilize the binomial distribution to calculate probability of observing "X or more" PETS at each peak pair (p-value)
	#	(4) Remove peak pairs achieving extraordinary significance (i.e., beats Bonferroni correction)
	#	(5) Repeat Steps 1-3 to re-derive probabilities after removing "most likely" interactions
    totalcombos = 0
	for (reps in (1:2))
    {
      #--------------- Distance Normalization ---------------#
    
      # determine borders to distance bins
      distanceborders = binmaker(putpairs$distances,binmethod="equalocc",numberbins=numofbins)
    
      # model IAB vs distance
      distance_IAB_model = model_chia(x=putpairs$distances,y=putpairs[,12],borders=distanceborders,yvals=TRUE)
      distance_IAB_model_file   = paste(outname ,".distance_IAB_model.",reps, ".text",sep="")
      write.table(distance_IAB_model,file=distance_IAB_model_file,quote = FALSE, sep = "\t",row.names = FALSE,col.names = TRUE)
      distance_IAB_spline =   smooth.spline(log10(distance_IAB_model[,1]),distance_IAB_model[,3],spar=.75)
    
      #--------------- Depth Normalization ---------------#
    
      # determine borders to depth bins
      depthborders = binmaker(putpairs$depths,binmethod="equalocc",numberbins=numofbins)
    
      # model IAB vs depth
      depth_IAB_model = model_chia(x=putpairs$depths,y=putpairs[,12],borders=depthborders,yvals=TRUE)
      depth_IAB_model_file   = paste(outname ,".depth_IAB_model.",reps, ".text",sep="")
      write.table(depth_IAB_model,file=depth_IAB_model_file,quote = FALSE, sep = "\t",row.names = FALSE,col.names = TRUE)
      depth_IAB_spline =   smooth.spline(log10(depth_IAB_model[,1]),depth_IAB_model[,3],spar=.75)
    
      # model Combos vs distance
      meanofx_dist  = rep(0,numofbins)
      sumofy_dist   = rep(0,numofbins)
      pvals_dist    = rep(0,numofbins)
      sumofx_dist   = rep(0,numofbins)
      countofx_dist = rep(0,numofbins)
      meanofx_depth  = rep(0,numofbins)
      sumofy_depth   = rep(0,numofbins)
      pvals_depth    = rep(0,numofbins)
      sumofx_depth   = rep(0,numofbins)
      countofx_depth = rep(0,numofbins)

      for (chrom in chromosomes)
      {
        # make combos
        combos = makecombos(chrom,outname,mindist=distancecutoff,maxdist=maxinteractingdist)
      
        # calculate distances
        combos$distance = abs( (combos[,2] + combos[,3] ) / 2 - (combos[,5] + combos[,6] ) / 2  )
      
        # calculate depths
        combos$depths = calcDepths(combos[,7:8],type="product")
      
        # model Combo vs distance
        distance_combo_model_chrom = model_chia(x=combos$distance,y=NA,borders=distanceborders,yvals=FALSE)
      
        sumofy_dist   = sumofy_dist   + distance_combo_model_chrom[,2]
        sumofx_dist   = sumofx_dist   + distance_combo_model_chrom[,4]
        countofx_dist = countofx_dist + distance_combo_model_chrom[,5]
       
        # model Combo vs depth
        depth_combo_model_chrom  = model_chia(x=combos$depths,y=NA,borders=depthborders,yvals=FALSE)
      
        sumofy_depth   = sumofy_depth   + depth_combo_model_chrom[,2]
        sumofx_depth   = sumofx_depth   + depth_combo_model_chrom[,4]
        countofx_depth = countofx_depth + depth_combo_model_chrom[,5]
      }
    
      # combine data from all chromosomes
      depth_combo_model    = cbind(sumofx_depth/countofx_depth, sumofy_depth, sumofy_depth/sum(sumofy_depth))
      distance_combo_model = cbind(sumofx_dist /countofx_dist,  sumofy_dist, sumofy_dist/sum(sumofy_dist))
    
      depth_combo_model_file   = paste(outname ,".depth_combo_model.",reps, ".text",sep="")
      write.table(depth_combo_model,file=depth_combo_model_file,quote = FALSE, sep = "\t",row.names = FALSE,col.names = TRUE)
    
      distance_combo_model_file   = paste(outname ,".distance_combo_model.",reps, ".text",sep="")
      write.table(distance_combo_model,file=distance_combo_model_file,quote = FALSE, sep = "\t",row.names = FALSE,col.names = TRUE)
    
      depth_combo_spline    =   smooth.spline(log10(depth_combo_model[,1]),depth_combo_model[,3],spar=.75)
      distance_combo_spline =   smooth.spline(log10(distance_combo_model[,1]),distance_combo_model[,3],spar=.75)

      if (reps == 2)
      {
        #--------------- Gather IAB data again ---------------#
      
        # gather all putative interactions
        putpairs = combineputativepairs(chromosomes,outname)
      
        # calculate interaction distances
        putpairs$distances = abs( (putpairs[,2] + putpairs[,3] ) / 2 - (putpairs[,5] + putpairs[,6] ) / 2  )
      
        # filter out putative interactions that don't fall into distance range
        putpairs = putpairs[which(putpairs$distances < maxinteractingdist & putpairs$distances > distancecutoff),]
      
        # calculate depths
        putpairs$depths = calcDepths(putpairs[,10:11],type="product")
      }

      #--------------- Score putative interactions ---------------#
    
  #    putpairstemfile       = paste(outname ,".putpairs",sep="")
  #    write.table(putpairs,file=putpairstemfile,quote = FALSE, sep = "\t",row.names = FALSE)
    
      # Assing the four probabilities
      putpairs$P_IAB_distance    = predict(distance_IAB_spline, log10(putpairs$distances))$y
      putpairs$P_combos_distance = predict(distance_combo_spline,log10(putpairs$distances))$y
      putpairs$P_IAB_depth       = predict(depth_IAB_spline,log10(putpairs$depths))$y
      putpairs$P_combos_depth    = predict(depth_combo_spline,log10(putpairs$depths))$y
    
  #     # fix negative values
  #     putpairs$P_IAB_distance[which(putpairs$P_IAB_distance <= 0)] = 
  #       min(putpairs$P_IAB_distance[which(putpairs$P_IAB_distance > 0)])
  #     putpairs$P_combos_distance[which(putpairs$P_combos_distance <= 0)] = 
  #       min(putpairs$P_combos_distance[which(putpairs$P_combos_distance > 0)])
  #     putpairs$P_IAB_depth[which(putpairs$P_IAB_depth <= 0)] = 
  #       min(putpairs$P_IAB_depth[which(putpairs$P_IAB_depth > 0)])
  #     putpairs$P_combos_depth[which(putpairs$P_combos_depth <= 0)] = 
  #       min(putpairs$P_combos_depth[which(putpairs$P_combos_depth > 0)])

      # cap values to min and max
      putpairs$P_IAB_distance[which(putpairs$P_IAB_distance <= min(distance_IAB_model[,3]))] = min(distance_IAB_model[,3])
      putpairs$P_IAB_distance[which(putpairs$P_IAB_distance >= max(distance_IAB_model[,3]))] = max(distance_IAB_model[,3]) 

      putpairs$P_combos_distance[which(putpairs$P_combos_distance <= min(distance_combo_model[,3]))] = min(distance_combo_model[,3])
      putpairs$P_combos_distance[which(putpairs$P_combos_distance >= max(distance_combo_model[,3]))] = max(distance_combo_model[,3])

      putpairs$P_IAB_depth[which(putpairs$P_IAB_depth <=  min(depth_IAB_model[,3]))] =  min(depth_IAB_model[,3])  
      putpairs$P_IAB_depth[which(putpairs$P_IAB_depth >=  max(depth_IAB_model[,3]))] =  max(depth_IAB_model[,3])
                                                                                          
      putpairs$P_combos_depth[which(putpairs$P_combos_depth <= min(depth_combo_model[,3]))] =  min(depth_combo_model[,3])   
      putpairs$P_combos_depth[which(putpairs$P_combos_depth >= max(depth_combo_model[,3]))] =  max(depth_combo_model[,3])   

      # calculate the binomial probability
      totalcombos =  sum(sumofy_dist)
      putpairs$p_binom           = (putpairs$P_IAB_distance * putpairs$P_IAB_depth) / 
        (putpairs$P_combos_distance * putpairs$P_combos_depth * totalcombos)
   
      # calculate the total IABs
      totalIAB = sum(distance_IAB_model[,2])

      # calculate the final interaction P values
      putpairs$P = apply(cbind(putpairs$V12,rep(totalIAB,nrow(putpairs)),putpairs$p_binom),1,calcP)

      if (reps == 1)
      {
        putpairs[which(putpairs$P < 1/totalcombos ),]$V12 = 0
      }
    }
	
	#Step 4.8: Correct for multiple hypothesis testing (default: total tests = total peak pairs)
    print ("correcting for multiple hypothesis testing")
    n=nrow(putpairs)
    if (MHT == "all")
    {
      n=totalcombos
    }

    putpairs$Q = p.adjust(putpairs$P,method=corrMethod,n=n)
	
	#Step 4.9: Provide pretty header names; discard unique/random identified given to each peak
    pairnames = paste("pair_",(1:nrow(putpairs)),sep="")
    putpairs = cbind(putpairs[,c(1,2,3,4,5,6)],pairnames,putpairs[,c(10,11,12,14,16,17,18,19,20,21,22)])
    names(putpairs) = c("chrom1","start1","end1","chrom2","start2","end2","name",
                        "peak1","peak2","PETs","distance",
                        "P_IAB_distance","P_combos_distance","P_IAB_depth","P_combos_depth",
                        "p_binom","P","Q")
	
	#Step 4.10: Filter out non-significant interation (FDR cutoff) and interactions supported by fewer then minPETs
	sig = putpairs[which(putpairs$Q < FDR & putpairs$PETs >= minPETS),]
	
	#Step 4.11: Output the data-frame to file (with all information or limited information)
    print ("writing output files")
	if (verboseoutput == TRUE)
	{
		if (reportallpairs == TRUE)
		{
			write.table(x=putpairs,file=allpairsfile,quote = FALSE, sep = "\t",row.names = FALSE)
		}
		write.table(x=sig,file=fdrpairsfile,quote = FALSE, sep = "\t",row.names = FALSE)
	}
	if (verboseoutput == FALSE)
	{
		if (reportallpairs == TRUE)
		{
			write.table(x=putpairs[,c("chrom1","start1","end1","chrom2","start2","end2","PETs","Q")],file=allpairsfile,quote = FALSE, sep = "\t",row.names = FALSE,col.names = FALSE)
		}
		write.table(x=sig[,c("chrom1","start1","end1","chrom2","start2","end2","PETs","Q")],file=fdrpairsfile,quote = FALSE, sep = "\t",row.names = FALSE,col.names = FALSE)
	} 
	
	#Step 4.12: If requested save all data-frames in R storage format
	if(saveRoutput)
	{
		print("saving all relevant data in R dataformat.")
		save(	putpairs,
				distance_IAB_model,distance_IAB_spline,
				distance_combo_model,distance_combo_spline,
				depth_IAB_model,depth_IAB_spline,
				depth_combo_model,depth_combo_spline,
				totalcombos,
				file=routputfile)
		
	}
	
	#Step 4.13: Save plots illustrating how the four probabilities were estimated
	print ("plotting results")

	pdf(modelspdf)
	par(mfrow=c(2,2))
	par(mgp=c(3,.3,0))
	plot(log10(distance_IAB_model[,1]),   distance_IAB_model[,3],pch=19,col="dodgerblue2",xlab="distance (bp)",ylab="IAB")
	lines(x=log10(distance_IAB_model[,1]),   predict(distance_IAB_spline,log10(distance_IAB_model[,1]))$y)   

	plot(log10(distance_combo_model[,1]), distance_combo_model[,3],  pch=19,col="firebrick2",xlab="distance (bp)",ylab="# combos")  
	lines(x=log10(distance_combo_model[,1]),   predict(distance_combo_spline,log10(distance_combo_model[,1]))$y)  

	plot(log10(depth_IAB_model[,1])   ,  depth_IAB_model[,3],    pch=19,col="dodgerblue2",xlab="depth (p1 * p2)",ylab="IAB") 
	lines(x=log10(depth_IAB_model[,1]),   predict(depth_IAB_spline,log10(depth_IAB_model[,1]))$y)

	plot(log10(depth_combo_model[,1]),   depth_combo_model[,3],  pch=19,col="firebrick2",xlab="depth (p1 * p2)",ylab="# combos") 
	lines(x=log10(depth_combo_model[,1]),   predict(depth_combo_spline,log10(depth_combo_model[,1]))$y)
	dev.off()
	
	#Step 4.14: Delete temporary files generated (i.e., chromsome specific files)
    print ("deleting temporary files")
    if (file.exists(distancefile)) file.remove(distancefile)
    for (chrom in originalchroms)
    {
      peaksizecount  = paste(outname,"." ,chrom, "_peaks.count.slopPeak",sep="")
      pairsbedpe     = paste(outname,"." ,chrom, ".pairs.bedpe",sep="")
      bedpefile      = paste(outname,"." ,chrom, ".bedpe",sep="")
      bedfile        = paste(outname,"." ,chrom, ".bed",sep="")
      overlapfile    = paste(outname,"." ,chrom, ".bedNpeak",sep="")
      if (file.exists(peaksizecount)) file.remove(peaksizecount)
      if (file.exists(pairsbedpe)) file.remove(pairsbedpe)
      if (file.exists(bedpefile)) file.remove(bedpefile)
      if (file.exists(bedfile)) file.remove(bedfile)
      if (file.exists(overlapfile)) file.remove(overlapfile)
    }
	
	return(c(nrow(putpairs),nrow(sig)))
}