
args=commandArgs(trailingOnly=T)

print(args)
metadata.f=args[1]
motif.f=args[2]
scoredRegion2peak.f=args[3]
motifScores.f=args[4]
psignal.f=args[5]
out=args[6]
matches.f=args[7] #set this to "NA" if you've already pre-filtered your SNPQTL matrix
totalPermsEmpDistr=100000 
motifMAFthreshold=as.numeric(as.character(args[8]))

#############################################################
#              FUNCTIONS USED HERE                          #
#############################################################
testing=function(){
metadata.f='/srv/gsfs0/projects/snyder/oursu/histoneQTL/motif_analysis/results/2015-01-13/metadata.75unrelatedYRI'
motif.f='/srv/gs1/projects/snyder/jzaugg/histoneQTL/motif_analysis/data/motifPWM/matchesPerTF/correctedByPouyaControls/motifFreqMatrices/motif.pouya.Motif.EBF1_known1scanThresh0'
scoredRegion2peak.f='/srv/gsfs0/projects/snyder/oursu/histoneQTL/motif_analysis/results/2015-01-13/Local/SNPQTLmatrix/motifMatched/SNPQTLmatrix.dhs.MotifMatch_EBF1_known1scanThresh0.bed.OverlapChIPseq.gz'
motifScores.f='/srv/gsfs0/projects/snyder/oursu/histoneQTL/motif_analysis/results/2015-01-13/Motif_analysis_results/Motif_Analysis_dhs/motif.pouya.Motif.EBF1_known1scanThresh0/Scores/motif.pouya.Motif.EBF1_known1scanThresh0_maternalScores.compiled.gz,/srv/gsfs0/projects/snyder/oursu/histoneQTL/motif_analysis/results/2015-01-13/Motif_analysis_results/Motif_Analysis_dhs/motif.pouya.Motif.EBF1_known1scanThresh0/Scores/motif.pouya.Motif.EBF1_known1scanThresh0_paternalScores.compiled.gz'
psignal.f='/srv/gsfs0/projects/snyder/oursu/histoneQTL/motif_analysis/results/2015-01-13//Data/CombinedPeakSignal.gz'
#/srv/gs1/projects/snyder/jzaugg/histoneQTL/hQTL_analysis/PeakData/hMat.norm.ALL.dhs.peer_lS_5.txt'
matches.f='/srv/gsfs0/projects/snyder/oursu/histoneQTL/motif_analysis/results/2015-01-13/Motif_Data/motifMatches/MotifMatch_EBF1_known1scanThresh0.bed.gz.OverlapChIPseq.gz'
motifMAFthreshold=0.05
scoredRegion2peak.f='/srv/gsfs0/projects/snyder/oursu/histoneQTL/motif_analysis/results/2015-01-13/Distal/SNPQTLmatrix/SNPQTLmatrix.LocalPeakIsdhs.gz'
}

getPeopleOfInterest=function(metadata.f){
        metadata=read.table(metadata.f,sep='\t')
        people=as.character(unique(metadata[,1]))
	return(people)
}
getMotifLength=function(motif.f){
        #Get max motif score
        motif=read.table(motif.f,skip=1)
        motif.norm=motif/rowSums(motif)
        motifLen=nrow(motif)
	return(motifLen)
}
readScoredRegion2Peak=function(scoredRegion2peak.f,motifLen){
	scoredRegion2peak=read.table(scoredRegion2peak.f,header=T)
        toremove=union(which(duplicated(paste(as.character(scoredRegion2peak$SNP),'_',as.character(scoredRegion2peak$gene),sep=''))),
		which(as.character(scoredRegion2peak$chr)=='chr'))
        if (length(toremove)>0){
           scoredRegion2peak=scoredRegion2peak[-toremove,]
        }
        #add name of the scored region
	chrs=as.character(scoredRegion2peak$chr)
	#### corrected snp positions!!
	CORRECTED=-1 #because snp coords were 1based, but the scored region comes from a bed file
	starts=as.numeric(as.character(scoredRegion2peak$snp.position))-motifLen+1+CORRECTED
	ends=as.character(as.numeric(as.character(scoredRegion2peak$snp.position))+motifLen+CORRECTED)
        scoredRegion.names=paste('chr',chrs,'_',as.character(starts),'_',as.character(ends),sep='')
	#build the scoredRegion2Peak data structure
        scoredRegion2peak=cbind(scoredRegion2peak,
				scoredRegionName=scoredRegion.names,
				start=starts,
				end=ends)
	return(scoredRegion2peak)
}
readMotifScores=function(motifScores.f){
	first=TRUE
        scoreFiles=strsplit(motifScores.f,',')[[1]]
        for (f in scoreFiles){
        	if (first==FALSE){
        		cur_motif=read.table(f,header=T)
                 rownames(cur_motif)=cur_motif[,1]
                 cur_motif=cur_motif[,-1]
                 motif.score.avg=motif.score.avg+cur_motif[rownames(motif.score.avg),colnames(motif.score.avg)]
        	}
            if (first==TRUE){
               motif.score.avg=read.table(f,header=T)
               rownames(motif.score.avg)=motif.score.avg[,1]
               motif.score.avg=motif.score.avg[,-1]
               first=FALSE
            }
           
        }
        motif.score.avg=motif.score.avg/(length(scoreFiles))
        rownames(motif.score.avg)=gsub(':','_',gsub('-','_',rownames(motif.score.avg)))
	return(motif.score.avg)
}
getVariableMotifs=function(motif.score.avg,motifMAFthreshold){
	motif.score.max=apply(motif.score.avg,1,max)
    motif.score.min=apply(motif.score.avg,1,min)
    motif.score.variable=which(as.numeric(as.character(motif.score.min))!=as.numeric(as.character(motif.score.max))) #### this and motif match define the studied TF motifs
	
	npeople=dim(motif.score.avg)[2]
	mode_counts=apply(motif.score.avg,1, function (x) max(table((x))))/npeople
	#filter by the motif maf
	
	return(rownames(motif.score.avg)[as.numeric(as.character(intersect(motif.score.variable,
		which(as.numeric(as.character(mode_counts))<=(1-motifMAFthreshold)))))])
}

getSNPregionsWithMotifMatch=function(matches.f,correlation_result){
        require(rtracklayer)
        require(GenomicRanges)

	#make range for matches
        m_data=read.table(matches.f) #matches bed file              #+1 b/c granges is 1based          #end included in interval, and 1based from bed already                              
        motif_matches=GRanges(seq=m_data[,1],ranges=IRanges(as.numeric(as.character(m_data[,2]))+1,end=as.numeric(as.character(m_data[,3]))))

        #make range for snps
	snpbed=data.frame(do.call(rbind,strsplit(as.character(correlation_result$scoredRegionName),'_')))
        snpbed.dupli=which(duplicated(snpbed))
        if (length(snpbed.dupli)>0){
           snpbed=snpbed[-snpbed.dupli,]
        }
	snps=GRanges(seq=snpbed[,1],ranges=IRanges(as.numeric(as.character(snpbed[,2]))+1,
                end=as.numeric(as.character(snpbed[,3])),names=paste(snpbed[,1],'_',snpbed[,2],'_',snpbed[,3],sep='')))

        #overlap them to get motif matches
        snp.vs.motif.match=as.data.frame(findOverlaps(query=motif_matches,subject=snps,type='within',select='all')) #motif must fall inside scored region
        duplicated_stuff=which(duplicated(snp.vs.motif.match[,2]))
        if (length(duplicated_stuff)>0){
           snp.vs.motif.match=snp.vs.motif.match[-which(duplicated(snp.vs.motif.match[,2])),] #remove times when multiple motifs in the same SNP
        }
	return(rownames(as.data.frame(snps)[as.numeric(as.character(snp.vs.motif.match[,2])),])) #the snp regions that contain a motif match
}
getMaxObservedMotif=function(motif.scores,scoredRegionNames){
	motif.maxs=apply(motif.scores,1,max)
	return(motif.maxs[as.character(scoredRegionNames)])
}
##################################################################

###############################
#     MAIN ACTION             #
###############################
	.libPaths('~/devtools/R-3.0.2/library') #### remove this soon.

	#Get people of interest
        #======================
	people=getPeopleOfInterest(metadata.f)

	#Get motif length
	#==============
	motifLen=getMotifLength(motif.f)
	print(paste('Got motif length: ',motifLen,sep=''))

	# Read in mapping from scored region to peak
	#===========================================
	scoredRegion2peak=readScoredRegion2Peak(scoredRegion2peak.f,motifLen)
	#snp_col=union(which(colnames(scoredRegion2peak)=='snp'),which(colnames(scoredRegion2peak)=='SNP'))[1]
	#colnames(scoredRegion2peak)[snp_col]='snp'
	print('Read scored regions to peaks')

	# Read in peak signal. 
	#=====================
	#Index is index present in the scoredRegion2peak.
	psignal=read.table(psignal.f,header=T)
	print(paste('Read in peak signal from ',psignal.f,sep=''))	

	# read in motif score averages
	#=============================
	motif.score.avg=readMotifScores(motifScores.f)
	print(paste('Read in motif score averages from: ',motifScores.f,sep=''))

	# keep only things pertaining to the individuals for whom we have data
	#====================================================================
	people=intersect(people,intersect(colnames(motif.score.avg),colnames(psignal)))
	motif.scores=motif.score.avg[,people]
	psignal.scores=psignal[,people]
	print('Subsetted people to: ')
	print(people)

	#Setup data strcuture for correlations
	#=====================================
	scoredRegionNames=as.character(scoredRegion2peak$scoredRegionName)
	localpeaks=rep('NA',length(scoredRegionNames))
	if ('localpeak' %in% colnames(scoredRegion2peak)) {
		localpeaks=scoredRegion2peak$localpeak
	}
	correlation_result=data.frame(scoredRegionName=scoredRegionNames,
					snp=scoredRegion2peak$SNP,
					peak=scoredRegion2peak$gene,
					snp.pos=scoredRegion2peak$snp.position,
					localpeak=localpeaks,
					spearmanCorr_MotifScoreVSpeak=NA,
					pearsonCorr_MotifScoreVSpeak=NA,
					maxObsMotif=NA,
					motifName=basename(motif.f))
	print('Setup correlation structure')

	#keep only SNP regions with motif matches and variable scores
	#============================================================
	#if (matches.f!="NA"){
	scoredRegionsToKeep=as.character(intersect(getVariableMotifs(motif.scores,motifMAFthreshold),getSNPregionsWithMotifMatch(matches.f,correlation_result)))
	#}
	correlation_result=correlation_result[which(as.character(correlation_result$scoredRegionName) %in% scoredRegionsToKeep),]
	correlation_result$maxObsMotif=getMaxObservedMotif(motif.scores,as.character(correlation_result$scoredRegionName))
	#Compute correlations
	#====================
	print('Starting motif correlations')
	print(paste(length(unique(correlation_result$snp)),' SNP instances tested',sep=''))	
	rownames(correlation_result)=paste(correlation_result$scoredRegionName,'_',correlation_result$peak,sep='')
	group_tests=split(1:nrow(correlation_result),rep(1:100,each=nrow(correlation_result)/100+1))
	require(reshape)
	for (group_id in names(group_tests)){
	    print(group_id)
	    if (length(group_tests[[group_id]])==0){
	       next
	    }
	    motifs=as.character(unique(correlation_result[group_tests[[group_id]],'scoredRegionName']))
	    peaks=as.character(unique(correlation_result[group_tests[[group_id]],'peak']))
	    cur_corr_sp=cor(t(motif.scores[motifs,]),t(psignal.scores[peaks,]),method='spearman',use="pairwise.complete.obs")
	    cur_corr_sp.long=melt(cur_corr_sp)
	    rownames(cur_corr_sp.long)=paste(cur_corr_sp.long[,1],'_',cur_corr_sp.long[,2],sep='')
	    correlation_result[group_tests[[group_id]],'spearmanCorr_MotifScoreVSpeak']=cur_corr_sp.long[rownames(correlation_result)[group_tests[[group_id]]],3]
	    
	    cur_corr_pe=cor(t(motif.scores[motifs,]),t(psignal.scores[peaks,]),method='pearson',use="pairwise.complete.obs")
            cur_corr_pe.long=melt(cur_corr_pe)
            rownames(cur_corr_pe.long)=paste(cur_corr_pe.long[,1],'_',cur_corr_pe.long[,2],sep='')
            correlation_result[group_tests[[group_id]],'pearsonCorr_MotifScoreVSpeak']=cur_corr_pe.long[rownames(correlation_result)[group_tests[[group_id]]],3]
	}

	#Now that motif correlations are done, do permutations to get significance.
	#FUNCTION FOR PERMUTATIONS
	compute_highres_p=function(i,true_corr,motif.score.avg,psignal,people,permNum){
		motif=as.character(true_corr[i,'scoredRegionName'])
        	peak=as.character(true_corr[i,'peak'])
        	motif.data=t(motif.score.avg[motif,people])
        	signal.data=t(psignal[peak,people])
		perms=sapply(1:permNum,function(z) {cor(sample(signal.data,replace=FALSE),motif.data,use="pairwise.complete.obs")})
		return(perms)
	}

	numGroups=100
	group_tests=split(1:nrow(correlation_result),rep(1:numGroups,each=nrow(correlation_result)/numGroups+1))
        first=T
	perSitePerm=max(50,ceiling(totalPermsEmpDistr/dim(correlation_result)[1]))
        for (group_id in names(group_tests)){
            print(group_id)
	    cur_perm=do.call(rbind,lapply(group_tests[[group_id]],compute_highres_p,
		true_corr=correlation_result,motif.score.avg=motif.score.avg,psignal=psignal,people=people,permNum=perSitePerm))
	    if (first==FALSE){
               full_perm=rbind(full_perm,cur_perm)
            }
            if (first==TRUE){
               full_perm=cur_perm
               first=FALSE
            }
        }
	write.table(full_perm,file=paste(out,'permCorrTable.txt',sep=''),quote=F,sep='\t',row.names=F,col.names=T)

	#Finally get significance measures.
	#Site-specific pvalue from site-specific permutations
	#Genome-wide p
	correlation_result=cbind(correlation_result,site_p.sp=NA,genomeWide_match.p.sp=NA)
	get_site_p=function(i,permTable,correlation_result){
		return(1-ecdf(abs(permTable[i,]))(abs(as.numeric(as.character(correlation_result[i,'spearmanCorr_MotifScoreVSpeak'])))))
	}
	correlation_result[,'site_p.sp']=as.vector(do.call(rbind,
		lapply(1:nrow(correlation_result),get_site_p,permTable=full_perm,correlation_result=correlation_result)))
	
	if (nrow(full_perm)*ncol(full_perm)>totalPermsEmpDistr){
	   perm_fixedSize=as.matrix(t(full_perm))[1:totalPermsEmpDistr]
	}
	if (nrow(full_perm)*ncol(full_perm)<=totalPermsEmpDistr){
           perm_fixedSize=as.matrix(t(full_perm))[1:(nrow(full_perm)*ncol(full_perm))]
	   perm_fixedSize=rep(perm_fixedSize,times=ceiling(totalPermsEmpDistr/length(perm_fixedSize)))[1:totalPermsEmpDistr]
        }

	correlation_result[,'genomeWide_match.p.sp']=(1-ecdf(abs(perm_fixedSize))(abs(as.numeric(as.character(correlation_result[,'spearmanCorr_MotifScoreVSpeak']))))+1/totalPermsEmpDistr)/(1+1/totalPermsEmpDistr)

	write.table(correlation_result,file=out,quote=F,sep='\t',row.names=F,col.names=T)
	print(paste('Wrote motif correlations to ',out,sep=''))
	print("======================")
	quit()


