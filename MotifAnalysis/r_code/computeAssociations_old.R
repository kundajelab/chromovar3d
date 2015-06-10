
args=commandArgs(trailingOnly=T)

print(args)
metadata.f=args[1]
motif.f=args[2]
scoredRegion2peak.f=args[3]
motifScores.f=args[4]
psignal.f=args[5]
out=args[6]
matches.f=args[7]
totalPermsEmpDistr=100000 
##motif_q=10


	.libPaths('~/devtools/R-3.0.2/library') #### remove this soon.

	#Get people of interest
	metadata=read.table(metadata.f,sep='\t')
	people=as.character(unique(metadata[,1]))

	print(people)

	print('Starting to compute correlations')
        #Get max motif score
        motif=read.table(motif.f,skip=1)
        motif.norm=motif/rowSums(motif)
	motifLen=nrow(motif)
        motifMax=sum(apply(motif.norm,1,max))-motifLen*log(0.25,base=10) #this is w
	print('Got motifmax')

	# Read in mapping from scored region to peak
	scoredRegion2peak=read.table(scoredRegion2peak.f,header=T)
	toremove=union(which(duplicated(paste(as.character(scoredRegion2peak$SNP),'_',as.character(scoredRegion2peak$gene),sep=''))),which(as.character(scoredRegion2peak$chr)=='chr'))
	if (length(toremove)>0){
	   scoredRegion2peak=scoredRegion2peak[-toremove,]
	}

	#add name of the scored region
	scoredRegion.names=paste('chr',as.character(scoredRegion2peak$chr),'_',as.character(as.numeric(as.character(scoredRegion2peak$snp.position))-motifLen+1),
	'_',as.character(as.numeric(as.character(scoredRegion2peak$snp.position))+motifLen),sep='')
	scoredRegion2peak=cbind(scoredRegion2peak,scoredRegionName=scoredRegion.names,
		start=as.numeric(as.character(scoredRegion2peak$snp.position))-motifLen+1,
		end=as.numeric(as.character(scoredRegion2peak$snp.position))+motifLen)
	print('Mapped scored regions to peaks')

	# Read in peak signal. Index is index present in the scoredRegion2peak.
	psignal=read.table(psignal.f,header=T)
	#Keep only people both wanted and in the dataset
	people=intersect(people,colnames(psignal))
	psignal=psignal[,people]
	print(paste('Read in peak signal from ',psignal.f,sep=''))	

	# read in motif score averages
	first=TRUE
	scoreFiles=strsplit(motifScores.f,',')[[1]]
	for (f in scoreFiles){
	    if (first==TRUE){
	       motif.score.avg=read.table(f,header=T)
	       rownames(motif.score.avg)=motif.score.avg[,1]
	       motif.score.avg=motif.score.avg[,-1]
	       first=FALSE
	    }
	    else {
	    	 cur_motif=read.table(f,header=T)
		 rownames(cur_motif)=cur_motif[,1]
		 cur_motif=cur_motif[,-1]
	     	 motif.score.avg=motif.score.avg+cur_motif[rownames(motif.score.avg),colnames(motif.score.avg)]
            }
	}
	motif.score.avg=motif.score.avg/(length(scoreFiles))
	motif.score.avg=motif.score.avg[,people]
	rownames(motif.score.avg)=gsub(':','_',gsub('-','_',rownames(motif.score.avg)))

	motif.score.max=apply(motif.score.avg,1,max)
        motif.score.min=apply(motif.score.avg,1,min)
        motif.score.variable=which(motif.score.min!=motif.score.max) #### this and motif match define the studied TF motifs

	#Setup data strcuture for correlations
	correlation_result=data.frame(scoredRegionName=as.character(scoredRegion2peak$scoredRegionName),snp=scoredRegion2peak$SNP,
	peak=scoredRegion2peak$gene,spearmanCorr_MotifScoreVSpeak=NA,pearsonCorr_MotifScoreVSpeak=NA,
	motifMaxInLog10=motifMax,maxObsMotif=NA,motifName=basename(motif.f))

	require(rtracklayer)
	require(GenomicRanges)
	m_data=read.table(matches.f)
	motif_matches=GRanges(seq=m_data[,1],ranges=IRanges(as.numeric(as.character(m_data[,2])),end=as.numeric(as.character(m_data[,3]))))
	snpbed=data.frame(do.call(rbind,strsplit(as.character(correlation_result$scoredRegionName),'_')))
	snpbed.dupli=which(duplicated(snpbed))
	if (length(snpbed.dupli)>0){
	   snpbed=snpbed[-snpbed.dupli,]
	}
	snps=GRanges(seq=snpbed[,1],ranges=IRanges(as.numeric(as.character(snpbed[,2])),
		end=as.numeric(as.character(snpbed[,3])),names=paste(snpbed[,1],'_',snpbed[,2],'_',snpbed[,3],sep='')))
	#overlap them to get motif matches
	snp.vs.motif.match=as.data.frame(findOverlaps(query=motif_matches,subject=snps,select='all'))
	duplicated_stuff=which(duplicated(snp.vs.motif.match[,2]))
	if (length(duplicated_stuff)>0){
	   snp.vs.motif.match=snp.vs.motif.match[-which(duplicated(snp.vs.motif.match[,2])),] #remove times when multiple motifs in the same SNP
	}
	snp.vs.motif.match=cbind(snp.vs.motif.match,snp_name=rownames(as.data.frame(snps)[snp.vs.motif.match[,2],]))
	matched.snp.idx=data.frame(snp.vs.motif.match)[,2]
	matched.snp=as.data.frame(snps)[matched.snp.idx,]
	matched.snp.name=paste(matched.snp[,1],'_',matched.snp[,2],'_',matched.snp[,3],sep='')
	matched_rows=intersect(matched.snp.name,as.character(correlation_result$scoredRegionName))
	correlation_result[which(correlation_result$scoredRegionName %in% matched_rows),'motifMatch']=TRUE
	correlation_result=cbind(correlation_result,motif_idx=NA)
	rownames(snp.vs.motif.match)=snp.vs.motif.match$snp_name
	correlation_result[,'motif_idx']=snp.vs.motif.match[as.character(correlation_result$scoredRegionName),1]

	#Keep only items with variable motif score and with a motif match.
	varScore.and.match=intersect(which(correlation_result$scoredRegionName %in% names(motif.score.variable)),
								   which(correlation_result$motifMatch==TRUE))
	correlation_result=correlation_result[varScore.and.match,]
	correlation_result[,'maxObsMotif']=motif.score.max[as.character(correlation_result$scoredRegionName)]
	### ## new! set motifMatch ==FALSE if max motif score <=0
	noMatch=which(as.numeric(as.character(correlation_result[,'maxObsMotif']))<=0)
	correlation_result[noMatch,'motifMatch']=FALSE
	true_matches=which(correlation_result$motifMatch==TRUE)
	correlation_result=correlation_result[true_matches,]
	#####take the top quantile scores
	####motif_quantile_value=quantile(as.numeric(as.character(correlation_result$maxObsMotif)),1-motif_q/100)
	####correlation_result=correlation_result[which(as.numeric(as.character(correlation_result$maxObsMotif))>=motif_quantile_value),]
	####print(paste("motif score threshold: ",motif_quantile_value,sep=''))

	motif.scores=motif.score.avg
	psignal.scores=psignal

	print('Starting motif correlations')
	print(paste(length(unique(correlation_result$motif_idx)),' motif instances tested',sep=''))	
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
	correlation_result[,'site_p.sp']=as.vector(do.call(rbind,lapply(1:nrow(correlation_result),get_site_p,permTable=full_perm,correlation_result=correlation_result)))
	
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


