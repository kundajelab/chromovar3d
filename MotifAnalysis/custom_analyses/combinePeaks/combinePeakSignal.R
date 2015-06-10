inTemplate='/srv/gs1/projects/snyder/jzaugg/histoneQTL/hQTL_analysis/PeakData/hMat.norm.ALL.HMARK.peer_lS_PEERS.txt'
out='/srv/gsfs0/projects/snyder/oursu/histoneQTL/motif_analysis/results/2015-01-13/Data/CombinedPeakSignal'

get_data=function(inTemplate,hmark,peers,columnsWanted){
	inFile=gsub('HMARK',hmark,gsub('PEERS',peers,inTemplate))
	data=read.table(inFile,header=TRUE)
	rownames(data)=paste(hmark,'_',as.character(rownames(data)),sep='')
	if (length(columnsWanted)>0){
		columns_needed=setdiff(columnsWanted,colnames(data))
		for (colName in columns_needed){
			data[,colName]=NA
		}
	}
	return(data)
}

#combine signal from all peaks
combined=get_data(inTemplate,'H3K4ME1','5',c())
columnsWeWant=colnames(combined)
combined=rbind(combined,get_data(inTemplate,'H3K4ME3','5',columnsWeWant)[,colnames(combined)])
combined=rbind(combined,get_data(inTemplate,'H3K27AC','5',columnsWeWant)[,colnames(combined)])
combined=rbind(combined,get_data(inTemplate,'dhs','5',columnsWeWant)[,colnames(combined)])
combined=rbind(combined,get_data(inTemplate,'RNA','10',columnsWeWant)[,colnames(combined)])

#write results
write.table(combined,
	file=out,sep='\t',quote=F,row.names=TRUE,col.names=TRUE)
system(paste('rm ',out,'.gz',sep=''))
system(paste('gzip ',out,sep=''))


