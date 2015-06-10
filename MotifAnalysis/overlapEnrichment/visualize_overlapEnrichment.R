read_enrichment_file=function(infile,SIG_THRESH,fillin,col_of_interest){
  data=read.table(infile,header=TRUE)
  data.enrich=data.frame(enrichment=data[,col_of_interest])
  rownames(data.enrich)=data$tf
  to_remove=which(as.numeric(as.character(data$BH))>SIG_THRESH)
  #keep only the significant ones
  if (length(to_remove)>0){
    data.enrich[to_remove,'enrichment']=fillin
  }
  return(data.enrich)
}

one_enrichment_plot=function(outpdf,data,top_value){
  print('doing enrichment plot')
  print(head(data))
  p300=which(grepl('EP300',as.character(data$tf)))
  if (length(p300)!=0){
    data=data[-p300,]
  }
  require(ggplot2)
  data[,'tf']=gsub('bed.OverlapChIPseq','',
                            gsub('MotifMatch_','',
                              gsub('MergedPeaks_ChIPseq_','',
                                gsub('correlatedMotifs.motif.pouya.Motif.','',
                                  gsub('scanThresh0','',data$tf)))))
    significance=(data$BH<=0.05)
    sig_vector=rep('Not significant',times=dim(data)[1])
    sig_vector[which(significance==TRUE)]='Significant'
    data=cbind(data,Significant=factor(sig_vector, levels=c('Significant','Not significant')),
      TFname=factor(data$tf,levels=data[order(data$enrichment),'tf']))
    pdf(outpdf,width=5,height=12)
    print(ggplot(data, aes(y=enrichment, x=TFname,col=Significant))+
      coord_flip()+geom_point()+scale_colour_manual(values=c("red","gray"))+geom_errorbar(aes(ymax = confLow,
      ymin=confHigh))+ylim(0,top_value)+theme_bw() + ylab('Enrichment of TF in QTL peaks') + xlab('TF')+theme(panel.border = element_blank(), 
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")))
    dev.off()
}

optimal_ordering=function(m,meth){
  require(cba)
  d <- dist(as.matrix(m),method=meth)
  hc <- hclust(d)
  co <- order.optimal(d, hc$merge)
  m.optimalRows=as.matrix(m)[co$order,]
  return(m.optimalRows)
}

order_by_column=function(m,column_name,decreasing_desired){
  return(m[order(m[,column_name],decreasing=decreasing_desired),])
}

heatmap_enrichments=function(data,out,meth,top_value){
  #data.optimal=optimal_ordering(data,meth)
  data.optimal=data
  require(pheatmap)
  pdf(out,width=10,height=20)
  pheatmap(as.matrix(data.optimal),cluster_rows=FALSE,cluster_cols=FALSE,fontsize=10,breaks=seq(from=1,to=top_value,by=(top_value-1)/20),
        color=colorRampPalette(c("white", "red",'black'))(n = 20),cellwidth=10,cellheight=10)#,legend_breaks=seq(from=1,to=top_value,by=1))
  #pheatmap(-log(as.matrix(data.optimal)),cluster_rows=FALSE,cluster_cols=FALSE,fontsize=5,
  #  cellwidth=5,cellheight=5,breaks=seq(from=1,to=600,by=1),
  #     color=colorRampPalette(c("gray", "red","black"))(n = 600))
 
  dev.off()
}

overlapEnrichment_distalQTL=function(){
  enrichfiles='/srv/gsfs0/projects/snyder/oursu/histoneQTL/motif_analysis/results/2015-05-30/OverlapEnrichment/ENRICHDIR/ENRICHPREF'
  enrichments=c('TFBS_overlap_','Motif_overlap_','MotifCorrelatedLocal_overlap_') #add in disrupted motif overlaps, and hQTL overlaps

  hmarks=c('H3K27AC','H3K4ME1','H3K4ME3')
  for (enrich in enrichments){
    if (enrich=='TFBS_overlap_'){
      top_value=3
    }
    if (enrich=='Motif_overlap_'){
      top_value=5
    }
    if (enrich=='MotifCorrelatedLocal_overlap_'){
      top_value=25
    }
    print(enrich)
    first=TRUE
    for (suffix in c('HMARK.QTLpeaks','LocalPeakIsHMARK.QTLpeaks_affectingDistalPeaks')){
      for (hmark in hmarks){
        f=gsub('ENRICHDIR',paste(enrich,'QTLpeaks0kb',sep=''),
                            gsub('ENRICHPREF',paste(enrich,hmark,'QTLpeaks0kb___.overlapEnrichIN',gsub('HMARK',hmark,suffix),sep=''),enrichfiles))
        cur_data=read_enrichment_file(f,0.05,1,'enrichment')
        cur_total=read.table(f,header=TRUE)
        #cur_data=cur_data/max(cur_data[,1])
        if (suffix=='HMARK.QTLpeaks'){
          addon='Local: HMARK'
        }
        if (suffix=='LocalPeakIsHMARK.QTLpeaks_affectingDistalPeaks'){
          addon='Distal: HMARK'
        }
        rownames(cur_data)=gsub('bed.OverlapChIPseq','',
                            gsub('MotifMatch_','',
                              gsub('MergedPeaks_ChIPseq_','',
                                gsub('correlatedMotifs.motif.pouya.Motif.','',
                                  gsub('scanThresh0','',rownames(cur_data))))))
        #condition=paste(gsub('_',' ',enrich),hmark,' ',addon,sep='')
        condition=gsub('HMARK',hmark,addon)
        cur_total=cbind(cur_total,condition=condition)
        one_enrichment_plot(paste(dirname(f),condition,'overlapEnrichmentPlot.pdf',sep=''),cur_total,top_value)
        if (first==FALSE){
          data=cbind(data,cur_data[rownames(data),])
          total=rbind(total,cur_total)
          colnames(data)[ncol(data)]=condition
        }
        if (first==TRUE){
         data=cur_data
         total=cur_total
          first=FALSE
          colnames(data)[1]=condition
        }
      }
    }

    ###### heatmap ###################################################
    #add in the k27AC again, to sort by it and its pvalue
    sortf=gsub('ENRICHDIR',paste(enrich,'QTLpeaks0kb',sep=''),
                            gsub('ENRICHPREF',paste(enrich,'H3K27AC','QTLpeaks0kb___.overlapEnrichIN',gsub('HMARK','H3K27AC','HMARK.QTLpeaks'),sep=''),enrichfiles))
    k27ac_data=read_enrichment_file(sortf,1.1,1,c('enrichment'))
    rownames(k27ac_data)=gsub('bed.OverlapChIPseq','',
                            gsub('MotifMatch_','',
                              gsub('MergedPeaks_ChIPseq_','',
                                gsub('correlatedMotifs.motif.pouya.Motif.','',
                                  gsub('scanThresh0','',rownames(k27ac_data))))))
    k27ac_sorted_rows=rownames(k27ac_data)[order(k27ac_data[,1],decreasing=TRUE)]
      
    data=data[k27ac_sorted_rows,]
    heatmap_enrichments(data,paste(dirname(f),'overlapEnrichmentHeatmap.pdf',sep=''),'euclidean',top_value)
    ####################################################################

  }
}

overlapEnrichment_distalQTL()

