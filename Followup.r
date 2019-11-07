## code useful in annotation data and plotting


library('biomaRt')
ensembl <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))
filters = listFilters(ensembl)

tmpfolder <- "../output/"
ph='y';peakI <- c(39,52);CHR <- '15'; 

Pcutoff = 4

diplo<- read.table(file=paste0(tmpfolder,"sigSNP_info_v3.txt"),header=TRUE,sep="\t",stringsAsFactors=FALSE, colClasses=c("sdp"="character"))

diplo <- diplo[which(diplo$gene_name!='unknown_intergenic'),c('pos','sdp','logP.merge','logP.founder','gene_name','cq','transcript_name')]
write.table(diplo,file=paste0(tmpfolder,"sigSNP_v4out.txt"),quote=FALSE,col.names = TRUE,row.names = FALSE,sep="\t")

#### New add gene name and other info #####
for (i in 1:dim(diplo)[1]){ diplo$gene_name[i] <-  strsplit(diplo$gene_name[i],split=',')[[1]][1] }

geneanno <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","mgi_symbol","description"),values=unique(diplo$gene_name),mart= ensembl)

diplo <- merge(diplo,geneanno,by.x='gene_name',by.y='ensembl_gene_id')
diplo <- diplo[order(-diplo$logP.merge),c(2,3,4,5,1,8,9,6,7)]

write.table(diplo,file=paste0(tmpfolder,'sigSNP_v5out.txt'),quote=FALSE,col.names = TRUE,row.names = FALSE,sep="\t")


diplo$posm<-diplo$pos/10^6
geno <- diplo[which(diplo$logP.merge>=Pcutoff),]
geno <- geno[which(geno$posm>=peakI[1] & geno$posm<=peakI[2]),]

diplo<- read.table(file=paste0(tmpfolder,"diplo_res_sample.txt"),header=TRUE,stringsAsFactors=FALSE)
diplo$posm<-diplo$pos/10^6
diplo <- diplo[order(diplo$posm),c('posm','efs')]
#diplo <- diplo[sort(sample(1:length(diplo[[1]]),1000)),c("posm","logP.founder")]
colnames(diplo)[2] <- 'logP.founder'

# order sdp by number of variants
StrDis_df <- aggregate(logP.merge ~ sdp, data=geno,length )
StrDis_df <- StrDis_df[order(StrDis_df$logP.merge),]
StrDis <- StrDis_df$sdp

library(RColorBrewer)
rincol <- brewer.pal(8, 'Dark2')[1:8];# FIND Eight color here
rincol <- c(rep('darkblue',(length(StrDis)-8)),rincol)
names(rincol) <- StrDis;
geno$col <- rincol[geno$sdp]

library(scales)

pdf(paste0(tmpfolder,ph,"_merge_",peakI[1],"_",peakI[2],".pdf"),width=18,height=12)
par(mfrow=c(2,1),mar=c(4,5,0.5,1.1))

plot(diplo$posm,diplo$logP.founder,xlim=peakI,ylim=c(0,max(geno$logP.merge,diplo$logP.founder)),
     col='red',ylab="-log10(Pvalue)",xlab=paste0("Mb on chr",CHR),type='l',lwd=2)
points(diplo$posm,diplo$logP.merge,pch=20,col='grey',cex=0.8)
points(geno$posm,geno$logP.merge,pch=20,col=geno$col)
#abline(v=peakP,lty=3,lwd=2)

blcknm <- length(StrDis)
plot(geno$posm,geno$posm,xlim=peakI,ylim=c(0,(blcknm-0.5)),type='n',yaxt = "n",xlab=paste0('Mb on chr',CHR),ylab='')
axis(2, at = seq(0.5, (blcknm-0.5) , by = 1), label = StrDis, las=1,cex.axis=1)
for (ii in c(1:floor(blcknm/2))){
  rect(peakI[1], (2*ii-1), peakI[2], (2*ii), col='gray', angle = 0,lwd=0)
}

for ( i in c(1:length(geno$posm))){
  if (geno$sdp[i] %in% StrDis){
    lineid=which(StrDis==geno$sdp[i])
    segments(geno$posm[i], (lineid-1) ,geno$posm[i],(lineid),col=alpha(geno$col[i]),lwd=2)
  }
}

dev.off()



