
gene_res_filter <- function(peakI,tmp.file.path=tmp.file.path){
GR <- data.frame()
for (i in c(peakI[1]:peakI[2])){
	tmp <-read.table(file=paste0(tmp.file.path,"geno_res",i,".txt"),header=TRUE)
	GR <- rbind(GR,tmp[tmp[[2]]>3,])
}
write.table(GR,file=paste0(tmp.file.path,"gene_res_filter.txt"),row.names=FALSE,col.names=TRUE,quote=FALSE,sep="\t")
}


diplo_res_filter <- function(peakI,tmp.file.path=tmp.file.path){
GR <- data.frame()
for (i in c(peakI[1]:peakI[2])){
	message(i)
        tmp <-read.table(file=paste0(tmp.file.path,"diplo_res",i,".txt"),header=TRUE,colClasses=c("sdp"="character"))

	tmpout <- data.frame()
        for (mypos in sort(unique(tmp$pos))){
		message(mypos)
		this <- tmp[tmp$pos==mypos,]	
		if (dim(this)[1]>1){
			out <- this[which.max(this$prob),]
			out$fs <- median(this$fs[sample(c(1:dim(this)[1]),size=20,prob=this$prob,replace=TRUE)])
			this <- out
		}	
		tmpout <- rbind(tmpout,this)
	}
	tmpout2 <- tmpout[tmpout$fs<2 | tmpout$fs < tmpout$efs,]
        tmpout2 <- tmpout2[sample(1:dim(tmpout2)[1],100),]

        tmpout <- tmpout[tmpout$fs>2 & tmpout$fs > tmpout$efs,]
	GR <- rbind(GR,tmpout,tmpout2)
	write.table(GR,file=paste0(tmp.file.path,"diplo_res_filter.txt"),row.names=FALSE,col.names=TRUE,quote=FALSE,sep="\t")
}
	write.table(GR,file=paste0(tmp.file.path,"diplo_res_sample.txt"),row.names=FALSE,col.names=TRUE,quote=FALSE,sep="\t")

	GR <- GR[GR$fs>2 & GR$fs > GR$efs,]
	GR <- GR[order(-GR$fs),]	
	write.table(GR,file=paste0(tmp.file.path,"diplo_res_filter.txt"),row.names=FALSE,col.names=TRUE,quote=FALSE,sep="\t")
}

diplo_res_sample <- function(peakI,tmp.file.path=tmp.file.path){

GR <- data.frame()
for (i in c(peakI[1]:peakI[2])){
        tmp <-read.table(file=paste0(tmp.file.path,"diplo_res",i,".txt"),header=TRUE)
        GR <- rbind(GR,tmp[sort(sample(nrow(tmp), 500)),c(1,5)])
}

write.table(GR,file=paste0(tmp.file.path,"diplo_res_sample.txt"),row.names=FALSE,col.names=TRUE,quote=FALSE,sep="\t")

}


