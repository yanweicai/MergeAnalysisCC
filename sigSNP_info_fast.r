sig_SNP_fast <- function(CHR,peakI,genodb,tmp.file.path){

indf <- read.table(file=paste0(tmp.file.path,"diplo_res_filter.txt"),header=TRUE,colClasses=c('sdp'='character'))
file_d=genodb
library(data.table)

indf_sub <- indf[,1:2]


ccgom <- data.table()
for (CC in c("A_J","C57BL6J","129S1_SvImJ","NOD_ShiLtJ","NZO_HlLtJ","CAST_EiJ","PWK_PhJ","WSB_EiJ")){
	message(CC)
        filepath<- paste0(file_d,"genotype/",CC,"/",CHR,".txt.tar.gz")
       	ccgt<-data.table(read.table(filepath,header=TRUE,skip=1,sep='\t'),stringsAsFactors=FALSE)
        setkey(ccgt, pos)
        ccgt <- ccgt[pos>peakI[1]*10^6 & pos < (peakI[2]+1)*10^6]
	ccgt <- ccgt[pos %in% indf_sub[[1]]]

	ccgt<-merge( data.frame(indf_sub,CC=CC),ccgt)
        ccgom <- rbind(ccgom,ccgt)
}
ccgom<-ccgom[prob>0.5]
genotype_founder_list<-ccgom[order(pos)]
	
foundertd <- data.frame()
for (mypos in indf_sub$pos){
	indfsubsub <- genotype_founder_list[pos==mypos]

	cq <- unique(c(as.character(indfsubsub$consequence_1), as.character(indfsubsub$consequence_2)))
	cq <- unique(unlist(strsplit(cq,'&')))
	if ('None' %in% cq){ 
		cq <- cq[-which(cq=='None')]
		cq <- c(cq,'unknown')
	}
	if ('reference' %in% cq){ cq <- cq[-which(cq=='reference')] }
	cq <- paste(cq,collapse=",")
	gene_name <- paste(unique(indfsubsub$gene_name),collapse=",")
	transcript_name <- paste(unique(indfsubsub$transcript_name),collapse=",")
		
	gene_name <- ifelse(gene_name=="","unknown_intergenic",gene_name)
	transcript_name <- ifelse(transcript_name=="","unknown",transcript_name)
	foundertd <- rbind(foundertd,data.frame(pos=mypos,gene_name,cq,transcript_name))
}
indf <- merge(indf,foundertd,by='pos')

colnames(indf)[4:5] <- c('logP.merge','logP.founder')
colnames(indf)[6:13] <-c("AJ.geno","B6.geno","X129.geno","NOD.geno","NZO.geno","CAST.geno","PWK.geno","WSB.geno")
indf <- indf[order(-indf$logP.merge),]
indf$sdp_c <- indf$sdp

sdpc2sdp <- function(x){
	out <- c()
	for (sdp in x){
		sdpv <- strsplit(sdp,split="")[[1]]
		sdpvu <- unique(sdpv)
		sdp_dic <- as.character(0:(length(sdpvu)-1))
		names(sdp_dic) <- sdpvu
		out <- c(out,paste0(sdp_dic[sdpv],collapse=''))
	}
	out
}

indf$sdp <- sdpc2sdp(indf$sdp_c)

write.table(indf,file=paste0(tmp.file.path,"sigSNP_info_v3.txt"),sep="\t",col.names=TRUE, row.names=FALSE,quote=FALSE)

}



