source("sigSNP_info_fast.r")

sig_SNP <- function(CHR,peakI,genodb,tmp.file.path){

indf <- read.table(file=paste0(tmp.file.path,"gene_res_filter.txt"),header=TRUE)[,1:2]
colnames(indf) <- c('pos','nlog10P')
indf$chr <- CHR
chr_list <- sort(unique(indf$chr))

file_d=genodb
library(data.table)

res<-data.frame()
for (CHR in chr_list){
	indf_sub <- indf[which(indf$chr==CHR),]
	# Get the strain distribution pattern

	ccgom <- data.table()
        for (CC in c("A_J","C57BL6J","129S1_SvImJ","NOD_ShiLtJ","NZO_HlLtJ","CAST_EiJ","PWK_PhJ","WSB_EiJ")){
                message(CC)
                filepath<- paste0(file_d,"genotype/",CC,"/",CHR,".txt.tar.gz")
                ccgt<-data.table(read.table(filepath,header=TRUE,skip=1,sep='\t'),stringsAsFactors=FALSE)
		setkey(ccgt, pos)
		ccgt <- ccgt[pos>peakI[1]*10^6 & pos < (peakI[2]+1)*10^6]
		
		#ccgt <- ccgt[which(ccgt$allele_1!='None' & ccgt$allele_2!='None'),]
		
		ccgt <- ccgt[pos %in% indf_sub[[1]]]
                
		ccgt<-merge( data.frame(indf_sub,CC=CC),ccgt)
                ccgom <- rbind(ccgom,ccgt)
        }
        ccgom<-ccgom[prob>0.5]
        genotype_founder_list<-ccgom[order(pos)]
	
	foundertd <- data.frame()
	for (mypos in indf_sub$pos){
		message(mypos)
		indfsubsub <- genotype_founder_list[pos==mypos & prob ]
		indfsubsubsub <- unique(indfsubsub[,1:8])
		#replace None to bls
		indfsubsub[allele_1=='None']$allele_1 <- indfsubsub[CC=='C57BL6J']$allele_1
		indfsubsub[allele_2=='None']$allele_2 <- indfsubsub[CC=='C57BL6J']$allele_2
		indfsubsub_orgin  <- indfsubsub
		indfsubsub <- unique(indfsubsub[,1:8])	

		if (dim(indfsubsub)[1]==8){
			geno_list <- paste(indfsubsub$allele_1, indfsubsub$allele_2,sep='/')
			ref <- indfsubsub$allele_1[2]
			# problems here
			allall <- sort(unique(   c(  as.character(indfsubsub$allele_1) , as.character(indfsubsub$allele_2))  )  )
			alt <- paste(allall[which(allall!=ref)],collapse=",")
			dp<-rep(2,8)
			allall <- sort(unique(geno_list))
			# dispattern #0/1 by default
			if (length(allall)==2){
				dp[which(geno_list==geno_list[2])]<-0
				dp[which(geno_list!=geno_list[2])]<-1
			}
			if (length(allall)==3){
				dp[which(geno_list==geno_list[2])]<-0
				dp[which(geno_list==geno_list[which(geno_list!=geno_list[2])[1]])]<-1	
			}
			cq <- unique(c(as.character(indfsubsub_orgin[consequence_1!='reference']$consequence_1), as.character(indfsubsub_orgin[consequence_2!='reference']$consequence_2)))
			cq <- paste(cq,collapse=",")
			gene_name <- paste(unique(indfsubsub_orgin$gene_name),collapse=",")
			transcript_name <- paste(unique(indfsubsub_orgin$transcript_name),collapse=",")
		
			gene_name <- ifelse(gene_name=="","unknown",gene_name)
			transcript_name <- ifelse(transcript_name=="","unknown",transcript_name)
	
			dps <- paste(dp,collapse="")

			geno_list <- paste(indfsubsubsub$allele_1, indfsubsubsub$allele_2,sep='/')
			foundertd <- rbind(foundertd,data.frame(CHROM=CHR,pos=mypos,REF=ref,ALT=alt,t(geno_list),sdp=dps,logP.merge=indfsubsub$nlog10P[1],gene_name,cq,transcript_name))
		}
	}

	res<-rbind(res,foundertd)
}

colnames(res)[5:12] <-c("AJ.geno","B6.geno","X129.geno","NOD.geno","NZO.geno","CAST.geno","PWK.geno","WSB.geno")
res <- res[order(-res$logP.merge),]

write.table(res,file=paste0(tmp.file.path,"sigSNP_info_v1.txt"),sep="\t",col.names=TRUE, row.names=FALSE,quote=FALSE)

nn <- 'None/None'
res2 <- res[which(res[5]!=nn & res[7]!=nn & res[8]!=nn & res[9]!=nn & res[10]!=nn & res[11]!=nn & res[12]!=nn),] 

write.table(res2,file=paste0(tmp.file.path,"sigSNP_info_v2.txt"),sep="\t",col.names=TRUE, row.names=FALSE,quote=FALSE)

}



