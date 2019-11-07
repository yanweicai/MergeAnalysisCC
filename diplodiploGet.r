

diplo.region.isvdb <- function (df, genodb, CHR , peakI , tmp.file.path='./'){ 
	STRbp <- peakI[1]*10^6
	ENDbp <- (peakI[2]+1)*10^6

	`%+=%` = function(e1,e2) eval.parent(substitute(e1 <- e1 + e2))
	
	library(data.table)
	# find collums of CCname indexing
	CCinx <- which(colnames(df)=='CC')
	CClist <- sort(unique(df[[CCinx]])) # Find collum of CCname
	# Split CCrix by spliting 'x' 
	IS.RIX <- 'x' %in% unlist(strsplit(as.character(CClist),''))

	if (IS.RIX){
		CCrix <- CClist
		CClist <- sort(unique(unlist(strsplit(CCrix,split='x'))))
		message("Generating genotype probability file for the region CCrix. (Step1 of 2)")
		CCout<-data.table()
		for (CC in CClist){
			message(CC)
			CCgs<-read.table(file=paste0(genodb,'diplotypeSampling/',CC,"/",CHR,'.txt.tar.gz'),skip=1,header=TRUE,sep="\t")
			CCgs <- data.table(CCgs)
			CCgs <- CCgs[pos>=STRbp & pos<=ENDbp]	
			CCout<-rbind(CCout,data.table(CC=CC,CCgs))
		}
		#write.table(CCout,file=paste0(tmp.file.path,"diplotype_tmp_dsampling.txt"),quote=FALSE,col.names=TRUE,row.names=FALSE,sep="\t")

		message("Calculating cross for CCrix. (Step 2 of 2)")
		#CCout <- data.table(read.table(file="genotype_2_94_110.txt",header=TRUE,sep="\t"))
		CCrixout<-data.frame()
		for (CCp in CCrix){
			message(CCp)
			CCps <- strsplit(CCp,split='x')[[1]]
			CC1 <- CCps[1]; CC2 <- CCps[2]
			CC1df <- as.data.frame(CCout[CC==CC1]); CC2df<- as.data.frame(CCout[CC==CC2])
			mergedf<-merge(CC1df,CC2df,by=c("variant_id","pos","gene_name"))
			mergedf2 <- cbind(mergedf[,c("variant_id","pos","founder.x","founder.y")],
					data.frame(prob = mergedf$prob.x * mergedf$prob.y),
		   			data.frame(gene_name=mergedf$gene_name))
			message(dim(mergedf2)[1])
			CCrixout<-rbind(CCrixout,data.frame(CCrix=CCp, mergedf2))
		}
		#write.table(CCrixout,file=paste0(tmp.file.path,"CCrix_geno_tmp.txt"),quote=FALSE,col.names=TRUE,row.names=FALSE,sep="\t")
		CCout <- unique(CCrixout[,1:6])
		colnames(CCout)[4:5] <- c('founder_1','founder_2')
		colnames(CCout)[1] <- 'CC'
	}else{
		CCout<-data.table()
		message("Generating diplotype probability file for the region.")
                for (CC in CClist){
                        message(CC)
                        CCgs<-read.table(file=paste0(genodb,'diplotype/',CC,"/",CHR,'.txt.tar.gz'),skip=1,header=TRUE,sep="\t")
                        CCgs <- data.table(CCgs)
                        CCgs <- CCgs[pos>=STRbp & pos<=ENDbp]
                        CCout<-rbind(CCout,data.table(CC=CC,CCgs))
                }
#                write.table(CCout,file=paste0(tmp.file.path, "diplotype_tmp.txt"),quote=FALSE,col.names=TRUE,row.names=FALSE,sep="\t")	
		CCout <- unique(CCout[,c(1:6)])
	}
	

	#CCout <- fread(file=paste0(tmp.file.path, "diplotype_tmp.txt"),header=TRUE)
	CCout$A=0;CCout$B=0;CCout$C=0;CCout$D=0;CCout$E=0;CCout$F=0;CCout$G=0;CCout$H=0;
	
	CCout[which(CCout$founder_1=='A_J'),'A'] %+=% 1
	CCout[which(CCout$founder_1=='C57BL6J'),'B'] %+=% 1
	CCout[which(CCout$founder_1=='129S1_SvImJ'),'C'] %+=% 1
	CCout[which(CCout$founder_1=='NOD_ShiLtJ'),'D'] %+=% 1
	CCout[which(CCout$founder_1=='NZO_HlLtJ'),'E'] %+=% 1
	CCout[which(CCout$founder_1=='CAST_EiJ'),'F'] %+=% 1
	CCout[which(CCout$founder_1=='PWK_PhJ'),'G'] %+=% 1
	CCout[which(CCout$founder_1=='WSB_EiJ'),'H'] %+=% 1
	
	CCout[which(CCout$founder_2=='A_J'),'A'] %+=% 1
        CCout[which(CCout$founder_2=='C57BL6J'),'B'] %+=% 1
        CCout[which(CCout$founder_2=='129S1_SvImJ'),'C'] %+=% 1
        CCout[which(CCout$founder_2=='NOD_ShiLtJ'),'D'] %+=% 1
        CCout[which(CCout$founder_2=='NZO_HlLtJ'),'E'] %+=% 1
        CCout[which(CCout$founder_2=='CAST_EiJ'),'F'] %+=% 1
        CCout[which(CCout$founder_2=='PWK_PhJ'),'G'] %+=% 1
        CCout[which(CCout$founder_2=='WSB_EiJ'),'H'] %+=% 1

	CCout$A <- CCout$A * CCout$prob;
        CCout$B <- CCout$B * CCout$prob;
        CCout$C <- CCout$C * CCout$prob;
        CCout$D <- CCout$D * CCout$prob;
        CCout$E <- CCout$E * CCout$prob;
        CCout$F <- CCout$F * CCout$prob;
        CCout$G <- CCout$G * CCout$prob;
        CCout$H <- CCout$H * CCout$prob;

	CCoutput <- CCout[,c('CC','variant_id','pos','A','B','C','D','E','F','G','H')]
	write.table(CCoutput,file=paste0(tmp.file.path, "diplotype_tmp.txt"),quote=FALSE,col.names=TRUE,row.names=FALSE,sep="\t")

}

diplo.founder.isvdb <- function(CHR,peakI,tmp.file.path){
	library(data.table)
	FD <- fread(file=paste0("/nas/depts/006/valdar-lab/PUBLIC/2018_1_isvdb/full_1504/","temp/founders/chr_",CHR,".txt"),header=TRUE,sep='\t')
	FD <- FD[pos >= (peakI[1]*10^6) & pos < (peakI[2]+1)*10^6]	
	FD <- unique(FD[,1:6])
	FD <- FD[order(pos),]
	write.table(FD,file=paste0(tmp.file.path, "founder_tmp.txt"),quote=FALSE,col.names=TRUE,row.names=FALSE,sep="\t")
}

