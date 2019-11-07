vec2sdp<- function(sdpvec){
        sdpvec<- c(unlist(sdpvec))
        names(sdpvec) <- c(0:7)
        for (thisi in 8:1){
                sdpvec[which(sdpvec==sdpvec[thisi])] <- names(sdpvec)[thisi]
        }
        names(sdpvec) <- c('A','B','C','D','E','F','G','H')
        sdpvec
}

sdpCon <- c('A','B','C','D','E','F','G','H');
names(sdpCon) <- c(0:7)
sdp2form<- function(sdp){
        mysdp <- sdpCon[sdp]
        names(mysdp) <- c('A','B','C','D','E','F','G','H')
        if (length(unique(mysdp))==1){
                fomula <- 'Pheno ~ 1'
        }else{
                fomula <- paste0('Pheno',' ~ ',paste(unique(mysdp)[-1],collapse=' + '))
        }
        as.formula(fomula)
}

diplotype_merge_mapping <- function(Phenofile,diplotype_file,MB,K=NULL,tmp.file.path="./"){

        for (nm in list.files("~/SNPmerge/flu/miqtl-master/R/", pattern = "\\.[RrSsQq]$")) {
                source(file.path("~/SNPmerge/flu/miqtl-master/R/", nm))
        }

diplotype_list <- fread(file=diplotype_file,header=TRUE,sep="\t",fill = TRUE)
diplotype_list <-diplotype_list[pos >= (MB[1]*10^6) & pos < (MB[length(MB)]+1)*10^6]

FD <- fread(file=paste0(tmp.file.path,"founder_tmp.txt"),header=TRUE,sep='\t')
FD <- FD[pos >= (MB[1]*10^6) & pos < (MB[length(MB)]+1)*10^6]

# Filter 

poslist <- sort(unique(diplotype_list$pos))
counti <- 0
outdf <-data.frame()
for( mypos in poslist ){
	counti <- counti + 1
	diplotype_list_sub<-diplotype_list[which(pos==mypos)]

	sdp_table_full <- FD[pos==mypos]
	# variant id check
	vtb <- table(sdp_table_full$variant_id)
        if (length(vtb)>=2){next} # remove duplicated information points
        sdp_table <- unique(sdp_table_full[variant_id==names(which(vtb==max(vtb)))[1],3:6])

        if (all(c(sdp_table$allele_1!='None',sdp_table$allele_2!='None'))){CompleteN = TRUE
        }else{CompleteN = FALSE }

        if (!CompleteN) { next }

	# get dosage explaination
	sdp_table_dose <- sdp_table
	SNPlist <- sort(unique(c(sdp_table$allele_1,sdp_table$allele_2)))
	for (snp in SNPlist){
		thisdose <- (ifelse(sdp_table_dose$allele_1==snp,1,0)+ifelse(sdp_table_dose$allele_2==snp,1,0))*sdp_table_dose$prob
		sdp_table_dose <- cbind(sdp_table_dose,data.frame(snp=thisdose))
	}
	colnames(sdp_table_dose)[5:dim(sdp_table_dose)[2]]<-SNPlist	
	sdp_table_dose <- aggregate( as.formula(paste0('cbind(',paste(SNPlist,collapse=','),')~strain')) ,sdp_table_dose,sum)
	sdp_matrix <- as.matrix(sdp_table_dose[,-1])
	rownames(sdp_matrix) <- sdp_table_dose$strain
	sdp_matrix <- sdp_matrix[c('A_J','C57BL6J','129S1_SvImJ','NOD_ShiLtJ','NZO_HlLtJ','CAST_EiJ','PWK_PhJ','WSB_EiJ'),]/2
	
	TBfounder <- function(founder,fid){
        	TB <- sdp_table[strain==founder];
        	# replace None to NA_ABCDEFGH
	        TB$allele_1[which(TB$allele_1=='None')] <- paste0('NA_',fid)
	        TB$allele_2[which(TB$allele_2=='None')] <- paste0('NA_',fid)
	        # get heterzypos
	        TB$allele='None/None'
	        for (thisi in 1:dim(TB)[1]){
	                TB$allele[thisi] = paste(sort(c(unlist(TB[thisi,2:3]))),collapse='/')
	        }

        	colnames(TB)[4:5] <- c(paste0(fid,'p'),fid)
        	as.data.frame(TB[,c(5,4)])
	}

	sdp_sum <- list(TBfounder('A_J','A'),TBfounder('C57BL6J','B'),
			TBfounder('129S1_SvImJ','C'),TBfounder('NOD_ShiLtJ','D'),
			TBfounder('NZO_HlLtJ','E'), TBfounder('CAST_EiJ','F'),
			TBfounder('PWK_PhJ','G'),TBfounder('WSB_EiJ','H'))	
	sdp_sum <- Reduce(function(...) merge(..., all=T),sdp_sum)
	sdp_sum$prob <- sdp_sum$Ap * sdp_sum$Bp * sdp_sum$Cp * sdp_sum$Dp * sdp_sum$Ep * sdp_sum$Fp * sdp_sum$Gp * sdp_sum$Hp
	sdp_sum <- sdp_sum[,c('A','B','C','D','E','F','G','H','prob')]
	sdp_sum <- sdp_sum[sdp_sum$prob>0.001,]

	mydf <- aggregate( cbind(A,B,C,D,E,F,G,H) ~ CC, diplotype_list_sub,sum)
	mydf[,c('A','B','C','D','E','F','G','H')] <- 2*mydf[,c('A','B','C','D','E','F','G','H')]/apply(as.matrix(mydf[,c('A','B','C','D','E','F','G','H')]),1,sum)
	mydf -> mydf_clean

	dosedf <- cbind(data.frame(CC=mydf$CC),as.data.frame(as.matrix(mydf[,-1]) %*% sdp_matrix))
	dosedf <- merge(Phenofile,dosedf,by='CC')	
	# check for MAF
	maf <- apply(dosedf[,-c(1:3)],2,sum)
	maf_min <- min(maf/sum(maf))
		
	mydf <- merge(Phenofile,mydf,by='CC')	

	if (!is.null(K)){
		CClist <- intersect(unique(mydf$CC),rownames(K))	
		Ksub <- K[match(CClist,rownames(K)),match(CClist,colnames(K))]
		mydf <- mydf[which(mydf$CC %in% CClist),]

		fit0 <- lmmbygls.replicates( Pheno ~ 1, data=mydf,pheno.id="SUBJECT.NAME",geno.id="CC", little.K=Ksub)	
		fit1 <- lmmbygls.replicates( Pheno ~ A + B + C + D + E + F + G, data=mydf,pheno.id="SUBJECT.NAME",geno.id="CC", little.K=Ksub)
		LHR<-2*(fit1$logLik - fit0$logLik)

		fit2 <- lmmbygls.replicates( as.formula( paste0('Pheno ~',paste(SNPlist[-1],collapse='+'))), data=dosedf,pheno.id="SUBJECT.NAME",geno.id="CC", little.K=Ksub)
		LHR2 <- 2*(fit2$logLik - fit0$logLik)
		df2 <- length(SNPlist)-1
	}else{
		fit0 <- lm(Pheno ~ 1,data=mydf,weights=NUM.OBS)
		fit1 <- lm(Pheno ~  A + B + C + D + E + F + G ,data=mydf, weights=NUM.OBS)
		LHR <- 2*(as.numeric(logLik(fit1))-as.numeric(logLik(fit0)))

	#	fit2 <- lm(as.formula( paste0('Pheno ~',paste(SNPlist[-1],collapse='+'))), data=mydf, weights=NUM.OBS)
	#	LHR2 <- 2*(as.numeric(logLik(fit2))-as.numeric(logLik(fit0)))
	#	df2 <- length(SNPlist)-1
	}

	df <- 7 - ifelse( all(mydf$A==0) ,1,0) -ifelse( all(mydf$B==0) ,1,0) - ifelse( all(mydf$C==0) ,1,0) - ifelse( all(mydf$D==0) ,1,0) - ifelse( all(mydf$E==0) ,1,0) - ifelse( all(mydf$F==0) ,1,0) - ifelse( all(mydf$G==0) ,1,0)
	fs <- -log10(1-pchisq(LHR, df))
	#fsd <- -log10(1-pchisq(LHR2, df2))
	rm(fit1);

	# sdp model
	for (i in 1:dim(sdp_sum)[1]){
		sdp_vec <- sdp_sum[i,1:8]
		sdp_prob <- sdp_sum[i,9]
		
		if (length(unique(as.character(sdp_vec)))==1){
			next
		}

		sdp <- vec2sdp(sdp_vec)
		sdpform <- sdp2form(sdp)

		thisdf <- mydf
		for (ff in unique(sdp)){
			fflist <- names(sdp)[which(sdp==ff)]
			if (length(fflist)>1){
				thisdf[,sdpCon[ff]]<- apply(thisdf[,names(sdp)[which(sdp==ff)]],1,sum)
			}
		}
		thisdf <- cbind(thisdf[,1:3], thisdf[,sdpCon[unique(sdp)]])

		
		if (!is.null(K)){
			fit1 <- lmmbygls.replicates(sdpform, data=thisdf,pheno.id="SUBJECT.NAME",geno.id="CC", little.K=Ksub)
			LHR<-2*(fit1$logLik - fit0$logLik)	
		}else{
			fit1 <- lm(sdpform, data=thisdf,weights=NUM.OBS)
			LHR <- 2*(as.numeric(logLik(fit1))-as.numeric(logLik(fit0)))
        	}
		df <- length(unique(sdp))-1
		thisfs <- -log10(1-pchisq(LHR, df))	
		thisdf <- data.frame(pos=mypos,sdp=paste0(sdp,collapse=''),prob=sdp_prob,fs=thisfs,efs=fs,maf=maf_min)
		outdf <- rbind(outdf,cbind(thisdf,sdp_vec))
	}
	
	if (!counti %% 1000) {
		print(counti)
		write.table(outdf,file=paste0(tmp.file.path,"diplo_res",MB,".txt"),col.names=TRUE,row.names=FALSE,quote=FALSE)
	}
}
	write.table(outdf,file=paste0(tmp.file.path,"diplo_res",MB,".txt"),col.names=TRUE,row.names=FALSE,quote=FALSE)
}

diplotype_timber_mapping <- function(Phenofile,diplotype_file,MB,K=NULL,tmp.file.path="./"){

prior.M.8allele <- list(model.type='fixed',M.IDs=c("0,1,2,3,4,5,6,7"), probs=c(1),hash.names=TRUE)

for (nm in list.files("~/R/TIMBR-master/R/", pattern = "\\.[RrSsQq]$")) {
	source(file.path("~/R/TIMBR-master/R/", nm))
}

diplotype_list <- fread(file=diplotype_file,header=TRUE,sep="\t",fill = TRUE)
diplotype_list <-diplotype_list[pos >= (MB[1]*10^6) & pos < (MB[length(MB)]+1)*10^6]

FD <- fread(file=paste0(tmp.file.path,"founder_tmp.txt"),header=TRUE,sep='\t')
FD <- FD[pos >= (MB[1]*10^6) & pos < (MB[length(MB)]+1)*10^6]

sdp_ref <- read.table(file="~/SNPmerge/Rscripts/src/sdp_table_crp_tmp.txt",header=TRUE,sep="\t",stringsAsFactors=FALSE)
sdp_ref$persudo <- sdp_ref$Freq + 1

TBfounder <- function(founder,fid){
        TB <- sdp_table[strain==founder];
        # replace None to NA_ABCDEFGH
        TB$allele_1[which(TB$allele_1=='None')] <- paste0('NA_',fid)
        TB$allele_2[which(TB$allele_2=='None')] <- paste0('NA_',fid)
        # get heterzypos
        TB$allele='None/None'
        for (thisi in 1:dim(TB)[1]){
                TB$allele[thisi] = paste(sort(c(unlist(TB[thisi,2:3]))),collapse='/')
        }

        colnames(TB)[4:5] <- c(paste0(fid,'p'),fid)
        as.data.frame(TB[,c(5,4)])
}

vec2sdp<- function(sdpvec){
        sdpvec<- c(unlist(sdpvec))
        names(sdpvec) <- c(0:7)
        for (thisi in 8:1){
                sdpvec[which(sdpvec==sdpvec[thisi])] <- names(sdpvec)[thisi]
        }
        names(sdpvec) <- c('A','B','C','D','E','F','G','H')
        sdpvec
}

poslist <- unique(diplotype_list$pos)
counti <- 0
outdf <-data.frame()
for( mypos in poslist ){
	message(mypos)
	counti <- counti + 1
	diplotype_list_sub<-diplotype_list[which(pos==mypos)]

	sdp_table_full <- FD[pos==mypos]
	# variant id check
	vtb <- table(sdp_table_full$variant_id)
	if (length(vtb)>=2){next} # remove duplicated information points
	sdp_table <- unique(sdp_table_full[variant_id==names(which(vtb==max(vtb)))[1],3:6])

	if (all(c(sdp_table$allele_1!='None',sdp_table$allele_2!='None'))){CompleteN = TRUE
	}else{CompleteN = FALSE }
	
	if (CompleteN) { next }

	sdp_sum <- list(TBfounder('A_J','A'),TBfounder('C57BL6J','B'),
			TBfounder('129S1_SvImJ','C'),TBfounder('NOD_ShiLtJ','D'),
			TBfounder('NZO_HlLtJ','E'), TBfounder('CAST_EiJ','F'),
			TBfounder('PWK_PhJ','G'),TBfounder('WSB_EiJ','H'))	
	sdp_sum <- Reduce(function(...) merge(..., all=T),sdp_sum)
	sdp_sum$prob <- sdp_sum$Ap * sdp_sum$Bp * sdp_sum$Cp * sdp_sum$Dp * sdp_sum$Ep * sdp_sum$Fp * sdp_sum$Gp * sdp_sum$Hp
	sdp_sum <- sdp_sum[,c('A','B','C','D','E','F','G','H','prob')]
	sdp_sum <- sdp_sum[sdp_sum$prob>0.001,]

	mydf <- aggregate( cbind(A,B,C,D,E,F,G,H) ~ CC, diplotype_list_sub,sum)
	mydf[,c('A','B','C','D','E','F','G','H')] <- 2*mydf[,c('A','B','C','D','E','F','G','H')]/apply(as.matrix(mydf[,c('A','B','C','D','E','F','G','H')]),1,sum)
	mydf <- merge(Phenofile,mydf,by='CC')	
	
	y <- as.matrix(mydf$Pheno); 
	rownames(y) <- mydf$CC;colnames(y) <- 'ph'
	P <- as.matrix(mydf[,c('A','B','C','D','E','F','G','H')]); 
	rownames(P) <- mydf$CC;
	A <- diag(8)
	colnames(A) <- c("A.J","C57BL.6J","m129.SvlmJ","NOD.LtJ","NZO.LtJ","CAST.EiJ","PWK.PhJ.2","WSB.EiJ")
	w <- c(mydf$NUM.OBS);

	prior.D <- list(P=P,A=A,fixed.diplo=FALSE)	

	# sdp model
	thesesdp <- data.frame()
	for (i in 1:dim(sdp_sum)[1]){
		sdp_vec <- as.character(sdp_sum[i,1:8])
		withinfo <- which(!grepl('NA',sdp_vec))
		ale_list <- unique(sdp_vec[withinfo])
		thissdp <- sdp_ref[,c('A','B','C','D','E','F','G','H','sdptim','persudo')]
		# X1 and Y1 should have same strain
		for (sgt in ale_list){
			sgt_ix <- which(sdp_vec==sgt)
			if (length(sgt_ix)>1){
				for (j in 2:length(sgt_ix)){
					thissdp <- thissdp[which(thissdp[,sgt_ix[j]]==thissdp[,sgt_ix[(j-1)]]), ]
				}
			}
		}
		# X2 and Y1 should have different strain
		if (length(ale_list)>1){
			for (k in 2:length(ale_list)){
				thissdp <- thissdp[which( thissdp[,which(sdp_vec==ale_list[k])[1]] !=
						          thissdp[,which(sdp_vec==ale_list[(k-1)])[1]] ),]
			}
		}
		thissdp$prob.local <- thissdp$persudo/sum(thissdp$persudo)
		thissdp$prob <- thissdp$prob.local * sdp_sum$prob[i]
		thesesdp <- rbind(thesesdp,thissdp)
	}

	sdp_sum <- aggregate(prob ~ sdptim,thesesdp,sum)
	prior.M <- list(model.type='list',M.IDs=sdp_sum$sdptim, probs=sdp_sum$prob,hash.names=TRUE)

	results.informed <- TIMBR(y, prior.D, prior.M,W=w,prior.v.b = 4 )
	bf1 <- results.informed$ln.ml - results.informed$ln.ml.null
	sdpp <- rev(sort(table(results.informed$post.M)))[1]/results.informed$samples

	#results.informed <- TIMBR(y, prior.D, prior.M.8allele,W=w,prior.v.b = 4 )
        #bf2 <- results.informed$ln.ml - results.informed$ln.ml.null
	bf2 <- 0

	rm(P);rm(prior.M)
	
	thisdf <- data.frame(pos=mypos,CompleteN=CompleteN,bf=bf1, bf8=bf2,sdp=names(sdpp[1]),prob=as.numeric(sdpp[1]))
	outdf <- rbind(outdf,cbind(thisdf))

	if (!counti %% 500) {
		write.table(outdf,file=paste0(tmp.file.path,"TIMBR_res_crp_",MB,".txt"),col.names=TRUE,row.names=FALSE,quote=FALSE)
	}
}
	write.table(outdf,file=paste0(tmp.file.path,"TIMBR_res_crp_",MB,".txt"),col.names=TRUE,row.names=FALSE,quote=FALSE)
outdf
}

