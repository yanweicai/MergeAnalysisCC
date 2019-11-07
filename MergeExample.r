source("diplotype_mapping_perMB.r")
source("diplodiploGet.r")
source("gene_res_combine_filter.r")
source("sigSNP_info.r")
library(data.table)
library(lme4)

genodb="/nas/depts/006/valdar-lab/PUBLIC/2018_1_isvdb/full_1504/"
# dont have CC078-CC082

# Define the region and path for tmp files and results
CHR <- '15'
peakI <- c(39,41) # 39,52
tmp.file.path='./output/'

# Data preparation (Strain mean), with toy dataset
# might need to do more about the dataset
mydf <- read.csv(file="at_data_modified.csv",stringsAsFactors=FALSE)
df <- merge( aggregate(y ~ strain,mydf,length),aggregate(y ~ strain,mydf,mean),by="strain")
colnames(df)<-c("CC","NUM.OBS","Pheno") # NUM.OBS number of observations serve as weight for regression
df<- df[1:51,] # remove CC078-CC082 as not in database

df -> Phenofile

##################
###### MAIN ######
## Pull out the required Genotype and Diplotype information of ISVdb
diplo.region.isvdb(Phenofile,genodb,CHR,peakI,tmp.file.path)
diplo.founder.isvdb(CHR,peakI,tmp.file.path)
## Run the tests, in per MB to save memory
## could chage to parallel 
for (MB in peakI[1]:peakI[2]){
    diplotype_merge_mapping(Phenofile,diplotype_file=paste0(tmp.file.path,"diplotype_tmp.txt"),MB,K=NULL,tmp.file.path=tmp.file.path)
}
## Filtering the results and generate output tables
diplo_res_filter(peakI,tmp.file.path)
sig_SNP_fast(CHR,peakI,genodb,tmp.file.path)

