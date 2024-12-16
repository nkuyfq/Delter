library(dplyr)

args=commandArgs(T)

res_dir = args[1]
homo_Q_thres = as.numeric(args[2])  #threshold for homopolymeric deletions
other_Q_thres = as.numeric(args[3])
depthres = 20

#res_dir="C:/Users/Lenovo/Desktop/Snakemake_Res/veillonella_rogosae"
# homo_Q_thres = 23
# other_Q_thres = 20.6

#setwd(res_dir)
VariantRes <- data.frame()
VariantInfo <- data.frame()
#=======================load result files====================
file <- paste0(res_dir,"/fq.Qscore.info.txt")
VariantRes <- read.table(file,sep="\t",header=TRUE)
file1 <- paste0(res_dir,"/variant.info.txt")
tmp1<- read.table(file1,sep="\t",header=TRUE)
tmp1$Pos <- paste0("pos",tmp1$POS+1)
VariantInfo <- tmp1

#=============add homo and non-homo to VariantRes============
VariantRes$Type <- NA
VariantRes$AF <- NA
for(i in 1:dim(VariantRes)[1]){
	#print(i)
	tmp <- VariantInfo %>% subset(Pos == VariantRes[i,"Pos"])
	if(dim(tmp)[1]>=1){
		VariantRes[i,"Type"] <- tmp[1,"Location"]
		VariantRes[i,"AF"] <- tmp[1,"AF"]
	}

}

#=================filter by sequencing depth=================
newidx <- c()
for(k in 1:(dim(VariantRes)[1]/2)){
	if(VariantRes[2*k-1,"Group1_Num"] >= depthres & VariantRes[2*k-1,"Group2_Num"] >= depthres & VariantRes[2*k,"Group1_Num"] >= depthres & VariantRes[2*k,"Group2_Num"] >= depthres){
		newidx <- c(newidx,c(2*k-1,2*k))
	}
}

if(length(newidx) >= 1){
	#=================identify artificial deletions=================
	#positions with flags marked as "FP" are potential artificial deletions
	tmpres <- VariantRes[newidx,]
	tmpres$flag <- "NA"
	for(k in 1:(dim(tmpres)[1]/2)){
		curthres <- NA
		if(tmpres[2*k-1,"Type"] == "Homo"){
			curthes <- homo_Q_thres
		}
		else {
			if(tmpres[2*k-1,"Type"] == "Non-Homo"){
				curthes <- other_Q_thres
			}
		}

		if(tmpres[2*k-1,"Group1_Mean"]<=curthes & tmpres[2*k,"Group1_Mean"]<=curthes){
			tmpres[2*k-1,"flag"] <- "FP"
			tmpres[2*k,"flag"] <- "FP"
		}
	}

	if(dim(tmpres)[1] >= 1){
		write.table(tmpres,file=paste0(res_dir,"/Qscore.filtered.txt"),sep="\t",quote=FALSE,row.names=FALSE)
	}
} else {
	print("No variations are retained after filtering!")
}

# #===============display the number of removed FP deletions=============
# #filtered FP info
# table(tmpres[tmpres$flag == "FP",]$Type)/2

# #===============display MuAF of removed FP homo-dels===================
# #filtered FP info
# tmp <- tmpres %>% subset(flag == "FP" & Type == "Homo") %>% select(AF)
# c(max(tmp$AF),min(tmp$AF))
