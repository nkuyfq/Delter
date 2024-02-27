library(dplyr)
library(ShortRead) #for fastq manipulation
library(stringr)
library(vegan)

args=commandArgs(T)

curpos = as.numeric(args[1])
dir = args[2]
outputfile =  args[3]


#setwd('D:/Data/20200304-nCov/basecalling/zilin_server/20221006_Twist/figures')


#define a function: to divide an int array into N bins
intoBins <- function(array,Binnum){
    sum = matrix(0,Binnum,1)
    count = matrix(0,Binnum,1)
    binned = matrix(0,Binnum,1)
    scal = length(array)/Binnum
    for (i in 1:length(array)){
        bin = as.integer(i/scal)
        sum[bin] = sum[bin]+array[i]
        count[bin] = count[bin]+1
    }
    for (bin in 1:Binnum){
        if(count[bin] != 0){
            binned[bin] = as.integer(sum[bin]/count[bin])
        }
    }
    return(binned)
}

#====================load Variant fastq==================

curvariant = c("target")
group1 = c("plus.del","minus.del")
group2 = c("plus.match","minus.match")

vs = c()
positions = c()
g1 = c()
g2 = c()
group1mean = c()
group2mean = c()
pvalues = c()

for(k in 1:length(group1)){
	curgroup1 = group1[k]
	curgroup2 = group2[k]
	print(curgroup1)
	DelSeqs=c()
	#load variant fq
	file = paste0("./",dir,"/",curvariant,".",curgroup1,".fq")
	rfq <- readFastq(file)
	#rfq
	DelSeqs = c(DelSeqs,rfq)
	names(DelSeqs) = curvariant
	#load control fq
	MatchSeqs=c()
	file = paste0("./",dir,"/",curvariant,".",curgroup2,".fq")
	rfq <- readFastq(file)
	#rfq
	MatchSeqs = c(MatchSeqs,rfq)
	names(MatchSeqs) = curvariant
	Delinfo <- data.frame(Variant = rep(curvariant,length(DelSeqs[[curvariant]]@quality)),
					 Pos = rep(curpos,length(DelSeqs[[curvariant]]@quality)),
					 Group = rep(curgroup1,length(DelSeqs[[curvariant]]@quality)))
	#then load variant quality score
	Delinfo$Quality <- rep("NaN",length(DelSeqs[[curvariant]]@quality))
	for (m in 1:length(DelSeqs[[curvariant]]@quality)){
		quality <- toString(DelSeqs[[curvariant]]@quality[[m]])
		Delinfo[m,"Quality"] <- quality
	}
	Matchinfo <- data.frame(Variant = rep(curvariant,length(MatchSeqs[[curvariant]]@quality)),
					 Pos = rep(curpos,length(MatchSeqs[[curvariant]]@quality)),
					 Group = rep(curgroup2,length(MatchSeqs[[curvariant]]@quality)))
	#then load control quality score
	Matchinfo$Quality <- rep("NaN",length(MatchSeqs[[curvariant]]@quality))
	for (m in 1:length(MatchSeqs[[curvariant]]@quality)){
		quality <- toString(MatchSeqs[[curvariant]]@quality[[m]])
		Matchinfo[m,"Quality"] <- quality
	}
	#now generate a new dataframe to combine the above info
	combinedinfo <- data.frame()
	variantQscore <- c()
	for (m in 1:length(DelSeqs[[curvariant]]@quality)){
		curquality <- Delinfo[m,"Quality"]
		qscore <- as.numeric(charToRaw(curquality))-33
		variantQscore <- c(variantQscore,qscore)
	}
	combinedinfo <- data.frame(Variant=rep(curvariant,length(variantQscore)),Pos=rep(curpos,length(variantQscore)),
							   Group=rep(curgroup1,length(variantQscore)),Qscore=variantQscore)
	controlQscore <- c()
	for (m in 1:length(MatchSeqs[[curvariant]]@quality)){
		curquality <- Matchinfo[m,"Quality"]
		qscore <- as.numeric(charToRaw(curquality))-33
		controlQscore <- c(controlQscore,qscore)
	}
	tmpinfo <- data.frame(Variant=rep(curvariant,length(controlQscore)),Pos=rep(curpos,length(controlQscore)),
						Group=rep(curgroup2,length(controlQscore)),Qscore=controlQscore)
	combinedinfo <- rbind(combinedinfo,tmpinfo)
	t <- wilcox.test(combinedinfo[combinedinfo$Group == curgroup1,]$Qscore,combinedinfo[combinedinfo$Group == curgroup2,]$Qscore)
	vs = c(vs,curvariant)
	positions = c(positions,curpos)
	g1 = c(g1,curgroup1)
	g2 = c(g2,curgroup2)
	group1mean = c(group1mean,mean(combinedinfo[combinedinfo$Group == curgroup1,]$Qscore))
	group2mean = c(group2mean,mean(combinedinfo[combinedinfo$Group == curgroup2,]$Qscore))
	pvalues = c(pvalues,t[[3]])
	#t[[3]]
}


finalinfo <- data.frame(Variant=vs,Pos=positions,Group1=g1,Group2=g2,Group1_Mean=group1mean,Group2_Mean=group2mean,Pvalues=pvalues)

write.table(finalinfo,file=outputfile,sep="\t",quote=FALSE,row.names=FALSE)
