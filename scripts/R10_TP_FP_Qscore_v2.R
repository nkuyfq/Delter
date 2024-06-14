library(dplyr)
library(ShortRead) #for fastq manipulation
library(stringr)
library(vegan)

args=commandArgs(T)
#curpos = as.numeric(args[1])
dir = args[1]
outputfile =  args[2]


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

#variants=c("Alpha","Beta","Delta","Epsilon","Gamma","Kappa","Iota","BA.1","BA.2")
variants = c("target")

#load Variant and corresponding Postions
VariantPos = data.frame()
vs = c()
pos = c()

for (i in 1:length(variants)){
	files = list.files(dir,pattern=paste0(variants[i],".pos\\d+.plus.del.fq"))
	for(j in 1:length(files)){
		print(files[j])
		tmp = str_extract(files[j],"pos\\d+")
		vs = c(vs,variants[i])
		pos = c(pos,tmp)
	}
}

VariantPos = data.frame(Variant=vs,Pos=pos)

#====================load Variant fastq==================
group1 = c("plus.del","minus.del")
group2 = c("plus.match","minus.match")

vs = c()
positions = c()
g1 = c()
g2 = c()
group1mean = c()
group2mean = c()
pvalues = c()
group1num = c()
group2num = c()


allQscoreinfo = data.frame()

for (i in 1:length(variants)){
	curvariant = variants[i]
	curposinfo <- VariantPos %>%
	              subset(Variant == curvariant)
	for(j in 1:length(curposinfo$Pos)){
		curpos = curposinfo[j,"Pos"]
		for(k in 1:length(group1)){
			curgroup1 = group1[k]
			curgroup2 = group2[k]
			DelSeqs=c()
			#load variant fq
			file = paste0(dir,"/",curvariant,".",curpos,".",curgroup1,".fq")
			rfq <- readFastq(file)
			#rfq
			DelSeqs = c(DelSeqs,rfq)
			group1num = c(group1num,length(rfq))
			names(DelSeqs) = curvariant
			#load control fq
			MatchSeqs=c()
			file = paste0(dir,"/",curvariant,".",curpos,".",curgroup2,".fq")
			rfq <- readFastq(file)
			#rfq
			MatchSeqs = c(MatchSeqs,rfq)
			group2num = c(group2num,length(rfq))
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
			#allQscoreinfo <- rbind(allQscoreinfo,combinedinfo)
			t <- wilcox.test(combinedinfo[combinedinfo$Group == curgroup1,]$Qscore,combinedinfo[combinedinfo$Group == curgroup2,]$Qscore)
			vs = c(vs,curvariant)
			positions = c(positions,curpos)
			g1 = c(g1,curgroup1)
			g2 = c(g2,curgroup2)
			group1mean = c(group1mean,mean(combinedinfo[combinedinfo$Group == curgroup1,]$Qscore))
			group2mean = c(group2mean,mean(combinedinfo[combinedinfo$Group == curgroup2,]$Qscore))
			pvalues = c(pvalues,t[[3]])
			
		}
	}
}

finalinfo <- data.frame(Variant=vs,Pos=positions,Group1=g1,Group2=g2,Group1_Num=group1num,Group2_Num=group2num,Group1_Mean=group1mean,Group2_Mean=group2mean,Pvalues=pvalues)

write.table(finalinfo,file=outputfile,sep="\t",quote=FALSE,row.names=FALSE)


