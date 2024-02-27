library(dplyr)
library(stringr)
library(vegan)

args=commandArgs(T)

subsample = as.numeric(args[1])
dir = args[2]
outputfile =  args[3]


# #setwd("D:/Data/snakemake-tutorial/data/test")
# subsample = 2000
# inputfile = "target.upstream_downstream.bases.singal_length.txt"
# outputfile =  "target.upstream_downstream.bases.singal_length.comres.test.txt"

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

#define a function: to sort a matrix via order function in increasing order and return the new index
sortidx <- function(id,m){
    rawidx = which(id == TRUE)
    tmp = as.data.frame(m[rawidx,])
    newidx = rawidx[do.call(order, tmp)]
    return(newidx)
}


#=======load file=======
files = list.files(dir,pattern="target.pos\\d+.upstream_downstream.\\d+bases\\+middle.\\d+bases")
df1 = data.frame()
for (i in 1:length(files)){
   file = files[i]
   file1 = paste0(dir,"/",file)
   if(i == 1){
      print(i)
      tmp1 = read.table(file1,sep="\t",header=TRUE)
      tmp1$Pos = paste0("Pos",tmp1$Position)
      df1 = tmp1
   }
   else {
      tmp1 = read.table(file1,sep="\t",header=TRUE)
      tmp1$Pos = paste0("Pos",tmp1$Position)
      df1 = rbind(df1,tmp1)
   }
}

variants = c("target")

FP_Pos = list()
FP_Pos[[1]] = unique(df1$Position)
names(FP_Pos) = variants

variantdfs = list()
variantdfs[[1]] = df1
names(variantdfs) = variants

#========================do some comparison work of FP positions========================
FP_readid_prefix = c()
FP_Positions = c()
FP_Group1 = c()
FP_Group2 = c()
FP_Group1_Num = c()
FP_Group2_Num = c()
FP_Group1_Mean = c()
FP_Group2_Mean = c()
FP_Pvalues = c()

for (j in 1:length(variants)){
    FP_lists = FP_Pos[[variants[j]]]
    FP_comparisons1 = c("plus_match","minus_match")
    FP_comparisons2 =c("plus_del","minus_del")
    for (i in 1:length(FP_lists)){
        for (k in 1:length(FP_comparisons1)){
            a <- variantdfs[[variants[j]]] %>%
                 subset(Pos == paste0("Pos",FP_lists[i]) & Group == FP_comparisons1[k])
            b <- variantdfs[[variants[j]]] %>%
                 subset(Pos == paste0("Pos",FP_lists[i]) & Group == FP_comparisons2[k])
            t <- wilcox.test(a$Group_Signal_Length,b$Group_Signal_Length)
            FP_readid_prefix = c(FP_readid_prefix,unique(a$readid_prefix))
            FP_Positions = c(FP_Positions,FP_lists[i])
            FP_Group1 = c(FP_Group1,FP_comparisons1[k])
            FP_Group2 = c(FP_Group2,FP_comparisons2[k])
            FP_Group1_Num = c(FP_Group1_Num,length(a$Group_Signal_Length))
            FP_Group2_Num = c(FP_Group2_Num,length(b$Group_Signal_Length))
            FP_Group1_Mean = c(FP_Group1_Mean,mean(a$Group_Signal_Length))
            FP_Group2_Mean = c(FP_Group2_Mean,mean(b$Group_Signal_Length))
            FP_Pvalues = c(FP_Pvalues,t[[3]])
        }
    }
}

Guppy_FP_df <- data.frame(readid_prefix=FP_readid_prefix,Pos=FP_Positions,Group1=FP_Group1,
                   Group2=FP_Group2,Group1_Num=FP_Group1_Num,Group2_Num=FP_Group2_Num,Group1_Mean=FP_Group1_Mean,Group2_Mean=FP_Group2_Mean,Pvalues=FP_Pvalues)

#========================do some cluster work of FP positions via kmeans,pca========================
Guppy_FP_df$MRPP2P = rep("NaN",length(Guppy_FP_df$Group1))
Guppy_FP_df$MRPP2A = rep("NaN",length(Guppy_FP_df$Group1))

for (k in 1:length(variants)){
    FP_lists = FP_Pos[[k]]
    FP_comparisons1 = c("plus_match","minus_match")
    FP_comparisons2 =c("plus_del","minus_del")
    for (j in 1:length(FP_lists)){
        target = paste0("Pos",FP_lists[j])
        print(target)
        a2 <- variantdfs[[variants[k]]] %>%
              subset(Pos == target & (Group == "minus_del" | Group == "plus_del"))
        b2 <- variantdfs[[variants[k]]] %>%
              subset(Pos == target & (Group == "minus_match" | Group == "plus_match"))
        #len2 <- min(min(a2$Group_Signal_Length),min(b2$Group_Signal_Length))
        len2 <- 60
        c2 <-rbind(a2,b2)
        #construct a matrix
        m2 <- matrix(0,length(c2$Group_Signal_Length),len2)
        for(i in 1:length(c2$Group_Signal_Length)){
            tmp = c2[i,"Group_Signals"]
            tmp1 = as.numeric(unlist(strsplit(tmp,split=",")))
            #divide tmp1 into len1 bins, referred to Ni's script
            m2[i,] = intoBins(tmp1[tmp1>0],len2)
            #or randomly subsample
            #m2[i,] = sample(tmp1,len2)
        }
        # #update 
        # m2=t(scale(t(m2),center = F,scale = T))

        #===visualize by heatmap===
        condition <- factor(c(rep("del",dim(a2)[1]),rep("match",dim(b2)[1])), levels = c("del","match"))
        rowData2 <- data.frame(row.names=rownames(m2), condition, subgroup = c(a2$Group,b2$Group), strain = c(a2$readid_prefix,b2$readid_prefix))
        
        #define the color mapping
        # col_fun = colorRamp2(c(-2, 0, 2), c("green", "white", "red"))
        #col<- colorRamp2(c(0, 200, 400, 600), c("white","#0d4be5","yellow","#e90304"))
        plusdel_idx2 <- c2$Group=="plus_del"
        plusmatch_idx2 <- c2$Group=="plus_match"
        minusdel_idx2 <- c2$Group=="minus_del"
        minusmatch_idx2 <- c2$Group=="minus_match"

        #only subsample at most 500 reads
        if(length(which(plusdel_idx2 == TRUE)) > subsample){
            tidx = which(plusdel_idx2 == TRUE)
            tidx = sample(tidx,subsample)
            plusdel_idx2[which(plusdel_idx2 == TRUE)] = FALSE 
            plusdel_idx2[tidx] = TRUE
        }

        if(length(which(plusmatch_idx2 == TRUE)) > subsample){
            tidx = which(plusmatch_idx2 == TRUE)
            tidx = sample(tidx,subsample)
            plusmatch_idx2[which(plusmatch_idx2 == TRUE)] = FALSE 
            plusmatch_idx2[tidx] = TRUE
        }

        if(length(which(minusdel_idx2 == TRUE)) > subsample){
            tidx = which(minusdel_idx2 == TRUE)
            tidx = sample(tidx,subsample)
            minusdel_idx2[which(minusdel_idx2 == TRUE)] = FALSE 
            minusdel_idx2[tidx] = TRUE
        }

        if(length(which(minusmatch_idx2 == TRUE)) > subsample){
            tidx = which(minusmatch_idx2 == TRUE)
            tidx = sample(tidx,subsample)
            minusmatch_idx2[which(minusmatch_idx2 == TRUE)] = FALSE 
            minusmatch_idx2[tidx] = TRUE
        }

        newplusdel_idx2 <- sortidx(plusdel_idx2,m2)
        newplusmatch_idx2 <- sortidx(plusmatch_idx2,m2)
        newminusdel_idx2 <- sortidx(minusdel_idx2,m2)
        newminusmatch_idx2 <- sortidx(minusmatch_idx2,m2)

        #first heatmap for plus strand data
        idx <- c(newplusdel_idx2,newplusmatch_idx2)
        #idx <- plusdel_idx2 | plusmatch_idx2

        #MRPP
        print("Plus MRPP")
        tmpmrpp <- mrpp(m2[idx,],rowData2$subgroup[idx],permutations = 999, distance = "bray", parallel = 8)

        Guppy_FP_df[Guppy_FP_df$readid_prefix == unique(variantdfs[[variants[k]]]$readid_prefix) & Guppy_FP_df$Pos == FP_lists[j] & Guppy_FP_df$Group2 == "plus_del", 
        c("MRPP2P","MRPP2A")] <- c(tmpmrpp$Pvalue,tmpmrpp$A)

        #then heatmap for minus strand data
        idx <- c(newminusdel_idx2,newminusmatch_idx2)
        #idx <- minusdel_idx2 | minusmatch_idx2

        #MRPP
        print("Minus MRPP")
        tmpmrpp <- mrpp(m2[idx,],rowData2$subgroup[idx],permutations = 999, distance = "bray", parallel = 8)

        Guppy_FP_df[Guppy_FP_df$readid_prefix == unique(variantdfs[[variants[k]]]$readid_prefix) & Guppy_FP_df$Pos == FP_lists[j] & Guppy_FP_df$Group2 == "minus_del", 
        c("MRPP2P","MRPP2A")] <- c(tmpmrpp$Pvalue,tmpmrpp$A)

    }
}

write.table(Guppy_FP_df,file=outputfile,sep="\t",quote=FALSE,row.names=FALSE)

