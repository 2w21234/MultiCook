library(data.table)
library(dplyr)
options(warn=-1)

args = commandArgs(trailingOnly=TRUE)

path= args[1] ######################### 이거 input으로 받아서 코드 수정
#print(path)
output='HIBAG_OUT.vcfh'

setwd(path)

files=list.files('./','*bagout')
genes=c()
#print(files)
for(i in 1:length(files)){
    file=files[i]
    genes[i]=strsplit(strsplit(file,'_')[[1]][3],'.bagout')[[1]]
}


total_probs=data.frame()


for(f in 1:length(files)){
    #print(files[f])
    gene=genes[f]
    a=fread(files[f])
    allele_pairs<-a[[1]]
    #rownames(a) <- a[[1]]
    a=a[,2:ncol(a)]
    rownames(a) <- allele_pairs
    
    k=1
    alleles=c()
    for(i in 1:length(allele_pairs)){
        a1=strsplit(allele_pairs[i], '/')[[1]][1]
        a2=strsplit(allele_pairs[i], '/')[[1]][2]
        #i=1
        if(!(a1 %in% alleles)){
            alleles[k]=a1
            k=k+1
        }
        if(!(a2 %in% alleles)){
            alleles[k]=a2
            k=k+1
        }
        
    }    

    new_alleles=c()
    for(l in 1:length(alleles)){
        tmp=paste0(strsplit(alleles[l],':')[[1]],collapse='')
        new_allele=paste0(c('HLA_',gene,'_',tmp,'_exon'), collapse='')
        new_alleles[l] <- new_allele
    }

    
    probs=data.frame(matrix(0, nrow=length(new_alleles), ncol=ncol(a)))
    rownames(probs)=new_alleles
    colnames(probs)=colnames(a)
    probs = as.data.frame(probs)
    
    for(j in 1:length(alleles)){
        #print(j)
        allele=alleles[j]
        k=1
        for(i in 1:length(allele_pairs)){
            ap=allele_pairs[i]
            ap=strsplit(ap, '/')[[1]]        
            if(allele %in% ap){
                if(k==1){
                    prob_sum <- a[i,]
                    k=k+1
                }else{
                    prob_sum <- prob_sum+a[i,]
                    k=k+1
                }
                probs[j,] <- prob_sum
            }   
        }
    }
    for(k in 1:ncol(probs)){
        x=probs[k]
        tmp=top_n(x,2)
        if(tmp[1,] < tmp[2,]*2){
            x <- x/2
            probs[k] <- x
        }
    }
    total_probs <- rbind(total_probs, probs)
    
}
total_probs_sub <- total_probs[!grepl('NA', rownames(total_probs)),]
total_probs_sub <-total_probs_sub[order(rownames(total_probs_sub)),]


res = data.frame(matrix(0,nrow=nrow(total_probs_sub), ncol=ncol(total_probs_sub)+9))
res[10:ncol(res)] <- total_probs_sub
res[,3] <- row.names(total_probs_sub)

colnames(res)[10:ncol(res)] <- colnames(total_probs_sub)
colnames(res)[0:9] <- 0

fwrite(res, output, sep='\t', row.names=F)

