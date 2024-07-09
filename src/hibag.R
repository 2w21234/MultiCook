
args <- commandArgs(trailingOnly = TRUE)

inputname <- args[1]
refname <- args[2]
outputname <-args[3]
gene <- args[4:length(args)]

print(inputname)
print(refname)
print(outputname)
print(gene)

zz <- file(paste0(outputname,"_temp.log"),open="a")
sink(zz,append=TRUE, type="message")

library(HIBAG)

model.list <- get(load(refname))
print(2)
print(inputname)
data <- hlaBED2Geno(bed.fn=paste0(inputname,".bed"),fam.fn=paste0(inputname,".fam"),bim.fn=paste0(inputname,".bim"))

print(1)
print(gene)
for(i in gene){
	print(i)
	hla.id <- i
	if(inherits(model.list[[hla.id]],"hlaAttrBagObj")){
		model <- hlaModelFromObj(model.list[[hla.id]])
		pred.guess <- predict(model,data,type="response+prob")
	        print('running')
		write.table(pred.guess$postprob, paste0(outputname,"_",i,".bagout"), quote=F, col.names=T, row.names=T)
	}
}

sink()
