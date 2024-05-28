args <- commandArgs(trailingOnly = TRUE)

inputname <- args[1]
refname <- args[2]
out <- args[3]
gene <- args[4]

zz <- file(paste0(out,"_temp.log"),open="a")
sink(zz,append=TRUE, type="message")

library(HIBAG)

model.list <- get(load(refname))
data <- hlaBED2Geno(bed.fn=paste0(inputname,".bed"),fam.fn=paste0(inputname,".fam"),bim.fn=paste0(inputname,".bim"),assembly="hg18")

	
model <- hlaModelFromObj(model.list)
pred.guess <- predict(model,data,type="response+prob")
	
write.table(pred.guess$postprob, paste0(out,"_",gene,".bagout"), quote=F, col.names=T, row.names=T)


