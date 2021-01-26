#R code
#Before running this script, combine the output file of individual sequencing library column by column.
#Name combined file as Combined_k100Mm_AllSenseTEm_counttable.txt
#The order of the column start with control libraries, and then the experimental conditions
#Here we have one control and two expermintal conditions. Each conditions has four timepoints, eich with three replicates

#select expr. TEs from htseq-count data
#htseq-count reports “mapped read pairs” when analysing pair-end sequencing data.
#the threshold of the pipeline is set to 10 reads, i.e. at least more than 5 read pairs
RE = read.table('Combined_k100Mm_AllSenseTEm_counttable.txt', header=T, row.names=1)
RE=RE[!grepl("__", row.names(RE)),]

RE2=RE[,1:3] #ctrl samples (T=0) are from column 1 to column 3
RE2$mean=rowMeans(RE2[,1:3])
RE2=subset(RE2, RE2$mean>5)
RE2$id=rownames(RE2)
RE2=RE2[,c(5,1:4)]
RE2$id=gsub("_ClassI;", "_", RE2$id)
RE2$id=gsub("_ClassII;", "_", RE2$id)
RE2$id=gsub(":", "_", RE2$id)
write.table(RE2, "k100Mm_AllSenseTEm_curatedGTF_exprCtrl05.txt", col.names=T, row.names=F, sep="\t", quote=F)

RE2=RE[,4:15] #mock samples are from column 4 to column 15
RE2$meanM01=rowMeans(RE2[,1:3])
RE2$meanM03=rowMeans(RE2[,4:6])
RE2$meanM06=rowMeans(RE2[,7:9])
RE2$meanM12=rowMeans(RE2[,10:12])
min_read = 5
RE2 <- RE2[apply(RE2[,13:16],1,function(x){max(x)}) > min_read,]
RE2$id=rownames(RE2)
RE2=RE2[,c(17,1:16)]
RE2$id=gsub("_ClassI;", "_", RE2$id)
RE2$id=gsub("_ClassII;", "_", RE2$id)
RE2$id=gsub(":", "_", RE2$id)
write.table(RE2, "k100Mm_AllSenseTEm_curatedGTF_exprMock05.txt", col.names=T, row.names=F, sep="\t", quote=F)

RE2=RE[,16:27] #yeast samples are from column 16 to column 27
RE2$meanY01=rowMeans(RE2[,1:3])
RE2$meanY03=rowMeans(RE2[,4:6])
RE2$meanY06=rowMeans(RE2[,7:9])
RE2$meanY12=rowMeans(RE2[,10:12])
min_read = 5
RE2 <- RE2[apply(RE2[,13:16],1,function(x){max(x)}) > min_read,]
RE2$id=rownames(RE2)
RE2=RE2[,c(17,1:16)]
RE2$id=gsub("_ClassI;", "_", RE2$id)
RE2$id=gsub("_ClassII;", "_", RE2$id)
RE2$id=gsub(":", "_", RE2$id)
write.table(RE2, "k100Mm_AllSenseTEm_curatedGTF_exprYeast05.txt", col.names=T, row.names=F, sep="\t", quote=F)
