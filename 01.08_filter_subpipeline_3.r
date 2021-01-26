#R code
#select expressed TEs from TEFingerprint data
dangler=read.table("AllRep_comp_intersect_V4.txt", header=F, sep="\t")
sense=read.table("BedCov_OverlapBP_sense_count.txt", header=T, sep="\t")

#ctrl
dangler2=dangler[,c(146,54:56)]
dangler2$mean=rowMeans(dangler2[,2:4])
dangler2=subset(dangler2, dangler2$mean>10)

colnames(dangler2)=c("id", "ctrl_a", "ctrl_b", "ctrl_c", "mean")
df=merge(dangler2, sense, by.x="id", by.y="TEm")
df$countMean=rowMeans(df[,11:13])
df2=subset(df, df$countMean>10)

df2=df2[,1:5]
dftest=df2[!duplicated(df2$id),]
write.table(df2, "AllRep_comp_intersect_exprCtrl_V4.txt", col.names=T, row.names=F, sep="\t", quote=F)

#mock
dangler2=dangler[,c(146,57:68)]
dangler2$meanM01=rowMeans(dangler2[,2:4])
dangler2$meanM03=rowMeans(dangler2[,5:7])
dangler2$meanM06=rowMeans(dangler2[,8:10])
dangler2$meanM12=rowMeans(dangler2[,11:13])
min_read = 10
dangler2 <- dangler2[apply(dangler2[,14:17],1,function(x){max(x)}) > min_read,]

colnames(dangler2)=c("id", "M01_a", "M01_b", "M01_c","M03_a", "M03_b", "M03_c","M06_a",
"M06_b", "M06_c","M12_a", "M12_b", "M12_c", "meanM01", "meanM03", "meanM06",
"meanM12")
df=merge(dangler2, sense, by.x="id", by.y="TEm")
df$countMean01=rowMeans(df[,26:28])
df$countMean03=rowMeans(df[,29:31])
df$countMean06=rowMeans(df[,32:34])
df$countMean12=rowMeans(df[,35:37])
min_read = 10
df2 <- df[apply(df[,62:65],1,function(x){max(x)}) > min_read,]

df2=df2[,1:17]
dftest=df2[!duplicated(df2$id),]
write.table(df2, "AllRep_comp_intersect_exprMock_V4.txt", col.names=T, row.names=F, sep="\t", quote=F)

#yeast
dangler2=dangler[,c(146,69:80)]
dangler2$meanY01=rowMeans(dangler2[,2:4])
dangler2$meanY03=rowMeans(dangler2[,5:7])
dangler2$meanY06=rowMeans(dangler2[,8:10])
dangler2$meanY12=rowMeans(dangler2[,11:13])
min_read = 10
dangler2 <- dangler2[apply(dangler2[,14:17],1,function(x){max(x)}) > min_read,]

colnames(dangler2)=c("id", "Y01_a", "Y01_b", "Y01_c","Y03_a", "Y03_b", "Y03_c","Y06_a",
"Y06_b", "Y06_c","Y12_a", "Y12_b", "Y12_c", "meanY01", "meanY03", "meanY06",
"meanY12")
df=merge(dangler2, sense, by.x="id", by.y="TEm")
df$countMean01=rowMeans(df[,38:40])
df$countMean03=rowMeans(df[,41:43])
df$countMean06=rowMeans(df[,44:46])
df$countMean12=rowMeans(df[,47:49])
min_read = 10
df2 <- df[apply(df[,62:65],1,function(x){max(x)}) > min_read,]

df2=df2[,1:17]
dftest=df2[!duplicated(df2$id),]
write.table(df2, "AllRep_comp_intersect_exprYeast_V4.txt", col.names=T, row.names=F, sep="\t", quote=F)
