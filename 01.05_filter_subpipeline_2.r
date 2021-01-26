#R code
df01=read.table("01_BedCov_OverlapBasePair_senseRead.txt", header=F, sep="\t")
df02=read.table("02_BedCov_OverlapBasePair_senseRead.txt", header=F, sep="\t")
df03=read.table("03_BedCov_OverlapBasePair_senseRead.txt", header=F, sep="\t")
df04=read.table("04_BedCov_OverlapBasePair_senseRead.txt", header=F, sep="\t")
df05=read.table("05_BedCov_OverlapBasePair_senseRead.txt", header=F, sep="\t")
df06=read.table("06_BedCov_OverlapBasePair_senseRead.txt", header=F, sep="\t")
df07=read.table("07_BedCov_OverlapBasePair_senseRead.txt", header=F, sep="\t")
df08=read.table("08_BedCov_OverlapBasePair_senseRead.txt", header=F, sep="\t")
df09=read.table("09_BedCov_OverlapBasePair_senseRead.txt", header=F, sep="\t")
df10=read.table("10_BedCov_OverlapBasePair_senseRead.txt", header=F, sep="\t")
df11=read.table("11_BedCov_OverlapBasePair_senseRead.txt", header=F, sep="\t")
df12=read.table("12_BedCov_OverlapBasePair_senseRead.txt", header=F, sep="\t")
df13=read.table("13_BedCov_OverlapBasePair_senseRead.txt", header=F, sep="\t")
df14=read.table("14_BedCov_OverlapBasePair_senseRead.txt", header=F, sep="\t")
df15=read.table("15_BedCov_OverlapBasePair_senseRead.txt", header=F, sep="\t")
df16=read.table("16_BedCov_OverlapBasePair_senseRead.txt", header=F, sep="\t")
df17=read.table("17_BedCov_OverlapBasePair_senseRead.txt", header=F, sep="\t")
df18=read.table("18_BedCov_OverlapBasePair_senseRead.txt", header=F, sep="\t")
df19=read.table("19_BedCov_OverlapBasePair_senseRead.txt", header=F, sep="\t")
df20=read.table("20_BedCov_OverlapBasePair_senseRead.txt", header=F, sep="\t")
df21=read.table("21_BedCov_OverlapBasePair_senseRead.txt", header=F, sep="\t")
df22=read.table("22_BedCov_OverlapBasePair_senseRead.txt", header=F, sep="\t")
df23=read.table("23_BedCov_OverlapBasePair_senseRead.txt", header=F, sep="\t")
df24=read.table("24_BedCov_OverlapBasePair_senseRead.txt", header=F, sep="\t")
df25=read.table("25_BedCov_OverlapBasePair_senseRead.txt", header=F, sep="\t")
df26=read.table("26_BedCov_OverlapBasePair_senseRead.txt", header=F, sep="\t")
df27=read.table("27_BedCov_OverlapBasePair_senseRead.txt", header=F, sep="\t")
df28=read.table("28_BedCov_OverlapBasePair_senseRead.txt", header=F, sep="\t")
df29=read.table("29_BedCov_OverlapBasePair_senseRead.txt", header=F, sep="\t")
df30=read.table("30_BedCov_OverlapBasePair_senseRead.txt", header=F, sep="\t")
df31=read.table("31_BedCov_OverlapBasePair_senseRead.txt", header=F, sep="\t")
df32=read.table("32_BedCov_OverlapBasePair_senseRead.txt", header=F, sep="\t")
df33=read.table("33_BedCov_OverlapBasePair_senseRead.txt", header=F, sep="\t")
df34=read.table("34_BedCov_OverlapBasePair_senseRead.txt", header=F, sep="\t")
df35=read.table("35_BedCov_OverlapBasePair_senseRead.txt", header=F, sep="\t")
df36=read.table("36_BedCov_OverlapBasePair_senseRead.txt", header=F, sep="\t")
df37=read.table("37_BedCov_OverlapBasePair_senseRead.txt", header=F, sep="\t")
df38=read.table("38_BedCov_OverlapBasePair_senseRead.txt", header=F, sep="\t")
df39=read.table("39_BedCov_OverlapBasePair_senseRead.txt", header=F, sep="\t")
df01=df01[order(df01$V9),]
df02=df02[order(df02$V9),]
df03=df03[order(df03$V9),]
df04=df04[order(df04$V9),]
df05=df05[order(df05$V9),]
df06=df06[order(df06$V9),]
df07=df07[order(df07$V9),]
df08=df08[order(df08$V9),]
df09=df09[order(df09$V9),]
df10=df10[order(df10$V9),]
df11=df11[order(df11$V9),]
df12=df12[order(df12$V9),]
df13=df13[order(df13$V9),]
df14=df14[order(df14$V9),]
df15=df15[order(df15$V9),]
df16=df16[order(df16$V9),]
df17=df17[order(df17$V9),]
df18=df18[order(df18$V9),]
df19=df19[order(df19$V9),]
df20=df20[order(df20$V9),]
df21=df21[order(df21$V9),]
df22=df22[order(df22$V9),]
df23=df23[order(df23$V9),]
df24=df24[order(df24$V9),]
df25=df25[order(df25$V9),]
df26=df26[order(df26$V9),]
df27=df27[order(df27$V9),]
df28=df28[order(df28$V9),]
df29=df29[order(df29$V9),]
df30=df30[order(df30$V9),]
df31=df31[order(df31$V9),]
df32=df32[order(df32$V9),]
df33=df33[order(df33$V9),]
df34=df34[order(df34$V9),]
df35=df35[order(df35$V9),]
df36=df36[order(df36$V9),]
df37=df37[order(df37$V9),]
df38=df38[order(df38$V9),]
df39=df39[order(df39$V9),]
data1=cbind(df01[,c(1,4,5,9,6,7)],df01[,10],df02[,10],df03[,10],df04[,10],df05[,10], df06[,10],df07[,10],df08[,10],df09[,10],df10[,10],df11[,10],df12[,10],df13[,10], df14[,10],df15[,10],df16[,10],df17[,10],df18[,10],df19[,10],df20[,10],df21[,10], df22[,10],df23[,10],df24[,10],df25[,10],df26[,10],df27[,10],df28[,10],df29[,10], df30[,10],df31[,10],df32[,10],df33[,10],df34[,10],df35[,10],df36[,10],df37[,10], df38[,10],df39[,10])
data2=cbind(df01[,c(1,4,5,9,6,7)],df01[,13],df02[,13],df03[,13],df04[,13],df05[,13], df06[,13],df07[,13],df08[,13],df09[,13],df10[,13],df11[,13],df12[,13],df13[,13],  df14[,13],df15[,13],df16[,13],df17[,13],df18[,13],df19[,13],df20[,13],df21[,13], df22[,13],df23[,13],df24[,13],df25[,13],df26[,13],df27[,13],df28[,13],df29[,13], df30[,13],df31[,13],df32[,13],df33[,13],df34[,13],df35[,13],df36[,13],df37[,13], df38[,13],df39[,13])
data3=cbind(df01[,c(1,4,5,9,6,7)],df01[,15],df02[,15],df03[,15],df04[,15],df05[,15], df06[,15],df07[,15],df08[,15],df09[,15],df10[,15],df11[,15],df12[,15],df13[,15], df14[,15],df15[,15],df16[,15],df17[,15],df18[,15],df19[,15],df20[,15],df21[,15], df22[,15],df23[,15],df24[,15],df25[,15],df26[,15],df27[,15],df28[,15],df29[,15], df30[,15],df31[,15],df32[,15],df33[,15],df34[,15],df35[,15],df36[,15],df37[,15], df38[,15],df39[,15])
write.table(data1, "tempdata1.txt", col.names=F, row.names=F, sep="\t", quote=F)
write.table(data2, "tempdata2.txt", col.names=F, row.names=F, sep="\t", quote=F)
write.table(data3, "tempdata3.txt", col.names=F, row.names=F, sep="\t", quote=F)
#then use text editor to replace ":" "_ClassI;" and "_ClassII;"  to "_" ; replace ";" to nothing; files to be modified are “tempdata1.txt”, “tempdata2.txt” and “tempdata3.txt”

data01=read.table(textConnection(gsub(" ", "\t", readLines("tempdata1.txt"))))
data02=read.table(textConnection(gsub(" ", "\t", readLines("tempdata2.txt"))))
data03=read.table(textConnection(gsub(" ", "\t", readLines("tempdata3.txt"))))
data01=data01[,c(1:3, 7, 12:52)]
data02=data02[,c(1:3, 7, 12:52)]
data03=data03[,c(1:3, 7, 12:52)]
name=c("chr","start","end","TEm","score","str","Ctrl_a","Ctrl_b","Ctrl_c","Mock01_a", "Mock01_b","Mock01_c","Mock03_a","Mock03_b","Mock03_c","Mock06_a","Mock06_b", "Mock06_c","Mock12_a","Mock12_b","Mock12_c","Yeast01_a","Yeast01_b","Yeast01_c", "Yeast03_a","Yeast03_b","Yeast03_c","Yeast06_a","Yeast06_b","Yeast06_c", "Yeast12_a","Yeast12_b","Yeast12_c","Botrytis01_a","Botrytis01_b","Botrytis01_c", "Botrytis03_a","Botrytis03_b","Botrytis03_c","Botrytis06_a","Botrytis06_b", "Botrytis06_c","Botrytis12_a","Botrytis12_b","Botrytis12_c")
colnames(data01)=name
colnames(data02)=name
colnames(data03)=name
write.table(data01, "BedCov_OverlapBP_sense_count.txt", col.names=T, row.names=F, sep="\t", quote=F)
write.table(data02, "BedCov_OverlapBP_sense_breadthCov.txt", col.names=T, row.names=F, sep="\t", quote=F)
write.table(data03, "BedCov_OverlapBP_sense_depth.txt", col.names=T, row.names=F, sep="\t", quote=F)

count=read.table("BedCov_OverlapBP_sense_count.txt", header=T, sep="\t")
depth=read.table("BedCov_OverlapBP_sense_depth.txt", header=T, sep="\t")

#ctrl
df=count
df$Cmean00=rowMeans(df[,7:9])
data1 <- subset(df, df$Cmean00 > 10)
df=depth
df$Cmean00=rowMeans(df[,7:9])
data2 <-subset(df, df$Cmean00 > 5)
ctrl=merge(data1, data2, by="TEm")
ctrl=ctrl[,-c(10:45,47:51,55:90)]
colnames(ctrl)=c("TEm","chr","start","end","score","str","count_Ctrl_a","count_Ctrl_b", "count_Ctrl_c","count_Cmean00","depth_Ctrl_a","depth_Ctrl_b","depth_Ctrl_c", "depth_Cmean00")
write.table(ctrl, "BedCov_OverlapBP_PassC.txt", col.names=T, row.names=F, sep="\t", quote=F)

#mock
df=count
df=df[,c(1:6,10:21)]
df$Mmean01=rowMeans(df[,7:9])
df$Mmean03=rowMeans(df[,10:12])
df$Mmean06=rowMeans(df[,13:15])
df$Mmean12=rowMeans(df[,16:18])
min_read = 10
data1 <- df[apply(df[,19:22],1,function(x){max(x)}) > min_read,]
df=depth
df=df[,c(1:6,10:21)]
df$Mmean01=rowMeans(df[,7:9])
df$Mmean03=rowMeans(df[,10:12])
df$Mmean06=rowMeans(df[,13:15])
df$Mmean12=rowMeans(df[,16:18])
min_depth = 5
data2 <- df[apply(df[,19:22],1,function(x){max(x)}) > min_depth,]
data=merge(data1, data2, by="TEm")
data=data[,-c(23:27)]
colnames(data)=c("TEm","chr","start","end","score","str","count_Mock01_a", "count_Mock01_b","count_Mock01_c","count_Mock03_a","count_Mock03_b", "count_Mock03_c","count_Mock06_a","count_Mock06_b","count_Mock06_c", "count_Mock12_a","count_Mock12_b","count_Mock12_c","count_Mmean01", "count_Mmean03","count_Mmean06","count_Mmean12","depth_Mock01_a", "depth_Mock01_b","depth_Mock01_c","depth_Mock03_a","depth_Mock03_b", "depth_Mock03_c","depth_Mock06_a","depth_Mock06_b","depth_Mock06_c", "depth_Mock12_a","depth_Mock12_b","depth_Mock12_c","depth_Mmean01", "depth_Mmean03","depth_Mmean06","depth_Mmean12")
write.table(data, "BedCov_OverlapBP_PassM.txt", col.names=T, row.names=F, sep="\t", quote=F)

#yeast
df=count
df=df[,c(1:6,22:33)]
df$Mmean01=rowMeans(df[,7:9])
df$Mmean03=rowMeans(df[,10:12])
df$Mmean06=rowMeans(df[,13:15])
df$Mmean12=rowMeans(df[,16:18])
min_read = 10
data1 <- df[apply(df[,19:22],1,function(x){max(x)}) > min_read,]
df=depth
df=df[,c(1:6,22:33)]
df$Mmean01=rowMeans(df[,7:9])
df$Mmean03=rowMeans(df[,10:12])
df$Mmean06=rowMeans(df[,13:15])
df$Mmean12=rowMeans(df[,16:18])
min_depth = 5
data2 <- df[apply(df[,19:22],1,function(x){max(x)}) > min_depth,]
data=merge(data1, data2, by="TEm")
data=data[,-c(23:27)]
colnames(data)=c("TEm","chr","start","end","score","str","count_Yeast01_a", "count_Yeast01_b","count_Yeast01_c","count_Yeast03_a","count_Yeast03_b", "count_Yeast03_c","count_Yeast06_a","count_Yeast06_b","count_Yeast06_c", "count_Yeast12_a","count_Yeast12_b","count_Yeast12_c","count_Ymean01", "count_Ymean03","count_Ymean06","count_Ymean12","depth_Yeast01_a", "depth_Yeast01_b","depth_Yeast01_c","depth_Yeast03_a","depth_Yeast03_b", "depth_Yeast03_c","depth_Yeast06_a","depth_Yeast06_b","depth_Yeast06_c", "depth_Yeast12_a","depth_Yeast12_b","depth_Yeast12_c","depth_Ymean01", "depth_Ymean03","depth_Ymean06","depth_Ymean12")
write.table(data, "BedCov_OverlapBP_PassY.txt", col.names=T, row.names=F, sep="\t", quote=F)
