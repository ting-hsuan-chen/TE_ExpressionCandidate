#R code
library(VennDiagram)

cal.overlap=function (x) {
    if (1 == length(x)) {
        overlap <- x
    }
    else if (2 == length(x)) {
        A <- x[[1]]
        B <- x[[2]]
        nab <- intersect(A, B)
        a1 = A[! A %in% nab]
        a2 = B[! B %in% nab]
        overlap <- list(a1 = a1, a2 = a2, a3 = nab)
    }
    else if (3 == length(x)) {
        A <- x[[1]]
        B <- x[[2]]
        C <- x[[3]]
        nab <- intersect(A, B)
        nbc <- intersect(B, C)
        nac <- intersect(A, C)
        nabc <- intersect(nab, C)
        a5 = nabc
        a2 = nab[! nab %in% a5]
        a4 = nac[-which(nac %in% a5)]
        a6 = nbc[-which(nbc %in% a5)]
        a1 = A[-which(A %in% c(a2, a4, a5))]
        a3 = B[-which(B %in% c(a2, a5, a6))]
        a7 = C[-which(C %in% c(a4, a5, a6))]
        overlap <- list(a5 = a5, a2 = a2, a4 = a4, a6 = a6, a1 = a1,
            a3 = a3, a7 = a7)
    }
    else if (4 == length(x)) {
        A <- x[[1]]
        B <- x[[2]]
        C <- x[[3]]
        D <- x[[4]]
        nab <- intersect(A, B)
        nbc <- intersect(B, C)
        nac <- intersect(A, C)
        nbd <- intersect(B, D)
        nad <- intersect(A, D)
        ncd <- intersect(C, D)
        nabc <- intersect(nab, C)
        nabd <- intersect(nab, D)
        nacd <- intersect(nac, D)
        nbcd <- intersect(nbc, D)
        nabcd <- intersect (nabc, D)
        a6  = nabcd
        a12 = nabc[-which(nabc %in% a6)]
        a11 = nabd[-which(nabd %in% a6)]
        a5  = nacd[-which(nacd %in% a6)]
        a7  = nbcd[-which(nbcd %in% a6)]
        a15 = nab[-which(nab %in% c(a6,a11,a12))]
        a4  = nac[-which(nac %in% c(a6,a5,a12))]
        a10 = nad[-which(nad %in% c(a6,a5,a11))]
        a13 = nbc[-which(nbc %in% c(a6,a7,a12))]
        a8  = nbd[-which(nbd %in% c(a6,a7,a11))]
        a2  = ncd[-which(ncd %in% c(a6,a5,a7))]
        a9  = A[-which(A %in% c(a4,a5,a6,a10,a11,a12,a15))]
        a14 = B[-which(B %in% c(a6,a7,a8,a11,a12,a13,a15))]
        a1  = C[-which(C %in% c(a2,a4,a5,a6,a7,a12,a13))]
        a3  = D[-which(D %in% c(a2,a5,a6,a7,a8,a10,a11))]
        overlap <- list(a6=a6, a12=a12, a11=a11, a5=a5, a7=a7, a15=a15,
        a4=a4, a10=a10, a13=a13, a8=a8, a2=a2, a9=a9, a14=a14, a1=a1, a3=a3)
    }
    else {
        flog.error("Invalid size of input object", name = "VennDiagramLogger")
        stop("Invalid size of input object")
    }
    return(overlap)
}

#ctrl
RE=read.table("k100Mm_AllSenseTEm_curatedGTF_exprCtrl05_V4.txt", header=T, sep="\t") #from sub-pipeline 1
df=read.table("BedCov_OverlapBP_PassC.txt", header=T, sep="\t") #from sub-pipeline 2
dangler=read.table("AllRep_comp_intersect_exprCtrl_V4.txt", header=T, sep="\t") #from sub-pipeline 3

RE=RE[!grepl("__", RE$id),]
RE=RE[order(RE$id),]
dangler=dangler[order(dangler$id),]
df=df[order(df$TEm),]

REuniq=RE[!duplicated(RE$id),]
dfuniq=df[!duplicated(df$TEm),]
dangleruniq=dangler[!duplicated(dangler$id),]

overlap<-cal.overlap(x=list("count&depth"=dfuniq[,1], "htseq-count"=REuniq[,1], "TEFingerprint"=dangleruniq[,1]))
capture.output(summary(overlap), file = "OverlapSummary3groups_exprCand_ctrl.txt")
capture.output(print(overlap), file = "OverlapPrint3groups_exprCand_ctrl.txt")
a2=data.frame(id=overlap$a2)
a4=data.frame(id=overlap$a4)
a5=data.frame(id=overlap$a5)
a1=data.frame(id=overlap$a1)
a3=data.frame(id=overlap$a3)
a6=data.frame(id=overlap$a6)
a7=data.frame(id=overlap$a7)

a2$C_tracking="countdepth_htseq"
a4$C_tracking="countdepth_TEF"
a5$C_tracking="countdepth_htseq_TEF"
a1$C_tracking="countdepthonly"
a3$C_tracking="htseqonly"
a6$C_tracking="htseqTEFonly"
a7$C_tracking="TEFonly"
allexpr=rbind(a1,a2,a3,a4,a5,a6,a7)
write.table(allexpr,"AllExpeCandidate_tracking_ctrl.txt", col.names=T, row.names=F, sep="\t", quote=F)

setwd("./analysis/ECstress_TEalignment/ExprCandidate")
df=read.table("AllExpeCandidate_tracking_ctrl.txt", header=T, sep="\t")
ref=read.table("../../Reference/Curated/AllTEsExpanded_curated_tags_V4.txt", header=T, sep="\t")
newdf=merge(ref, df, by="id")
write.table(newdf,"AllExpeCandidate_ctrl_tag_new.txt", col.names=T, row.names=F, sep="\t", quote=F)

#mock
RE=read.table("k100Mm_AllSenseTEm_curatedGTF_exprMock05_V4. txt", header=T, sep="\t") #from sub-pipeline 1
dangler=read.table("/AllRep_comp_intersect_exprMock_V4.txt", header=T, sep="\t") #from sub-pipeline 2
df=read.table("BedCov_OverlapBP_PassM.txt", header=T, sep="\t") #from sub-pipeline 3

RE=RE[!grepl("__", RE$id),]
RE=RE[order(RE$id),]
dangler=dangler[order(dangler$id),]
df=df[order(df$TEm),]

dfuniq=df[!duplicated(df$TEm),]
REuniq=RE[!duplicated(RE$id),]
dangleruniq=dangler[!duplicated(dangler$id),]

overlap<-cal.overlap(x=list("count&depth"=dfuniq[,1], "htseq-count"=REuniq[,1], "TEFingerprint"=dangleruniq[,1]))
capture.output(summary(overlap), file = "OverlapSummary3groups_exprCand_mock.txt")
capture.output(print(overlap), file = "OverlapPrint3groups_exprCand_mock.txt")
a2=data.frame(id=overlap$a2)
a4=data.frame(id=overlap$a4)
a5=data.frame(id=overlap$a5)
a1=data.frame(id=overlap$a1)
a3=data.frame(id=overlap$a3)
a6=data.frame(id=overlap$a6)
a7=data.frame(id=overlap$a7)

a2$C_tracking="countdepth_htseq"
a4$C_tracking="countdepth_TEF"
a5$C_tracking="countdepth_htseq_TEF"
expr1=rbind(a2,a4,a5)
a1$C_tracking="countdepthonly"
a3$C_tracking="htseqonly"
a6$C_tracking="htseqTEFonly"
a7$C_tracking="TEFonly"
allexpr=rbind(a1,a2,a3,a4,a5,a6,a7)
write.table(allexpr,"AllExpeCandidate_tracking_mock.txt", col.names=T, row.names=F, sep="\t", quote=F)

setwd("./analysis/ECstress_TEalignment/ExprCandidate")
df=read.table("AllExpeCandidate_tracking_mock.txt", header=T, sep="\t")
ref=read.table("../../Reference/Curated/AllTEsExpanded_curated_tags_V4.txt", header=T, sep="\t")
newdf=merge(ref, df, by="id")
write.table(newdf,"AllExpeCandidate_mock_tag_new.txt", col.names=T, row.names=F, sep="\t", quote=F)


#yeast
RE=read.table("k100Mm_AllSenseTEm_curatedGTF_exprYeast05_V4.txt", header=T, sep="\t") #from sub-pipeline 1
dangler=read.table("AllRep_comp_intersect_exprYeast_V4.txt", header=T, sep="\t") #from sub-pipeline 2
df=read.table("BedCov_OverlapBP_PassY.txt", header=T, sep="\t") #from sub-pipeline 3

RE=RE[!grepl("__", RE$id),]
RE=RE[order(RE$id),]
dangler=dangler[order(dangler$id),]
df=df[order(df$TEm),]

dfuniq=df[!duplicated(df$TEm),]
REuniq=RE[!duplicated(RE$id),]
dangleruniq=dangler[!duplicated(dangler$id),]

overlap<-cal.overlap(x=list("count&depth"=dfuniq[,1], "htseq-count"=REuniq[,1], "TEFingerprint"=dangleruniq[,1]))
capture.output(summary(overlap), file = "OverlapSummary3groups_exprCand_yeast.txt")
capture.output(print(overlap), file = "OverlapPrint3groups_exprCand_yeast.txt")
a2=data.frame(id=overlap$a2)
a4=data.frame(id=overlap$a4)
a5=data.frame(id=overlap$a5)
a1=data.frame(id=overlap$a1)
a3=data.frame(id=overlap$a3)
a6=data.frame(id=overlap$a6)
a7=data.frame(id=overlap$a7)

a2$C_tracking="countdepth_htseq"
a4$C_tracking="countdepth_TEF"
a5$C_tracking="countdepth_htseq_TEF"
expr1=rbind(a2,a4,a5)
a1$C_tracking="countdepthonly"
a3$C_tracking="htseqonly"
a6$C_tracking="htseqTEFonly"
a7$C_tracking="TEFonly"
allexpr=rbind(a1,a2,a3,a4,a5,a6,a7)
write.table(allexpr,"AllExpeCandidate_tracking_yeast.txt", col.names=T, row.names=F, sep="\t", quote=F)

setwd("./analysis/ECstress_TEalignment/ExprCandidate")
df=read.table("AllExpeCandidate_tracking_yeast.txt", header=T, sep="\t")
ref=read.table(" ../../Reference/Curated/AllTEsExpanded_curated_tags_V4.txt", header=T, sep="\t")
newdf=merge(ref, df, by="id")
write.table(newdf,"AllExpeCandidate_yeast_tag_new.txt", col.names=T, row.names=F, sep="\t", quote=F)
