args <- commandArgs(trailingOnly = TRUE)
### using reference fasta, TDNA fasta files, Mosaic REF_TDNA, Mosaic length  as input and output mosaic diagram in pdf files and mosaic sequence in fasta files
library(ShortRead)
##agrs[1] = reference fasta files
##args[2] = TDNA fasta files
##args[3] = Mosaic REF_TDNA files
##args[4] = Output Mosaic pdf  files
##args[5] = Output Mosaic fasta files
##args[6] = Sample names
##args[7] = Mosaic length
Sample_name =""
if(length(args)>5){
        Sample_name =  args[6]
}

Mosaic_length = 1000
if(length(args)>6){
        Mosaic_length = as.integer( args[7])
}

library(ShortRead)
genome = readDNAStringSet(args[1], "fasta")
TDNA = readDNAStringSet(args[2], "fasta")
TDNA [[1]]= toupper(TDNA[[1]])
names(TDNA)[1] = gsub(" .*","", names(TDNA)[1])
Mosaic = read.table(args[3],header = TRUE,sep="\t")
ID = data.frame(table(Mosaic[,10]))
##Type1, RefPos_TDNAPos, Type2, RefPos_TDNANeg
##Type3, TDNAPos_RefPos, Type4, TDNANeg_RefPos
##ID1 Mosaic on the left side
##ID2 Mosaic on the right side
ID1 =data.frame(table(Mosaic[grepl(":RefPos_",Mosaic[,9]),10]))
tmp = gsub(":.*$","", ID1[,1])
Pos_REF = as.numeric(gsub(".*__","", tmp))
Chr_REF =gsub("REF__|__\\d+$","", tmp)
Pos_TDNA =as.numeric(gsub(".*__","", ID1[,1]))
ID1 = data.frame(ID1,Chr_REF, Pos_REF, Pos_TDNA)
ID2 =data.frame(table(Mosaic[grepl("_RefPos$",Mosaic[,9]),10]))
tmp = gsub(":.*$","", ID2[,1])
Pos_REF = as.numeric(gsub(".*__","", ID2[,1]))
Chr_REF =gsub(".*:|REF__|__\\d+$","", ID2[,1])
Pos_TDNA =as.numeric(gsub(".*__","", tmp))
ID2 = data.frame(ID2,Chr_REF, Pos_REF, Pos_TDNA)
if(dim(ID1)[1]>1){
	ID1 = ID1[order(ID1[,3], ID1[,4], ID1[,5]),]
}
if(dim(ID2)[1]>1){
	ID2 = ID2[order(ID2[,3], ID2[,4], ID2[,5]),]
}
##divided left mosaic side and right mosaic side into two groups
insert_l = character(0)
insert_r = character(0)
m = 1
n = 1
if(dim(ID1)[1]==0){
	for(n in 1:length(ID2[,1])){
		insert_l = append(insert_l,  "")
		insert_r = append(insert_r, as.character(ID2[n,1]) )
	}
} else if (dim(ID2)[1]==0){
        for(m in 1:length(ID1[,1])){
                insert_l = append(insert_l,  as.character(ID1[m,1]))
                insert_r = append(insert_r, "" )
        }
} else while((m<= length(ID1[,1])) |( n<= length(ID2[,1]))){
    if((m> length(ID1[,1])) & (n> length(ID2[,1]))){
        break
    }else if((m> length(ID1[,1]) )& (n<= length(ID2[,1]))){
        insert_l = append(insert_l,  "")
        insert_r = append(insert_r, as.character(ID2[n,1]))
        n = n+1
    } else if((m<= length(ID1[,1]) )& (n> length(ID2[,1]))){

        insert_l = append(insert_l, as.character(ID1[m,1]))
        insert_r = append(insert_r,  "")
        m= m+1
    }else{

        if((ID1[m,3]==ID2[n,3]) & (abs(ID1[m,4]- ID2[n,4]) < 10000) ){
            insert_l = append(insert_l, as.character(ID1[m,1]))
            insert_r = append(insert_r, as.character(ID2[n,1]))
            m= m+1
            n = n+1
        } else if((n+1<= length(ID2[,1])) & (ID1[m,3]==ID2[n+1,3]) & (abs( ID1[m,4]- ID2[n+1,4]) < 10000)  ){
            insert_l = append(insert_l, "")
            insert_r = append(insert_r, as.character(ID2[n,1]))
            insert_l = append(insert_l, as.character(ID1[m,1]))
            insert_r = append(insert_r, as.character(ID2[n+1,1]))
            m = m+1
            n = n+2

        }else if((m+1<= length(ID1[,1]) )& (ID1[m+1,3]==ID2[n,3]) & (abs(ID1[m+1,4]- ID2[n,4] ) < 10000)  ){
            insert_l = append(insert_l, as.character(ID1[m,1]))
            insert_r = append(insert_r, "")
            insert_l = append(insert_l, as.character(ID1[m+1,1]))
            insert_r = append(insert_r, as.character(ID2[n,1]))
            m = m+2
            n = n+1
	 }else if((length(ID1[,1])-m) > (length(ID2[,1])-n)){
                insert_l = append(insert_l, as.character(ID1[m,1]))
                insert_r = append(insert_r, "")
                m = m+1
        } else if((length(ID2[,1])-n) >= (length(ID1[,1])-m)){
                insert_l = append(insert_l, "")
                insert_r = append(insert_r, as.character(ID2[n,1]))
                n = n+1
        }
        }

}

index = c("01","02","03","04","05","06","07","08","09","10","11","12","13","14","15","16","17","18","19","10","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42","43","44","45","46","47","48","49","50")
pdf(file = args[4], 8.5,11)
for (num in 1:length(insert_l)){

    A1 = Mosaic[Mosaic[,10] %in% insert_l[num],]
    A2 = Mosaic[Mosaic[,10] %in% insert_r[num],]
    num_read = max(length(A1[,1]), length(A2[,1]))
    par(mfrow = c(3,1))
    read_l = 150
    if(dim(A1)[1]>0)
    {
        read_l = nchar(A1[1,2])
    }
    if(dim(A2)[1]>0)
    {
        read_l = nchar(A2[1,2])
    }
    plot(c(1, read_l* 5), c(1, (3+num_read)*40), type= "n", xlab = "", ylab = "",axes=FALSE, main = paste(Sample_name,"_Mosaic_Site_",num,sep=""))
    rect(read_l* 2, 30, read_l* 3, 70,border = "purple",density=10)
    if(dim(A1)[1]>0){
        rect(1, 30, read_l, 70,col="blue",border = "blue")
        rect(read_l, 30, read_l* 2, 70,col="purple",border = "purple")
        legend("top" ,legend = c("REF","T-DNA"),col= c("blue","purple"),fill = c("blue","purple"))
        for( i in 1:length(A1[,1]))
        {
            rect(read_l-nchar(A1[i,4]),120+ (i-1)* 40,read_l, 140+ (i-1)* 40,col="blue",border = "blue" )
            rect(read_l,120+ (i-1)* 40,read_l + nchar(A1[i,5]), 140+  (i-1)* 40,col="purple",border = "purple")
        }
        label_REF = paste("REF:", ID1[ID1[,1] == insert_l[num],3][1],"-",ID1[ID1[,1] == insert_l[num],4][1],sep="" )
        label_TDNA = paste("TDNA:", ID1[ID1[,1] == insert_l[num],5][1],sep="" )
        text(read_l/2,90, label_REF,col="blue")
        text(read_l/2,1,  label_TDNA,col="purple")
        if(A1[1,9]=="type1:RefPos_TDNAPos"){
            arrows(read_l, 10,read_l+60, 10, length = 0.1, angle = 5,code = 2)
        }
        if(A1[1,9]=="type2:RefPos_TDNANeg"){
            arrows(read_l+60, 10,read_l, 10, length = 0.1, angle = 5,code = 2)
        }
    }

    if(dim(A2)[1]>0){
        rect(read_l*3, 30, read_l*4, 70,col="purple",border = "purple")
        rect(read_l*4, 30, read_l*5, 70,col="blue",border = "blue")
        for( i in 1:length(A2[,1])){
            rect(read_l*4-nchar(A2[i,5]),120+ (i-1)* 40,read_l*4, 140+ (i-1)* 40 ,col="purple",border = "purple")
            rect(read_l*4,120+ (i-1)* 40,read_l*4+nchar(A2[i,4]), 140+  (i-1)* 40,col="blue",border = "blue")
        }
        label_REF = paste("REF:", ID2[ID2[,1] == insert_r[num],3][1],"-",ID2[ID2[,1] == insert_r[num],4][1],sep="" )
        label_TDNA = paste("TDNA:", ID2[ID2[,1] == insert_r[num],5][1],sep="" )

        text(read_l*4,90,  label_REF,col="blue")
        text(read_l*4,1,  label_TDNA ,col="purple")
        if(A2[1,9]=="type3:TDNAPos_RefPos"){
            arrows(read_l*3 + 50, 10,read_l*3 + 110, 10, length = 0.1, angle = 5,code = 2)
        }
        if(A2[1,9]=="type4:TDNANeg_RefPos"){
            arrows(read_l*3 + 110, 10,read_l*3 + 50, 10, length = 0.1, angle = 5,code = 2)
        }
    }
    if(dim(A1)[1]>0){
        plot(c(1, read_l*5), c(1, (1+length(A1[,2]))*60), type= "n", xlab = "", ylab = "",axes=FALSE,main = "Read sequences: 5'--sequence at reference----sequence at TDNA--3'")
        cc =0.6
        for( i in 1:length(A1[,1])){
            seq1 =paste("Read",index[i], ": 5'--",A1[i,4],"----",sep="")
            seq2 = paste(A1[i,5],"--3'",sep="")
            if(A1[i,9]=="type2:RefPos_TDNANeg")
            {
                seq2 = paste(as.character(reverseComplement(DNAString(A1[i,5]))),"--3'",sep="")
            }
            text(20,50 + i* 60, seq1,col="blue",cex =cc ,pos = 4)
            text(20,30 + i* 60, seq2,col="purple",cex = cc,pos = 4)
        }
    }
    if(dim(A2)[1]>0){
        plot(c(1, read_l*5), c(1, (1+length(A2[,2]))*60), type= "n", xlab = "", ylab = "",axes=FALSE,main ="Read sequences: 5'--sequence at TDNA----sequence at reference--3'")
        cc = 0.6
        for( i in 1:length(A2[,1])){
            seq1 =paste("Read",index[i], ": 5'--",A2[i,5],"----",sep="")
            seq2 = paste(A2[i,4],"--3'",sep="")
            if(A2[i,9]=="type4:TDNANeg_RefPos")
            {
                seq1 =paste("Read",index[i], ": 5'--",as.character(reverseComplement(DNAString(A2[i,5]))),"----",sep="")
            }
            text(20,50 + i* 60, seq1,col="purple",cex = cc,pos = 4)
            text(20,30 + i* 60, seq2,col="blue",cex = cc,pos = 4)
        }
    }
}

dev.off()

Seq_all = BStringSet()
for (num in 1:length(insert_l)){
    A1 = Mosaic[Mosaic[,10] %in% insert_l[num],]
    A2 = Mosaic[Mosaic[,10] %in% insert_r[num],]
    id1 =ID1[ID1[,1] %in% insert_l[num],]
    id2 =ID2[ID2[,1] %in% insert_r[num],]
    if(dim(A1)[1]>0){
        A1seq_Ref = as.character(Views(genome[names(genome)==as.character(id1[1,3])][[1]],start = max(1,id1[1,4]-Mosaic_length+1), end= id1[1,4]))
        A1seq_TDNA =as.character(Views(TDNA[[1]],start = id1[1,5], end= min(id1[1,5]+ Mosaic_length-1,width(TDNA))))
        if(A1[1,9]=="type2:RefPos_TDNANeg"){

            A1seq_TDNA =as.character(reverseComplement(Views(TDNA[[1]],start = max(1,id1[1,5]-Mosaic_length+1), end=id1[1,5])))
        }
        A1seq = BStringSet(paste(A1seq_Ref, A1seq_TDNA ,sep="----"))
        names(A1seq) =  paste(Sample_name,"::",paste("Mosaic_Site_",num,sep=""),"::REF",Mosaic_length, "--TDNA",Mosaic_length,"::", id1[1,1],sep="")
        Seq_all = c(Seq_all, A1seq)
    }
    if(dim(A2)[1]>0){
        A2seq_TDNA =as.character(Views(TDNA[[1]],start = max(1,id2[1,5]-Mosaic_length+1), end=id2[1,5]))
        if(A2[1,9]=="type4:TDNANeg_RefPos"){
            A2seq_TDNA =as.character(reverseComplement(Views(TDNA[[1]],start = id2[1,5], end= min(id2[1,5]+ Mosaic_length-1 ,width(TDNA)))))
        }
        A2seq_Ref = as.character(Views(genome[names(genome)==as.character(id2[1,3])][[1]],start = id2[1,4], end= min(id2[1,4]+ Mosaic_length-1,width(genome[names(genome)==as.character(id2[1,3])]))))
        A2seq = BStringSet(paste( A2seq_TDNA ,A2seq_Ref,sep="----"))
        names(A2seq) =paste(Sample_name,"::",paste("Mosaic_Site_",num,sep=""),"::TDNA",Mosaic_length, "--REF",Mosaic_length,"::", id2[1,1],sep="")
        Seq_all = c(Seq_all, A2seq)
    }
}

writeXStringSet(Seq_all ,args[5],compress=FALSE)


