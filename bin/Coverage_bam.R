args <- commandArgs(trailingOnly = TRUE)
### using  Mosaic REF_TDNA, Reference alignment files  as input and output coverage at the mosaic sites in pdf files 
##args[1] = Mosaic REF_TDNA files
##args[2] = Reference Alignment files
##args[3] = Output coverage pdf  files
##args[4] = folder to put middle Bam files at mosaic region with flanking 2000bp and coverage data
##args[5] = Sample_name
Sample_name = ""
if(length(args)>4){
	Sample_name  = args[5]
}
Mosaic = read.table(args[1],header = TRUE,sep="\t")
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
} else while(m<= length(ID1[,1]) | n<= length(ID2[,1])){
    if(m> length(ID1[,1]) & n> length(ID2[,1])){
        break
    }else if(m> length(ID1[,1]) & n<= length(ID2[,1])){
        insert_l = append(insert_l,  "")
        insert_r = append(insert_r, as.character(ID2[n,1]))
        n = n+1
    } else if(m<= length(ID1[,1]) & n> length(ID2[,1])){

        insert_l = append(insert_l, as.character(ID1[m,1]))
        insert_r = append(insert_r,  "")
        m= m+1
    }else{

        if(ID1[m,3]==ID2[n,3] & abs(ID1[m,4]- ID2[n,4]) < 10000 ){
            insert_l = append(insert_l, as.character(ID1[m,1]))
            insert_r = append(insert_r, as.character(ID2[n,1]))
            m= m+1
            n = n+1
        } else if(n+1<= length(ID2[,1]) & ID1[,3]==ID2[n+1,3] & abs( ID1[m,4]- ID2[n+1,4]) < 10000  ){
            insert_l = append(insert_l, "")
            insert_r = append(insert_r, as.character(ID2[n,1]))
            insert_l = append(insert_l, as.character(ID1[m,1]))
            insert_r = append(insert_r, as.character(ID2[n+1,1]))
            m = m+1
            n = n+2

        }else if(m+1<= length(ID1[,1]) & ID1[m+1,3]==ID2[n,3] & abs(ID1[m+1,4]- ID2[n,4] ) < 10000  ){
            insert_l = append(insert_l, as.character(ID1[m,1]))
            insert_r = append(insert_r, "")
            insert_l = append(insert_l, as.character(ID1[m+1,1]))
            insert_r = append(insert_r, as.character(ID2[n,1]))
            m = m+2
            n = n+1
	}
	else if((length(ID1[,1])-m) > (length(ID2[,1])-n)){
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


pdf(file = args[3], 8.5,11)
par(mfrow = c(4,1))
for (num in 1:length(insert_l)){
    A1 = Mosaic[Mosaic[,10] %in% insert_l[num],]
    A2 = Mosaic[Mosaic[,10] %in% insert_r[num],]
    id1 =ID1[ID1[,1] %in% insert_l[num],]
    id2 =ID2[ID2[,1] %in% insert_r[num],]
    region1 = ""
    region2 = ""
    if(dim(id1)[1]>0 & dim(id2)[1]>0){
        dis = abs(id2[1,4]-id1[1,4])
        start =min(id2[1,4],id1[1,4])-dis
        end =max(id2[1,4],id1[1,4])+ dis
        region1 = paste(id1[1,3],":", start,"-",end,sep="")
        region2 = paste(id1[1,3],":", start - 2000,"-",end + 2000,sep="")
    }
    if(dim(id2)[1]==0 ){
        region1 = paste(id1[1,3],":", id1[1,4]-nchar(A1[1,2]),"-", id1[1,4]+ nchar(A1[1,2]),sep="")
        region2 = paste(id1[1,3],":", id1[1,4]-nchar(A1[1,2])- 2000,"-", id1[1,4]+ nchar(A1[1,2])+ 2000,sep="")
    }
    if(dim(id1)[1] == 0 ){
        region1 = paste(id2[1,3],":", id2[1,4]-nchar(A2[1,2]),"-", id2[1,4]+ nchar(A2[1,2]),sep="")
        region2 = paste(id2[1,3],":", id2[1,4]-nchar(A2[1,2])- 2000,"-", id2[1,4]+ nchar(A2[1,2])+ 2000,sep="")
    }
    cmd1 = paste("samtools view -bhS ", args[2], " ", region2, " -o ",args[4], "/",args[5], "_Mosaic_",num,".bam", sep="")
    system(cmd1)
    cmd2 = paste("samtools index ", args[4],"/",args[5], "_Mosaic_",num,".bam",sep="")
    system(cmd2)
    cmd3 = paste("python /home/shaofangli/SOFT/T_LOC/bin/Coverage.py --bam ", args[4], "/",args[5], "_Mosaic_",num,".bam", " --Region ", region1 ," --output ", args[4],"/",args[5], "_Mosaic_",num,"_coverage.txt",sep="")
    system(cmd3)
    cov = read.table(paste(args[4],"/",args[5], "_Mosaic_",num,"_coverage.txt",sep=""))
    barplot(cov[,2],ylim = c(0, as.integer(max(cov[,2]) * 2/100+1) * 100),las =1,main = Sample_name,names.arg = cov[,1],xlab = gsub(":.*","",region1),ylab = "Read number")
}
dev.off()
