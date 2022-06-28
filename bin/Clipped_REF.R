args <- commandArgs(trailingOnly = TRUE)
### using reference fasta, TDNA fasta files, Reference clipped sam files  and control_ID as input
library(ShortRead)
##agrs[1] = reference fasta files
##args[2] = TDNA fasta files
##args[3] = REF clipped sam files
##args[4] = REF control ID files
##args[5] = output Refence_TDNA mosaic file
##args[6] = output Ref_Ref  mosaic file
##args[7] = output All Ref  mosaic file
##args[8] = archor length default is 20
##args[9] = minimal required mosaic reads numbers, default is 2
##args[10] = minimal required mosaic reads numbers, default is 4
##args[11] = maximal allowed sequence insert length between mosaice reads, default = 30
anchor = 20
if(length(args)>7){
	anchor = as.integer( args[8])
}

lim = 2
if(length(args)>8){
        lim = as.integer( args[9])
}

lim2 = 4
if(length(args)>9){
        lim2 = as.integer( args[10])
}

insert = 30
if(length(args)>10){
        insert = as.integer( args[11])
}


##step1, get rid of read alignment with deletions and mismatches
##step2, get rid of clipped reads with clipping on both sides
##step3, get rid of clipped reads with anchor with shorter length and clipped from the arbitrary start and end site of TDNA
##step4: seprate reads into two type: type1, mismatch-TDNAmatch, type2, TDNAmatch-mismatch
##step4, get rid of clipped sites with minimal number of clipped reads 
##step5, search for reads with matches from TDNA or reference seqeunces
##step6, get rid of reads fully matches to the reference genome
##step7, output the mosaic reads in three diffrent files
genome = readDNAStringSet(args[1], "fasta")
TDNA = readDNAStringSet(args[2], "fasta")
TDNA [[1]]= toupper(TDNA[[1]])
names(TDNA)[1] = gsub(" .*","", names(TDNA)[1]) 
c_id = read.table(args[4])
c_id = c_id[,1]
a=read.delim(file=args[3],header=F,col.names=1:25)
a1 = a[!duplicated(a[,10]),]
a1 = a1[!grepl("I|D",a1[,6]),]
a1 = a1[!grepl("A|T|C|G",a1[,13]),]
type = rep("type1",length(a1[,1]))
a1 = data.frame(a1[,c(1,3,6,10)],start = a1[,4],end = a1[,4]-1 + as.numeric( gsub("\\d+S|M","",as.character(a1[,6]))),type)
colnames(a1)[1:4] = c("read_id","Chr","CIGAR","read_seq")
a1[grepl("^\\d+S",a1[,3]) & grepl("\\d+M$",a1[,3]),7] = "type1"
a1[grepl("^\\d+M",a1[,3]) & grepl("\\d+S$",a1[,3]),7] ="type2"
a1[grepl("^\\d+S",a1[,3]) & grepl("\\d+S$",a1[,3]),7] ="type3"
a1 = a1[a1[,7]!="type3",]

B1 = a1[a1[,7]=="type1",]
B1 = data.frame(B1, anchor = as.numeric( gsub("S.*","",as.character(B1[,3]))),seq = character(length(B1[,1])),matchREF_P = numeric(length(B1[,1])),matchREF_N = numeric(length(B1[,1])),matchTDNA_P = numeric(length(B1[,1])),matchTDNA_N = numeric(length(B1[,1])), matchREF_Mismatch_P = numeric(length(B1[,1])),matchREF_Mismatch_N = numeric(length(B1[,1])),matchTDNA_Mismatch_P = numeric(length(B1[,1])),matchTDNA_Mismatch_N = numeric(length(B1[,1])),matchREF_Insert_P = numeric(length(B1[,1])),matchREF_Insert_N = numeric(length(B1[,1])),matchTDNA_Insert_P = numeric(length(B1[,1])),matchTDNA_Insert_N = numeric(length(B1[,1])), ID = paste(B1[,2],B1[,5], "type1",sep="_"),REF_match = rep("false", length(B1[,1])))
B1 = B1[!B1[,22]%in%c_id,]
B1[,9] = substr(as.character(B1[,4]),1,B1[,8])
B1 = B1[order(B1[,5],-B1[,8]),]
B1 = B1[B1[,8]> anchor & B1[,8]< nchar(B1[1,4]) -anchor,]
Bt1 = table(B1[,22])
Bt1 =Bt1[Bt1>=lim]
B1 = B1[B1[,22] %in% names(Bt1),]

C1 = a1[a1[,7]=="type2",]
C1 = data.frame(C1, anchor = as.numeric( gsub("\\d+M|S","",as.character(C1[,3]))),seq = character(length(C1[,1])), matchREF_P = numeric(length(C1[,1])),matchREF_N = numeric(length(C1[,1])),matchTDNA_P = numeric(length(C1[,1])),matchTDNA_N = numeric(length(C1[,1])), matchREF_Mismatch_P = numeric(length(C1[,1])),matchREF_Mismatch_N = numeric(length(C1[,1])),matchTDNA_Mismatch_P = numeric(length(C1[,1])),matchTDNA_Mismatch_N = numeric(length(C1[,1])),matchREF_Insert_P = numeric(length(C1[,1])),matchREF_Insert_N = numeric(length(C1[,1])),matchTDNA_Insert_P = numeric(length(C1[,1])),matchTDNA_Insert_N = numeric(length(C1[,1])),ID = paste(C1[,2], C1[,6], "type2",sep="_"), REF_match = rep("false", length(C1[,1])))
C1 = C1[!C1[,22]%in%c_id,]
C1[,9] = substr(as.character(C1[,4]),nchar(C1[,4])-C1[,8]+1,nchar(C1[,4]))
C1 = C1[order(C1[,6],-C1[,8]),]
C1 = C1[C1[,8]> anchor & C1[,8]< nchar(C1[1,4]) -anchor,]
Ct1 = table(C1[,22])
Ct1 =Ct1[Ct1>= lim]
C1 = C1[C1[,22] %in% names(Ct1),]


MATCH_REF_mosaic = function(genome, TDNA, clip_read){
        s = as.character(clip_read[9])
        ss = as.character(reverseComplement(DNAString(s)))
    clip_read[12] =countPattern(s, TDNA[[1]])
    clip_read[13] =countPattern(ss, TDNA[[1]])
    search1 = vcountPattern(s, genome)
    search2 = vcountPattern(ss, genome)
    clip_read[10] =sum(search1)
    clip_read[11] =sum(search2)
    if(clip_read[12]>0){
        pp = sum(vcountPattern(as.character(clip_read[4]), TDNA, max.mismatch =1,with.indels=TRUE))
        if(pp>0){
                clip_read[23] = paste(gsub("false","", clip_read[23]),"perfectPos_TDNAGenome",sep=";")
        }
    }
    if(clip_read[13]>0){
        pp = sum(vcountPattern(as.character(reverseComplement(DNAString(as.character(clip_read[4])))), TDNA, max.mismatch =1,with.indels=TRUE))
        if(pp>0){
                clip_read[23] = paste(gsub("false","", clip_read[23]),"perfectNeg_TDNAGenome",sep=";")
        }
    }
        if(clip_read[23]!="false"){
                 return(clip_read)
        }
    if(clip_read[10]>0){
            tt = data.frame(names(genome),search1)
            tt =tt[tt[,2]>0,]
            match_REF = data.frame(unlist(vmatchPattern(s, genome[names(genome)%in% tt[tt[,2]>0,1]])))
            ID = paste(match_REF[1,4], match_REF[1,1],match_REF[1,2],sep="__")
            match_REF_row = 2
            while(match_REF_row<= length(match_REF[,1])){
                    ID = paste(ID, paste(match_REF[match_REF_row ,4], match_REF[match_REF_row ,1],match_REF[match_REF_row ,2],sep="__"),sep=";")
                    match_REF_row = match_REF_row +1
            }
            clip_read[10] = ID
            clip_read[23] = paste(gsub("false","", clip_read[23]),"REF_POS",sep=";")
    }

     if(clip_read[11]>0 ){
            tt = data.frame(names(genome),search2)
            tt =tt[tt[,2]>0,]
            match_REF = data.frame(unlist(vmatchPattern(ss, genome[names(genome)%in% tt[tt[,2]>0,1]])))
            ID = paste(match_REF[1,4], match_REF[1,1],match_REF[1,2],sep="__")
            match_REF_row = 2
            while(match_REF_row<= length(match_REF[,1])){
                    ID = paste(ID, paste(match_REF[match_REF_row ,4], match_REF[match_REF_row ,1],match_REF[match_REF_row ,2],sep="__"),sep=";")
                    match_REF_row = match_REF_row +1
            }
            clip_read[11] = ID
            clip_read[23] = paste(gsub("false","", clip_read[23]),"REF_NEG",sep=";")
    }

    if(clip_read[12]>0){
            match_TDNA = data.frame(unlist(vmatchPattern(s, TDNA[1])))
            ID = paste(match_TDNA[1,4], match_TDNA[1,1],match_TDNA[1,2],sep="__")
            match_TDNA_row = 2
            while(match_TDNA_row<= length(match_TDNA[,1])){
                    ID = paste(ID, paste(match_TDNA[match_TDNA_row ,4], match_TDNA[match_TDNA_row ,1],match_TDNA[match_TDNA_row ,2],sep="__"),sep=";")
                    match_TDNA_row = match_TDNA_row +1
            }
            clip_read[12] = ID
            clip_read[23] = paste(gsub("false","", clip_read[23]),"TDNA_POS",sep=";")

    }
    if(clip_read[13]>0){
            match_TDNA = data.frame(unlist(vmatchPattern(ss, TDNA[1])))
            ID = paste(match_TDNA[1,4], match_TDNA[1,1],match_TDNA[1,2],sep="__")
            match_TDNA_row = 2
            while(match_TDNA_row<= length(match_TDNA[,1])){
                    ID = paste(ID, paste(match_TDNA[match_TDNA_row ,4], match_TDNA[match_TDNA_row ,1],match_TDNA[match_TDNA_row ,2],sep="__"),sep=";")
                    match_TDNA_row = match_TDNA_row +1
            }
            clip_read[13] = ID
            clip_read[23] = paste(gsub("false","", clip_read[23]),"TDNA_NEG",sep=";")

    }
    if(clip_read[23]!="false"){
                 return(clip_read)
        }
        max_mis = min(5,  as.integer(nchar(s)/10))
        clip_read[16] =countPattern(s, TDNA[[1]],max.mismatch = max_mis)
        clip_read[17] =countPattern(ss, TDNA[[1]],max.mismatch = max_mis)
        search3 = vcountPattern(s, genome,max.mismatch = max_mis)
        search4 = vcountPattern(ss, genome,max.mismatch = max_mis)
        clip_read[14] =sum(search3)
        clip_read[15] =sum(search4)
        if(clip_read[16]>0){
                pp = sum(vcountPattern(as.character(clip_read[4]), TDNA, max.mismatch = max_mis))
                if(pp>0){
                        clip_read[23] = paste(gsub("false","", clip_read[23]),"MismatchPos_TDNAGenome",sep="")
                }
        }
        if(clip_read[17]>0){
                pp = sum(vcountPattern(as.character(reverseComplement(DNAString(as.character(clip_read[4])))), TDNA,  max.mismatch = max_mis))
                if(pp>0){
                        clip_read[23] = paste(gsub("false","", clip_read[23]),"MismatchNeg_TDNAGenome",sep="")
                }
        }
        if(clip_read[23]!="false"){
                 return(clip_read)
        }
        if(clip_read[14]>0){
            tt = data.frame(names(genome),search3)
            tt =tt[tt[,2]>0,]
            match_REF = data.frame(unlist(vmatchPattern(s, genome[names(genome)%in% tt[tt[,2]>0,1]], max.mismatch = max_mis)))
            ID = paste(match_REF[1,4], match_REF[1,1],match_REF[1,2],sep="__")
            match_REF_row = 2
            while(match_REF_row<= length(match_REF[,1])){
                    ID = paste(ID, paste(match_REF[match_REF_row ,4], match_REF[match_REF_row ,1],match_REF[match_REF_row ,2],sep="__"),sep=";")
                    match_REF_row = match_REF_row +1
            }
            clip_read[14] = ID
            clip_read[23] = paste(gsub("false","", clip_read[23]),"mismatch_REF_POS",sep=";")
    }

     if(clip_read[15]>0 ){
            tt = data.frame(names(genome),search4)
            tt =tt[tt[,2]>0,]
            match_REF = data.frame(unlist(vmatchPattern(ss, genome[names(genome)%in% tt[tt[,2]>0,1]],max.mismatch = max_mis)))
            ID = paste(match_REF[1,4], match_REF[1,1],match_REF[1,2],sep="__")
            match_REF_row = 2
            while(match_REF_row<= length(match_REF[,1])){
                    ID = paste(ID, paste(match_REF[match_REF_row ,4], match_REF[match_REF_row ,1],match_REF[match_REF_row ,2],sep="__"),sep=";")
                    match_REF_row = match_REF_row +1
            }
            clip_read[15] = ID
            clip_read[23] = paste(gsub("false","", clip_read[23]),"mismatch_REF_NEG",sep=";")
    }

    if(clip_read[16]>0){
            match_TDNA = data.frame(unlist(vmatchPattern(s, TDNA[1],max.mismatch = max_mis)))
            ID = paste(match_TDNA[1,4], match_TDNA[1,1],match_TDNA[1,2],sep="__")
            match_TDNA_row = 2
            while(match_TDNA_row<= length(match_TDNA[,1])){
                    ID = paste(ID, paste(match_TDNA[match_TDNA_row ,4], match_TDNA[match_TDNA_row ,1],match_TDNA[match_TDNA_row ,2],sep="__"),sep=";")
                    match_TDNA_row = match_TDNA_row +1
            }
            clip_read[16] = ID
            clip_read[23] = paste(gsub("false","", clip_read[23]),"mismatch_TDNA_POS",sep=";")

    }
    if(clip_read[17]>0){
            match_TDNA = data.frame(unlist(vmatchPattern(ss, TDNA[1],max.mismatch = max_mis)))
            ID = paste(match_TDNA[1,4], match_TDNA[1,1],match_TDNA[1,2],sep="__")
            match_TDNA_row = 2
            while(match_TDNA_row<= length(match_TDNA[,1])){
                    ID = paste(ID, paste(match_TDNA[match_TDNA_row ,4], match_TDNA[match_TDNA_row ,1],match_TDNA[match_TDNA_row ,2],sep="__"),sep=";")
                    match_TDNA_row = match_TDNA_row +1
            }
            clip_read[17] = ID
            clip_read[23] = paste(gsub("false","", clip_read[23]),"mismatch_TDNA_NEG",sep=";")

    }
    if(clip_read[23]!="false"){
                 return(clip_read)
        }

      if( as.integer(clip_read[8])> (anchor + 6) ){
	if(as.integer(clip_read[8])<= (anchor + insert)){
              insert = as.integer(clip_read[8]) - anchor
      	}
	
        sss =substr(s, 1, nchar(s)- insert)
        if(clip_read[7]=="type2"){
        sss = substr(s, insert +1, nchar(s))
        }
        ssss = as.character(reverseComplement(DNAString(sss)))
        clip_read[20] =countPattern(sss, TDNA[[1]])
        clip_read[21] =countPattern(ssss, TDNA[[1]])
        search5 = vcountPattern(sss, genome)
        search6 = vcountPattern(ssss, genome)
        clip_read[18] =sum(search5)
        clip_read[19] =sum(search6)
        if(clip_read[18]>0){
                    tt = data.frame(names(genome),search5)
                    tt =tt[tt[,2]>0,]
                    match_REF = data.frame(unlist(vmatchPattern(sss, genome[names(genome)%in% tt[tt[,2]>0,1]])))
                    ID = paste(match_REF[1,4], match_REF[1,1],match_REF[1,2],sep="__")
                    match_REF_row = 2
                    while(match_REF_row<= length(match_REF[,1])){
                            ID = paste(ID, paste(match_REF[match_REF_row ,4], match_REF[match_REF_row ,1],match_REF[match_REF_row ,2],sep="__"),sep=";")
                            match_REF_row = match_REF_row +1
                    }
                    clip_read[18] = ID
                    clip_read[23] = paste(gsub("false","", clip_read[23]),paste("insert_", insert,"_REF_POS",sep=""),sep=";")
            }

             if(clip_read[19]>0 ){
                    tt = data.frame(names(genome),search6)
                    tt =tt[tt[,2]>0,]
                    match_REF = data.frame(unlist(vmatchPattern(ssss, genome[names(genome)%in% tt[tt[,2]>0,1]])))
                    ID = paste(match_REF[1,4], match_REF[1,1],match_REF[1,2],sep="__")
                    match_REF_row = 2
                    while(match_REF_row<= length(match_REF[,1])){
                            ID = paste(ID, paste(match_REF[match_REF_row ,4], match_REF[match_REF_row ,1],match_REF[match_REF_row ,2],sep="__"),sep=";")
                            match_REF_row = match_REF_row +1
                    }
                    clip_read[19] = ID
                    clip_read[23] = paste(gsub("false","", clip_read[23]),paste("insert_", insert,"_REF_NEG",sep=""),sep=";")
            }

            if(clip_read[20]>0){
                    match_TDNA = data.frame(unlist(vmatchPattern(sss, TDNA[1])))
                    ID = paste(match_TDNA[1,4], match_TDNA[1,1],match_TDNA[1,2],sep="__")
                    match_TDNA_row = 2
                    while(match_TDNA_row<= length(match_TDNA[,1])){
                            ID = paste(ID, paste(match_TDNA[match_TDNA_row ,4], match_TDNA[match_TDNA_row ,1],match_TDNA[match_TDNA_row ,2],sep="__"),sep=";")
                            match_TDNA_row = match_TDNA_row +1
                    }
                    clip_read[20] = ID
                    clip_read[23] = paste(gsub("false","", clip_read[23]),paste("insert_", insert,"_TDNA_POS",sep=""),sep=";")

            }
            if(clip_read[21]>0){
                    match_TDNA = data.frame(unlist(vmatchPattern(ssss, TDNA[1])))
                    ID = paste(match_TDNA[1,4], match_TDNA[1,1],match_TDNA[1,2],sep="__")
                    match_TDNA_row = 2
                    while(match_TDNA_row<= length(match_TDNA[,1])){
                            ID = paste(ID, paste(match_TDNA[match_TDNA_row ,4], match_TDNA[match_TDNA_row ,1],match_TDNA[match_TDNA_row ,2],sep="__"),sep=";")
                            match_TDNA_row = match_TDNA_row +1
                    }
                    clip_read[21] = ID
                    clip_read[23] = paste(gsub("false","", clip_read[23]),paste("insert_", insert,"_TDNA_NEG",sep=""),sep=";")

            }
    }
   return(clip_read)

}

##get the first round of search
re_B = list()
for ( i in 1:length(Bt1)){
        tmp = B1[B1[,22] %in% names(Bt1)[i],]
        tmp_list = list(tmp,tmp[tmp[,1]==1000,])
        index = 1
        while(dim(tmp_list[[index]])[1]>=1){
                tmp1 = data.frame(tmp_list[[index]],log = rep(0,length(tmp_list[[index]][,1])))
                seq_t = DNAString(as.character(tmp1[1,9]))

                for (index2 in 1:length(tmp1[,1])){
                        tmp1[index2,24] = countPattern(tmp1[index2,9], seq_t, max.mismatch =1,with.indels=TRUE)
                }
                tmp_list[[index]] = tmp1[tmp1[,24]>0,1:23]
                tmp_list[[index+1]] = tmp1[tmp1[,24]==0,1:23]
                index = index +1

        }
        tmp_list = tmp_list[1:(index-1)]
        re_B = c(re_B,tmp_list)

}


for ( i in 1:length(Ct1)){
        tmp = C1[C1[,22] %in% names(Ct1)[i],]
        tmp_list = list(tmp,tmp[tmp[,1]==1000,])
        index = 1
        while(dim(tmp_list[[index]])[1]>=1){
                tmp1 = data.frame(tmp_list[[index]],log = rep(0,length(tmp_list[[index]][,1])))
                seq_t = DNAString(as.character(tmp1[1,9]))

                for (index2 in 1:length(tmp1[,1])){
                        tmp1[index2,24] = countPattern(tmp1[index2,9], seq_t, max.mismatch =1,with.indels=TRUE)
                }
                tmp_list[[index]] = tmp1[tmp1[,24]>0,1:23]
                tmp_list[[index+1]] = tmp1[tmp1[,24]==0,1:23]
                index = index +1

        }
        tmp_list = tmp_list[1:(index-1)]
        re_B = c(re_B,tmp_list)

}


for ( i in 1:length(re_B)){
        len = dim(re_B[[i]])[1]
        if(len<lim){
                re_B[[i]][1,23] ="SingleRead"
        }
        else{
        re_B[[i]][1,] = MATCH_REF_mosaic (genome, TDNA, as.vector(re_B[[i]][1,]))
        re_B[[i]][2:len,10:23] = re_B[[i]][1,10:23]
	 if(re_B[[i]][1,23]=="false"){

 		re_B[[i]][2,] = MATCH_REF_mosaic (genome, TDNA, as.vector(re_B[[i]][2,]))
		re_B[[i]][1,10:23] = re_B[[i]][2,10:23]
                re_B[[i]][2:len,10:23] = re_B[[i]][1,10:23]
                }
        }
}


all = B1[B1[,1]==1000,]
for ( i in 1:length(re_B)){
	len = dim(re_B[[i]])[1]
	if(len>=lim){
        all = rbind(all, re_B[[i]])
        print(dim(re_B[[i]])[1])
	}
}

all = all[!grepl("Genome|SingleRead",all[,23]),]
write.table(all, args[7], sep="\t",quote= FALSE, row.names = FALSE,col.names = TRUE)
ALL = all[!grepl("Genome|SingleRead|false",all[,23]),]
ALL_t = table(ALL[,22])
ALL_t = ALL_t[ALL_t>=lim]
ALL = ALL[ALL[,22] %in% names(ALL_t),]
##For TDNA_Ref sequences, there are four types. REF reads is alway from positive strand, TDNA can be either in positive or negative strand
##Type1, RefPos_TDNAPos, Type2, RefPos_TDNANeg
##Type3, TDNAPos_RefPos, Type4, TDNANeg_RefPos
ALL_REF = data.frame(all[all[,1]==1000,c(1:6,8)],seqREFPos = character(0), seqTDNAPos = character(0),TDNA_ID = character(0), genome_ID = character(0),match_type =character(0), ID_type= character(0), ID_all=character(0))
##For REFREF sequences, there are two types
##Type1, REFPos_REFPos, type2:REFNeg_REFPos, type3:REFPos_REFPos,type4, REFPos_REFNeg
ALL_REFREF = data.frame(all[all[,1]==1000,c(1:6,8)],REF_seqPOS1 = character(0),REF_seqPOS2 = character(0),REF_ID1 = character(0),REF_ID2 = character(0),match_type = character(0),ID_type= character(0), ID_all=character(0))
for( i in 1:length(ALL[,1])){

        if(ALL[i,7]=="type1"){
            ##seq1 is the part not match in the genome and seq2 is the part match in the genome
        seq1 = as.character(ALL[i,9])
        seq1RC =as.character(reverseComplement(DNAString(seq1)))
        seq2 = substr(as.character(ALL[i,4]), ALL[i,8] +1, nchar(as.character(ALL[i,4])))
        seq2RC =as.character(reverseComplement(DNAString(seq2)))
        ##for type 1, clipped read in the front that was not mapped in the first round
        if(grepl(";TDNA_POS", ALL[i,23])){
            id = gsub("__\\d+__","__", ALL[i,12])
            ##id, end
          re = data.frame(ALL[i,c(1:6,8)],seqREFPos = seq2, seqTDNAPos = seq1 ,TDNA_ID = id,genome_ID  = ALL[i,22] ,match_type = "perfect",ID_type= "type3:TDNAPos_RefPos", ID_all= paste(id,":REF__",ALL[i,2], "__", ALL[i,5],sep=""))
          ALL_REF= rbind(ALL_REF,re)
        }
        if(grepl(";mismatch_TDNA_POS", ALL[i,23])){
            id = gsub("__\\d+__","__", ALL[i,16])
            re = data.frame(ALL[i,c(1:6,8)],seqREFPos = seq2, seqTDNAPos = seq1 ,TDNA_ID = id,genome_ID  = ALL[i,22] ,match_type = "mismatch_TDNA",ID_type= "type3:TDNAPos_RefPos", ID_all= paste(id,":REF__",ALL[i,2], "__", ALL[i,5],sep=""))
            ALL_REF= rbind(ALL_REF,re)
        }
        if(grepl(";insert_\\w+_TDNA_POS", ALL[i,23])){
            id = gsub("__\\d+__","__", ALL[i,20])
          re = data.frame(ALL[i,c(1:6,8)],seqREFPos = seq2, seqTDNAPos = seq1 ,TDNA_ID = id,genome_ID  = ALL[i,22] ,match_type = "insertion_TDNA3end",ID_type= "type3:TDNAPos_RefPos", ID_all= paste(id,":REF__",ALL[i,2], "__", ALL[i,5],sep=""))
          ALL_REF= rbind(ALL_REF,re)
        }

        if(grepl(";TDNA_NEG", ALL[i,23])){
            id = gsub("__\\d+$|__\\d+;",";", ALL[i,13])
            id = gsub(";$","", id)
            ##id, start
            re = data.frame(ALL[i,c(1:6,8)],seqREFPos = seq2, seqTDNAPos = seq1RC ,TDNA_ID = id,genome_ID  = ALL[i,22] ,match_type = "perfect",ID_type= "type4:TDNANeg_RefPos", ID_all= paste(id,":REF__",ALL[i,2], "__", ALL[i,5],sep=""))
            ALL_REF = rbind(ALL_REF,re)
        }
        if(grepl(";mismatch_TDNA_NEG", ALL[i,23])){
            id = gsub("__\\d+$|__\\d+;",";", ALL[i,17])
            id = gsub(";$","", id)
            re = data.frame(ALL[i,c(1:6,8)],seqREFPos = seq2, seqTDNAPos = seq1RC ,TDNA_ID = id,genome_ID  = ALL[i,22] ,match_type = "mismatch_TDNA",ID_type= "type4:TDNANeg_RefPos", ID_all= paste(id,":REF__",ALL[i,2], "__", ALL[i,5],sep=""))
            ALL_REF = rbind(ALL_REF,re)
        }
        if(grepl(";insert_\\w+_TDNA_NEG", ALL[i,23])){
            id = gsub("__\\d+$|__\\d+;",";", ALL[i,21])
            id = gsub(";$","", id)
          re = data.frame(ALL[i,c(1:6,8)],seqREFPos = seq2, seqTDNAPos = seq1RC ,TDNA_ID = id,genome_ID  = ALL[i,22] ,match_type = "insertion_TDNA5end",ID_type= "type4:TDNANeg_RefPos", ID_all= paste(id,":REF__",ALL[i,2], "__", ALL[i,5],sep=""))
          ALL_REF = rbind(ALL_REF,re)
        }

                if(grepl(";REF_POS", ALL[i,23])){
                  id = gsub("__\\d+__","__", ALL[i,10])
                  re = data.frame(ALL[i,c(1:6,8)],REF_seqPOS1 = seq1, REF_seqPOS2 = seq2,REF_ID1= id, REF_ID2 = ALL[i,22] ,match_type = "perfect", ID_type= "type1:REFPos_REFPos", ID_all= paste(id,":REF__",ALL[i,2], "__", ALL[i,5],sep=""))
                  ALL_REFREF = rbind(ALL_REFREF,re)
                }
                if(grepl(";mismatch_REF_POS", ALL[i,23])){
                  id = gsub("__\\d+__","__", ALL[i,14])
                  re = data.frame(ALL[i,c(1:6,8)],REF_seqPOS1 = seq1, REF_seqPOS2 = seq2,REF_ID1= id, REF_ID2 = ALL[i,22] ,match_type = "mismatch", ID_type= "type1:REFPos_REFPos", ID_all= paste(id,":REF__",ALL[i,2], "__",  ALL[i,5],sep=""))
                  ALL_REFREF = rbind(ALL_REFREF,re)
                }
                if(grepl(";insert_\\w+_REF_POS", ALL[i,23])){
                  id = gsub("__\\d+__","__", ALL[i,18])
                  re = data.frame(ALL[i,c(1:6,8)],REF_seqPOS1 = seq1, REF_seqPOS2 = seq2,REF_ID1= id, REF_ID2 = ALL[i,22] ,match_type = "insert_REF3end", ID_type= "type1:REFPos_REFPos", ID_all= paste(id,":REF__",ALL[i,2], "__",  ALL[i,5],sep=""))
                  ALL_REFREF = rbind(ALL_REFREF,re)
                }
                if(grepl(";REF_NEG", ALL[i,23])){
                  id = gsub("__\\d+$|__\\d+;",";", ALL[i,11])
                  id = gsub(";$","", id)
                  re = data.frame(ALL[i,c(1:6,8)],REF_seqPOS1 = seq1RC, REF_seqPOS2 = seq2,REF_ID1= id, REF_ID2 = ALL[i,22] ,match_type = "perfect", ID_type= "type2:REFNeg_REFPos", ID_all= paste(id,":REF__",ALL[i,2],"__",  ALL[i,5],sep=""))
                  ALL_REFREF = rbind(ALL_REFREF,re)
                }

                if(grepl(";mismatch_REF_NEG", ALL[i,23])){
                    id = gsub("__\\d+$|__\\d+;",";", ALL[i,15])
                    id = gsub(";$","", id)
                  re = data.frame(ALL[i,c(1:6,8)],REF_seqPOS1 = seq1RC, REF_seqPOS2 = seq2,REF_ID1= id, REF_ID2 = ALL[i,22] ,match_type = "mismatch", ID_type= "type2:REFNeg_REFPos", ID_all= paste(id,":REF__",ALL[i,2], "__", ALL[i,5],sep=""))
                  ALL_REFREF = rbind(ALL_REFREF,re)
                }

                if(grepl(";insert_\\w+_REF_NEG", ALL[i,23])){
                    id = gsub("__\\d+$|__\\d+;",";", ALL[i,19])
                    id = gsub(";$","", id)
                  re = data.frame(ALL[i,c(1:6,8)],REF_seqPOS1 = seq1RC, REF_seqPOS2 = seq2,REF_ID1= id, REF_ID2 = ALL[i,22] ,match_type = "insert_REF5end", ID_type= "type2:REFNeg_REFPos", ID_all= paste(id,":REF__",ALL[i,2], "__", ALL[i,5],sep=""))
                  ALL_REFREF = rbind(ALL_REFREF,re)
                }
        }
        if(ALL[i,7]=="type2"){
            ##seq1 is the part not match in the genome and seq2 is the part match in the genome
            ##for type 2, clipped read in the end that was not mapped in the first round
                   seq1 = as.character(ALL[i,9])
                   seq1RC =as.character(reverseComplement(DNAString(seq1)))
                   seq2 = substr(as.character(ALL[i,4]), 1,nchar(as.character(ALL[i,4]))- ALL[i,8] )
                   seq2RC =as.character(reverseComplement(DNAString(seq2)))
                   if(grepl(";TDNA_POS", ALL[i,23])){
                     id = gsub("__\\d+$|__\\d+;",";", ALL[i,12])
                     id = gsub(";$","", id)
                     ##id, start
                     re = data.frame(ALL[i,c(1:6,8)],seqREFPos = seq2, seqTDNAPos = seq1 ,TDNA_ID = id,genome_ID  = ALL[i,22] ,match_type = "perfect",ID_type= "type1:RefPos_TDNAPos", ID_all= paste("REF__",ALL[i,2], "__", ALL[i,6],":",id, sep=""))

                     ALL_REF = rbind(ALL_REF,re)
                   }
                   if(grepl(";mismatch_TDNA_POS", ALL[i,23])){
                       id = gsub("__\\d+$|__\\d+;",";", ALL[i,16])
                       id = gsub(";$","", id)
                       re = data.frame(ALL[i,c(1:6,8)],seqREFPos = seq2, seqTDNAPos = seq1 ,TDNA_ID = id,genome_ID  = ALL[i,22] ,match_type = "mismatch",ID_type= "type1:RefPos_TDNAPos", ID_all= paste("REF__",ALL[i,2], "__", ALL[i,6],":",id, sep=""))
                     ALL_REF = rbind(ALL_REF,re)
                   }
                   if(grepl(";insert_\\w+_TDNA_POS", ALL[i,23])){
                       id = gsub("__\\d+$|__\\d+;",";", ALL[i,20])
                       id = gsub(";$","", id)
                       re = data.frame(ALL[i,c(1:6,8)],seqREFPos = seq2, seqTDNAPos = seq1 ,TDNA_ID = id,genome_ID  = ALL[i,22] ,match_type = "insertion_TDNA5end",ID_type= "type1:RefPos_TDNAPos", ID_all= paste("REF__",ALL[i,2], "__", ALL[i,6],":",id, sep=""))
                     ALL_REF = rbind(ALL_REF,re)
                   }

                   if(grepl(";TDNA_NEG", ALL[i,23])){
                       id = gsub("__\\d+__","__", ALL[i,13])
                       ##id, end
                       re = data.frame(ALL[i,c(1:6,8)],seqREFPos = seq2, seqTDNAPos = seq1RC ,TDNA_ID = id,genome_ID  = ALL[i,22] ,match_type = "perfect",ID_type= "type2:RefPos_TDNANeg", ID_all= paste("REF__",ALL[i,2], "__", ALL[i,6],":",id, sep=""))
                     ALL_REF = rbind(ALL_REF,re)
                   }
                   if(grepl(";mismatch_TDNA_NEG", ALL[i,23])){
                       id = gsub("__\\d+__","__", ALL[i,17])
                       re = data.frame(ALL[i,c(1:6,8)],seqREFPos = seq2, seqTDNAPos = seq1RC ,TDNA_ID = id,genome_ID  = ALL[i,22] ,match_type = "mismatch",ID_type= "type2:RefPos_TDNANeg", ID_all= paste("REF__",ALL[i,2], "__", ALL[i,6],":",id, sep=""))
                     ALL_REF = rbind(ALL_REF,re)
                   }
                   if(grepl(";insert_\\w+_TDNA_NEG", ALL[i,23])){
                       id = gsub("__\\d+__","__", ALL[i,21])
                       re = data.frame(ALL[i,c(1:6,8)],seqREFPos = seq2, seqTDNAPos = seq1RC ,TDNA_ID = id,genome_ID  = ALL[i,22] ,match_type = "insertion_TDNA3end",ID_type= "type2:RefPos_TDNANeg", ID_all= paste("REF__",ALL[i,2], "__", ALL[i,6],":",id, sep=""))
                     ALL_REF = rbind(ALL_REF,re)
                   }


                   if(grepl(";REF_POS", ALL[i,23])){
                                           id = gsub("__\\d+$|__\\d+;",";", ALL[i,10])
                                           id = gsub(";$","", id)
                                           re = data.frame(ALL[i,c(1:6,8)],REF_seqPOS1 = seq2, REF_seqPOS2 = seq1,REF_ID1 = ALL[i,22], REF_ID2 = id ,match_type = "perfect", ID_type= "type3:REFPos_REFPos", ID_all= paste("REF__",ALL[i,2],"__", ALL[i,6],":", id,sep=""))
                                       ALL_REFREF = rbind(ALL_REFREF,re)
                                       }
                                       if(grepl(";mismatch_REF_POS", ALL[i,23])){
                                           id = gsub("__\\d+$|__\\d+;",";", ALL[i,14])
                                           id = gsub(";$","", id)
                                           re = data.frame(ALL[i,c(1:6,8)],REF_seqPOS1 = seq2, REF_seqPOS2 = seq1,REF_ID1 = ALL[i,22], REF_ID2 = id ,match_type = "mismatch", ID_type= "type3:REFPos_REFPos", ID_all= paste("REF__",ALL[i,2],"__", ALL[i,6],":", id,sep=""))
                                         ALL_REFREF = rbind(ALL_REFREF,re)
                                       }
                                       if(grepl(";insert_\\w+_REF_POS", ALL[i,23])){
                                           id = gsub("__\\d+$|__\\d+;",";", ALL[i,18])
                                           id = gsub(";$","", id)
                                           re = data.frame(ALL[i,c(1:6,8)],REF_seqPOS1 = seq2, REF_seqPOS2 = seq1,REF_ID1 = ALL[i,22], REF_ID2 = id ,match_type = "insert_REF5end", ID_type= "type3:REFPos_REFPos", ID_all= paste("REF__",ALL[i,2],"__", ALL[i,6],":", id,sep=""))
                                         ALL_REFREF = rbind(ALL_REFREF,re)
                                       }

                                        if(grepl(";REF_NEG", ALL[i,23])){
                                            id = gsub("__\\d+__","__", ALL[i,11])
                                            re = data.frame(ALL[i,c(1:6,8)],REF_seqPOS1 = seq2, REF_seqPOS2 = seq1RC,REF_ID1 = ALL[i,22], REF_ID2 = id ,match_type = "perfect", ID_type= "type4:REFPos_REFNeg", ID_all= paste("REF__",ALL[i,2],"__", ALL[i,6],":", id,sep=""))
                                         ALL_REFREF = rbind(ALL_REFREF,re)
                                       }

                                       if(grepl(";mismatch_REF_NEG", ALL[i,23])){
                                           id = gsub("__\\d+__","__", ALL[i,15])
                                           re = data.frame(ALL[i,c(1:6,8)],REF_seqPOS1 = seq2, REF_seqPOS2 = seq1RC,REF_ID1 = ALL[i,22], REF_ID2 = id ,match_type = "mismatch", ID_type= "type4:REFPos_REFNeg", ID_all= paste("REF__",ALL[i,2],"__", ALL[i,6],":", id,sep=""))
                                         ALL_REFREF = rbind(ALL_REFREF,re)
                                       }

                                       if(grepl(";insert_\\w+_REF_NEG", ALL[i,23])){
                                           id = gsub("__\\d+__","__", ALL[i,19])
                                           re = data.frame(ALL[i,c(1:6,8)],REF_seqPOS1 = seq2, REF_seqPOS2 = seq1RC,REF_ID1 = ALL[i,22], REF_ID2 = id ,match_type = "insert_REF3end", ID_type= "type4:REFPos_REFNeg", ID_all= paste("REF__",ALL[i,2],"__", ALL[i,6],":", id,sep=""))
                                         ALL_REFREF = rbind(ALL_REFREF,re)
                                       }
                            }
}

ALL_REF = ALL_REF[order(ALL_REF[,10], ALL_REF[,11]),]
ALL_REF_t = table(ALL_REF[,14])
ALL_REF_t = ALL_REF_t[ALL_REF_t>=lim]
ALL_REF = ALL_REF[ALL_REF[,14] %in% names(ALL_REF_t),]
write.table(ALL_REF[,-c(2,3,5,6)], args[5], sep="\t",quote= FALSE, row.names = FALSE,col.names = TRUE)
ALL_REFREF = ALL_REFREF[order(ALL_REFREF[,10], ALL_REFREF[,11]),]
ALL_REFREF_t = table(ALL_REFREF[,14])
ALL_REFREF_t = ALL_REFREF_t[ALL_REFREF_t>=lim2]
ALL_REFREF = ALL_REFREF[ALL_REFREF[,14] %in% names(ALL_REFREF_t),]
if(dim(ALL_REFREF)[1]>0){
        write.table(ALL_REFREF[,-c(2,3,5,6)], args[6], sep="\t",quote= FALSE, row.names = FALSE,col.names = TRUE)
}


