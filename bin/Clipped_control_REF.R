args <- commandArgs(trailingOnly = TRUE)
##using one or multiple control reference clipped sample files to output the clipped ID
library(ShortRead)
##agrs[1] = one of multipe control clipped sam files seperated by comma
##args[2] = output control ID file
##args[3] = the anchor length
files = unlist(strsplit(args[1],","))
print (files)
ID = character()
for (i in 1:length(files)){
	a=read.delim(file=files[i],header=F,col.names=1:20)
	a1 = a[!duplicated(a[,10]),]
	a1 = a1[!grepl("I|D",a1[,6]),]
	a1 = a1[!grepl("A|T|C|G",a1[,13]),]
	type = rep("type1",length(a1[,1]))
	a1 = data.frame(a1[,c(1,3,6,10)],start = a1[,4],end = a1[,4]-1 + as.numeric( gsub("\\d+S|M","",as.character(a1[,6]))))
	colnames(a1)[1:4] = c("read_id","Chr","CIGAR","read_seq")
	a1 = data.frame(a1,type)
	a1[grepl("^\\d+S",a1[,3]) & grepl("\\d+M$",a1[,3]),7] = "type1"
	a1[grepl("^\\d+M",a1[,3]) & grepl("\\d+S$",a1[,3]),7] ="type2"
	#a1[grepl("^\\d+S",a1[,3]) & grepl("\\d+S$",a1[,3]),7] ="type3"
	#a1 = a1[a1[,7]!="type3",]
	C1 = a1[a1[,7]=="type2",]
	C1 = data.frame(C1,anchor = as.numeric( gsub("\\d+M|S","",as.character(C1[,3]))),ID = paste(C1[,2], C1[,6],C1[,7],sep="_"))
	C1 = C1[C1[,8]> anchor & C1[,8]< nchar(C1[1,4]) -anchor,]
	
	B1 =  a1[a1[,7]=="type1",]
	B1 = data.frame(B1,anchor = as.numeric( gsub("S.*","",as.character(B1[,3]))),ID = paste(B1[,2], B1[,5],B1[,7],sep="_"))
	B1 = B1[B1[,8]> anchor & B1[,8]< nchar(B1[1,4]) -anchor,]
	ID = unique(c(ID, as.character(B1[,11]), as.character(C1[,11])))
}

write.table(ID, args[2], quote = FALSE, row.names = FALSE, col.names = FALSE)

