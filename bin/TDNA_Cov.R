args <- commandArgs(trailingOnly = TRUE)
### output the coverage of T-DNA
##args[1] = Reference Alignment file
##args[2] = TDNA Alignment file
##args[3] = Genome Seqeunce file
##args[4] = TDNA Sequence file
##args[5] = sample names
##args[6] = Output coverage pdf  files

library(ShortRead)
genome = readDNAStringSet(args[3],"fasta")
TDNA = readDNAStringSet(args[4],"fasta")
g_align= readGAlignments(args[1])
t_align= readGAlignments(args[2])

gg = data.frame(coverage(g_align))

g_fold = sum(gg[,3])/sum(width(genome))
tt = data.frame(coverage(t_align))
pdf(file = args[6])
plot(tt[,3],type = "l",main = args[5], las = 1, ylab ="TDNA coverage", xlab = "TDNA reads",ylim = c(0, min(max(tt[,3]),10*g_fold)))
abline(h = g_fold,col="red")
abline(h = 2*g_fold,col="blue")
abline(h = g_fold/2,col="blue")
dev.off()
