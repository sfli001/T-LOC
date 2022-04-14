#!/bin/python
import getopt,copy,re,os,sys,logging,time,datetime;
import pysam,os.path
options, args = getopt.getopt(sys.argv[1:], 'o:',['bam=','Region=','output='])
bam ='';
Region = '';
output ='';
for opt, arg in options:
	if opt in ('-o','--output'):
		output = arg
	elif opt in ('--bam'):
                bam = arg
	elif opt in ('--Region'):
                Region = arg
	elif opt in ('--Start'):
		Start = int(arg)
        elif opt in ('--End'):
                End = int(arg)
if (not output or not bam or not Region):
	print "Not enough parameters!"
	print "Program : ", sys.argv[0]
	print "          A python program to output the coverage given sorted.bam files."
	print "Usage :", sys.argv[0], " --bam: the sorted bam file;"
        print "Usage :", sys.argv[0], " --Region: the region need to calculated the coverage, in the format of Chr:start-end;"
	print "Usage :", sys.argv[0], " --output: col;"
	print datetime.datetime.now()
	print "Author  : Shaofang Li"
	print "Contact : shaofangli2021@hotmail.com"
	sys.exit()
##START,END, 1-based coordiate
Chr = re.sub(":.*","",Region)
Start = int(re.sub(".*:|-\\d+","",Region))
End =  int(re.sub(".*-","",Region))
fw = open(output,"w")
cov = dict()
for p in range(Start,End+1):
    cov[p]=0
fr = pysam.Samfile(bam, "rb")
for iter in fr:
    pos = iter.get_reference_positions()
    for p in pos:
        if(p+1 in cov):
            cov[p+1] = cov[p+1] +1 

for p in range(Start,End+1):
    fw.write("%s\t%s\n" % (p, cov[p]))

fr.close()
fw.close()
