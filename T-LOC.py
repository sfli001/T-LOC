#!/bin/python
import getopt,copy,re,os,sys,logging,time,datetime;
options, args = getopt.getopt(sys.argv[1:], 'o:',['fastq=','bamR=','bamT=','genome=','TDNA=','genome_Bwa=','TDNA_Bwa=','anchor=','output=','resume=', 'read_min_TDNA=','read_min_REF=','control_bamR=','insert=','control_fastq=', 'control_ID=', 'Mosaic_length=','Sample_name='])
fastq = ""
bamR=""
control_bamR=""
control_fastq=""
control_ID=""
bamT=""
genome=""
genome_Bwa=""
TDNA=""
TDNA_Bwa=""
anchor = 20
output = "./"
resume = "true"
read_min_TDNA = 2
read_min_REF = 4
insert = 30
clipped_samR = ""
clipped_samT = ""
Mosaic_length = 500
Sample_name = 'ss'
for opt, arg in options:
        if opt in ('--fastq'):
                fastq = arg
        elif opt in ('--control_bamR'):
                control_bamR = arg
        elif opt in ('--control_fastq'):
                control_fastq = arg
        elif opt in ('--control_ID'):
                control_ID = arg
        elif opt in ('--bamR'):
                bamR = arg
        elif opt in ('--bamT'):
                bamT = arg
        elif opt in ('--genome'):
                genome = arg
        elif opt in ('--TDNA'):
                TDNA = arg
        elif opt in ('--genome_Bwa'):
                genome_Bwa = arg
        elif opt in ('--TDNA_Bwa'):
                TDNA_Bwa = arg
        elif opt in ('--anchor'):
                anchor = int(arg)
        elif opt in ('--insert'):
                insert = int(arg)
        elif opt in ('--read_min_TDNA'):
                read_min_TDNA = int(arg)
        elif opt in ('--read_min_REF'):
                read_min_REF = int(arg)
        elif opt in ('--Mosaic_length'):
                Mosaic_length = int(arg)
        elif opt in ('--resume'):
                resume  = arg
        elif opt in ('--Sample_name'):
                Sample_name  = arg
        elif opt in ('-o','--output'):
                output = arg
run = "true"
if(not fastq):
        if( not bamR ):
                print "Please provide either fastq reads or reference alignment files"
                run ="false"
        if( not bamT ):
                print "Please provide either fastq reads or TDNA alignment files"
                run ="false"
if(not genome or not TDNA or not output):
        print "please provide the reference genome sequence and TDNA sequences and fastq reads"
        run = "false"

if (run =="false"):
        print "Not enough parameters!"
        print "Program : ", sys.argv[0]
        print "          A python program to locate the T_DNA insertion site.\n"
        print "Usage :", sys.argv[0], " Required Parameters;"
        print "Usage :", sys.argv[0], " --fastq: the sequence fastq files, paired reads were separated by comma;"
        print "Usage :", sys.argv[0], " --bamR: sorted bam files of alignment to reference genome;"
        print "Usage :", sys.argv[0], " --bamT: sorted bam files of TDNA alignment;"
        print "Usage :", sys.argv[0], " --genome: fasta format of given genome sequences;"
        print "Usage :", sys.argv[0], " --TDNA: fasta format of given T-DNA sequences;"
        print "Usage :", sys.argv[0], " Optional  Parameters;"
        print "Usage :", sys.argv[0], " --control_bamR: a series of bam files from plants without T-DNA insertion, different bam files were separated by comma;"
        print "Usage :", sys.argv[0], " --control_fastq: a series of fastq files from plants without T-DNA insertion, different fastq files were sperated by colon, paired reads were separated by comma;"
        print "Usage :", sys.argv[0], " --control_ID: the control ID files;"
        print "Usage :", sys.argv[0], " --genome_Bwa: bwa index of given reference genome;"
        print "Usage :", sys.argv[0], " --TDNA_Bwa: bwa index of given T-DNA sequences;"
        print "Usage :", sys.argv[0], " --anchor: The minimal anchor length required for the clipped reads. The default value is 20;"
        print "Usage :", sys.argv[0], " --insert: The maximal insertion sequence length at the TDNA insertion sites. The default value is 3;"
        print "Usage :", sys.argv[0], " --read_min_REF: The minimal number of mosaic reads bwteen two reference positions. The default value is 4. It can be set to any integer larger or equal to 2. If the program run slow, please increase the number.;"
        print "Usage :", sys.argv[0], " --read_min_TDNA: The minimal number of mosaic reads required for each TDNA insertion sites as well as mosaic reads between two TDNA positions. The default value is 2. It can be set to any integer larger or equal to 2. ;"
        print "Usage :", sys.argv[0], " --Mosaic_length: the length of output sequences of each mosaic side;"
        print "Usage :", sys.argv[0], " --Sample_name: sample names will be used to label output files,defalut is ss;"
        print "Usage :", sys.argv[0], " --resume: Whether to resume previous run. The default is true;"
        print "Usage :", sys.argv[0], " -o/--output: The output directory. The default is current directory;"
        print datetime.datetime.now()
        print "Author  : Shaofang Li"
        print "Contact : shaofangli2021@hotmail.com"
        sys.exit()

def listToString(ss):
  Str = '';
  for a in ss:
    Str += a+' ';
  return Str;

if (not os.path.exists(output)):
        os.system("mkdir %s" % output)

### setting up the logging format 
logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(message)s',
                    filename=output+'/log.T_LOC' + str(datetime.datetime.now())+'.txt' ,
                    filemode='w')
##### Getting Start Time ######
if(resume =="true"):
        logging.debug('Resume the program with [%s]\n', listToString(sys.argv));
else:
        logging.debug('Start the program with [%s]\n', listToString(sys.argv));
startTime = time.time();

#full program would include fastq files, Reference genome files, 
#get path of the main program
path = os.path.abspath(os.path.dirname(__file__));
##get the path  of the python  programs
bin_path = "%s/bin" % path
REF_path ="%s/REF" % output 
TDNA_path ="%s/TDNA" % output
Result_path = "%s/Result" % output
Control_path = "%s/Control" % output
if (not os.path.exists(REF_path)):
    os.system("mkdir %s" % REF_path)
if (not os.path.exists(TDNA_path)):
    os.system("mkdir %s" % TDNA_path)
if (not os.path.exists(Result_path)):
    os.system("mkdir %s" % Result_path)
if (not os.path.exists(Control_path)):
    os.system("mkdir %s" % Control_path)
##do the bwa alignment to the reference genome
if(not bamR):
        logging.debug("#######################Run fastq alignment with fastq files#######################\n")
        if(resume=="true" and os.path.exists("%s/Bwa_genome/Done_Bwa_index_genome.txt" %output)):
            genome_Bwa = "%s/Bwa_genome/%s" %(output, re.sub(".*\\/","",genome))
        if(not genome_Bwa):
                cmd1 = "mkdir %s/Bwa_genome" % output
                cmd2 = "cp %s %s/Bwa_genome/" %(genome, output)
                cmd3 = "bwa index  %s/Bwa_genome/%s" %(output, re.sub(".*\\/","",genome))
                logging.debug(cmd1)
                logging.debug(cmd2)
                logging.debug(cmd3)
                os.system(cmd1)
                os.system(cmd2)
                os.system(cmd3)
                genome_Bwa ="%s/Bwa_genome/%s" %(output, re.sub(".*\\/","",genome))
                logging.debug("Generate the bwa index for the reference genome\n")
                fw = open("%s/Bwa_genome/Done_Bwa_index_genome.txt" %output,"w")
                fw.close()
                logging.debug("Output the file for finishing the indexing the genome\n")
        if(resume=="true" and os.path.exists("%s/Reference_%s.sorted.bam.bai" %(REF_path, Sample_name))):
            bamR = "%s/Reference_%s.sorted.bam" % (REF_path,Sample_name) 
        if(not bamR):
            cmd_bwa1 = "bwa mem %s %s > %s/Reference_%s.sam 2> %s/Reference_bwa_%s.txt" % (genome_Bwa, re.sub(","," ", fastq), REF_path,Sample_name,REF_path,Sample_name)
            bamR = "%s/Reference_%s.sorted.bam" % (REF_path,Sample_name)
            logging.debug(cmd_bwa1)
            os.system(cmd_bwa1)
            logging.debug("using bwa mem  to do the alignment of reference genome\n")
            cmd1 = "samtools view -bhST %s %s/Reference_%s.sam -o %s/Reference_%s.bam" % (genome, REF_path,Sample_name, REF_path,Sample_name)
            cmd2 = "samtools sort -o %s/Reference_%s.sorted.bam %s/Reference_%s.bam"% (REF_path,Sample_name, REF_path,Sample_name)
            cmd3 = "samtools index %s" % (bamR)
            logging.debug(cmd1)
            logging.debug(cmd2)
            logging.debug(cmd3)
            os.system(cmd1)
            os.system(cmd2)
            os.system(cmd3)
            logging.debug("using samtools to generate sorted.bam files \n")
clipped_samR = "%s/Reference_clipped_%s.sam" % (REF_path,Sample_name)
cmd = "samtools view " + bamR + "|awk '$6~/[S]/' - > " + clipped_samR
logging.debug(cmd)
os.system(cmd)
logging.debug("output reference clipped reads\n")
##do the bwa alignment to the TDNA
if(not bamT ):
    logging.debug("#######################Run fastq alignment with fastq files#######################\n")
    fq = fastq.split(",")
    if(not TDNA_Bwa):
            cmd1 = "mkdir %s/Bwa_TDNA" % output
            cmd2 = "cp %s %s/Bwa_TDNA/" %(TDNA, output)
            cmd3 = "bwa index  %s/Bwa_TDNA/%s" %(output, re.sub(".*\\/","",TDNA))
            logging.debug(cmd1)
            logging.debug(cmd2)
            logging.debug(cmd3)
            os.system(cmd1)
            os.system(cmd2)
            os.system(cmd3)
            TDNA_Bwa ="%s/Bwa_TDNA/%s" %(output, re.sub(".*\\/","",TDNA))
            logging.debug("Generate the bwa index for the T-DNA sequence\n")
    if(resume=="true" and os.path.exists("%s/TDNA_%s.sorted.bam.bai" %(TDNA_path, Sample_name))):
            bamT = "%s/TDNA_%s.sorted.bam" % (TDNA_path,Sample_name)
    if(not bamT):
        cmd_bwa = "bwa mem %s %s > %s/TDNA_%s.sam 2> %s/TDNA_bwa_%s.txt" % (TDNA_Bwa,re.sub(","," ", fastq), TDNA_path,Sample_name,TDNA_path,Sample_name)
        bamT = "%s/TDNA_%s.sorted.bam" % (TDNA_path,Sample_name)
        logging.debug(cmd_bwa)
        os.system(cmd_bwa)
        logging.debug("using bwa mem files to do the alignment of T-DNA\n")
        cmd1 = "samtools view -hS -F4   %s/TDNA_%s.sam |samtools view -bhST  %s -   -o %s/TDNA_%s.bam" %(TDNA_path,Sample_name,TDNA, TDNA_path,Sample_name)
        cmd2 = "samtools sort -o %s %s/TDNA_%s.bam"% (bamT, TDNA_path,Sample_name)
        cmd3 = "samtools index %s" % (bamT)
        logging.debug(cmd1)
        logging.debug(cmd2)
        logging.debug(cmd3)
        os.system(cmd1)
        os.system(cmd2)
        os.system(cmd3)
        logging.debug("using samtools to generate sorted.bam files \n")
clipped_samT = "%s/TDNA_clipped_%s.sam" % (TDNA_path,Sample_name)
cmd = "samtools view " + bamT + "|awk '$6~/[S]/' - > " + clipped_samT
logging.debug(cmd)
os.system(cmd)
logging.debug("output TDNA clipped reads\n")
if(not control_ID):
    logging.debug("#######################Generated the control ID  files or skip the calculation on reference clipped files #######################\n")
if(not control_ID and control_bamR):
    C_bamR = control_bamR.split(",")
    clipped_controlR = ""
    for i in range(0, len(C_bamR)):
        cmd = "samtools view " + C_bamR[i] + "|awk '$6~/[S]/' - > " + Control_path + "/Control_clipped_" + str(i+1)+ ".sam"
        logging.debug(cmd)
        os.system(cmd)
        clipped_controlR = clipped_controlR + "," + Control_path + "/Control_clipped_" + str(i+1) + ".sam"
    clipped_controlR = re.sub(",","", clipped_controlR)
    control_ID = Control_path + "/Control_id.txt"
    cmd = "Rscript %s/Clipped_control_REF.R %s %s %s" % (bin_path, clipped_controlR, control_ID, anchor )
    os.system(cmd)
    logging.debug("output ID from control\n")
if( not control_ID and control_fastq):
    C_fastq = control_fastq.split(":")
    clipped_controlR = ""
    for i in range(0, len(C_fastq)):
        cmd_bwa1 = "bwa mem  %s %s > %s/Reference_control_%s.sam 2> %s/Reference_control_%s_bwa.txt" % (genome_Bwa, re.sub(","," ", C_fastq[i]), Control_path,i+1,Control_path,i+1)
        cmd = "samtools view " + Control_path + "/Reference_control_"+  str(i +1) + ".sam" + "|awk '$6~/[S]/' - > " + Control_path + "/Control_clipped_" + str(i+1)+ ".sam"
        logging.debug(cmd_bwa1)
        os.system(cmd_bwa1)
        logging.debug(cmd)
        os.system(cmd)
        clipped_controlR = clipped_controlR + "," + Control_path + "/Control_clipped_" + str(i+1)+ ".sam"
        logging.debug(clipped_controlR)
    clipped_controlR = re.sub(",","", clipped_controlR)
    control_ID = Control_path + "/Control_id.txt"
    cmd = "Rscript %s/Clipped_control_REF.R %s %s %s" % (bin_path, clipped_controlR, control_ID, anchor )
    os.system(cmd)
    logging.debug("output ID from control\n")

if (control_ID):
    cmd = "Rscript %s/Clipped_REF.R %s %s %s %s %s/%s_mosaic_REF_TDNA.txt  %s/%s_mosaic_REF_REF.txt  %s/%s_mosaic_all.txt %s %s %s %s" % (bin_path,genome, TDNA, clipped_samR, control_ID, REF_path, Sample_name, REF_path,Sample_name,REF_path,Sample_name, anchor, read_min_TDNA, read_min_REF, insert)
    logging.debug(cmd)
    os.system(cmd)
    logging.debug("output mosaic read from clipped reads in reference\n")

cmd = "Rscript %s/TDNA_Cov.R %s %s %s %s %s %s/%s_TDNA_Cov.pdf " % (bin_path,bamR, bamT, genome, TDNA,Sample_name,Result_path, Sample_name)
logging.debug(cmd)
os.system(cmd)
logging.debug("Produce the TDNA coverage\n")


cmd = "Rscript %s/Clipped_TDNA.R %s %s %s  %s/%s_mosaic_REF_TDNA.txt  %s/%s_mosaic_TDNA_TDNA.txt  %s/%s_mosaic_all.txt %s %s %s" % (bin_path,genome, TDNA, clipped_samT,TDNA_path, Sample_name,TDNA_path,Sample_name,TDNA_path,Sample_name, anchor, read_min_TDNA,  insert)
logging.debug(cmd)
os.system(cmd)
logging.debug("output mosaic read from clipped reads in TDNA\n")

if(os.path.exists("%s/%s_mosaic_REF_TDNA.txt" % (REF_path, Sample_name) )):
    cmd = "Rscript %s/Mosaic_output.R  %s %s %s/%s_mosaic_REF_TDNA.txt  %s/Mosaic_%s_FromREF.pdf %s/Mosaic_%s_FromREF.fa %s %s" % (bin_path,genome, TDNA, REF_path,Sample_name,Result_path,Sample_name, Result_path,Sample_name,Sample_name,  Mosaic_length)
    logging.debug(cmd)
    os.system(cmd)
    logging.debug("output Mosaic diagram and fasta files from clipped REF reads\n")
    cmd = "Rscript %s/Coverage_bam.R  %s/%s_mosaic_REF_TDNA.txt  %s %s/Cov_%s_FromREF.pdf %s %s %s" % (bin_path,REF_path,Sample_name,bamR, Result_path,Sample_name, REF_path, Sample_name,bin_path)
    logging.debug(cmd)
    os.system(cmd)
    logging.debug("output coverage for mosaic sites get from clipped REF reads\n")    

cmd = "Rscript %s/Mosaic_output.R  %s %s %s/%s_mosaic_REF_TDNA.txt  %s/Mosaic_%s_FromTDNA.pdf %s/Mosaic_%s_FromTDNA.fa %s %s" % (bin_path,genome, TDNA, TDNA_path,Sample_name,Result_path,Sample_name, Result_path,Sample_name, Sample_name,Mosaic_length)
logging.debug(cmd)
os.system(cmd)
logging.debug("output Mosaic diagram and fasta files from clipped TDNA reads\n")
cmd = "Rscript %s/Coverage_bam.R  %s/%s_mosaic_REF_TDNA.txt  %s %s/Cov_%s_FromTDNA.pdf %s %s %s" % (bin_path,TDNA_path,Sample_name,bamR, Result_path,Sample_name, TDNA_path,Sample_name,bin_path)
logging.debug(cmd)
os.system(cmd)
logging.debug("output coverage for mosaic sites get from clipped TDNA reads\n")   


logging.debug("Finish runing T_LOC\n")
