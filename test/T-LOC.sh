#python ../T_LOC.py 
#Please provide either fastq reads or reference alignment files
#Please provide either fastq reads or TDNA alignment files
#please provide the reference genome sequence and TDNA sequences and fastq reads
#Not enough parameters!
#Program :  ../T_LOC.py
#          A python program to locate the T_DNA insertion site.
#
#Usage : ../T_LOC.py  Required Parameters;
#Usage : ../T_LOC.py  --fastq: the sequence fastq files, paired reads were separated by comma;
#Usage : ../T_LOC.py  --bamR: bam files of alignment to reference genome;
#Usage : ../T_LOC.py  --bamT: bam files of TDNA alignment;
#Usage : ../T_LOC.py  --genome: fasta format of given genome sequences;
#Usage : ../T_LOC.py  --TDNA: fasta format of given T-DNA sequences;
#Usage : ../T_LOC.py  Optional  Parameters;
#Usage : ../T_LOC.py  --control_bamR: a series of bam files from plants without T-DNA insertion, different bam files were separated by comma;
#Usage : ../T_LOC.py  --control_fastq: a series of fastq files from plants without T-DNA insertion, different fastq files were sperated by comma, paired reads were separated by colon;
#Usage : ../T_LOC.py  --control_ID: the control ID files;
#Usage : ../T_LOC.py  --genome_Bwa: bwa index of given genome sequences;
#Usage : ../T_LOC.py  --TDNA_Bwa: bwa index of given T-DNA sequences;
#Usage : ../T_LOC.py  --LB: The left border sequence position;
#Usage : ../T_LOC.py  --RB: The  right border seuqnce position;
#Usage : ../T_LOC.py  --anchor: the anchor length, with default value of 20 ;
#Usage : ../T_LOC.py  --Mosaic_length: the length of output sequences of each mosaic side;
#Usage : ../T_LOC.py  --resume: Whether to resume previous run. The default is false;
#Usage : ../T_LOC.py  -o/--output: The output directory. The default is current directory;
python ../T_LOC.py   --genome /home/shaofangli/genome/rice/KitaakeChrCM/KitaakeChrCM.fa --TDNA /home/Store/T_LOC/TDNAFasta/V46bAG.fa --fastq /home/Store/OFF-target/DNA/Bwa_ChrCM/46bAG_s2_1.fq.gz,/home/Store/OFF-target/DNA/Bwa_ChrCM/46bAG_s2_2.fq.gz --output test3 --control_ID /home/Store/T_LOC/REF/control_id.txt
python ../T_LOC.py   --genome /home/shaofangli/genome/rice/KitaakeChrCM/KitaakeChrCM.fa --TDNA /home/Store/T_LOC/TDNAFasta/V46bAG.fa --fastq /home/Store/OFF-target/DNA/Bwa_ChrCM/46bAG_s2_1.fq.gz,/home/Store/OFF-target/DNA/Bwa_ChrCM/46bAG_s2_2.fq.gz --output test1
python ../T_LOC.py   --genome /home/shaofangli/genome/rice/KitaakeChrCM/KitaakeChrCM.fa --TDNA /home/Store/T_LOC/TDNAFasta/V46bAG.fa --fastq /home/Store/OFF-target/DNA/Bwa_ChrCM/46bAG_s2_1.fq.gz,/home/Store/OFF-target/DNA/Bwa_ChrCM/46bAG_s2_2.fq.gz --output test2 --control_fastq /home/Store/OFF-target/DNA/Bwa_ChrCM/kit001_1.fq.gz,/home/Store/OFF-target/DNA/Bwa_ChrCM/kit001_2.fq.gz


python ../T_LOC.py   --genome /home/shaofangli/genome/rice/KitaakeChrCM/KitaakeChrCM.fa --TDNA /home/Store/T_LOC/TDNAFasta/V46bAG.fa --fastq /home/Store/OFF-target/DNA/Bwa_ChrCM/46bAG_s2_1.fq.gz,/home/Store/OFF-target/DNA/Bwa_ChrCM/46bAG_s2_2.fq.gz --output test4 --control_fastq /home/Store/OFF-target/DNA/Bwa_ChrCM/kit001_1.fq.gz,/home/Store/OFF-target/DNA/Bwa_ChrCM/kit001_2.fq.gz:/home/Store/OFF-target/DNA/Bwa_ChrCM/kit002_1.fq.gz,/home/Store/OFF-target/DNA/Bwa_ChrCM/kit002_2.fq.gz

python ../T_LOC.py   --genome /home/shaofangli/genome/rice/KitaakeChrCM/KitaakeChrCM.fa --TDNA /home/Store/T_LOC/TDNAFasta/V46bAG.fa --fastq /home/Store/OFF-target/DNA/Bwa_ChrCM/46bAG_s2_1.fq.gz,/home/Store/OFF-target/DNA/Bwa_ChrCM/46bAG_s2_2.fq.gz --output test5 --control_bamR /home/Store/OFF-target/DNA/kit001.sorted.bam,/home/Store/OFF-target/DNA/kit002.sorted.bam
