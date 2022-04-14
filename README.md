## T-LOC: A T-DNA insertion sites Locator
A program using whole-genome sequencing data to localize T-DNA insertion site in plant genome and output the diagram of T-DNA insertion site, mosaic sequence reads, flanking sequences of T-DNA insertion sites as well as the coverage of  TDNA insertion sites flanking region.
The program were based on the mosaic reads aligned to both reference and TDNA. Meanwhile, the program also output mosaic reads aligned to two difference position at the reference genomes as well as two difference positions at the TDNA.

Requirements
------------
1. Install [Python 2.7.x](https://www.python.org/downloads)
2. Install [bwa](http://bio-bwa.sourceforge.net/)
3. Install [R](https://www.r-project.org)
4. Install [ShortRead](http://www.bioconductor.org/packages/release/bioc/html/ShortRead.html)

Installation
------------
The source code can be directly called from Python.

Example Usage
--------------------------------
Run T-LOC with: sequenced fastq files with T-DNA insertion, genome fasta files,  TDNA fasta files, 

        python T-LOC.py  --fastq samples_1.fq.gz,samples_2.fq.gz --genome Reference.fa  --TDNA TDNA.fa --ouput outputfolder

Run T-LOC with: bam files alignment to both reference and TDNA, genome fasta files,  TDNA fasta files,

        python T-LOC.py  --bamR samples_Ref.bam --bamT samples_TDNA.bam --genome Reference.fa  --TDNA TDNA.fa --ouput outputfolder

Run T-LOC with: sequenced fastq files with T-DNA insertion,  sequenced fastq files without T-DNA insertion as control, genome fasta files,  TDNA fasta files

        python T-LOC.py  --fastq samples_1.fq.gz,samples_2.fq.gz --genome Reference.fa  --TDNA TDNA.fa --control_fastq control_s1_1.fq.gz,control_s1_2.fq.gz:control_s2_1.fq.gz,control_s2_2.fq.gz --ouput outputfolder

Run T-LOC with: sequenced fastq files with T-DNA insertion,  bam files alignment to reference without T-DNA insertion as control, genome fasta files,  TDNA fasta files

        python T-LOC.py  --fastq samples_1.fq.gz,samples_2.fq.gz --genome Reference.fa  --TDNA TDNA.fa --control_bamR control_s1.sam,control_s2.sam --ouput outputfolder


Required Parameters
------------
	--genome:
		The fasta files of the reference genome
	--TDNA:
		The fasta files of the TDNA  genome
	--fastq: 
		s1_1.fq.gz[,s1_2.fq.gz]. The raw sequencing reads in fastq format from plant with TDNA insertions that are used to align to the reference genome and TDNA genomes. It can be omitted if both bamR and bamT were provided.
	--bamR:
		Sorted bam file alignment to the reference genome from sequencing reads on plant with TDNA insertions. It can be omitted if fastq files  were provided.
	--bamT:
                Sorted Bam file alignment to the TDNA genome from sequencing reads on plant with TDNA insertions. It can be omitted if fastq files  were provided.

Optional Parameters
------------	
	--control_fastq:
		c1_1.fq.gz[,c1_2.fq.gz][:c2_1.fq.gz[,c2_2.fq.gz],...]. The raw sequencing reads in fastq format from plant without TDNA insertions that are used to align to the reference genome to produce control clipped read ID
	--control_bamR:
		c1.bam[,c2.bam,...]. Bam/sam file alignment to the reference genome from sequencing reads on plant without TDNA insertions that are used to produce control clipped read ID
	--control_ID:
		A txt file that are clipped read ID from plant without TDNA insertions. If control_fastq, control_bamR and control_ID were all missed, TDNA insertion sites were localized through TDNA clipped reads only.
	--anchor:
		The minimal anchor length required for the clipped reads. The default value is 20.
	--insert:
		The maximal insertion sequence length at the TDNA insertion sites. The default value is 30.
	--genome_Bwa: 
		bwa index of given reference genome
	--TDNA_Bwa:    
                bwa index of given TDNA sequence
	--read_min_TDNA:
		The minimal number of mosaic reads required for each TDNA insertion sites as well as mosaic reads between two TDNA positions. The default value is 2. It can be set to any integer larger or equal to 2. 
	--read_min_REF:
		The minimal number of mosaic reads bwteen two reference positions. The default value is 4. It can be set to any integer larger or equal to 2. If the program run slow, please increase the number.
	--Sample_name:
		sample names will be used to label output files,defalut is ss
	--resume:
		 Whether to resume previous run. The default is true
	--Mosaic_length: 
		the length of output sequences of each mosaic side
	--ouptut, -o
                The output folder. The default is current directory
 


Output list
--------------------------------
Result

	The folder contains the final output three types of files: 
	MosaicFromTDNA.pdf
		It contains the TDNA insertion diagram and mosaic reads at both sides for each of the TDNA insertion sites calucated from TDNA clipped alignment reads.	i
	MosaicFromREF.pdf
       		It contains the TDNA insertion diagram and mosaic reads at both sides for each of the TDNA insertion sites calucated from reference clipped alignment reads.
	MosaicFromTDNA.fasta
		It contains the flanking reference sequences as well as the TDNA sequence for each of the TDNA insertion sites calucated from TDNA clipped alignment reads.
	MosaicFromREF.fasta
       		It contains the flanking reference sequences as well as the TDNA sequence for each of the TDNA insertion sites calucated from reference clipped alignment reads. 
	CovREF.pdf
		The coverage of reference genome for each of these flanking regions of T-DNA insertion sites.

log.T-LOC

	Log file for running T-LOC pipeline

Bwa_genome

	The folder contains the bwa index files of reference genome

TDNA_genome

	The folder contains the bwa index files of TDNA sequences

REF

	The folder contains the middle files using reference clipped reads

TDNA

	The folder contains the middle files using TDNA clipped reads

Control

	The folder contains the middle files to generate control clip ID from plants without TDNA insertions



Copyright and License Information
---------------------------------
Copyright (C) 2021, Shaofang Li

Authors: Shaofang Li

A patent has been applied for this program. For commercial users, please contact me at shaofangli2021@hotmail.com.
This program is free for academic users. 


