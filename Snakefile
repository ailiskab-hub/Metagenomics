rule all:
  input:
	  "group1_1_R1_paired.fastq.gz",
	  "group1_2_R1_paired.fastq.gz",
	  "group1_3_R1_paired.fastq.gz",
	  "group1_4_R1_paired.fastq.gz",
	  "group1_5_R1_paired.fastq.gz",
	  "group1_6_R1_paired.fastq.gz",
	  "group2_1_R1_paired.fastq.gz",
	  "group2_2_R1_paired.fastq.gz",
	  "group2_3_R1_paired.fastq.gz",
	  "group2_4_R1_paired.fastq.gz",
	  "group2_5_R1_paired.fastq.gz",
	  "group2_6_R1_paired.fastq.gz"
	  

rule trim:
	input:
		forv="{seq}_R1_001.fastq.gz",
		rev="{seq}_R2_001.fastq.gz"
	output:
		"{seq}_R1_paired.fastq.gz",
		"{seq}_R1_unpaired.fastq.gz",
		"{seq}_R2_paired.fastq.gz",
		"{seq}_R2_unpaired.fastq.gz",
	shell:
		"trimmomatic PE {input.forv} {input.rev} {output} LEADING:22 TRAILING:22 SLIDINGWINDOW:4:24 MINLEN:36"
