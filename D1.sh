#!/bin/bash
#SBATCH --job-name=RNA_seq_D1
#SBATCH --time=15:00:00
#SBATCH --partition=largemem
#SBATCH --account=zanglab

mkdir /scratch/jw4xtu/hisat_RNA
mkdir /scratch/jw4xtu/hisat_RNA/sam_file
mkdir /scratch/jw4xtu/hisat_RNA/bam_file
mkdir /scratch/jw4xtu/hisat_RNA/gft_file
cd /scratch/jw4xtu/hisat_RNA

###quality control with FASTQC###
module load fastqc
fastqc -t 6 -o qc_result reads/LNCaP_DHT_RNA_rep_1_1.fastq.gz
fastqc -t 6 -o qc_result reads/LNCaP_DHT_RNA_rep_1_2.fastq.gz

###trimming and quality check again###
mkdir /scratch/jw4xtu/hisat_RNA/reads/trimmed
module load trimmomatic/0.39
java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE -threads 12 reads/LNCaP_DHT_RNA_rep1_1.fastq.gz reads/LNCaP_DHT_RNA_rep1_2.fastq.gz reads/trimmed/DHT_1_1_paired.fastq.gz reads/trimmed/DHT_1_1_unpaired.fastq.gz reads/trimmed/DHT_1_2_paired.fastq.gz reads/trimmed/DHT_1_2_unpaired.fastq.gz LEADING:3 TRAILING:3 HEADCROP:13
mkdir /scratch/jw4xtu/hisat_RNA/reads/trimmed/qc_result
module load fastqc
fastqc -t 6 -o reads/trimmed/qc_result reads/trimmed/DHT_1_1_paired.fastq.gz
fastqc -t 6 -o reads/trimmed/qc_result reads/trimmed/DHT_1_2_paired.fastq.gz

###creating index (downloading from hisat2 official website)###
wget https://genome-idx.s3.amazonaws.com/hisat/hg38_genome.tar.gz
tar -zxvf hg38_genome.tar.gz 

###alignment with hisat2###
module load gcc/9.2.0
module load hisat2/2.1.0
cd hg38
hisat2 -x genome -1 /scratch/jw4xtu/hisat_RNA/reads/trimmed/DHT_1_1_paired.fastq.gz -2 /scratch/jw4xtu/hisat_RNA/reads/trimmed/DHT_1_2_paired.fastq.gz -S /scratch/jw4xtu/hisat_RNA/sam_file/DHT_1.sam

###converting sam file to bam file and sorting and indexing with samtools###
cd /scratch/jw4xtu/hisat_RNA
module load samtools/1.10
samtools view -bS sam_file/DHT_1.sam > bam_file/DHT_1.bam
samtools sort bam_file/DHT_1.bam -o bam_file/DHT_1.sorted.bam
samtools index bam_file/DHT_1.sorted.bam

###quality check again with fastqc###
module load fastqc
mkdir /scratch/jw4xtu/hisat_RNA/hg38/bam_file/qc_result
fastqc -t 6 -o bam_file/qc_result bam_file/DHT_1.sorted.bam

###summarizing read counts with stringtie and featurecounts###
module load gcc stringtie
stringtie bam_file/DHT_1.sorted.bam -o gtf_file/DHT1.gtf -p 5 -G /scratch/jw4xtu/hiat_RNA/gencode.v42.basic.annotation.gtf




