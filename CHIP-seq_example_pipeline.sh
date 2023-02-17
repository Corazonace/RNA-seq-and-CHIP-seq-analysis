#!/bin/bash
#SBATCH --job-name=DHT_Veh_1
#SBATCH --time=15:00:00
#SBATCH --partition=largemem
#SBATCH --mem=60GB         
#SBATCH --account=zanglab

mkdir /scratch/jw4xtu/bowtie_CHIP
mkdir /scratch/jw4xtu/bowtie_CHIP/sam_file
mkdir /scratch/jw4xtu/bowtie_CHIP/bam_file
mkdir /scratch/jw4xtu/bowtie_CHIP/peak_calling_result
cd /scratch/jw4xtu/bowtie_CHIP

###creating index (downloaded from Bowtie2 official website)###
wget https://genome-idx.s3.amazonaws.com/bt/GRCh38_noalt_as.zip
unzip GRCh38_noalt_as.zip

###alignment###
cd GRCh38_noalt_as
bowtie2 -p 5 -x GRCh38_noalt_as -U reads/LNCaP_DHT_AR_1.fastq.gz -S sam_file/DHT_AR_1.sam
bowtie2 -p 5 -x GRCh38_noalt_as -U reads/LNCaP_Veh_AR_1.fastq.gz -S sam_file/Veh_AR_1.sam

###converting sam file to bam file and sorting and indexing###
cd /scratch/jw4xtu/bowtie_CHIP
module load gcc samtools
samtools view -S -b sam_file/DHT_AR_1.sam > bam_file/DHT_AR_1.bam
samtools sort bam_file/DHT_AR_1.bam -o bam_file/DHT_AR_1.sort.bam
samtools index bam_file/DHT_AR_1.sort.bam
samtools view -S -b sam_file/Veh_AR_1.sam > bam_file/Veh_AR_1.bam
samtools sort bam_file/Veh_AR_1.bam -o bam_file/Veh_AR_1.sort.bam
samtools index bam_file/Veh_AR_1.sort.bam

###quality check again###
mkdir /scratch/jw4xtu/bowtie_CHIP/bam_file/qc_result
module load fastqc
fastqc -t 6 -o qc_result bam_file/DHT_AR_1.sort.bam
fastqc -t 6 -o qc_result bam_file/Veh_AR_1.sort.bam

###removing PCR duplication###
cd /scratch/jw4xtu/bowtie_CHIP/bam_file
module load picard
java -jar $EBROOTPICARD/picard.jar MarkDuplicates I=DHT_AR_1.sort.bam O=DHT_1_rmdup.bam M=rmdup_metrics_D1.txt REMOVE_DUPLICATES=TRUE ASSUME_SORTED=TRUE
java -jar $EBROOTPICARD/picard.jar MarkDuplicates I=Veh_AR_1.sort.bam O=Veh_1_rmdup.bam M=rmdup_metrics_V1.txt REMOVE_DUPLICATES=TRUE ASSUME_SORTED=TRUE

####peak_calling with q-value 0.01###
cd /scratch/jw4xtu/bowtie_CHIP
module load macs2
macs2 callpeak -t DHT_1_rmdup.bam -c Veh_1_rmdup.bam -f BAM -g hs -B -q 0.01 -o peak_calling_result