#!/bin/bash

# Bash script to process chromatin experiments
# A. M. Chakrabarti
# 11th November 2017

#SBATCH --job-name=cvalid
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=12:00:00
#SBATCH --mem-per-cpu=4GB
#SBATCH --workdir=/home/camp/chakraa2/working/nobby/projects/crispr-indels/validation
#SBATCH --output=log/validation-%A-%a.out
#SBATCH --array=1-48

# Load modules
ml BBMap
ml SAMtools

# Get parameters
FILELIST=$1
pair1=`sed -n "${SLURM_ARRAY_TASK_ID}p" $FILELIST`
pair2=`echo $pair1 | sed 's/R1/R2/'`
out=${pair1%%_*}

echo $pair1
echo $pair2
echo $out
echo

# Run steps

# Merge
bbmerge.sh usejni=t k=150 qtrim=t trimq=10 adapter=AGATCGGAAGAGC in1=fastq/$pair1 in2=fastq/$pair2 out=merged/$out.merged.fastq.gz outu1=merged/unmerged/unmerged-$pair1 outu2=merged/unmerged/unmerged-$pair2

# Map
bbmap.sh -Xmx23G path=/home/camp/chakraa2/working/nobby/genomes/index/bbmap_hg19 in=merged/$out.merged.fastq.gz out=mapped/sam1.3/$out.bam sam=1.3 usejni=t maxindel=5000 trimreaddescriptions=t local=t pigz=t showprogress=10000 statsfile=mapped/log/$out.sam1.3.metrics
bbmap.sh -Xmx23G path=/home/camp/chakraa2/working/nobby/genomes/index/bbmap_hg19 in=merged/$out.merged.fastq.gz out=mapped/sam1.4/$out.bam sam=1.4 usejni=t maxindel=5000 trimreaddescriptions=t local=t pigz=t showprogress=10000 statsfile=mapped/log/$out.sam1.4.metrics

# Sort
samtools view -q 38 -hu mapped/sam1.3/$out.bam | sambamba sort -t 8 -o mapped/sam1.3/$out.q38.bam /dev/stdin
samtools view -q 38 -hu mapped/sam1.4/$out.bam | sambamba sort -t 8 -o mapped/sam1.4/$out.q38.bam /dev/stdin

# Metrics
input=$((`zcat fastq/$pair1 | wc -l ` / 4))
merged=$((`zcat merged/$out.merged.fastq.gz | wc -l ` / 4))
mapped=`samtools view -F 4 -c mapped/sam1.3/$out.bam`
qcmapped=`samtools view -c mapped/sam1.3/$out.q38.bam`

echo -e $out'\t'$input'\t'$merged'\t'$mapped'\t'$qcmapped >> metrics.tsv

