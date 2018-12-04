# Align CRISPR data (array) - Phase 2
# A. M. Chakrabarti
# 26th November 2016

setwd("/SAN/luscombelab/general/nobby/CRISPR")
# setwd("~/Scaffidi/metadata")
library(data.table)

# =============================================================================
# Align
# =============================================================================

Align <- function(read1, read2, stem) {
  
  # 1. Align to hg19 setting max indel size of 2000
  cmd <- paste0("bbmap.sh -Xmx23G path=/SAN/luscombelab/general/nobby/CRISPR/ref/hg19 in1=", read1, ".1.fq.gz in2=", read2, ".2.fq.gz out=", stem, ".tmp.bam maxindel=2000 trimreaddescriptions=t local=t showprogress=1000000 bhist=", stem, "_bhist.txt qhist=", stem, "_qhist.txt aqhist=", stem, "_aqhist.txt lhist=", stem, "_lhist.txt ihist=", stem, "_ihist.txt ehist=", stem, "_ehist.txt qahist=", stem, "_qahist.txt indelhist=", stem, "_indelhist.txt mhist=", stem, "_mhist.txt gchist=", stem, "_gchist.txt idhist=", stem, "_idhist.txt scafstats=", stem, "_scafstats.txt statsfile=", stem, ".metrics")
  print(cmd)
  system(cmd)

  # 2. Sort by coordinate for Picard MarkDuplicates
  cmd <- paste0("sambamba sort -t 8 -o ", stem, ".bam ", stem, ".tmp.bam")
  print(cmd)
  system(cmd)

  # 3. Mark duplicates (but do not remove at this stage, to keep aligned marked BAM file in case adjust filtering later)
  cmd <- paste0("java -jar /SAN/luscombelab/general/resources/bin/picard-2.1.0/picard.jar MarkDuplicates I=", stem, ".bam O=", stem, ".indel.bam M=", stem, ".duplicates.txt")
  print(cmd)
  system(cmd)

  # 4. Remove previous BAM files if the indel.bam has been created
  if(file.exists(paste0(stem, ".indel.bam"))) {
  	file.remove(paste0(stem, ".bam"))
  	file.remove(paste0(stem, ".bam.bai"))
  	file.remove(paste0(stem, ".tmp.bam"))
	}

  # 5. Filter properly paired reads that are not duplicates with a MAPQ of 38 or higher for downstream analysis (could also filter by pulldown BED file to reduce size of later BAM file)
  cmd <- paste0("samtools view -q 38 -f 2 -F 1024 -hu ", stem, ".indel.bam | sambamba sort -t 8 -o ", stem, ".indel.q38.bam /dev/stdin")
  print(cmd)
  system(cmd)

}

metadata <- fread("crispr_sequencing_metadata.tsv")

i <- as.numeric(Sys.getenv("SGE_TASK_ID")) # Get the job id from the environment

print(paste0("Processing: ", i))

r1 <- file.path("phase1", metadata$id[i])
r2 <- file.path("phase1", metadata$id[i])
st <- file.path("phase2", metadata$id[i])

ptm <- proc.time()
Align(read1 = r1, read2 = r2, stem = st)
proc.time() - ptm