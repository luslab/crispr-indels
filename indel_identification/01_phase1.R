# Align CRISPR data (array) - Phase 1
# A. M. Chakrabarti
# 26th November 2016

setwd("/SAN/luscombelab/general/nobby/CRISPR")
# setwd("~/Scaffidi/metadata")
library(data.table)

# =============================================================================
# Align
# =============================================================================

Align <- function(read1, read2, stem) {
  
  cmd <- paste0("bbmap.sh -Xmx23G path=/SAN/luscombelab/general/nobby/CRISPR/ref/hg19 in1=", read1, " in2=", read2, " out=", stem, ".bam maxindel=0 strictmaxindel=t trimreaddescriptions=t local=t showprogress=1000000 bhist=", stem, "_bhist.txt qhist=", stem, "_qhist.txt aqhist=", stem, "_aqhist.txt lhist=", stem, "_lhist.txt ihist=", stem, "_ihist.txt ehist=", stem, "_ehist.txt qahist=", stem, "_qahist.txt indelhist=", stem, "_indelhist.txt mhist=", stem, "_mhist.txt gchist=", stem, "_gchist.txt idhist=", stem, "_idhist.txt scafstats=", stem, "_scafstats.txt statsfile=", stem, ".metrics")
  print(cmd)
  system(cmd)

  cmd <- paste0("samtools view -hu -F 2 ", stem, ".bam | sambamba sort -t 8 -N -o ", stem, ".notproper.bam /dev/stdin")
  print(cmd)
  system(cmd)

  cmd <- paste0("bedtools bamtofastq -i ", stem, ".notproper.bam -fq ", stem, ".1.fq -fq2 ", stem, ".2.fq")
  print(cmd)
  system(cmd)

  cmd <- paste0("pigz ", stem, ".1.fq ", stem, ".2.fq")
  print(cmd)
  system(cmd)

}

metadata <- fread("crispr_sequencing_metadata.tsv")

i <- as.numeric(Sys.getenv("SGE_TASK_ID")) # Get the job id from the environment

print(paste0("Preprocessing: ", i))

r1 <- file.path("raw", metadata$read1[i])
r2 <- file.path("raw", metadata$read2[i])
st <- file.path("phase1", metadata$id[i])

ptm <- proc.time()
Align(read1 = r1, read2 = r2, stem = st)
proc.time() - ptm