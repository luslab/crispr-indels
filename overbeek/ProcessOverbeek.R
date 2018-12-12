# Script to process Overbeek spacer files
# A. M. Chakrabarti
# Last updated: 14th September 2018

library(data.table)
library(GenomicAlignments)
library(stringr)
library(BSgenome.Hsapiens.UCSC.hg19)

# ml R
# ml SAMtools
# sbatch -N 1 --array 1-10 -J ovb -o logs/%A_%a.log --wrap="Rscript --vanilla ProcessOverbeek.R"
# sbatch -N 1 --array 1-4607%300 -J ovb -o logs/%A_%a.log --wrap="Rscript --vanilla ProcessOverbeek.R"

# setwd("~/Scaffidi/overbeek")

ptm <- proc.time()

# Get array id
i <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

# Merge spacer data with SRA data
# sra <- fread("/Volumes/lab-luscomben/working/nobby/projects/crispr-indels/overbeek/SraRunTable.txt")
sra <- fread("SraRunTable.txt")
sra <- sra[!grepl("MTSS", Sample_Name)] # select spacer experiments
sra <- sra[!grepl("dna_pkcs", Sample_Name)] # remove DNA repair experiment
sra <- sra[!grepl("Lenti", Sample_Name)] # Remove lentiviral
sra <- sra[!grepl("HSC", Sample_Name)] # Remove HSC

sra[treatment != "WT", replicate := sapply(strsplit(Library_Name, "_"), "[[", 6)]
sra[treatment == "WT", replicate := "R1"]

sra[, `:=` (genome = sapply(strsplit(chromosome_loc, "_"), "[[", 1),
            coord = sapply(strsplit(chromosome_loc, "_"), "[[", 2))]
# sra <- merge(sra, spacer, by.x = "coord", by.y = "Genomic.location.of.spacer.(hg19)", all.x = TRUE)
# stopifnot(all(!is.na(sra$Spacer)))
# sra[is.na(Spacer)] # But the coords are correct in the SRA annotation...

sra[, `:=` (chr = sapply(strsplit(coord, ":"), "[[", 1),
            start = as.numeric(sapply(strsplit(sapply(strsplit(coord, ":"), "[[", 2), "-"), "[[", 1)),
            end = as.numeric(sapply(strsplit(sapply(strsplit(coord, ":"), "[[", 2), "-"), "[[", 2)))]
sra[, sgrna := getSeq(Hsapiens, names = chr, start = start, end = end, as.character = TRUE)]
stopifnot(all(nchar(sra$sgrna) == 23))
sra[, `:=` (l_pam = str_sub(sgrna, 1, 2),
               r_pam = str_sub(sgrna, 22, 23))]
sra[r_pam == "GG" & l_pam != "CC", strand := "+"]
sra[l_pam == "CC" & r_pam != "GG", strand := "-"]

# 3 could go either way - checked and these are +ve strand
# unique(sra[is.na(strand), .(chr, start, end, coord, sgrna)])
sra[is.na(strand), strand := "+"]

stopifnot(all(!is.na(sra$strand)))
st <- data.table(openxlsx::read.xlsx("1-s2.0-S1097276516303252-mmc2.xlsx", startRow = 4)) # need to check before flipping sequence for it to match with supplementary table
stopifnot(all(sra$sgrna %in% st$Spacer.sequence))

# Now flip sequence for negative strand
# sra[, sgrna := getSeq(Hsapiens, names = chr, start = start, end = end, strand = strand, as.character = TRUE)]
sra[strand == "-", sgrna := as.character(reverseComplement(DNAString(sgrna))), by = Run]
stopifnot(all(str_sub(sra$sgrna, 22, 23) == "GG"))

sra[strand == "+", cleavage := end - 5]
sra[strand == "-", cleavage := start + 5]

# Download BAM
bamfile <- paste0(sra$Run[i], ".bam")
cmd <- paste0("sam-dump ", sra$Run[i], " | samtools view -hu - | sambamba sort -t 4 -o ", bamfile, " /dev/stdin")
message(cmd)
system(cmd)

# Cleavage site
cleavage.gr <- with(sra, GRanges(seqnames = chr,
                                 ranges = IRanges(start = cleavage, width = 1),
                                 strand = strand))[i]
# Get indels
GetIndelsOverbeek <- function(bam.file, cleavage.site, focus) {
  
  param <- Rsamtools::ScanBamParam(what = c("seq"))
  
  bam <- readGAlignments(bam.file, use.names = TRUE, param = param)

  # ==========
  # Deletions
  # ==========
  
  # First fet GRanges of deletions for each pair
  del <- bam[grepl("D", cigar(bam))]
  
  if(length(del) > 0) {
  
    del.ir <- cigarRangesAlongReferenceSpace(cigar(del), ops = "D", pos = start(del))
    names(del.ir) <- names(del)
    del.ir <- unlist(del.ir)
    del.gr <- GRanges(seqnames = Rle(seqnames(del)[match(names(del.ir), names(del))]),
                            ranges = del.ir,
                            strand = Rle("*"),
                            indel = "Deletion",
                            size = width(del.ir),
                            read = names(del.ir))
    # stopifnot(length(del.gr) == length(del))
    
    # Then apply filtering based on proximity to cleavage site
    cleavage.focus <- resize(cleavage.site, width = width(cleavage.site) + (2 * focus), fix = "center")
    correct.del.gr <- subsetByOverlaps(del.gr, cleavage.focus, type = "any", ignore.strand = TRUE)
  
  } else {
    
    correct.del.gr <- GRanges()
    
  }
  
  # ==========
  # Insertions
  # ==========
  
  # First get GRanges of insertions for each pair
  ins <- bam[grepl("I", cigar(bam))]
  
  if(length(ins) > 0) {
  
    ins.ir <- cigarRangesAlongReferenceSpace(cigar(ins), ops = "I", pos = start(ins))
    names(ins.ir) <- names(ins)
    ins.ir <- unlist(ins.ir)
    ins.width <- width(unlist(cigarRangesAlongQuerySpace(cigar(ins), ops = "I"))) # get width from query space, as width in reference space is 0
    
    # Get insertion nucleotide
    nt.ins.irl <- cigarRangesAlongQuerySpace(cigar(ins), ops = "I")
    names(nt.ins.irl) <- names(ins)
    nt.ins.ir <- unlist(nt.ins.irl)
    nt.ins.seq <- mcols(ins)$seq[match(names(nt.ins.ir), names(ins))]
    names(nt.ins.seq) <- names(ins)[match(names(nt.ins.ir), names(ins))]
    nt.ins.nt <- as.character(subseq(nt.ins.seq, start = start(nt.ins.ir), width = width(nt.ins.ir)))
    
    ins.gr <- GRanges(seqnames = Rle(seqnames(ins)[match(names(ins.ir), names(ins))]),
                            ranges = IRanges(start = start(ins.ir), width = 1),
                            strand = Rle("*"),
                            indel = "Insertion",
                            size = ins.width,
                            read = names(ins.ir),
                            nt = nt.ins.nt)
    # stopifnot(length(ins.gr) == length(ins))
    
    # Then apply filtering based on proximity to cleavage site (have kept this separate in case want different focusses for insertion and deletion)
    cleavage.focus <- resize(cleavage.site, width = width(cleavage.site) + (2 * focus), fix = "center")
    correct.ins.gr <- subsetByOverlaps(resize(ins.gr, 1), cleavage.focus, type = "any", ignore.strand = TRUE)
  
  } else {
    
    correct.ins.gr <- GRanges()
    
  }
  
  # ==========
  # Indels
  # ==========
  
  if(length(correct.del.gr) > 0) correct.del.gr$nt <- as.character(NA)
  
  if(length(correct.del.gr) > 0 | length(correct.ins.gr) > 0) {
    
    indels.gr <- c(correct.del.gr, correct.ins.gr)
    
    # Finally ensure same read does not appear twice (i.e. if same indel covered by both reads of a pair)
    names(indels.gr) <- NULL
    indels.gr <- indels.gr[!duplicated(indels.gr$read)]
    
    # Annotate with read count
    # if(length(indels.gr) > 0) indels.gr$total_reads <- length(bam)
    
    # Write out to read length table
    read.length.dt <- data.table(total_reads = length(bam),
                                 exp = sra[i]$Run,
                                 cell = sra[i]$isolate,
                                 time = sra[i]$treatment,
                                 replicate = sra[i]$replicate,
                                 sgrna_chr = sra[i]$chr,
                                 sgrna_start = sra[i]$start,
                                 sgrna_end = sra[i]$end,
                                 sgrna_strand = sra[i]$strand,
                                 sgrna = sra[i]$sgrna,
                                 cleavage = sra[i]$cleavage)
    write.table(read.length.dt, "spacer_read_length.tsv", append = TRUE, quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")
    
    return(indels.gr)
    
  } else {
    read.length.dt <- data.table(total_reads = length(bam),
                                 exp = sra[i]$Run,
                                 cell = sra[i]$isolate,
                                 time = sra[i]$treatment,
                                 replicate = sra[i]$replicate,
                                 sgrna_chr = sra[i]$chr,
                                 sgrna_start = sra[i]$start,
                                 sgrna_end = sra[i]$end,
                                 sgrna_strand = sra[i]$strand,
                                 sgrna = sra[i]$sgrna,
                                 cleavage = sra[i]$cleavage)
    write.table(read.length.dt, "spacer_read_length.tsv", append = TRUE, quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")    
    return(GRanges())
    
  }
  
}

indels.gr <- GetIndelsOverbeek(bam.file = bamfile, cleavage.site = cleavage.gr, focus = 5)

# Annotate
indels.dt <- as.data.table(indels.gr)
indels.dt[, `:=` (exp = sra[i]$Run,
                  cell = sra[i]$isolate,
                  time = sra[i]$treatment,
                  replicate = sra[i]$replicate,
                  sgrna_chr = sra[i]$chr,
                  sgrna_start = sra[i]$start,
                  sgrna_end = sra[i]$end,
                  sgrna_strand = sra[i]$strand,
                  sgrna = sra[i]$sgrna,
                  cleavage = sra[i]$cleavage)]

if(!dir.exists(paste0(sra[i]$isolate, "/", sra[i]$treatment))) dir.create(paste0(sra[i]$isolate, "/", sra[i]$treatment), recursive = TRUE)
write.table(indels.dt, paste0(sra[i]$isolate, "/", sra[i]$treatment, "/", sra[i]$Run, ".indels.tsv"), sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

# Tidy up
cmd <- paste0("pigz ", sra[i]$isolate, "/", sra[i]$treatment, "/", sra[i]$Run, ".indels.tsv")
system(cmd)

file.remove(bamfile)
file.remove(paste0(bamfile, ".bai"))

ptm <- proc.time() - ptm
print(paste("Finished in", round(ptm[[3]]/60), "minutes"))