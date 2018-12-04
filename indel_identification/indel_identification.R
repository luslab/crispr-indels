# Script to identify indels in screening data
# A. M. Chakrabarti
# Last updated: 16th October 2018

library(data.table)
library(GenomicAlignments)
library(stringr)
library(BSgenome.Hsapiens.UCSC.hg19)
library(rtracklayer)

setwd("~/Scaffidi")

source("~/Github/crispr-indels-dev/R/functions.R")

# ==========
# Get HepG2 indels
# ==========

bam.file <- "phase2/NGS-13068.indel.q38.bam"

param <- Rsamtools::ScanBamParam(what = c("seq"))

bam <- readGAlignmentPairs(bam.file, use.names = TRUE, param = param)

first <- GenomicAlignments::first(bam, real.strand = TRUE) # Get each pair of a read # needed to add :: else uses wrong function
second <- GenomicAlignments::second(bam, real.strand = TRUE)

# First get GRanges of deletions for each pair
first.del <- first[grepl("D", cigar(first))]
first.del.ir <- cigarRangesAlongReferenceSpace(cigar(first.del), ops = "D", pos = start(first.del))
names(first.del.ir) <- names(first.del)
first.del.ir <- unlist(first.del.ir)
first.del.gr <- GRanges(seqnames = Rle(seqnames(first.del)[match(names(first.del.ir), names(first.del))]),
                        ranges = first.del.ir,
                        strand = Rle("*"),
                        indel = "Deletion",
                        size = width(first.del.ir),
                        read = names(first.del.ir))
stopifnot(length(first.del.gr) == length(first.del))

second.del <- second[grepl("D", cigar(second))]
second.del.ir <- cigarRangesAlongReferenceSpace(cigar(second.del), ops = "D", pos = start(second.del))
names(second.del.ir) <- names(second.del)
second.del.ir <- unlist(second.del.ir)
second.del.gr <- GRanges(seqnames = Rle(seqnames(second.del)[match(names(second.del.ir), names(second.del))]),
                         ranges = second.del.ir,
                         strand = Rle("*"),
                         indel = "Deletion",
                         size = width(second.del.ir),
                         read = names(second.del.ir))
stopifnot(length(second.del.gr) == length(second.del))

wildtype.deletions.gr <- c(first.del.gr, second.del.gr)
names(wildtype.deletions.gr) <- NULL
wildtype.deletions.gr <- unique(wildtype.deletions.gr) # 179511
wildtype.deletions.gr <- wildtype.deletions.gr[width(wildtype.deletions.gr) < 2000] # 179508

# First get GRanges of insertions for each pair
first.ins <- first[grepl("I", cigar(first))]
first.ins.ir <- cigarRangesAlongReferenceSpace(cigar(first.ins), ops = "I", pos = start(first.ins))
names(first.ins.ir) <- names(first.ins)
first.ins.ir <- unlist(first.ins.ir)
first.ins.width <- width(unlist(cigarRangesAlongQuerySpace(cigar(first.ins), ops = "I"))) # get width from query space, as width in reference space is 0

# Get insertion nucleotide
ins.ir <- unlist(cigarRangesAlongQuerySpace(cigar(first.ins), ops = "I"))
ins.seq <- mcols(first.ins)$seq
ins.nt <- as.character(subseq(ins.seq, start = start(ins.ir), width = width(ins.ir)))

first.ins.gr <- GRanges(seqnames = Rle(seqnames(first.ins)[match(names(first.ins.ir), names(first.ins))]),
                        ranges = IRanges(start = start(first.ins.ir), width = 1),
                        strand = Rle("*"),
                        indel = "Insertion",
                        size = first.ins.width,
                        read = names(first.ins.ir),
                        nt = ins.nt)
stopifnot(length(first.ins.gr) == length(first.ins))

second.ins <- second[grepl("I", cigar(second))]
second.ins.ir <- cigarRangesAlongReferenceSpace(cigar(second.ins), ops = "I", pos = start(second.ins))
names(second.ins.ir) <- names(second.ins)
second.ins.ir <- unlist(second.ins.ir)
second.ins.width <- width(unlist(cigarRangesAlongQuerySpace(cigar(second.ins), ops = "I"))) # get width from query space, as width in reference space is 0

# Get insertion nucleotide
ins.ir <- unlist(cigarRangesAlongQuerySpace(cigar(second.ins), ops = "I"))
ins.seq <- mcols(second.ins)$seq
ins.nt <- as.character(subseq(ins.seq, start = start(ins.ir), width = width(ins.ir)))

second.ins.gr <- GRanges(seqnames = Rle(seqnames(second.ins)[match(names(second.ins.ir), names(second.ins))]),
                         ranges = IRanges(start = start(second.ins.ir), width = 1),
                         strand = Rle("*"),
                         indel = "Insertion",
                         size = second.ins.width,
                         read = names(second.ins.ir),
                         nt = ins.nt)
stopifnot(length(second.ins.gr) == length(second.ins))

wildtype.insertion.gr <- c(first.ins.gr, second.ins.gr)
names(wildtype.insertion.gr) <- NULL
wildtype.insertion.gr <- unique(wildtype.insertion.gr)

hepg2_deletions <- wildtype.deletions.gr
hepg2_insertions <- wildtype.insertion.gr

saveRDS(hepg2_deletions, "~/Scaffidi/revisions/hepg2_deletions.rds")
saveRDS(hepg2_insertions, "~/Scaffidi/revisions/hepg2_insertions.rds")

# ==========
# Identify indels in each experiment
# ==========

# Define pulldown region
pulldown.gr <- import.bed("sgrna/probed_regions.bed")

# 100s
pool <- import.bed("sgrna/Pool1.bed")
cleavage.gr <- resize(resize(pool, width = 6, fix = "end"), width = 1, fix = "start")
indel.1 <- GetIndels(bam.file = "phase2/NGS-13075.indel.q38.bam",
                     pulldown.gr = pulldown.gr,
                     cleavage.site = cleavage.gr,
                     focus = 5,
                     exp = "NGS-13075_100-1")

pool <- import.bed("sgrna/Pool2.bed")
cleavage.gr <- resize(resize(pool, width = 6, fix = "end"), width = 1, fix = "start")
indel.2 <- GetIndels(bam.file = "phase2/NGS-13076.indel.q38.bam",
                     pulldown.gr = pulldown.gr,
                     cleavage.site = cleavage.gr,
                     focus = 5,
                     exp = "NGS-13076_100-2")

pool <- import.bed("sgrna/Pool3.bed")
cleavage.gr <- resize(resize(pool, width = 6, fix = "end"), width = 1, fix = "start")
indel.3 <- GetIndels(bam.file = "phase2/NGS-13077.indel.q38.bam",
                     pulldown.gr = pulldown.gr,
                     cleavage.site = cleavage.gr,
                     focus = 5,
                     exp = "NGS-13077_100-3")

# 450-5
pool5 <- import.bed("sgrna/450_pool_5.bed")
cleavage.gr <- resize(resize(pool5, width = 6, fix = "end"), width = 1, fix = "start")
indel.5.T <- GetIndels(bam.file = "phase2/NGS-13069.indel.q38.bam",
                       pulldown = pulldown.gr,
                       cleavage.site = cleavage.gr,
                       focus = 5,
                       exp = "NGS-13069_450-5T")

indel.5.P <- GetIndels(bam.file = "phase2/NGS-13070.indel.q38.bam",
                       pulldown = pulldown.gr,
                       cleavage.site = cleavage.gr,
                       focus = 5,
                       exp = "NGS-13070_450-5P")

# 450-6
pool6 <- import.bed("sgrna/450_pool_6.bed")
cleavage.gr <- resize(resize(pool6, width = 6, fix = "end"), width = 1, fix = "start")
indel.6.T <- GetIndels(bam.file = "phase2/NGS-13072.indel.q38.bam",
                       pulldown = pulldown.gr,
                       cleavage.site = cleavage.gr,
                       focus = 5,
                       exp = "NGS-13072_450-6T")

indel.6.P <- GetIndels(bam.file = "phase2/NGS-13071.indel.q38.bam",
                       pulldown = pulldown.gr,
                       cleavage.site = cleavage.gr,
                       focus = 5,
                       exp = "NGS-13071_450-6P")

# 450-7
pool7 <- import.bed("sgrna/450_pool_7.bed")
cleavage.gr <- resize(resize(pool7, width = 6, fix = "end"), width = 1, fix = "start")
indel.7.T <- GetIndels(bam.file = "phase2/NGS-13074.indel.q38.bam",
                       pulldown = pulldown.gr,
                       cleavage.site = cleavage.gr,
                       focus = 5,
                       exp = "NGS-13074_450-7T")

indel.7.P <- GetIndels(bam.file = "phase2/NGS-13073.indel.q38.bam",
                       pulldown = pulldown.gr,
                       cleavage.site = cleavage.gr,
                       focus = 5,
                       exp = "NGS-13073_450-7P")

all.indels.gr <- c(indel.1, indel.2, indel.3,
                   indel.5.T, indel.5.P,
                   indel.6.T, indel.6.P,
                   indel.7.T, indel.7.P)
all.indels.gr[all.indels.gr$indel == "Insertion"]$size <- -all.indels.gr[all.indels.gr$indel == "Insertion"]$size

all.indels.gr <- all.indels.gr[all.indels.gr$sgrna != "NCOA6.7"]

saveRDS(all.indels.gr, "revisions/all.indels.revised.gr.rds")

# ==========
# Create summary table
# ==========

# Load sgRNA data
sgrna.dt <- fread("~/Dropbox (Lab)/CRISPR-indels/ref/sgrna.tsv")
setnames(sgrna.dt, "id", "sgrna")
targets <- import.bed("~/Scaffidi/sgrna/sgRNA_PAM_hg19_curated.bed")

selected.targets <- c(import.bed("sgrna/Pool1.bed")$name,
                      import.bed("sgrna/Pool2.bed")$name,
                      import.bed("sgrna/Pool3.bed")$name,
                      import.bed("sgrna/450_pool_5.bed")$name,
                      import.bed("sgrna/450_pool_6.bed")$name,
                      import.bed("sgrna/450_pool_7.bed")$name)

selected.targets <- sort(unique(selected.targets))

# ==========
# Get precision
# ==========
indels.dt <- as.data.table(all.indels.gr)
indels.dt[, indel_count := .N, by = .(seqnames, start, end, nt, sgrna)]
indels.dt[, total_indels := .N, by = .(sgrna)]
indels.dt[, indel_frequency := indel_count/total_indels]
indels.dt[, max_indel_frequency := max(indel_frequency), by = sgrna]

# Assign indel precision
indels.dt[max_indel_frequency > 0 & max_indel_frequency <= 0.25, group := "I"]
indels.dt[max_indel_frequency > 0.25 & max_indel_frequency <= 0.5, group := "M"]
indels.dt[max_indel_frequency > 0.5 & max_indel_frequency <= 1, group := "P"]

sgrna.threshold <- unique(indels.dt[total_indels >= 10]$sgrna)

# Print out table of max frequency and group
targets.summary.dt <- unique(indels.dt[total_indels >= 10, .(sgrna, max_indel_frequency, total_indels, group)])
fwrite(targets.summary.dt, "~/Dropbox (Lab)/CRISPR-indels/resubmission/data/649_targets_with_predictability.tsv", sep = "\t", quote = F)

targets.summary.dt <- fread("~/Dropbox (Lab)/CRISPR-indels/resubmission/data/649_targets_with_predictability.tsv")
sgrna.dt <- fread("~/Dropbox (Lab)/CRISPR-indels/ref/sgrna.tsv")
setnames(sgrna.dt, "id", "sgrna")
setkey(sgrna.dt, "sgrna")
setkey(targets.summary.dt, "sgrna")

fwrite(targets.summary.dt, "~/Dropbox (Lab)/CRISPR-indels/resubmission/data/649_targets_with_predictability_and_sgrna_metadata.tsv", sep = "\t", quote = F)

targets.summary.dt <- targets.summary.dt[sgrna.dt]


all.indels.gr <- readRDS("~/Scaffidi/revisions/all.indels.revised.gr.rds")

indels.dt <- as.data.table(all.indels.gr)
indels.dt[, indel_count := .N, by = .(seqnames, start, end, nt, sgrna)]
indels.dt[, total_indels := .N, by = .(sgrna)]
indels.dt[, indel_frequency := indel_count/total_indels]
indels.dt[, max_indel_frequency := max(indel_frequency), by = sgrna]

# Assign indel precision
indels.dt[max_indel_frequency > 0 & max_indel_frequency <= 0.25, group := "I"]
indels.dt[max_indel_frequency > 0.25 & max_indel_frequency <= 0.5, group := "M"]
indels.dt[max_indel_frequency > 0.5 & max_indel_frequency <= 1, group := "P"]

sgrna.threshold <- unique(indels.dt[total_indels >= 10]$sgrna)

# Print out table of max frequency and 
targets.summary.dt <- unique(indels.dt[total_indels >= 10, .(sgrna, max_indel_frequency, total_indels, group)])
setnames(targets.summary.dt, c("sgRNA/Target", "Commonest indel frequency", "Total number of indels", "Precision group"))
setorder(targets.summary.dt, -`Commonest indel frequency`, `sgRNA/Target`)
targets.summary.dt$`Commonest indel frequency` <- round(targets.summary.dt$`Commonest indel frequency`, 2)

test <- fread(file.path(output.dir, "649_targets_precision.tsv"))

test2 <- merge(test, targets.summary.dt, by.y = "sgrna", by.x = "sgRNA/Target")

all(test2$`Total number of indels` == test2$total_indels)
all(test2$`Commonest indel frequency` == round(test2$max_indel_frequency, 2))
all(test2$`Precision group` == test2$group)
