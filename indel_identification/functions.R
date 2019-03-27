# Functions for CRISPR-indels
# A. M. Chakrabarti
# Last updated: 11th December 2017

# ====================================================================================================
# Get indels from BAM files (for Phase 2)
# ====================================================================================================

GetIndels <- function(bam.file, pulldown.gr, cleavage.site, focus, exp) {
  
  param <- Rsamtools::ScanBamParam(what = c("seq"))
  
  bam <- readGAlignmentPairs(bam.file, use.names = TRUE, param = param)
  pulldown.bam <- subsetByOverlaps(bam, pulldown.gr, ignore.strand = TRUE, type = "any") # Get only reads that are in the pulldown readion
  
  first <- GenomicAlignments::first(pulldown.bam, real.strand = TRUE) # Get each pair of a read # needed to add :: else uses wrong function
  second <- GenomicAlignments::second(pulldown.bam, real.strand = TRUE)
  
  # ==========
  # Deletions
  # ==========
  
  # First fet GRanges of deletions for each pair
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
  
  del.gr <- c(first.del.gr, second.del.gr)
  
  # Then apply filtering based on proximity to cleavage site
  cleavage.focus <- resize(cleavage.site, width = width(cleavage.site) + (2 * focus), fix = "center")
  correct.del.gr <- subsetByOverlaps(del.gr, cleavage.focus, type = "any", ignore.strand = TRUE)
  
  # The apply filtering to ensure deletion not present in HepG2, using an exact match
  ol <- findOverlaps(correct.del.gr, hepg2_deletions, type = "equal")
  correct.del.gr <- correct.del.gr[!1:length(correct.del.gr) %in% queryHits(ol)]
  
  # ==========
  # Insertions
  # ==========
  
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
  
  ins.gr <- c(first.ins.gr, second.ins.gr)
  
  # Then apply filtering based on proximity to cleavage site (have kept this separate in case want different focusses for insertion and deletion)
  cleavage.focus <- resize(cleavage.site, width = width(cleavage.site) + (2 * focus), fix = "center")
  correct.ins.gr <- subsetByOverlaps(resize(ins.gr, 1), cleavage.focus, type = "any", ignore.strand = TRUE)
  
  # The apply filtering to ensure an insertion not present in HepG2 (not nucleotide specific)
  if(length(correct.ins.gr) != 0) {
    ol <- findOverlaps(correct.ins.gr, hepg2_insertions, type = "equal")
    correct.ins.gr <- correct.ins.gr[!1:length(correct.ins.gr) %in% queryHits(ol)]
  }
  
  # ==========
  # Indels
  # ==========
  
  correct.del.gr$nt <- as.character(NA)
  indels.gr <- c(correct.del.gr, correct.ins.gr)
  
  # Finally ensure same read does not appear twice (i.e. if same indel covered by both reads of a pair)
  names(indels.gr) <- NULL
  indels.gr <- indels.gr[!duplicated(indels.gr$read)]

  # ==========
  # Annotate with sgRNA
  # ==========
  
  ol <- findOverlaps(indels.gr, cleavage.focus, type = "any", ignore.strand = TRUE)
  indels.gr$sgrna <- as.character(NA)
  indels.gr[queryHits(ol)]$sgrna <- cleavage.focus[subjectHits(ol)]$name
  
  stopifnot(!is.na(indels.gr$sgrna))
  indels.gr$read <- NULL
  indels.gr$exp <- exp # add experiment
  
  return(indels.gr)
  
}

# ====================================================================================================
# Function to normalise CrispRVariants list
# ====================================================================================================

NormaliseCV <- function(CV.list, normalisation.method = c("Raw", "DESeq", "TMM", "Quantile", "Total"), quant = 0.9) {
  
  norm.vc.list <- lapply(1:length(CV.list), function(i) {
    
  # Get variant counts
  message(names(CV.list[i]))
  cv <- CV.list[[i]]
  vc <- variantCounts(cv, include.nonindel = FALSE)
  
  # Remove failed experiment if present
  if(names(cv.list)[i] == "SMARCD2.1") vc <- vc[, colnames(vc) != "SETD2_KO_RepA"] 
  
  # Get size factors
  if(normalisation.method == "DESeq") {
    norm.factors <- DESeq2::estimateSizeFactorsForMatrix(vc)
  } else if(normalisation.method == "TMM") {
    norm.factors <- edgeR::calcNormFactors(vc, method = "TMM")
  } else if(normalisation.method == "Quantile") {
    norm.factors <- EBSeq::QuantileNorm(vc, quant)
  } else if(normalisation.method == "Total") {
    norm.factors <- CalculateTotalNorm(vc)
  } else if(normalisation.method == "Raw") {
    norm.factors <- 1
  } else {
    stop("Normalisation method not supported.")
  }
    
  # Divide the two by row
  norm.vc <- sweep(vc, MARGIN = 2, norm.factors, `/`)
  
  return(norm.vc)
  
  })
  
  return(norm.vc.list)
  
}

# ====================================================================================================
# Function to get closest indel from CrispRVariants indel code
# ====================================================================================================

GetClosestIndelSize <- function(indel_code) {

  # Get each indel
  indels <- str_split(indel_code, ",")[[1]]

  # Get the absolute distance from the cleavage site
  indel <- str_split(indels, pattern = ":")
  
  distance_from_cleavage <- abs(sapply(indel, function(x) {
    
    # Need to add for deletion, but just take position for insertion
    if(grepl("D", x[2])) {
      as.numeric(x[1]) + as.numeric(gsub("I|D", "", x[2]))
    } else if(grepl("I", x[2])) {
      as.numeric(x[1])
    }
  
  }))

  # Select the closest one
  # In case of tie select first
  indels[distance_from_cleavage == min(distance_from_cleavage)][1]

}

# ====================================================================================================
# Function to calculate total normalisation factor
# ====================================================================================================

CalculateTotalNorm <- function(mat) {
  
  norm_factor <- colSums(mat)
    return(norm_factor/mean(norm_factor))
  
}