#' Generate Extended Gene Intervals for Coverage Mapping (Optimized Version)
#'
#' This function generates extended gene intervals from a GRanges object
#' containing gene annotations, properly identifying the TSS based on strand.
#' The function implements proper gene boundary checking to prevent regions
#' from extending past neighboring gene TSSs or gene bodies.
#'
#' @param genes A GRanges object containing gene annotations with at least 'symbol' in mcols
#' @param extendUpstream Numeric vector of length 2 specifying min and max bp to extend upstream of TSS
#' @param extendDownstream Numeric vector of length 2 specifying min and max bp to extend downstream of TSS
#' @param useGeneBoundaries Logical indicating whether to respect gene boundaries
#' @param includeChr Character vector of chromosome names to include
#' @param outputFile Path to save the output BED file (NULL to return GRanges without saving)
#' @param includeWeights Logical indicating whether to include gene weights in output
#' @param geneScaleFactor Numeric scaling factor for gene weights (based on gene length)
#'
#' @return GRanges object with extended gene intervals
#' @export
generateExtendedGeneIntervals <- function(
    genes,
    extendUpstream = c(1000, 100000),
    extendDownstream = c(1000, 100000),
    useGeneBoundaries = TRUE,
    includeChr = paste0("chr", c(1:22, "X")),
    outputFile = NULL,
    includeWeights = TRUE,
    geneScaleFactor = 5,
    cores = 1
) {
  # Load required libraries
  suppressPackageStartupMessages({
    library(GenomicRanges)
    library(rtracklayer)
    library(dplyr)
    library(pbmcapply)
  })

  # Start timing
  startTime <- Sys.time()
  message("Starting gene interval extension at ", format(startTime))

  # Keep only standard chromosomes in seqinfo
  stdChrs <- includeChr
  message("Keeping only standard chromosomes: ", paste(head(stdChrs, 5), collapse=", "),
          ifelse(length(stdChrs) > 5, "...", ""))

  # Filter genes to only include standard chromosomes
  genes <- genes[as.character(seqnames(genes)) %in% stdChrs]

  # Reset seqlevels to only keep standard chromosomes
  seqlevels(genes) <- stdChrs

  # Sort chromosomes in standard order
  genes <- sortSeqlevels(genes)
  genes <- sort(genes)

  # Filter out genes with missing symbols
  genes <- genes[!is.na(mcols(genes)$symbol)]

  message("Processing ", length(genes), " genes across ",
          length(unique(seqnames(genes))), " chromosomes")

  # Create gene regions and identify TSS
  genes$geneStart <- start(genes)
  genes$geneEnd <- end(genes)

  # Identify TSS based on strand
  genes$tss <- ifelse(
    as.character(strand(genes)) == "+",
    genes$geneStart,
    genes$geneEnd
  )

  # Calculate gene weights based on gene length
  if (includeWeights) {
    geneWidth <- width(genes)
    m <- 1 / geneWidth
    genes$geneWeight <- 1 + m * (geneScaleFactor - 1) / (max(m) - min(m))
  }

  # Get chromosome lengths from seqinfo
  chromLengths <- seqlengths(genes)

  # Create a copy for extending
  extendedGeneRegions <- genes

  # We'll implement gene boundary checking if requested
  if (useGeneBoundaries) {
    message("Respecting gene boundaries when extending regions...")

    # Split genes by chromosome for more efficient processing
    genesByChr <- split(genes, as.character(seqnames(genes)))

    # Process each chromosome in parallel if possible
    extendedRegions <- pbmclapply(names(genesByChr), function(chr) {
      chrStartTime <- Sys.time()
      message("Processing chromosome ", chr, "...")

      # Get genes on this chromosome
      chrGenes <- genesByChr[[chr]]
      chrLength <- as.numeric(chromLengths[chr])

      # Create the maximally extended intervals (without considering gene boundaries)
      tssPositions <- chrGenes$tss
      isPlus <- as.character(strand(chrGenes)) == "+"

      # Calculate preliminary extension limits
      upMax <- ifelse(isPlus,
                      pmax(1, tssPositions - max(extendUpstream)),
                      pmin(chrLength, tssPositions + max(extendUpstream)))

      downMax <- ifelse(isPlus,
                        pmin(chrLength, tssPositions + max(extendDownstream)),
                        pmax(1, tssPositions - max(extendDownstream)))

      # For each gene, create temporary GRanges for the maximally extended region
      maxExtendedRegions <- GRanges(
        seqnames = chr,
        ranges = IRanges(
          start = ifelse(isPlus, upMax, downMax),
          end = ifelse(isPlus, downMax, upMax)
        ),
        strand = strand(chrGenes),
        tss = tssPositions,
        geneIdx = seq_along(chrGenes)
      )

      # Create a GRanges of gene bodies
      geneBodies <- GRanges(
        seqnames = chr,
        ranges = IRanges(
          start = chrGenes$geneStart,
          end = chrGenes$geneEnd
        ),
        strand = strand(chrGenes),
        geneIdx = seq_along(chrGenes)
      )

      # Find all overlaps between extended regions and gene bodies
      overlaps <- findOverlaps(maxExtendedRegions, geneBodies, ignore.strand=TRUE)

      # Filter out self-overlaps
      overlaps <- overlaps[queryHits(overlaps) != subjectHits(overlaps)]

      # If there are overlaps, we need to adjust the extension limits
      if (length(overlaps) > 0) {
        # Prepare data structure to track updated extension limits
        updatedExtensions <- data.frame(
          geneIdx = seq_along(chrGenes),
          upstream = upMax,
          downstream = downMax,
          stringsAsFactors = FALSE
        )

        # Process each overlap
        for (i in seq_along(overlaps)) {
          qHit <- queryHits(overlaps)[i]
          sHit <- subjectHits(overlaps)[i]

          # Get the original gene and the overlapping gene
          gIdx <- maxExtendedRegions$geneIdx[qHit]
          oIdx <- geneBodies$geneIdx[sHit]

          # Skip if it's a self-overlap
          if (gIdx == oIdx) next

          # Get the gene strand
          isGPlus <- isPlus[gIdx]

          # Get gene and overlap coordinates
          gTss <- tssPositions[gIdx]
          oStart <- geneBodies@ranges@start[sHit]
          oEnd <- geneBodies@ranges@start[sHit] + geneBodies@ranges@width[sHit] - 1
          oTss <- chrGenes$tss[oIdx]

          # Adjust extension limits based on strand and overlap
          if (isGPlus) {
            # For + strand genes
            # If overlapping gene is to the left (upstream)
            if (oEnd < gTss && oEnd > updatedExtensions$upstream[gIdx]) {
              updatedExtensions$upstream[gIdx] <- oEnd + 1
            }
            # If overlapping gene is to the right (downstream)
            if (oStart > gTss && oStart < updatedExtensions$downstream[gIdx]) {
              updatedExtensions$downstream[gIdx] <- oStart - 1
            }
          } else {
            # For - strand genes
            # If overlapping gene is to the right (upstream)
            if (oStart > gTss && oStart < updatedExtensions$upstream[gIdx]) {
              updatedExtensions$upstream[gIdx] <- oStart - 1
            }
            # If overlapping gene is to the left (downstream)
            if (oEnd < gTss && oEnd > updatedExtensions$downstream[gIdx]) {
              updatedExtensions$downstream[gIdx] <- oEnd + 1
            }
          }
        }

        # Apply updated extension limits
        for (i in seq_along(chrGenes)) {
          if (isPlus[i]) {
            start(chrGenes)[i] <- updatedExtensions$upstream[i]
            end(chrGenes)[i] <- updatedExtensions$downstream[i]
          } else {
            start(chrGenes)[i] <- updatedExtensions$downstream[i]
            end(chrGenes)[i] <- updatedExtensions$upstream[i]
          }
        }
      } else {
        # If no overlaps, just use the maximally extended regions
        for (i in seq_along(chrGenes)) {
          if (isPlus[i]) {
            start(chrGenes)[i] <- upMax[i]
            end(chrGenes)[i] <- downMax[i]
          } else {
            start(chrGenes)[i] <- downMax[i]
            end(chrGenes)[i] <- upMax[i]
          }
        }
      }

      # Ensure start < end for all ranges
      needSwap <- start(chrGenes) > end(chrGenes)
      if (any(needSwap)) {
        idx <- which(needSwap)
        for (i in idx) {
          temp <- start(chrGenes)[i]
          start(chrGenes)[i] <- end(chrGenes)[i]
          end(chrGenes)[i] <- temp
        }
      }

      message("Completed ", chr, " in ", round(difftime(Sys.time(), chrStartTime, units="secs"), 2), " seconds")
      return(chrGenes)
    }, mc.cores = cores)

    # Combine all chromosomes back together
    extendedGeneRegions <- do.call(c, extendedRegions)

  } else {
    # Simple extension without respecting gene boundaries
    # For + strand genes
    plusIdx <- which(as.character(strand(genes)) == "+")
    if (length(plusIdx) > 0) {
      start(extendedGeneRegions)[plusIdx] <- pmax(1, genes$tss[plusIdx] - max(extendUpstream))
      end(extendedGeneRegions)[plusIdx] <- pmin(
        genes$tss[plusIdx] + max(extendDownstream),
        chromLengths[as.character(seqnames(genes)[plusIdx])]
      )
    }

    # For - strand genes
    minusIdx <- which(as.character(strand(genes)) == "-")
    if (length(minusIdx) > 0) {
      start(extendedGeneRegions)[minusIdx] <- pmax(1, genes$tss[minusIdx] - max(extendDownstream))
      end(extendedGeneRegions)[minusIdx] <- pmin(
        genes$tss[minusIdx] + max(extendUpstream),
        chromLengths[as.character(seqnames(genes)[minusIdx])]
      )
    }
  }

  # Add TSS to metadata
  mcols(extendedGeneRegions)$tss <- genes$tss

  # Explicitly validate and trim any remaining out-of-bound ranges
  message("Final check for out-of-bound ranges...")
  tryCatch({
    extendedGeneRegions <- GenomicRanges::trim(extendedGeneRegions)
    message("Trimming completed successfully.")
  }, error = function(e) {
    message("Warning: Error during final trimming - ", e$message)
  })

  # Final sorting
  extendedGeneRegions <- sortSeqlevels(extendedGeneRegions)
  extendedGeneRegions <- sort(extendedGeneRegions)

  # Save to BED file if requested
  if (!is.null(outputFile)) {
    message("Writing output to ", outputFile)
    bedDF <- data.frame(
      seqnames = seqnames(extendedGeneRegions),
      start = start(extendedGeneRegions) - 1,  # BED format is 0-based
      end = end(extendedGeneRegions),
      name = mcols(extendedGeneRegions)$symbol,
      score = ifelse(includeWeights, round(mcols(extendedGeneRegions)$geneWeight * 1000), 1000),
      strand = strand(extendedGeneRegions),
      tss = mcols(extendedGeneRegions)$tss,
      geneStart = mcols(extendedGeneRegions)$geneStart,
      geneEnd = mcols(extendedGeneRegions)$geneEnd
    )

    write.table(
      bedDF,
      file = outputFile,
      quote = FALSE,
      sep = "\t",
      row.names = FALSE,
      col.names = FALSE
    )
    message("Extended gene intervals written to: ", outputFile)
  }

  # Report execution time
  endTime <- Sys.time()
  timeTaken <- difftime(endTime, startTime, units="secs")
  message("Function completed in ", round(timeTaken, 2), " seconds")

  return(extendedGeneRegions)
}
