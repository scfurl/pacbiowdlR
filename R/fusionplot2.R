# Simplified Interactive Fusion Gene Visualization
#
# This R script provides functions to visualize gene fusions in an interactive manner.
# It is a simplified version of a more complex fusion visualization tool.
#
# Author: Claude
# Date: March 13, 2025
# License: MIT
#
# The package provides the following main functions:
# - visualize_fusion(): Creates a visualization of a gene fusion
# - run_fusion_examples(): Demonstrates the visualization with common fusion examples
#
#' Determine the type of gene fusion
#'
#' @param contig1 Chromosome of gene1
#' @param contig2 Chromosome of gene2
#' @param direction1 Direction of gene1 in fusion
#' @param direction2 Direction of gene2 in fusion
#' @param breakpoint1 Genomic coordinate of breakpoint in gene1
#' @param breakpoint2 Genomic coordinate of breakpoint in gene2
#' @return Character string describing the fusion type
#' @export
determine_fusion_type <- function(contig1, contig2, direction1, direction2, breakpoint1, breakpoint2) {
  if (contig1 != contig2) {
    return("Translocation")
  } else if (direction1 == direction2) {
    return("Inversion")
  } else if ((direction1 == "downstream") == (breakpoint1 < breakpoint2)) {
    return("Deletion")
  } else {
    return("Duplication")
  }
}

#' Visualize a gene fusion
#'
#' @description
#' Creates a simplified visualization of a fusion between two genes, showing
#' exon structures, breakpoints, and the resulting fusion arrangement.
#'
#' @param gene1 Character. Name of the 5' gene in the fusion.
#' @param gene2 Character. Name of the 3' gene in the fusion.
#' @param breakpoint1 Numeric. Genomic coordinate of the breakpoint in gene1.
#' @param breakpoint2 Numeric. Genomic coordinate of the breakpoint in gene2.
#' @param contig1 Character. Chromosome or contig name for gene1.
#' @param contig2 Character. Chromosome or contig name for gene2.
#' @param exons1 Data frame. Optional exon structure for gene1. Should contain
#'        columns 'start' and 'end' with genomic coordinates. Optional columns are
#'        'type' (e.g., "exon", "CDS") and 'exon_number'.
#' @param exons2 Data frame. Optional exon structure for gene2. Same format as exons1.
#' @param direction1 Character. Direction of gene1 in the fusion: "downstream" or "upstream".
#' @param direction2 Character. Direction of gene2 in the fusion: "downstream" or "upstream".
#' @param color1 Character. Color for gene1 representation (hex code or R color name).
#' @param color2 Character. Color for gene2 representation (hex code or R color name).
#'
#' @return A plot visualizing the gene fusion.
#'
#' @examples
#' # Basic example with BCR-ABL1 fusion
#' visualize_fusion(
#'   gene1 = "BCR",
#'   gene2 = "ABL1",
#'   breakpoint1 = 23632610,
#'   breakpoint2 = 133729451,
#'   contig1 = "chr22",
#'   contig2 = "chr9"
#' )
#'
#' # Example with custom exon structure
#' exons1 <- data.frame(
#'   start = c(23632000, 23632500, 23633000),
#'   end = c(23632400, 23632800, 23633500),
#'   type = c("exon", "exon", "CDS"),
#'   exon_number = c("1", "2", "3")
#' )
#'
#' exons2 <- data.frame(
#'   start = c(133728000, 133729000, 133730000),
#'   end = c(133728800, 133729800, 133730800),
#'   type = c("exon", "CDS", "CDS"),
#'   exon_number = c("1", "2", "3")
#' )
#'
#' visualize_fusion(
#'   gene1 = "BCR",
#'   gene2 = "ABL1",
#'   breakpoint1 = 23632610,
#'   breakpoint2 = 133729451,
#'   contig1 = "chr22",
#'   contig2 = "chr9",
#'   exons1 = exons1,
#'   exons2 = exons2
#' )
visualize_fusion <- function(gene1, gene2,
                             breakpoint1, breakpoint2,
                             contig1, contig2,
                             exons1, exons2,
                             direction1 = "downstream",
                             direction2 = "upstream",
                             color1 = "#e5a5a5",
                             color2 = "#a7c4e5") {

  # Define dark colors for gene strands
  getDarkColor <- function(color) {
    rgb(
      min(255, max(0, col2rgb(color)["red",] - 100)),
      min(255, max(0, col2rgb(color)["green",] - 100)),
      min(255, max(0, col2rgb(color)["blue",] - 100)),
      maxColorValue = 255
    )
  }

  darkColor1 <- getDarkColor(color1)
  darkColor2 <- getDarkColor(color2)

  # Create a new plot
  par(mar = c(2, 2, 2, 2))
  plot(0, 0, type = "n", xlim = c(0, 1), ylim = c(0, 1),
       xlab = "", ylab = "", xaxt = "n", yaxt = "n")
  title(main = paste0(gene1, "-", gene2, " Fusion"))

  # Vertical coordinates of plot elements
  yExons <- 0.7      # Position for exon drawing
  yGeneNames <- 0.6  # Position for gene names
  yFusion <- 0.4     # Position for fusion representation

  # Function to draw exons
  drawExon <- function(left, right, y, color, title = "", type = "exon") {
    exonHeight <- 0.03

    if (type == "CDS") {
      # Draw coding regions as thicker bars
      rect(left, y + exonHeight, right, y + exonHeight/2 - 0.001, col = color, border = NA)
      rect(left, y - exonHeight, right, y - exonHeight/2 + 0.001, col = color, border = NA)
      # Draw border
      lines(c(left, left, right, right), c(y + exonHeight/2, y + exonHeight, y + exonHeight, y + exonHeight/2),
            col = getDarkColor(color), lend = 2)
      lines(c(left, left, right, right), c(y - exonHeight/2, y - exonHeight, y - exonHeight, y - exonHeight/2),
            col = getDarkColor(color), lend = 2)
    } else if (type == "exon") {
      rect(left, y + exonHeight/2, right, y - exonHeight/2, col = color, border = getDarkColor(color))
      # Add exon label if provided
      if (title != "") {
        text((left + right)/2, y, title, cex = 0.8)
      }
    }
  }

  # Function to draw strand with arrows
  drawStrand <- function(left, right, y, color, strand) {
    if (strand %in% c("+", "-")) {
      # Draw strand
      lines(c(left + 0.001, right - 0.001), c(y, y), col = color, lwd = 2)
      lines(c(left + 0.001, right - 0.001), c(y, y), col = rgb(1, 1, 1, 0.1), lwd = 1)
      # Indicate orientation with arrows
      if (right - left > 0.01) {
        for (i in seq(left + 0.01, right - 0.01, by = 0.02)) {
          arrows(i, y, i + 0.01 * ifelse(strand == "+", 1, -1), y,
                 col = color, length = 0.05, lwd = 2, angle = 30)
        }
      }
    }
  }

  # Plot gene 1
  gene1_left <- 0.05
  gene1_right <- 0.45

  lines(c(gene1_left, gene1_right), c(yExons, yExons), col = darkColor1)
  drawStrand(gene1_left, gene1_right, yExons, darkColor1, "+")

  # Draw example exons for gene 1
  if (missing(exons1)) {
    # Default example exons
    exon_positions <- seq(gene1_left + 0.03, gene1_right - 0.03, length.out = 5)
    exon_widths <- rep(0.03, 5)

    for (i in 1:5) {
      drawExon(exon_positions[i] - exon_widths[i]/2, exon_positions[i] + exon_widths[i]/2,
               yExons, color1, paste("E", i, sep=""), "exon")
    }
  } else {
    # Use provided exons
    total_width <- gene1_right - gene1_left
    for (i in 1:nrow(exons1)) {
      ex_left <- gene1_left + (exons1$start[i] - min(exons1$start)) /
        (max(exons1$end) - min(exons1$start)) * total_width
      ex_right <- gene1_left + (exons1$end[i] - min(exons1$start)) /
        (max(exons1$end) - min(exons1$start)) * total_width

      drawExon(ex_left, ex_right, yExons, color1,
               ifelse("exon_number" %in% colnames(exons1), exons1$exon_number[i], ""),
               ifelse("type" %in% colnames(exons1), exons1$type[i], "exon"))
    }
  }

  # Plot gene 2
  gene2_left <- 0.55
  gene2_right <- 0.95

  lines(c(gene2_left, gene2_right), c(yExons, yExons), col = darkColor2)
  drawStrand(gene2_left, gene2_right, yExons, darkColor2, "+")

  # Draw example exons for gene 2
  if (missing(exons2)) {
    # Default example exons
    exon_positions <- seq(gene2_left + 0.03, gene2_right - 0.03, length.out = 4)
    exon_widths <- rep(0.04, 4)

    for (i in 1:4) {
      drawExon(exon_positions[i] - exon_widths[i]/2, exon_positions[i] + exon_widths[i]/2,
               yExons, color2, paste("E", i, sep=""), "exon")
    }
  } else {
    # Use provided exons
    total_width <- gene2_right - gene2_left
    for (i in 1:nrow(exons2)) {
      ex_left <- gene2_left + (exons2$start[i] - min(exons2$start)) /
        (max(exons2$end) - min(exons2$start)) * total_width
      ex_right <- gene2_left + (exons2$end[i] - min(exons2$start)) /
        (max(exons2$end) - min(exons2$start)) * total_width

      drawExon(ex_left, ex_right, yExons, color2,
               ifelse("exon_number" %in% colnames(exons2), exons2$exon_number[i], ""),
               ifelse("type" %in% colnames(exons2), exons2$type[i], "exon"))
    }
  }

  # Calculate breakpoint positions on the plot
  bp1_pos <- gene1_left + (gene1_right - gene1_left) * 0.7  # Example position
  bp2_pos <- gene2_left + (gene2_right - gene2_left) * 0.3  # Example position

  if (!missing(breakpoint1) && !missing(exons1)) {
    # Calculate relative position within gene1
    bp1_pos <- gene1_left + (breakpoint1 - min(exons1$start)) /
      (max(exons1$end) - min(exons1$start)) * (gene1_right - gene1_left)
  }

  if (!missing(breakpoint2) && !missing(exons2)) {
    # Calculate relative position within gene2
    bp2_pos <- gene2_left + (breakpoint2 - min(exons2$start)) /
      (max(exons2$end) - min(exons2$start)) * (gene2_right - gene2_left)
  }

  # Draw gene names
  text((gene1_left + gene1_right)/2, yGeneNames, gene1, font = 2)
  text((gene2_left + gene2_right)/2, yGeneNames, gene2, font = 2)

  # Draw breakpoint labels
  text(bp1_pos, yExons + 0.08, paste0(contig1, ":", breakpoint1), cex = 0.8)
  text(bp2_pos, yExons + 0.08, paste0(contig2, ":", breakpoint2), cex = 0.8)

  # Draw fusion representation
  # Gene 1 part
  if (direction1 == "downstream") {
    lines(c(0.25, bp1_pos), c(yFusion, yFusion), col = darkColor1, lwd = 2)
    drawStrand(0.25, bp1_pos, yFusion, darkColor1, "+")
  } else { # upstream
    lines(c(0.25, bp1_pos), c(yFusion, yFusion), col = darkColor1, lwd = 2)
    drawStrand(0.25, bp1_pos, yFusion, darkColor1, "-")
  }

  # Gene 2 part
  if (direction2 == "downstream") {
    lines(c(bp2_pos, 0.75), c(yFusion, yFusion), col = darkColor2, lwd = 2)
    drawStrand(bp2_pos, 0.75, yFusion, darkColor2, "+")
  } else { # upstream
    lines(c(bp2_pos, 0.75), c(yFusion, yFusion), col = darkColor2, lwd = 2)
    drawStrand(bp2_pos, 0.75, yFusion, darkColor2, "-")
  }

  # Connect breakpoints with dotted line
  lines(c(bp1_pos, bp1_pos), c(yExons, yFusion), col = "red", lty = 2)
  lines(c(bp2_pos, bp2_pos), c(yExons, yFusion), col = "red", lty = 2)

  # Draw fusion junction
  points(c(bp1_pos, bp2_pos), rep(yFusion, 2), pch = 19, col = "red", cex = 1)

  # Determine and add fusion type annotation
  fusion_type <- determine_fusion_type(contig1, contig2, direction1, direction2, breakpoint1, breakpoint2)

  text(0.5, 0.2, paste("Fusion Type:", fusion_type), cex = 1, font = 2)

  # Add legend
  legend("bottomleft",
         legend = c(gene1, gene2),
         fill = c(color1, color2),
         border = c(darkColor1, darkColor2),
         cex = 0.8)
}

# Interactive examples

#' @title Run an interactive fusion visualization demo
#' @description Demonstrates the fusion visualization with different examples
#' @export
run_fusion_examples <- function() {
  # Example 1: BCR-ABL1 fusion
  cat("Example 1: BCR-ABL1 fusion (Philadelphia chromosome)\n")
  visualize_fusion(
    gene1 = "BCR",
    gene2 = "ABL1",
    breakpoint1 = 23632610,
    breakpoint2 = 133729451,
    contig1 = "chr22",
    contig2 = "chr9"
  )
  cat("Press [Enter] to see the next example...")
  invisible(readline())

  # Example 2: EWSR1-FLI1 fusion (Ewing sarcoma)
  cat("Example 2: EWSR1-FLI1 fusion (Ewing sarcoma)\n")
  visualize_fusion(
    gene1 = "EWSR1",
    gene2 = "FLI1",
    breakpoint1 = 29685752,
    breakpoint2 = 128675241,
    contig1 = "chr22",
    contig2 = "chr11",
    direction1 = "downstream",
    direction2 = "upstream",
    color1 = "#ff9e81",
    color2 = "#81b1ff"
  )
  cat("Press [Enter] to see the next example...")
  invisible(readline())

  # Example 3: Custom exon structure
  cat("Example 3: BCR-ABL1 fusion with detailed exon structure\n")
  exons1 <- data.frame(
    start = c(23632000, 23632500, 23633000, 23634000),
    end = c(23632400, 23632800, 23633500, 23634500),
    type = c("exon", "exon", "CDS", "CDS"),
    exon_number = c("1", "2", "3", "4")
  )

  exons2 <- data.frame(
    start = c(133728000, 133729000, 133730000, 133731000),
    end = c(133728800, 133729800, 133730800, 133731800),
    type = c("exon", "exon", "CDS", "CDS"),
    exon_number = c("1", "2", "3", "4")
  )

  visualize_fusion(
    gene1 = "BCR",
    gene2 = "ABL1",
    breakpoint1 = 23632610,
    breakpoint2 = 133729451,
    contig1 = "chr22",
    contig2 = "chr9",
    direction1 = "downstream",
    direction2 = "upstream",
    exons1 = exons1,
    exons2 = exons2
  )

  cat("Demo completed. You can now try your own fusions!\n")
}

# To run the demo, uncomment the next line:
# run_fusion_examples()
