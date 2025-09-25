#' Read an Iso-Seq *summary.txt* report into tidy tibbles
#'
#' `read_isoseq_summary()` parses the plain-text summary produced by SQANTI3,
#' TALON, or similar long-read annotation tools and converts every metrics
#' block into a tidy tibble.  The function returns a named list that downstream
#' helpers such as [`plot_isoseq_overview()`] can consume directly.
#'
#' @param path `character(1)`
#'   File path to the *summary.txt* report.
#'
#' @return A named `list` of six tibbles
#'   \describe{
#'     \item{`headline`}{Headline counts (*Input transcripts*, *Unique
#'       transcripts*, *Unique genes*).}
#'     \item{`counts_isoform`}{Isoform-level classification counts (with `%`
#'       column if present).}
#'     \item{`counts_read`}{Read-level classification counts (with `%`
#'       column if present).}
#'     \item{`junction_df`}{Splice-junction type breakdown.}
#'     \item{`gene_df`}{Known *vs.* novel gene counts.}
#'     \item{`rt_df`}{RT-switching artefact counts.}
#'   }
#'
#' @details
#' Sections are identified by their header lines (e.g. *“Classifications, by
#' isoform”*).  Everything between two headers is captured and parsed with a
#' small regular-expression helper, making the function resilient to extra
#' blank lines or horizontal rulers (`------`).  Percentages in parentheses are
#' kept when present; otherwise the `pct` column is filled with `NA`.
#'
#' @importFrom stringr str_trim str_replace str_extract str_detect
#' @importFrom dplyr mutate select
#' @importFrom tibble tibble
#'
#' @examples
#' \dontrun{
#' s <- read_isoseq_summary("LL1_S1.summary.txt")
#' str(s, max.level = 1)
#' }
#' @export
read_isoseq_summary <- function(path) {
  stopifnot(length(path) == 1L, file.exists(path))

  ## helper: trim & integer-coerce ------------------------------------------------
  trim   <- function(x) str_trim(x, side = "both")
  to_int <- function(x) as.integer(gsub(",", "", x, fixed = TRUE))

  ## helper: parse one block → tibble<item, n, pct> ------------------------------
  block_tbl <- function(lines) {
    tibble(raw = trim(lines)) |>
      mutate(
        item = str_replace(raw, ":.*", ""),
        n    = to_int(str_extract(raw, "(?<=: )\\d[\\d,]*")),
        pct  = as.numeric(str_extract(raw, "(?<=\\()[0-9.]+(?=%\\))"))
      ) |>
      select(-raw)
  }

  ## ── 1  read & pre-clean -------------------------------------------------------
  L <- readLines(path, warn = FALSE) |> trim()
  L <- L[L != "" & !str_detect(L, "^[-]+$")]      # drop blanks & rulers

  ## ── 2  headline block ---------------------------------------------------------
  headline <- block_tbl(
    L[str_detect(L, "^(Input transcripts|Unique transcripts|Unique genes)")]
  ) |>
    select(metric = item, value = n)

  ## ── 3  locate section headers -------------------------------------------------
  idx <- list(
    iso  = which(str_detect(L, "^Classifications, by isoform")),
    read = which(str_detect(L, "^Classifications, by read")),
    junc = which(str_detect(L, "^Junctions")),
    gene = which(str_detect(L, "^Genes$")),
    rt   = which(str_detect(L, "^RT Switching"))
  )

  slice_between <- function(start, end) L[(start + 1L):(end - 1L)]

  ## ── 4  parse each block -------------------------------------------------------
  counts_isoform <- slice_between(idx$iso,  idx$read) |> block_tbl() |>
    select(class = item, n, pct)

  counts_read    <- slice_between(idx$read, idx$junc) |> block_tbl() |>
    select(class = item, n, pct)

  junction_df    <- slice_between(idx$junc, idx$gene) |> block_tbl() |>
    select(type = item, n)

  gene_df        <- slice_between(idx$gene, idx$rt)   |> block_tbl() |>
    select(status = item, n, pct)

  rt_df          <- L[(idx$rt + 1L):length(L)]        |> block_tbl() |>
    select(measure = item, n)

  ## ── 5  return -----------------------------------------------------------------
  list(
    headline       = headline,
    counts_isoform = counts_isoform,
    counts_read    = counts_read,
    junction_df    = junction_df,
    gene_df        = gene_df,
    rt_df          = rt_df
  )
}



#' Plot a multi-panel graphical overview of an Iso-Seq summary
#'
#' Produces a compact dashboard with headline counts, stacked‐bar
#' classifications, splice-junction composition, gene novelty, and
#' RT-switching counts.  Accepts either a *summary.txt* path (parsed with
#' [\code{read_isoseq_summary()}]) or the already-parsed list returned by that
#' function.
#'
#' @importFrom ggplot2 ggplot aes geom_text geom_col geom_hline position_stack
#'   scale_y_continuous scale_fill_manual labs theme_void theme_minimal
#'   theme element_text ggsave scale_fill_brewer
#' @importFrom dplyr bind_rows mutate group_by ungroup
#' @importFrom magrittr "%>%"
#' @importFrom forcats fct_reorder
#' @importFrom scales percent comma percent_format
#' @importFrom RColorBrewer brewer.pal
#' @import patchwork
#'
#' @param x Either a character path to *summary.txt* **or** the list produced
#'   by [\code{read_isoseq_summary()}].
#' @param save_to Optional file name; if supplied, the figure is written via
#'   \code{ggplot2::ggsave()}.  Any extension recognised by *ggplot2* is
#'   allowed (e.g. \code{"overview.pdf"}).
#' @param width,height Device size (in inches) used when \code{save_to} is not
#'   \code{NULL}.  Defaults: 11 × 8.5.
#'
#' @return A \pkg{patchwork} / \pkg{ggplot2} object invisibly.
#' @export
plot_isoseq_overview <- function(x,
                                 save_to = NULL,
                                 width   = 11,
                                 height  = 8.5) {

  ## ── 1  Input handling ──────────────────────────────────────────────────────
  if (is.character(x) && length(x) == 1L) {
    data <- read_isoseq_summary(x)
  } else if (is.list(x) &&
             all(c("headline", "counts_isoform", "counts_read",
                   "junction_df", "gene_df", "rt_df") %in% names(x))) {
    data <- x
  } else {
    stop("`x` must be a file path or the list returned by read_isoseq_summary().")
  }
  class_colors <- brewer.pal(9, name = "Pastel1")
  names <- c("Genic intron", "Other", "Antisense",  "Genic genomic", "Novel in catalog", "Novel not in catalog",  "Full splice match", "Intergenic", "Incomplete splice match")
  names_padded <- as.character(sapply(names, function(x) sprintf("%-*s", 24, x)))
  names(class_colors) <- names_padded
  with(data, {

    percent_lab <- function(n, total) percent(n / total, accuracy = 0.1)

    ## 1a  Headline strip ----
    p_head <- ggplot(headline,
                     aes(metric, value,
                         label = paste(metric, comma(value), sep = " : "))) +
      geom_text(size = 4, fontface = "bold", hjust = 0) +
      theme_void() +
      labs(title = "Iso-Seq summary metrics")

    ## 1b  Stacked classification bars ----
    class_df <- bind_rows(
      mutate(counts_isoform, level = "Isoform"),
      mutate(counts_read,    level = "Read")
    ) %>%
      group_by(level) %>%
      mutate(pct   = n / sum(n),
             class = fct_reorder(class, n, sum)) %>%
      ungroup()

    p_class <- ggplot(class_df, aes(level, pct, fill = class)) +
      geom_col(width = .7) +
      geom_text(aes(label = percent_lab(n, sum(n))),
                position = position_stack(vjust = 0.5), size = 2) +
      scale_y_continuous(labels = percent_format()) +
      scale_fill_manual(values = class_colors) +
      labs(x = NULL, y = "Proportion", fill = "Classification",
           title = "Transcript classifications") +
      theme_minimal(base_size = 11) +
      theme(legend.position = "right")

    ## 1c  Splice-junction types ----
    p_junc <- ggplot(junction_df,
                     aes(fct_reorder(type, n), n, fill = type)) +
      geom_col(width = .7, show.legend = FALSE) +
      geom_text(aes(label = comma(n)), vjust = -0.4, size = 2) +
      scale_fill_brewer(palette = "Set3") +
      labs(x = NULL, y = "Count", title = "Splice-junction types") +
      theme_minimal(base_size = 11) +
      theme(axis.text.x = element_text(angle = 25, hjust = 1))

    ## 1d  Gene novelty ----
    p_gene <- ggplot(gene_df, aes(status, n, fill = status)) +
      geom_col(width = .5, show.legend = FALSE) +
      geom_text(aes(label = comma(n)), vjust = -0.4, size = 2) +
      scale_fill_brewer(palette = "Set2") +
      labs(x = NULL, y = "Genes", title = "Genes discovered") +
      theme_minimal(base_size = 11)

    ## 1e  RT-switching artefacts ----
    p_rt <- ggplot(rt_df,
                   aes(fct_reorder(measure, n), n, fill = measure)) +
      geom_col(show.legend = FALSE) +
      geom_text(aes(label = comma(n)), vjust = -0.4, size = 2) +
      scale_fill_brewer(palette = "Set2") +
      labs(x = NULL, y = "Count", title = "RT-switching artefacts") +
      theme_minimal(base_size = 11) +
      theme(axis.text.x = element_text(angle = 20, hjust = 1))

    ## ── 2  Assemble layout ───────────────────────────────────────────────────
    # fig <- p_head /
    #   (p_class | (p_junc / (p_gene | p_rt))) +
    fig <-(p_class | (p_junc / (p_gene | p_rt))) +
      plot_layout(heights = c(1, 4), guides = "collect") &
      theme(plot.title = element_text(size = 10, face = "bold"))

    ## ── 3  Save if requested ─────────────────────────────────────────────────
    if (!is.null(save_to))
      ggsave(save_to, fig, width = width, height = height)

    fig
  })
}


#' Fit asymptotic models to an Iso-Seq saturation (rarefaction) curve
#'
#' `fit_isoseq_saturation()` reads a tab-separated *`*.saturation.txt`* file
#' produced by PacBio’s Iso-Seq pipeline, fits three common asymptotic models,
#' and returns the fitted objects, asymptote estimates, saturation fractions,
#' an AIC model comparison, and ready-made `ggplot2` overlays.
#'
#' @param file   `character(1)`
#'   Path to the saturation table.  The file **must** contain a column named
#'   `reads` (number of FLNC reads subsampled) followed by one or more feature
#'   columns (e.g. `unique_genes`, `unique_isoforms`).
#' @param feature `character` or `integer` (default `2`)
#'   Which column to model.  Supply either the *column name* or its position.
#'
#' @return (invisibly) a `list` with components
#'   \describe{
#'     \item{`data`}{Tibble-like `data.frame` with columns `reads`, `value`.}
#'     \item{`fits`}{Named list with `exp`, `mm`, and `log` model objects.}
#'     \item{`asymptote`}{Named numeric vector of plateau estimates.}
#'     \item{`saturation`}{Observed / asymptote fraction for each model.}
#'     \item{`aic`}{`stats::AIC()` table comparing the three fits.}
#'     \item{`plots`}{List of three `ggplot` overlays (one per model).}
#'   }
#'
#' @examples
#' \dontrun{
#' res <- fit_isoseq_saturation("LL1_S1.saturation.txt")
#' res$asymptote           # plateau estimates
#' res$plots$mm            # Michaelis–Menten overlay
#' }
#'
#' @importFrom utils read.table
#' @importFrom ggplot2 ggplot aes geom_point stat_function geom_hline annotate
#'   theme_bw labs geom_line
#' @importFrom stats nls SSasymp nls.control AIC coef median predict
#' @importFrom drc drm LL.4
#' @export
fit_isoseq_saturation <- function(file, feature = 2) {

  ## ── 1  read & pre-process ---------------------------------------------------
  dat <- utils::read.table(file, sep = "\t", header = TRUE, check.names = FALSE)

  if (is.character(feature)) {
    stopifnot(feature %in% colnames(dat))
    dat <- dat[, c("reads", feature)]
  } else {
    dat <- dat[, c(1, feature)]
  }
  colnames(dat) <- c("reads", "value")

  x <- dat$reads
  y <- dat$value

  ## ── 2  model fits -----------------------------------------------------------
  ## 2a  negative exponential  (self-starting)
  fit_exp <- stats::nls(y ~ stats::SSasymp(x, Asym, R0, lrc),
                        algorithm = "plinear")
  A_exp   <- stats::coef(fit_exp)["Asym"]

  ## 2b  Michaelis–Menten hyperbola
  mm <- function(x, Smax, Km) Smax * x / (Km + x)
  start <- list(Smax = max(y) * 1.2, Km = stats::median(x))
  fit_mm <- stats::nls(y ~ mm(x, Smax, Km), start = start,
                       control = stats::nls.control(maxiter = 2000,
                                                    warnOnly = TRUE))
  A_mm   <- stats::coef(fit_mm)["Smax"]

  ## 2c  4-parameter logistic
  fit_log <- drc::drm(y ~ x,
                      fct = drc::LL.4(names = c("Slope", "Lower",
                                                "Upper", "ED50")))
  A_log   <- stats::coef(fit_log)["Upper:(Intercept)"]

  ## ── 3  saturation fractions -------------------------------------------------
  sat_mm  <- y[length(y)] / A_mm
  sat_exp <- y[length(y)] / A_exp
  sat_log <- y[length(y)] / A_log

  ## ── 4  plotting helpers -----------------------------------------------------
  p_mm <- ggplot2::ggplot(dat, ggplot2::aes(reads, value)) +
    ggplot2::geom_point(size = 1.2, colour = "grey40") +
    ggplot2::stat_function(fun = function(x)
      stats::predict(fit_mm, list(x = x)),
      linewidth = 1, color = "steelblue") +
    ggplot2::geom_hline(yintercept = A_mm, linetype = 2) +
    ggplot2::annotate("text", x = max(x) * 0.6, y = A_mm,
                      label = paste0("Asymptote ≈ ", round(A_mm)),
                      vjust = -0.5) +
    ggplot2::theme_bw() +
    ggplot2::labs(title = "Rarefaction (Michaelis–Menten)",
                  x = "Reads", y = "Feature")

  p_exp <- ggplot2::ggplot(dat, ggplot2::aes(reads, value)) +
    ggplot2::geom_point(size = 1.2, colour = "grey40") +
    ggplot2::stat_function(fun = function(x)
      stats::predict(fit_exp, list(x = x)),
      linewidth = 1, color = "steelblue") +
    ggplot2::geom_hline(yintercept = A_exp, linetype = 2) +
    ggplot2::annotate("text", x = max(x) * 0.6, y = A_exp,
                      label = paste0("Asymptote ≈ ", round(A_exp)),
                      vjust = -0.5) +
    ggplot2::theme_bw() +
    ggplot2::labs(title = "Rarefaction (negative exponential)",
                  x = "Reads", y = "Feature")

  p_log <- ggplot2::ggplot(dat, ggplot2::aes(reads, value)) +
    ggplot2::geom_point(size = 1.2, colour = "grey40") +
    ggplot2::geom_line(size = 1, colour = "steelblue",
                       data = data.frame(reads = x,
                                         value = stats::predict(
                                           fit_log,
                                           newdata = data.frame(x = x)))) +
    ggplot2::geom_hline(yintercept = A_log, linetype = 2) +
    ggplot2::annotate("text", x = max(x) * 0.6, y = A_log,
                      label = paste0("Asymptote ≈ ", round(A_log)),
                      vjust = -0.5) +
    ggplot2::theme_bw() +
    ggplot2::labs(title = "Rarefaction (4-parameter logistic)",
                  x = "Reads", y = "Feature")

  ## ── 5  assemble & return ----------------------------------------------------
  res <- list(
    data       = dat,
    fits       = list(exp = fit_exp, mm = fit_mm, log = fit_log),
    asymptote  = c(exp = A_exp, mm = A_mm, log = A_log),
    saturation = c(exp = sat_exp, mm = sat_mm, log = sat_log),
    aic        = stats::AIC(fit_exp, fit_mm, fit_log),
    plots      = list(mm = p_mm, exp = p_exp, log = p_log)
  )

  res
}
