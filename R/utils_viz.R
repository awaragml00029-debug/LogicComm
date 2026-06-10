# R/utils_viz.R
#
# LogicComm publication visual system: brand palette, colorblind-safe scales,
# a polished theme, and a figure export helper. Every scale here is designed to
# stay distinguishable under deuteranopia/protanopia and to survive greyscale
# printing (the brand hues differ in luminance, not only in hue).

#' LogicComm Brand Colors
#'
#' Named list of the LogicComm brand colors used across all package figures.
#' Use these for bespoke annotations so custom layers match the package style.
#'
#' @return A named list of hex colors.
#' @export
logiccomm_brand <- list(
  primary   = "#1F6F78",  # deep teal (signature)
  secondary = "#33415C",  # slate indigo
  accent    = "#E8A33D",  # amber
  highlight = "#D1495B",  # coral
  positive  = "#1F6F78",  # teal  (up / high)
  negative  = "#B45A00",  # burnt orange (down / low)
  neutral   = "#8A8F99",  # grey
  ink       = "#23262B",  # near-black text
  paper     = "#FFFFFF"
)

# Internal palettes. Qualitative is Okabe-Ito-inspired but led by the brand
# hues; every fixed-meaning palette avoids red/green confusion.
.logiccomm_palettes <- list(
  qualitative = c("#1F6F78", "#E8A33D", "#33415C", "#D1495B", "#5AA9A3",
                  "#9B6FB0", "#7F8C8D", "#C9A227", "#2E86AB", "#B05E2A"),
  roles = c("Sender" = "#D1495B", "Receiver" = "#1F6F78",
            "Mediator" = "#E8A33D", "Influencer" = "#33415C"),
  balance = c(low = "#1F6F78", mid = "#F4F4F2", high = "#D1495B"),
  significance = c("Up" = "#D1495B", "Down" = "#2E86AB", "n.s." = "grey78"),
  # Ordered evidence tiers: darker brand teal = stronger; "broad / non-specific"
  # is a warm-neutral "set aside" colour; weakest is a light grey.
  tier = c(
    "Tier 1: strong, specific, supported" = "#10403B",
    "Tier 2: strong candidate"            = "#1F6F78",
    "Tier 3: emerging candidate"          = "#5AA9A3",
    "Tier 3: broad / non-specific"        = "#C9A227",
    "Tier 4: weak / context-dependent"    = "#BFC4CB"
  )
)

#' LogicComm Publication Theme
#'
#' A clean, export-ready ggplot2 theme used by every LogicComm figure.
#'
#' @param base_size Base font size in points (default 12).
#' @param base_family Base font family ("" uses the device default).
#' @param grid Which major gridlines to keep: "xy", "x", "y", or "none".
#' @param legend Legend position passed to \code{ggplot2::theme()}.
#' @return A ggplot2 theme object.
#' @export
theme_logiccomm <- function(base_size = 12, base_family = "",
                            grid = c("xy", "x", "y", "none"),
                            legend = "right") {
  grid <- match.arg(grid)
  th <- ggplot2::theme_minimal(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      text = ggplot2::element_text(color = logiccomm_brand$ink),
      plot.title = ggplot2::element_text(face = "bold", size = ggplot2::rel(1.25),
                                         margin = ggplot2::margin(b = 6)),
      plot.subtitle = ggplot2::element_text(color = "grey35", size = ggplot2::rel(0.95),
                                            margin = ggplot2::margin(b = 10)),
      plot.caption = ggplot2::element_text(color = "grey45", size = ggplot2::rel(0.78),
                                           hjust = 0, margin = ggplot2::margin(t = 8)),
      axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 5)),
      axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 5)),
      axis.text = ggplot2::element_text(color = "grey25"),
      panel.grid.major = ggplot2::element_line(color = "grey92", linewidth = 0.35),
      panel.grid.minor = ggplot2::element_blank(),
      legend.title = ggplot2::element_text(face = "bold", size = ggplot2::rel(0.9)),
      legend.text = ggplot2::element_text(size = ggplot2::rel(0.85)),
      legend.key.size = grid::unit(0.9, "lines"),
      legend.background = ggplot2::element_blank(),
      legend.position = legend,
      strip.background = ggplot2::element_rect(fill = "grey94", color = NA),
      strip.text = ggplot2::element_text(face = "bold", size = ggplot2::rel(0.85),
                                         margin = ggplot2::margin(4, 4, 4, 4)),
      plot.title.position = "plot",
      plot.caption.position = "plot"
    )
  if (grid == "x")    th <- th + ggplot2::theme(panel.grid.major.y = ggplot2::element_blank())
  if (grid == "y")    th <- th + ggplot2::theme(panel.grid.major.x = ggplot2::element_blank())
  if (grid == "none") th <- th + ggplot2::theme(panel.grid.major = ggplot2::element_blank())
  th
}

# ---- Colour / fill scales -------------------------------------------------

#' LogicComm Colour and Fill Scales
#'
#' Brand, colorblind-safe scales for LogicComm figures: discrete (\code{_d},
#' qualitative brand palette), continuous sequential (\code{_c}, the brand
#' "mako" ramp), diverging (teal-amber through near-white), and a fixed ordered
#' scale for \code{evidence_tier}.
#'
#' @param ... Passed to the underlying ggplot2 scale.
#' @param na.value Colour for missing values.
#' @param option,begin,end,direction Viridis arguments for the sequential scale.
#' @param midpoint Midpoint for the diverging scale.
#' @param name Legend title.
#' @return A ggplot2 scale.
#' @name logiccomm_scales
NULL

#' @rdname logiccomm_scales
#' @export
scale_color_logiccomm_d <- function(..., na.value = "grey80") {
  ggplot2::scale_color_manual(..., values = unname(.logiccomm_palettes$qualitative),
                              na.value = na.value)
}
#' @rdname logiccomm_scales
#' @export
scale_fill_logiccomm_d <- function(..., na.value = "grey80") {
  ggplot2::scale_fill_manual(..., values = unname(.logiccomm_palettes$qualitative),
                             na.value = na.value)
}
#' @rdname logiccomm_scales
#' @export
scale_color_logiccomm_c <- function(..., option = "mako", begin = 0.08, end = 0.94,
                                    direction = -1) {
  ggplot2::scale_color_viridis_c(..., option = option, begin = begin, end = end,
                                 direction = direction)
}
#' @rdname logiccomm_scales
#' @export
scale_fill_logiccomm_c <- function(..., option = "mako", begin = 0.08, end = 0.94,
                                   direction = -1) {
  ggplot2::scale_fill_viridis_c(..., option = option, begin = begin, end = end,
                                direction = direction)
}
#' @rdname logiccomm_scales
#' @export
scale_color_logiccomm_diverging <- function(..., midpoint = 0) {
  ggplot2::scale_color_gradient2(..., low = logiccomm_brand$primary, mid = "#F4F4F2",
                                 high = logiccomm_brand$accent, midpoint = midpoint)
}
#' @rdname logiccomm_scales
#' @export
scale_fill_logiccomm_diverging <- function(..., midpoint = 0) {
  ggplot2::scale_fill_gradient2(..., low = logiccomm_brand$primary, mid = "#F4F4F2",
                                high = logiccomm_brand$accent, midpoint = midpoint)
}
#' @rdname logiccomm_scales
#' @export
scale_color_tier <- function(..., name = "Evidence tier") {
  ggplot2::scale_color_manual(..., name = name, values = .logiccomm_palettes$tier,
                              na.value = "grey85")
}
#' @rdname logiccomm_scales
#' @export
scale_fill_tier <- function(..., name = "Evidence tier") {
  ggplot2::scale_fill_manual(..., name = name, values = .logiccomm_palettes$tier,
                             na.value = "grey85")
}

# Order an evidence_tier vector as an ordered factor (strong -> weak) so legends
# and colours read in the right order regardless of input ordering.
.tier_factor <- function(x) {
  lv <- names(.logiccomm_palettes$tier)
  extra <- setdiff(unique(as.character(x[!is.na(x)])), lv)
  factor(as.character(x), levels = c(lv, extra))
}

# ---- Export helper --------------------------------------------------------

#' Save a LogicComm Figure for Publication
#'
#' Thin \code{ggplot2::ggsave()} wrapper with journal-friendly defaults: physical
#' widths in millimetres (single-column 89 mm, double-column 183 mm), high DPI,
#' and a white background. Use a \code{.pdf} or \code{.svg} file for vector output.
#'
#' For \code{.pdf} output the \code{cairo_pdf} device is used by default (when
#' available) so the Unicode glyphs in axis labels -- the arrow in
#' "Sender \\u2192 receiver" and the ellipsis in truncated names -- embed cleanly.
#' The base \code{pdf()} device drops those glyphs with a "conversion failure"
#' warning, so prefer this helper (or pass \code{device = grDevices::cairo_pdf})
#' for vector export.
#'
#' @param plot A ggplot (or patchwork) object.
#' @param file Output path; the extension selects the device (pdf/svg/png/tiff).
#' @param width One of "single" (89 mm), "double" (183 mm), or "custom".
#' @param width_mm Custom width in mm (required when \code{width = "custom"}).
#' @param height_mm Figure height in mm.
#' @param dpi Resolution for raster devices.
#' @param device Optional device override passed to \code{ggsave()}; defaults to
#'   \code{cairo_pdf} for \code{.pdf} when cairo is available.
#' @param ... Passed to \code{ggplot2::ggsave()}.
#' @return The output path, invisibly.
#' @export
save_logiccomm_figure <- function(plot, file,
                                   width = c("single", "double", "custom"),
                                   width_mm = NULL, height_mm = 110, dpi = 600,
                                   device = NULL, ...) {
  if (is.character(width)) {
    width <- match.arg(width)
    width_mm <- switch(width, single = 89, double = 183, custom = width_mm)
  }
  if (is.null(width_mm) || !is.finite(width_mm)) {
    stop("Provide a finite 'width_mm' when width = 'custom'.")
  }
  # Unicode-safe vector PDF: the base pdf() device cannot render the arrow/ellipsis
  # glyphs LogicComm uses in labels, so default to cairo_pdf when exporting a PDF.
  if (is.null(device) && identical(tolower(tools::file_ext(file)), "pdf") &&
      isTRUE(capabilities("cairo"))) {
    device <- grDevices::cairo_pdf
  }
  ggplot2::ggsave(filename = file, plot = plot, width = width_mm, height = height_mm,
                  units = "mm", dpi = dpi, bg = logiccomm_brand$paper, device = device, ...)
  invisible(file)
}
