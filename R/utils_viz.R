# R/utils_viz.R

#' Premium LogicComm Visual Theme
#'
#' @param base_size Base font size.
#' @param base_family Base font family.
#' @return A ggplot2 theme object.
#' @export
theme_logiccomm <- function(base_size = 11, base_family = "") {
  ggplot2::theme_minimal(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = ggplot2::rel(1.3), margin = ggplot2::margin(b = 10)),
      plot.subtitle = ggplot2::element_text(color = "grey30", size = ggplot2::rel(1.0), margin = ggplot2::margin(b = 10)),
      axis.title = ggplot2::element_text(face = "bold"),
      panel.grid.minor = ggplot2::element_blank(),
      legend.title = ggplot2::element_text(face = "bold"),
      legend.background = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "grey95", color = "white"),
      strip.text = ggplot2::element_text(face = "bold", size = ggplot2::rel(0.9))
    )
}

#' LogicComm Color Palettes
#' @keywords internal
.logiccomm_palettes <- list(
  range = c("juxtacrine" = "#E41A1C", "paracrine" = "#377EB8", "distal/endocrine" = "#4DAF4A", "unknown" = "grey80"),
  roles = c("Sender" = "#b2182b", "Receiver" = "#2166ac", "Mediator" = "#ef8a62", "Influencer" = "#67a9cf"),
  balance = c("low" = "#2166ac", "mid" = "#f7f7f7", "high" = "#b2182b")
)
