## ---- echo = FALSE------------------------------------------------------------
embed_png <- function(path, dpi = NULL) {
  meta <- attr(png::readPNG(path, native = TRUE, info = TRUE), "info")
  if (!is.null(dpi)) meta$dpi <- rep(dpi, 2)
  knitr::asis_output(paste0(
    "<img src='", path, "'",
    " width=", round(meta$dim[1] / (meta$dpi[1] / 96)),
    " height=", round(meta$dim[2] / (meta$dpi[2] / 96)),
    " />"
  ))
}
knitr::opts_chunk$set(comment = "#>", collapse = TRUE)

## ---- echo = FALSE------------------------------------------------------------
embed_png("Fiji-toolbar.png")

## ---- echo = FALSE------------------------------------------------------------
embed_png("open_hillshade_2012.png")

## ---- echo = FALSE------------------------------------------------------------
embed_png("open_hillshade_2017.png")

## ---- echo = FALSE------------------------------------------------------------
embed_png("bunwarpj_toolbar.png")

## ---- echo = FALSE------------------------------------------------------------
embed_png("landmark_all_hillshade_2012.png")

## ---- echo = FALSE------------------------------------------------------------
embed_png("landmark_all_hillshade_2017.png")

## ---- echo = FALSE------------------------------------------------------------
embed_png("bunwarpj_io_menu.png")

## ---- echo = FALSE------------------------------------------------------------
embed_png("bunwarpj_settings.png")

## ---- echo = FALSE------------------------------------------------------------
embed_png("convert_direct_to_raw_select_elastic.png")

## ---- echo = FALSE------------------------------------------------------------
embed_png("convert_direct_to_raw_output_filename.png")

