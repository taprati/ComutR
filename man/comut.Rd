% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/comut.R
\name{comut}
\alias{comut}
\title{Comut Plot}
\usage{
comut(
  data,
  metadata,
  variant_colors,
  show_variant_legend = TRUE,
  col_maps,
  features_of_interest,
  text_annotation,
  text_annotation_col = "white",
  text_annotation_size = 8,
  barplot_data,
  grob = FALSE,
  body_width,
  body_height,
  body_border = FALSE,
  cell_height = 0.25,
  cell_width = 0.25,
  legend_side = "right",
  legend_fontsize = 10,
  bar_anno_height = 0.75,
  anno_fontsize = 8,
  add_borders = FALSE,
  show_barcodes = TRUE,
  id_order,
  ids,
  ...
)
}
\arguments{
\item{data}{df with Tumor_Sample_Barcode, Hugo_Symbol, Variant_Classification, and optional column for text annotations}

\item{metadata}{df with Tumor_Sample_Barcode, and metadata columns}

\item{variant_colors}{colorscheme for alteration types, named list}

\item{show_variant_legend}{whether the legend of variant types should be include. Default is TRUE.}

\item{col_maps}{named list of color maps. Names should match columns of metadata}

\item{features_of_interest}{optional vector of genes of interest}

\item{text_annotation}{How to annotate comut plot squares if desired, must be a column in data}

\item{text_annotation_col}{Color of text annotations, default is white.}

\item{text_annotation_size}{Size of text annotations, default is 8.}

\item{barplot_data}{named list of named lists. Each sub list should contain data, colors, and legend params for plots}

\item{grob}{whether to return grob object instead of plotting. Useful for other frameworks.}

\item{body_width}{width of the heatmap body in inches.}

\item{body_height}{height of the heatmap body in inches.}

\item{body_border}{boolean to add border to comut body.}

\item{cell_height}{height of the heatmap cells. Also applies to annotations.}

\item{cell_width}{width of the heatmap cells.}

\item{legend_side}{side of heatmap to put the legend on. Default is right}

\item{legend_fontsize}{Font size for legends. Default is 10}

\item{bar_anno_height}{Height of each barplot annotation in inches.}

\item{anno_fontsize}{Font size for annotations. Default is 8}

\item{add_borders}{whether padding should be added around each box in the comut body.}

\item{show_barcodes}{whether the sample ids should be shown in the plot}

\item{id_order}{optional vector with order of Tumor_Sample_Barcodes}

\item{ids}{optional vector of Tumor_Sample_Barcodes to show}

\item{...}{variable arguments to pass to Heatmap call}
}
\value{
Comut Plot
}
\description{
Comut Plot
}
\examples{

input_maf <- data.frame(
 Tumor_Sample_Barcode = c("1", "1", "1", "2", "3", "4", "4"),
 Hugo_Symbol = c("A", "B", "C", "C", "A", "A", "B"),
 Variant_Classification = c(
   "Missense_Mutation", "Nonsense_Mutation",
   "In_Frame_Del", "In_Frame_Del", "Missense_Mutation",
   "Nonsense_Mutation", "Nonsense_Mutation"))

ComutR::comut(data = input_maf)

}
