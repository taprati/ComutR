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
  variant_scheme,
  show_variant_legend = TRUE,
  col_maps,
  features_of_interest,
  text_annotation,
  barplot_data,
  grob = FALSE,
  body_width = 6,
  body_height = 5,
  add_borders = FALSE,
  show_barcodes = TRUE,
  id_order,
  ids
)
}
\arguments{
\item{data}{df with Tumor_Sample_Barcode, Hugo_Symbol, Variant_Classification, and optional column for text annotations}

\item{metadata}{df with Tumor_Sample_Barcode, and metadata columns}

\item{variant_colors}{colorscheme for alteration types, named list}

\item{variant_scheme}{naming scheme to convert alteration types. If provided should match variant_colors names}

\item{show_variant_legend}{whether the legend of variant types should be include. Default is TRUE.}

\item{col_maps}{named list of color maps. Names should match columns of metadata}

\item{features_of_interest}{optional vector of genes of interest}

\item{text_annotation}{How to annotate comut plot squares if desired, must be a column in data}

\item{barplot_data}{named list of named lists. Each sub list should contain data, colors, and legend params for plots}

\item{grob}{whether to return grob object instead of plotting. Useful for other frameworks.}

\item{body_width}{width of the heatmap body in inches. Defaults is 6}

\item{body_height}{height of the heatmap body in inches. Default is 5}

\item{add_borders}{whether padding should be added around each box in the comut body.}

\item{show_barcodes}{whether the sample ids should be shown in the plot}

\item{id_order}{optional vector with order of Tumor_Sample_Barcodes}

\item{ids}{optional vector of Tumor_Sample_Barcodes to show}
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
