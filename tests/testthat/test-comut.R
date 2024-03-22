# Testing Comut

# Setup input data
dat <- tibble::tribble(
  ~ID, ~GeneA, ~GeneB, ~GeneC,
  1, 1, 1, 1,
  2, 0, 0, 0,
  3, 1, 0, 0,
  4, 1, 1, 0
)
input_maf <- tibble::tribble(
  ~Tumor_Sample_Barcode, ~Hugo_Symbol, ~Variant_Classification,
  "1", "A", "Missense_Mutation",
  "1", "B", "Nonsense_Mutation",
  "1", "C", "In_Frame_Del",
  "2", "C", "In_Frame_Del",
  #"3", "A", "Missense_Mutation",
  "4", "A", "Nonsense_Mutation",
  "4", "B", "Nonsense_Mutation"
  )
  meta <- tibble::tribble(
  ~Tumor_Sample_Barcode, ~Gender, ~Biopsy,
  "1", "Male", "Primary",
  "2", "Female", "Metastasis",
  "3", "Male", "Metastasis",
  #"4", "Female", "Primary"
)
signatures <- tibble::tribble(
  ~Tumor_Sample_Barcode, ~Aging, ~Other,
  "1", 0.5, 0.5,
  "2", 0.25, 0.75,
  "3", 0.75, 0.25,
  "4", 1, 0
)
total_snp_data <- tibble::tribble(
  ~Tumor_Sample_Barcode, ~n_snps,
  "1", 10,
  "2", 4,
  "3", 2,
  "4", 7
)
ids <- as.character(c(1,2,3,4))

# Color scheme creation
signature_colors <- c(
  "Aging" = "tan",
  "Other" = "orange"
)
gender_colors <- c(
  "Male" = "lightblue",
  "Female" = "pink"
)
biopsy_type_colors <- c(
  "Primary" = "lightgreen",
  "Metastasis" = "maroon"
)
color_maps <- list(
  "Gender" = gender_colors,
  "Biopsy" = biopsy_type_colors
)

# Barplot specifications
barplot_data <- list(
  "Mutation Signatures" = list(
    "data" = signatures,
    "colors" = signature_colors,
    "legend" = TRUE
    ),
  "N SNPs" = list(
    "data" = total_snp_data,
    "colors" = c("n_snps" = "cadetblue3"),
    "legend" = FALSE
    )
)

# Test
test_that("Comut plots generate", {
          expect_warning(
            comut(
              data = input_maf,
              metadata = meta,
              col_maps = color_maps,
              barplot_data = barplot_data
            ), "Following annotation names are duplicated: value")
          })

test_that("Comut plots reorder", {
          expect_warning(
            comut(
              data = input_maf,
              metadata = meta,
              col_maps = color_maps,
              id_order = c("1", "2", "3", "4"),
              barplot_data = barplot_data
            ), "Following annotation names are duplicated: value")
          })

test_that("Comut plots reorder", {
          expect_warning(
            comut(
              data = input_maf,
              metadata = meta,
              col_maps = color_maps,
              id_order = c("1", "2", "3", "4"),
              show_barcodes = FALSE,
              add_borders = TRUE,
              barplot_data = barplot_data
            ), "Following annotation names are duplicated: value")
          })
