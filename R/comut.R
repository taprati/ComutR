#' Comut Plot
#'
#' @param data df with Tumor_Sample_Barcode, Hugo_Symbol, Variant_Classification, and optional column for text annotations
#' @param metadata df with Tumor_Sample_Barcode, and metadata columns
#' @param variant_colors colorscheme for alteration types, named list
#' @param show_variant_legend whether the legend of variant types should be include. Default is TRUE.
#' @param col_maps named list of color maps. Names should match columns of metadata
#' @param features_of_interest optional vector of genes of interest
#' @param text_annotation How to annotate comut plot squares if desired, must be a column in data
#' @param text_annotation_size Size of text annotations, default is 8.
#' @param text_annotation_col Color of text annotations, default is white.
#' @param barplot_data named list of named lists. Each sub list should contain data, colors, and legend params for plots
#' @param grob whether to return grob object instead of plotting. Useful for other frameworks.
#' @param body_width width of the heatmap body in inches.
#' @param body_height height of the heatmap body in inches.
#' @param body_border boolean to add border to comut body.
#' @param cell_width width of the heatmap cells.
#' @param cell_height height of the heatmap cells. Also applies to annotations.
#' @param legend_side side of heatmap to put the legend on. Default is right
#' @param legend_fontsize Font size for legends. Default is 10
#' @param anno_fontsize Font size for annotations. Default is 8
#' @param add_borders whether padding should be added around each box in the comut body.
#' @param show_barcodes whether the sample ids should be shown in the plot
#' @param ids optional vector of Tumor_Sample_Barcodes to show
#' @param id_order optional vector with order of Tumor_Sample_Barcodes
#' @param ... variable arguments to pass to Heatmap call
#'
#' @return Comut Plot
#' @export
#'
#' @examples
#'
#'input_maf <- data.frame(
#'  Tumor_Sample_Barcode = c("1", "1", "1", "2", "3", "4", "4"),
#'  Hugo_Symbol = c("A", "B", "C", "C", "A", "A", "B"),
#'  Variant_Classification = c(
#'    "Missense_Mutation", "Nonsense_Mutation",
#'    "In_Frame_Del", "In_Frame_Del", "Missense_Mutation",
#'    "Nonsense_Mutation", "Nonsense_Mutation"))
#'
#' ComutR::comut(data = input_maf)
#'
comut <-
  function(data,
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
           anno_fontsize = 8,
           add_borders = FALSE,
           show_barcodes = TRUE,
           id_order,
           ids, ...) {
    # Parse parameters
    if (missing(data)) stop("data is required!")
    if (missing(metadata)) { metadata = NULL }
    if (missing(col_maps)) { col_maps = NULL }
    if (missing(id_order)) { id_order = NULL }
    if (missing(body_width)) { body_width = NULL }
    if (missing(body_height)) { body_height = NULL }
    if (missing(features_of_interest)) { features_of_interest = NULL }
    if (missing(text_annotation)) { text_annotation = "none" }
    if (missing(barplot_data)) { barplot_data = NULL }
    if (missing(ids)) {
      # message("No ids specified. Defaulting to all ids!")
      # Default to the union of all ids given across all data
      bp_ids <- c()
      meta_ids <- c()
      if (!is.null(barplot_data)) {
        for (b in names(barplot_data)) {
          temp_ids <- barplot_data[[b]][["data"]]$Tumor_Sample_Barcode
          bp_ids <- c(bp_ids, temp_ids)
        }
      }
      if (!is.null(metadata)) {
        meta_ids <- unique(metadata$Tumor_Sample_Barcode)
      }
      maf_ids <- unique(data$Tumor_Sample_Barcode)
      ids <- union(maf_ids, union(meta_ids, bp_ids))
    }
    # Check for matching col_maps and metadata features
    if (
      (is.null(col_maps) & !is.null(metadata)) |
      (!is.null(col_maps) & is.null(metadata))
    ) {
      stop("If metadata provided, col_map needs to be defined! \n
            If col_map is defined, metadata should be provided!")
    } else {
      feats <- colnames(metadata)[colnames(metadata) != "Tumor_Sample_Barcode"]
      if (!all(feats %in% names(col_maps))) {
        stop("Some metadata features are not in the color maps!")
      }
    }

    # Filter the data to patients of interest, and format alterations
    filtered_data <- data %>%
      dplyr::filter(Tumor_Sample_Barcode %in% ids)

    # Filter to features of interest if needed
    if (!is.null(features_of_interest)) {
      filtered_data <- filtered_data %>%
        dplyr::filter(Hugo_Symbol %in% features_of_interest)

      # Create dummy rows for genes not mutated at all
      dummy_genes <- features_of_interest[!features_of_interest %in% filtered_data$Hugo_Symbol]
      dummy_matrix <-
        matrix(
          data = NA,
          ncol = length(ids),
          nrow = length(dummy_genes)
        )
      rownames(dummy_matrix) <- dummy_genes
      mutation_dummy_matrix <-
        matrix(
          data = 0,
          ncol = length(ids),
          nrow = length(dummy_genes)
        )
      rownames(mutation_dummy_matrix) <- dummy_genes
    } else {
      dummy_matrix <- NULL
      mutation_dummy_matrix <- NULL
    }

    # Set variant color map or check provided one
    # By default assigns colors in order of variant frequency
    if (missing(variant_colors)) {
      colors <- c("#24796c", "#e58606", "#5d69b1", "#52bca3", "#99c945",
                  "#cc61b0", "#daa51b", "#2f8ac4", "#764e9f", "#ed645a",
                  "#cc3a8e")
      var_class_values <- names(sort(table(filtered_data[["Variant_Classification"]]), decreasing = TRUE))

      if (length(var_class_values) > length(colors)) {
        stop("More variants than colors! Custom variant_colors needed!")
      }
      variant_colors <- setNames(colors[1:length(var_class_values)], var_class_values)
    } else {
      var_class_values <- names(sort(table(filtered_data[["Variant_Classification"]]), decreasing = TRUE))
      if (!all(var_class_values %in% names(variant_colors))) {
        stop("Provided variant_colors does not contain every variant classification!")
      }
    }
    # Add a color for multiple variants if needed
    if (any(table(filtered_data[c("Tumor_Sample_Barcode", "Variant_Classification")]) > 2)) {
      variant_colors["Multiple"] <- "#a5aa99"
    }

    # Create dummy columns for ids with no alterations in genes
    ids_with_no_alterations <- ids[!ids %in% unique(filtered_data$Tumor_Sample_Barcode)]
    if (length(ids_with_no_alterations) != 0) {
      present_features <- unique(filtered_data$Hugo_Symbol)
      dummy_cols <- matrix(
        data = NA,
        ncol = length(ids_with_no_alterations),
        nrow = length(present_features)
      )
      colnames(dummy_cols) <- ids_with_no_alterations
      mutation_dummy_cols <- matrix(
        data = 0,
        ncol = length(ids_with_no_alterations),
        nrow = length(present_features)
      )
      colnames(mutation_dummy_cols) <- ids_with_no_alterations
    } else {
      dummy_cols <- NULL
      mutation_dummy_cols <- NULL
    }

    # Create various matrices for text annotations
    if (text_annotation != "none") {
      text_matrix <- filtered_data %>%
        dplyr::select(Tumor_Sample_Barcode, Hugo_Symbol, tidyselect::all_of(text_annotation)) %>%
        tidyr::pivot_wider(
          names_from = Hugo_Symbol,
          values_from = text_annotation,
          values_fill = NA,
          values_fn = function(value) {
            paste(as.character(value), collapse = "\n")
          }
        ) %>%
        tibble::column_to_rownames(var = "Tumor_Sample_Barcode") %>%
        as.matrix() %>% t() %>%
        cbind(dummy_cols) %>%
        rbind(dummy_matrix)
    } else {
      text_matrix <- NULL
    }

    # Create main variant and alteration matrix
    variant_matrix <- filtered_data %>%
      dplyr::select(Tumor_Sample_Barcode, Hugo_Symbol, Variant_Classification) %>%
      tidyr::pivot_wider(
        names_from = Hugo_Symbol,
        values_from = Variant_Classification,
        values_fill = NA,
        values_fn = function(value) {
          paste(value, collapse = ";")
        }
      ) %>%
      tibble::column_to_rownames(var = "Tumor_Sample_Barcode") %>% as.matrix() %>% t() %>%
      cbind(dummy_cols) %>%
      rbind(dummy_matrix)

    alteration_matrix <- filtered_data %>%
      dplyr::select(Tumor_Sample_Barcode, Hugo_Symbol) %>%
      dplyr::mutate(has_mut = 1) %>%
      tidyr::pivot_wider(
        names_from = Hugo_Symbol,
        values_from = has_mut,
        values_fill = as.integer(0),
        values_fn = function(value) {
          return(as.integer(1))
        }
      ) %>%
      tibble::column_to_rownames(var = "Tumor_Sample_Barcode") %>%
      as.matrix() %>%
      t() %>%
      cbind(mutation_dummy_cols) %>%
      rbind(mutation_dummy_matrix)

    # Adjust column order if desired
    if (!is.null(id_order)) {
      alteration_matrix <- alteration_matrix[ , id_order]
      variant_matrix <- variant_matrix[ , id_order]
      # Adjust text matrix on the fly
      if (!is.null(text_matrix)) {
        text_matrix <- text_matrix[ , id_order]
      }
    }

    # Get heatmap order
    heatmap_gene_order <- alteration_matrix %>%
      t() %>%
      colSums() %>%
      sort(decreasing = TRUE) %>%
      names()

    # Initialize lists for legends and annotations
    all_annotations <- NULL
    all_lgds <- list()

    ### Build barplots
    # Filter barplot data to patients of interest and match orders
    for (b in names(barplot_data)) {
      barplot_data[[b]][["data"]] <- barplot_data[[b]][["data"]] %>%
        dplyr::filter(Tumor_Sample_Barcode %in% ids) %>%
        dplyr::arrange(match(Tumor_Sample_Barcode,
                      colnames(alteration_matrix))) %>%
        tibble::column_to_rownames(var = "Tumor_Sample_Barcode")
    }
    # Make heatmap barplots
    for (feature in names(barplot_data)) {
      dat <- barplot_data[[feature]][["data"]]

      bar_anno <- ComplexHeatmap::HeatmapAnnotation(
        name = feature,
        value = ComplexHeatmap::anno_barplot(
          dat,
          bar_width = 1,
          border = FALSE,
          gp = grid::gpar(fill = barplot_data[[feature]][["colors"]][colnames(dat)],
                    fontsize = anno_fontsize,
                    col = "white")
        ),
        # Split multi-word labels across lines
        annotation_label = paste(stringr::str_replace(feature, " ", "\n")),
        annotation_name_gp = grid::gpar(fontsize = anno_fontsize),
        annotation_height = grid::unit(0.75, "inches"),
        gp = grid::gpar(fontsize = anno_fontsize, col = "white"),
        # gap = grid::unit(5, "mm"),
        which = "column",
        annotation_name_side = "left",
        border = FALSE,
        show_legend = FALSE
      )
      # Rename the annotation to avoid collisions
      names(bar_anno) <- feature
      # Add barplots to annotation list
      all_annotations <- c(bar_anno, all_annotations)
      if (barplot_data[[feature]][["legend"]]) {
        bar_lgd <- ComplexHeatmap::Legend(
          labels = names(barplot_data[[feature]][["colors"]]),
          title = paste(feature),
          grid_width = grid::unit(0.5, "cm"),
          grid_height = grid::unit(0.5, "cm"),
          legend_gp = grid::gpar(fill = barplot_data[[feature]][["colors"]],
                           fontsize = legend_fontsize),
          labels_gp = grid::gpar(fontsize = legend_fontsize),
          title_gp = grid::gpar(fontsize = legend_fontsize, fontface = "bold")
        )
        all_lgds[[feature]] <- bar_lgd
      }
    }
    # Reverse legends to match annotation ordering visually
    all_lgds <- rev(all_lgds)

    # Add gap between bar annotations
    if (!is.null(all_annotations)) {
      all_annotations@gap = grid::unit(rep(2, length(all_annotations)), "mm")
    }

    ### Metadata annotations
    # Filter metadata and format
    if (!is.null(metadata)) {
      # Create placeholder metadata for missing patients
      ids_with_no_meta <- ids[!ids %in% unique(metadata$Tumor_Sample_Barcode)]
      if (length(ids_with_no_meta) != 0) {
        meta_cols <- colnames(metadata)
        missing_meta <- data.frame(
          Tumor_Sample_Barcode = ids_with_no_meta,
          matrix(
            ncol = length(meta_cols) - 1,
            nrow = length(ids_with_no_meta)
          )
        )
        colnames(missing_meta) <- meta_cols
      } else {
        missing_meta <- NULL
      }
      annodata <- metadata %>%
        rbind(missing_meta) %>%
        dplyr::filter(Tumor_Sample_Barcode %in% ids) %>%
        dplyr::arrange(match(Tumor_Sample_Barcode,
                      colnames(alteration_matrix))) %>%
        tibble::column_to_rownames(var = "Tumor_Sample_Barcode")

      # Filter colorschemes to values in data
      for (c in names(col_maps)) {
        if (c %in% colnames(annodata)) {
          values_in_data <- annodata %>%
            dplyr::pull(paste(c)) %>% unique() %>% sort()
          col_maps[[c]] <- col_maps[[c]][values_in_data]
        }
      }

      # Create legends and annotations
      # TODO: check for barplot data and handle errors
      meta_anno <-
        ComplexHeatmap::HeatmapAnnotation(
          df = annodata,
          gp = grid::gpar(fontsize = anno_fontsize, col = "white", lwd = 1),
          which = "column",
          col = col_maps,
          annotation_name_side = "left",
          annotation_name_gp = grid::gpar(fontsize = anno_fontsize),
          simple_anno_size = grid::unit(cell_height, "inches"),
          border = FALSE,
          gap = grid::unit(2, "points"),
          show_legend = FALSE
        )
      if (is.null(all_annotations)) { # If no barplot annotations, the data types will cause errors
        all_annotations <- meta_anno
      } else {
        all_annotations <- c(all_annotations, meta_anno)
      }

      # Create legends
      features <- colnames(annodata)
      for (feature in features) {
        all_lgds[[feature]] <-
          ComplexHeatmap::Legend(
            labels = names(col_maps[[feature]]),
            title = paste(feature),
            grid_width = grid::unit(0.5, "cm"),
            grid_height = grid::unit(0.5, "cm"),
            legend_gp = grid::gpar(fill = col_maps[[feature]], fontsize = legend_fontsize),
            labels_gp = grid::gpar(fontsize = legend_fontsize),
            title_gp = grid::gpar(fontsize = legend_fontsize, fontface = "bold")
          )
      }
    }

    # Add legend of alteration type if desired
    if (show_variant_legend) {
      all_lgds[["Alteration Type"]] <- ComplexHeatmap::Legend(
        labels = names(variant_colors),
        title = "Alteration Type",
        # by_row = FALSE,
        grid_width = grid::unit(0.5, "cm"),
        grid_height = grid::unit(0.5, "cm"),
        legend_gp = grid::gpar(fill = variant_colors, fontsize = legend_fontsize),
        labels_gp = grid::gpar(fill = variant_colors, fontsize = legend_fontsize),
        title_gp = grid::gpar(fill = variant_colors, fontsize = legend_fontsize, fontface = "bold")
      )
    }
    # Define cell dimensions
    # Set so cell height is the same as annotation heights by default
    if (is.null(body_width) & is.null(body_height)) {
      body_width = ncol(alteration_matrix) * grid::unit(cell_width, "inches")
      body_height = nrow(alteration_matrix) * grid::unit(cell_height, "inches")
    } else {
      body_width = grid::unit(body_width, "inches")
      body_height = grid::unit(body_height, "inches")
    }

    # Hack to set the base colormap
    # If the variants are homogeneous, the colormap needs to be a single item
    if (length(unique(as.vector(alteration_matrix))) == 1) {
      col <- c("white")
    } else {
      col <- c("white", "white")
    }

    # Create heatmap
    comut <- ComplexHeatmap::Heatmap(
      matrix = alteration_matrix,
      top_annotation = all_annotations,
      show_column_dend = FALSE,
      show_heatmap_legend = FALSE,
      show_column_names = show_barcodes,
      width = body_width,
      height = body_height,
      border = body_border,
      column_names_gp = grid::gpar(fontsize = 10),
      row_names_gp = grid::gpar(fontsize = 10),
      row_names_side = "left",
      row_order = heatmap_gene_order,
      column_order = id_order,
      col = col,
      cell_fun = function(j, i, x, y, width, height, fill) {
        if (alteration_matrix[i, j] != 0) {
          alteration <- sort(stringr::str_split_1(variant_matrix[i, j], ";"))
          color <- variant_colors[alteration]
          if (length(alteration) == 1) {
            grid::grid.rect(x, y, width, height,
              just = "center",
              gp = grid::gpar(
                col = "white",
                fill = color,
                lwd = 1
              )
            )
          } else if (length(alteration) == 2) {
            first <- alteration[[1]]
            second <- alteration[[2]]
            first_color <- variant_colors[first]
            second_color <- variant_colors[second]
            grid::grid.polygon(
              grid::unit.c(x + 0.5 * width, x + 0.5 * width, x - 0.5 * width),
              grid::unit.c(y + 0.5 * height, y - 0.5 * height, y + 0.5 * height),
              gp = grid::gpar(
                fill = first_color,
                col = "white",
                lwd = 1
              )
            )
            grid::grid.polygon(
              grid::unit.c(x - 0.5 * width, x - 0.5 * width, x + 0.5 * width),
              grid::unit.c(y - 0.5 * height, y + 0.5 * height, y - 0.5 * height),
              gp = grid::gpar(
                fill = second_color,
                col = "white",
                lwd = 1
              )
            )
          } else if (length(alteration) > 2) {
            grid::grid.rect(x, y, width, height,
              just = "center",
              gp = grid::gpar(
                col = "white",
                fill = variant_colors["Multiple"],
                lwd = 1
              )
            )
          }
          if (add_borders == TRUE) {
            grid::grid.rect(x, y, width, height,
              just = "center",
              gp = grid::gpar(
                col = "white",
                lwd = 5,
                fill = "transparent"
              )
            )
          }
          if (length(alteration) %in% c(1, 2)) {
            if (text_annotation != "none") {
              grid::grid.text(sprintf("%s", text_matrix[i, j]),
                        x, y,
                        gp = grid::gpar(fontsize = text_annotation_size,
                                        col = text_annotation_col))
            }
          }
        }
      },
      ...
    )

    # Return grob object or draw plot
    if (grob) {
      p <- grid::grid.grabExpr(
        ComplexHeatmap::draw(comut,
                             annotation_legend_side = legend_side,
                             annotation_legend_list = all_lgds))
      return(p)
    } else {
      ComplexHeatmap::draw(comut,
                           annotation_legend_side = legend_side,
                           annotation_legend_list = all_lgds)
    }
  }
