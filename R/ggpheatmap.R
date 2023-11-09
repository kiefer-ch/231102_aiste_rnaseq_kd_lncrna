#' Make nice heatmaps
#'
#' `ggpheatmap()` makes nice heatmaps similar to pheatmap but using ggplot.
#'
#' @param df A tible with one row_id column and one column per sample.
#' @param row_id Name of the column holding the row ids.
#' @param n_clust Number of clusters for y axis.
#' @param cluster_rows Should rows be clustered?
#' @param cluster_cols Should columns be clustered?
#' @param annotation_col A tibble with two columns, "col_id" and "col_annotation".
#' @param base_size Base size for text.
#' @param row_palette Name of a color brewer palette.
#' @param plot Should a ggplot object or a dataframe be returned?
#' @param label One in c("none", "rows").
#'
#' @return A heatmap with dendrograms.
#'
#' @export
ggpheatmap <- function(df, row_id = "gene_id", n_clust = 2, cluster_rows = TRUE,
    cluster_cols = TRUE, annotation_col, base_size = 12, row_palette = "Dark2",
    plot = TRUE, label = "none", label_size = 8) {


    cluster_matrix <- df %>%
        as.data.frame() %>%
        tibble::column_to_rownames(row_id) %>%
        data.matrix()


    cluster_df <- cluster_matrix %>%
        t() %>%
        scale(center = TRUE, scale = TRUE) %>%
        t() %>%
        dplyr::as_tibble(rownames = "row_id") %>%
        tidyr::gather("col_id", "z", -row_id) %>%
        dplyr::mutate(col_id = as.factor(col_id)) %>%
        dplyr::mutate(row_id = as.factor(row_id))


    if (isTRUE(cluster_rows)) {

        clust_row <- cluster_matrix %>%
            t() %>%
            scale(center = TRUE, scale = TRUE) %>%
            t() %>%
            dist("euclidean") %>%
            hclust("complete")


        # create dendrogram for rows
        dend_row <- as.dendrogram(clust_row)

        cluster_df$row_id <- factor(cluster_df$row_id, levels = labels(dend_row))

        cluster_ids_row <- cutree(clust_row, k = n_clust)[labels(dend_row)] %>%
            tibble::enframe("row_id", "cluster_id") %>%
            dplyr::mutate(row_id = factor(row_id, levels = row_id))

        cluster_lookup <- 1:length(unique(cluster_ids_row$cluster_id))
        names(cluster_lookup) <- unique(cluster_ids_row$cluster_id)

        # reorder cluster labels to fit numerical order in plot
        cluster_ids_row <- cluster_ids_row %>%
            dplyr::mutate(cluster_id = cluster_lookup[as.character(cluster_id)]) %>%
            dplyr::mutate(cluster_id = as.factor(cluster_id))

        dend_row <- dend_row %>%
            ggtree::ggtree(lwd = .1, ladderize = FALSE) +
            ggplot2::theme_void(base_size = base_size, base_family = "Helvetica") +
            scale_y_reverse(expand = expansion(0, 0.6))

        cluster_df <- cluster_df %>%
            dplyr::left_join(cluster_ids_row, by = "row_id")
    }


    # row labeling
    label_row <- cluster_ids_row %>%
        ggplot2::ggplot(ggplot2::aes(1, row_id)) +
        ggplot2::geom_tile(aes(fill = cluster_id)) +
        ggplot2::scale_y_discrete(limits = rev) +
        ggplot2::scale_fill_brewer(type = "qual", palette = row_palette, name = "Cluster") +
        ggplot2::theme_void(base_size = base_size, base_family = "Helvetica")


    if (cluster_cols) {

        clust_col <- cluster_matrix %>%
            scale(center = TRUE, scale = TRUE) %>%
            t() %>%
            dist("euclidean") %>%
            hclust("complete")


        dend_col <- as.dendrogram(clust_col)


        cluster_df$col_id <- factor(cluster_df$col_id, levels = labels(dend_col))

        dend_col <- dend_col %>%
            ggtree::ggtree(lwd = .1, ladderize = FALSE) +
            ggplot2::scale_x_reverse() +
            ggplot2::coord_flip() +
            ggplot2::theme_void(base_size = base_size, base_family = "Helvetica")

    }


    if (plot) {

        # column labeling
        label_col <- annotation_col %>%
            dplyr::mutate(col_id = factor(col_id, levels = levels(cluster_df$col_id))) %>%
            ggplot2::ggplot(ggplot2::aes(x = col_id, y = 1)) +
            ggplot2::geom_tile(aes(fill = col_annotation)) +
            ggthemes::scale_fill_colorblind(name = NULL) +
            ggplot2::theme_void(base_size = base_size, base_family = "Helvetica")


        # actual heatmap
        hm <- cluster_df %>%
            ggplot2::ggplot(ggplot2::aes(col_id, row_id)) +
            ggplot2::geom_tile(aes(fill = z)) +
            ggplot2::scale_y_discrete(limits = rev) +
            ggplot2::scale_fill_distiller(type = "div", palette = "RdYlBu")


        # labeling of y axis
        if (label == "none") {

            hm <- hm +
                ggplot2::theme_void()

        } else if (label == "rows") {

            hm <- hm +
                ggplot2::theme_classic(base_size = base_size, base_family = "Helvetica") +
                theme(
                    axis.title.y = ggplot2::element_blank(),
                    axis.ticks.y = ggplot2::element_blank(),
                    axis.text.y = element_text(size = label_size),
                    axis.title.x = ggplot2::element_blank(),
                    axis.ticks.x = ggplot2::element_blank(),
                    axis.text.x = ggplot2::element_blank(),
                    axis.line = ggplot2::element_blank(),
                    plot.margin = ggplot2::unit(c(0, 0, 0, 0), "lines")
                )

        }


        if (cluster_cols) {

            pl <- patchwork::plot_spacer() + patchwork::plot_spacer() + dend_col +
                patchwork::plot_spacer() + patchwork::plot_spacer() + label_col +
                dend_row + label_row + hm +
                patchwork::plot_layout(ncol = 3,
                    widths = c(.2, .1, 1),
                    heights = c(.1, .1, 2),
                    guides = 'collect')

        } else {

            pl <- patchwork::plot_spacer() + patchwork::plot_spacer() + label_col +
                dend_row + label_row + hm +
                patchwork::plot_layout(ncol = 3,
                    widths = c(.2, .1, 1),
                    heights = c(.1, 2),
                    guides = 'collect')

        }

        return(pl)

    } else {

        return(cluster_df)

    }

}

