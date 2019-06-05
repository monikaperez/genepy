#' Make volcano plot
#'
#' @param data data frame containing stats
#' @param effect_var effect size
#' @param p_var p value
#' @param q_var q value--when specified, defaults to highlighting points that pass q_thresh
#' @param q_thresh q value threshold for highlighting points
#' @param label_var column containing labels for data points
#' @param n_labeled if rank_by = 'effect', n_labeled points per side; if rank_by = 'pval', n_labeled points in total
#' @param label_size size of points
#' @param label_bool logical column indicating which points should be labeled
#' @param rank_by data type used to rank data points when labeling
#' @param ggrepel_type choose whether ggrepel's geom_text_repel or geom_label_repel should be used
#' @param color_var logical/categorical column for coloring points
#' @param color_values vector assigning categories from color_var to custom colors
#'
#' @importFrom magrittr "%>%"
#' @importFrom magrittr "%<>%"
#'
#' @export
#'
make_volcano <- function(data, effect_var, p_var, q_var = NULL, q_thresh = 0.1, label_var = NULL,
                         n_labeled = 10, label_size = 3, label_bool = NULL, rank_by = c('effect', 'pval'),
                         ggrepel_type = c('text', 'label'), color_var = NULL, color_values = NULL) {

  # set label for colors in the legend
  guide_title <- color_var
  # log 10 transform the p values
  transformed_p_var <- paste0('-log10(', p_var, ')')
  data[[transformed_p_var]] <- -log10(data[[p_var]])

  # if user has specified q values but no variable to color by, color points that pass q_thresh
  if (is.null(color_var) & !is.null(q_var)) {
    color_var <- 'internal_sig'
    data[[color_var]] <- data[[q_var]] < q_thresh
    guide_title <- sprintf('FDR < %.3f', q_thresh)
  }

  if (is.null(color_var)) {
    p <- data %>%
      ggplot(aes_string(effect_var, transformed_p_var)) +
      geom_point(color = '#333333', alpha = 0.7)
  } else {
    if (is.null(color_values)) { # user has not specified exact colors, create default color schemes
      if (is.logical(data[[color_var]])) color_values <- c(`TRUE` = '#BF2026', `FALSE` = '#4d4d4d')
      else color_values <- RColorBrewer::brewer.pal(length(unique(data[[color_var]])), 'Dark2')
    }
    # plot layers one by one
    layering <- sort(unique(as.character(data[[color_var]])))
    p <- ggplot(data, aes_string(effect_var, transformed_p_var, color = color_var))
    for (cur_layer in layering) p <- p + geom_point(data = data[data[[color_var]] == cur_layer,], alpha = 0.7)
    p <- p +
      scale_color_manual(values = color_values) +
      ggplot2::guides(color = guide_legend(title = guide_title))
  }

  if (!is.null(label_var)) { # user has specified column to label points
    if (is.null(label_bool)) { # default to labeling top 10 data points on each side by effect size
      label_bool <- 'to_label' # define points to label with this column
      if (rank_by[1] == 'effect') {
        left <- rank(data[[effect_var]], ties.method = 'random')
        right <- rank(-data[[effect_var]], ties.method = 'random')
        data[[label_bool]] <- (left <= n_labeled) | (right <= n_labeled)
      } else if (rank_by[1] == 'pval') {
        data[[label_bool]] <- rank(data[[p_var]], ties.method = 'random') <= n_labeled
      }
    }
    if (ggrepel_type[1] == 'text') p <- p +
        ggrepel::geom_text_repel(data = data[data[[label_bool]],], aes_string(label = label_var), size = label_size, show.legend = F)
    else if (ggrepel_type[1] == 'label') p <- p +
      ggrepel::geom_label_repel(data = data[data[[label_bool]],], aes_string(label = label_var), size = label_size, show.legend = F)
  }
  return(p + theme_Publication())
}



#' Make boxplot
#'
#' @param df: data frame
#' @param group: string name of variable to group by
#' @param yvar: string name of y variable
#' @param min_per_group: min samples per group [2]
#' @param group_order: how order groups. ('increasing' [default], 'decreasing', 'none')
#' @param show_points: whether or not to show individual data points [TRUE]
#'
#' @export
#'
make_group_boxplot <- function(df, group, yvar, min_per_group = 2,
                               group_order = 'increasing', show_points = TRUE) {
  stopifnot(group_order %in% c('increasing', 'decreasing', 'none'))
  gdf <- plyr::ddply(df, group, function(dd) {
    data.frame(avg = mean(dd[[yvar]], na.rm=T),
               n = sum(!is.na(dd[yvar])))
  })
  if (group_order != 'none') {
    gord <- gdf %>%
      dplyr::filter(n >= min_per_group)
    if (group_order == 'increasing') {
      gord %<>% dplyr::arrange(avg)
    } else if (group_order == 'decreasing') {
      gord %<>% dplyr::arrange(desc(avg))
    }
    gord %<>% .[[group]]
    df[[group]] <- factor(df[[group]], levels = gord)
  }
  df <- df[!is.na(df[[group]]),,drop=FALSE]

  if (show_points) {
    outlier_shape <- NA
  } else {
    outlier_shape <- 19
  }

  g <- ggplot2::ggplot(df, ggplot2::aes_string(group, yvar)) +
    ggplot2::geom_boxplot(outlier.shape = outlier_shape) +
    ggplot2::theme(axis.text.x = element_text(angle = 90, hjust = 1))
  if (show_points) {
    g <- g + ggplot2::geom_point(alpha = 0.75)
  }
  return(g)
}



#' Scatterplot comparing two variables with gene labels based on input subset of genes
#'
#' @param df
#' @param xlab
#' @param ylab
#' @param gene_field
#' @param top_genes
#' @param xvar
#' @param yvar
#' @param color_annot
#' @param lab_size
#' @param show_stats
#' @param make_plotly
#'
#' @importFrom magrittr "%>%"
#' @importFrom magrittr "%<>%"
#'
#' @export
#'
make_scatter_plot <- function(df, xvar, yvar, xlab = NULL, ylab = NULL,
                              gene_field = 'Gene', top_genes = 25,
                              lab_size = 2,
                              color_annot = NULL, show_stats = TRUE,
                              make_plotly = FALSE) {
  df$x1 <- df[[xvar]]; df$y1 <- df[[yvar]]
  lab_genes <- c(
    df %>%
      dplyr::arrange(desc(abs(x1))) %>%
      head(top_genes) %>%
      .[[gene_field]],
    df %>%
      dplyr::arrange(desc(abs(y1))) %>%
      head(top_genes) %>%
      .[[gene_field]]
  ) %>% unique()
  df$Gene <- df[[gene_field]]
  if (!is.null(color_annot)) {
    cur_aes <- ggplot2::aes_string(color = color_annot)
  } else {
    cur_aes <- NULL
  }
  g <- ggplot2::ggplot(df, ggplot2::aes_string(x=xvar, y=yvar, text = gene_field)) +
    ggplot2::geom_point(alpha = 0.75, size = 1, mapping = cur_aes) +
    ggplot2::geom_vline(xintercept = 0, linetype = 'dashed') +
    ggplot2::geom_hline(yintercept = 0, linetype = 'dashed') +
    ggrepel::geom_text_repel(data = df %>%
                               dplyr::filter(Gene %in% lab_genes),
                             ggplot2::aes(label = Gene), size = lab_size) +
    theme_Publication()
  if (!is.null(xlab)) {
    g <- g + ggplot2::xlab(xlab)
  }
  if (!is.null(ylab)) {
    g <- g + ggplot2::ylab(ylab)
  }
  if (!is.null(color_annot)) {
    g <- g + scale_colour_Publication()
  }
  if (show_stats) {
    g <- g + ggpubr::stat_cor()
  }
  if (make_plotly) {
    g <- plotly::ggplotly(g)
  }
  return(g)
}




#' Set publication theme on ggplot object
#'
#' @param base_size
#' @param base_family
#'
#' @return
#' @export
#'
#' @examples
#' https://rpubs.com/Koundy/71792
theme_Publication <- function(base_size=12, base_family="Helvetica") {
  library(grid)
  library(ggthemes)
  library(ggplot2)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(),
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.text = element_text(size = rel(1.2)),
            legend.direction = "horizontal",
            legend.key.size= unit(0.3, "cm"),
            legend.spacing = unit(0, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(5,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))

}

#' Set publication fill scheme
#'
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
scale_fill_Publication <- function(...){
  library(scales)
  discrete_scale("fill","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
}

#' Set publication color scheme
#'
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
scale_color_Publication <- function(...){
  library(scales)
  discrete_scale("colour","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
}


