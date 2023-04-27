library(compiler)
library(ggplot2)
compiler::setCompilerOptions(suppressAll = TRUE)
monocle_theme_opts = function ()
{
  ggplot2::theme(strip.background = ggplot2::element_rect(colour = "white",
                                        fill = "white")) + ggplot2::theme(panel.border = ggplot2::element_blank()) +
    ggplot2::theme(axis.line.x = ggplot2::element_line(size = 0.25, color = "black")) +
    ggplot2::theme(axis.line.y = ggplot2::element_line(size = 0.25, color = "black")) +
   # theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) +
    #theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()) +
    ggplot2::theme(panel.background = ggplot2::element_rect(fill = "white")) +
    ggplot2::theme(legend.key = ggplot2::element_blank())
}


table.ramp <- function(n, mid = 0.5, sill = 0.5, base = 1, height = 1)
{
    x <- seq(0, 1, length.out = n)
    y <- rep(0, length(x))
    sill.min <- max(c(1, round((n - 1) * (mid - sill / 2)) + 1))
    sill.max <- min(c(n, round((n - 1) * (mid + sill / 2)) + 1))
    y[sill.min:sill.max] <- 1
    base.min <- round((n - 1) * (mid - base / 2)) + 1
    base.max <- round((n - 1) * (mid + base / 2)) + 1
    xi <- base.min:sill.min
    yi <- seq(0, 1, length.out = length(xi))
    i <- which(xi > 0 & xi <= n)
    y[xi[i]] <- yi[i]
    xi <- sill.max:base.max
    yi <- seq(1, 0, length.out = length(xi))
    i <- which(xi > 0 & xi <= n)
    y[xi[i]] <- yi[i]
    height * y
}





#' Plot a curve of a single gene over pseudotime on different branches.
#'
#' @param df_prepare A list, unfitted gene expression data for two branches.
#' @param df1 A dataframe, gene expression data for fate 1.
#' @param df2 A dataframe, gene expression data for fate 2.
#' @param gene Gene name.
#' @param color_by Color by
#' @param cell_size Cell size
#' @param min_expr Min expr
#'
#' @return ggplot
#' @export
#' @import dplyr

plot_SingleGene_trend = function(df_prepare,df1,df2,gene,color_by = 'State',cell_size = 1.5,min_expr = NULL){

  cds1 = df_prepare[[1]][,c('State',gene)]
  cds1$Pseudotime = df1$Pseudotime
  cds1$exp = df1[,gene]  # Obtain fitted gene expression data
  cds1$branch = rep('fate1',nrow(cds1))
  cds1$cell_id = rownames(df_prepare[[1]])

  cds2 = df_prepare[[2]][,c('State',gene)]
  cds2$Pseudotime = df2$Pseudotime
  cds2$exp = df2[,gene]  # Obtain fitted gene expression data
  cds2$branch = rep('fate2',nrow(cds2))
  cds2$cell_id = rownames(df_prepare[[2]])

  cds = dplyr::bind_rows(cds1,cds2)
  colnames(cds)[2] = 'expression'
  cds$expression = log1p(cds$expression)
  cds$exp = log1p(cds$exp)


  q <- ggplot2::ggplot(ggplot2::aes(cds$Pseudotime, cds$expression), data = cds)

  if (is.null(color_by) == FALSE) {
    q <- q + ggplot2::geom_point(ggplot2::aes_string(color = color_by), size = I(cell_size))
  }
  q <- q + ggplot2::geom_line(ggplot2::aes_string(x = "Pseudotime", y = "exp",
                                linetype = "branch"), data = cds, size = 1)
  q <- q + ggplot2::ylab("Log(expression + 1)") + ggplot2::xlab("Pseudotime")
  q <- q + monocle_theme_opts()
  q + ggplot2::expand_limits(y = min_expr)+
    ggplot2::theme(panel.grid.major = ggplot2::element_line(linetype = 2, color = 'gray'))
}


#' Draw expression curves for two genes.
#'
#' @param df_prepare A list, unfitted gene expression data for two branches.
#' @param df1 A dataframe, gene expression data for fate 1.
#' @param df2 A dataframe, gene expression data for fate 2.
#' @param double_gene Gene vector
#' @param color_by Color by
#' @param cell_size cell size
#' @param min_expr Min expr
#'
#' @return ggplot
#' @export
plot_DoubleGenes_trend = function(df_prepare,df1,df2,double_gene,color_by = 'State',cell_size = 1.5,min_expr = NULL){
  # branch 1
  tmp1 = data.frame(df_prepare[[1]][,c('State')])
  colnames(tmp1)[1] = 'State'
  tmp1 = dplyr::mutate(tmp1,
                branch=rep('fate1',nrow(tmp1)),

                g1_full_model_expectation = df1[,double_gene[1]],  # Fitted value
                g2_full_model_expectation = df1[,double_gene[2]],

                g1_expression = df_prepare[[1]][,c(double_gene[1])],  # Unfitted values
                g2_expression = df_prepare[[1]][,c(double_gene[2])]
  )

  # branch 2
  tmp2 = data.frame(df_prepare[[2]][,c('State')])
  colnames(tmp2)[1] = 'State'
  tmp2 = dplyr::mutate(tmp2,
                branch=rep('fate2',nrow(tmp2)),

                g1_full_model_expectation = df2[,double_gene[1]],  # Fitted value
                g2_full_model_expectation = df2[,double_gene[2]],

                g1_expression = df_prepare[[2]][,c(double_gene[1])],  # Unfitted values
                g2_expression = df_prepare[[2]][,c(double_gene[2])]
  )

  cds = rbind(tmp1,tmp2)
  cds$g1_expression = log1p(cds$g1_expression)
  cds$g2_expression = log1p(cds$g2_expression)
  cds$g1_full_model_expectation = log1p(cds$g1_full_model_expectation)
  cds$g2_full_model_expectation = log1p(cds$g2_full_model_expectation)

  q <- ggplot2::ggplot(ggplot2::aes(x=cds$g1_expression, y=cds$g2_expression), data = cds)
  if (is.null(color_by) == FALSE) {
    q <- q + ggplot2::geom_point(ggplot2::aes_string(color = color_by), size = I(cell_size))
  }

  q <- q + ggplot2::geom_path(ggplot2::aes_string(x = 'g1_full_model_expectation', y = 'g2_full_model_expectation',
                                linetype = "branch"), data = cds, size = 1)

  q <- q + ggplot2::ylab(double_gene[2]) + ggplot2::xlab(double_gene[1])
  q <- q + monocle_theme_opts()
  q + ggplot2::expand_limits(y = min_expr)+
    ggplot2::theme(panel.grid.major = ggplot2::element_line(linetype = 2, color = 'gray'))

}

#' Draw a heat map of gene expression.
#'
#' @param df1 A dataframe, gene expression data for fate 1.
#' @param df2 A dataframe, gene expression data for fate 2.
#' @param G1_1 Gene vector.
#' @param branch_labels A vector.
#' @param branch_colors A vector
#' @param use_gene_short_name T/F
#' @param show_rownames T/F
#' @param hclust_method hclust method
#' @param num_clusters number of clusters
#' @param scale_max scale max
#' @param scale_min scale min
#'
#' @return ggplot
#' @export
plot_Genes_Heatmap = function(
    df1,df2,
    G1_1,
    branch_labels = c("fate1", "fate2"),
    branch_colors = c('#979797', '#F05662', '#7990C8'),
    use_gene_short_name = TRUE,
    show_rownames = T,
    hclust_method = "ward.D2",
    num_clusters = 6,
    scale_max=3,
    scale_min=-3){


  # branch 1
  tmp1 = df1[order(df1$Pseudotime,decreasing = TRUE),] # Reverse order
  # tmp1 = tmp1[,c(G1_1,G2_1)]
  tmp1 = tmp1[,G1_1]

  tmp1 = tmp1[as.integer(seq(1,nrow(tmp1),length.out=100)),] # sampling
  # branch 2
  tmp2 = df2[order(df2$Pseudotime),] # positive sequence
  # tmp2 = tmp2[,c(G1_1,G2_1)]
  tmp2 = tmp2[,G1_1]
  tmp2 = tmp2[as.integer(seq(1,nrow(tmp2),length.out=100)),]

  cds = rbind(tmp1,tmp2)
  cds = t(cds)

  cds=cds[!apply(cds, 1, stats::sd)==0,]
  cds=Matrix::t(scale(Matrix::t(cds),center=TRUE))
  cds=cds[is.na(row.names(cds)) == FALSE,]
  cds[is.nan(cds)] = 0
  cds[cds>scale_max] = scale_max
  cds[cds<scale_min] = scale_min

  col_gap_ind <- 101
  heatmap_matrix_ori <- cds
  cds <- cds[is.finite(cds[, 1]) & is.finite(cds[, col_gap_ind]), ]


  row_dist <- stats::as.dist((1 - stats::cor(Matrix::t(cds)))/2)
  row_dist[is.na(row_dist)] <- 1

  exp_rng <- range(cds) #bks is based on the expression range
  bks <- seq(exp_rng[1] - 0.1, exp_rng[2] + 0.1, by=0.1)

  ph = pheatmap::pheatmap(cds,
                color=c(grDevices::colorRampPalette(colors = c("blue","white"))(length(bks)/2),grDevices::colorRampPalette(colors = c("white","red"))(length(bks)/2)),
                breaks=bks,cluster_cols=FALSE,cluster_rows=FALSE,show_rownames=F,show_colnames=F,useRaster = T,silent=TRUE)
  colnames(cds) <- c(1:ncol(cds))
  annotation_col <- data.frame(row.names = c(1:ncol(cds)), "Cell Type" = c(rep(branch_labels[1], 70),
                                                                           rep("Pre-branch",  2 * 30),
                                                                           rep(branch_labels[2], 70)))
  colnames(annotation_col) <- "Cell Type"
  names(branch_colors) <- c("Pre-branch", branch_labels[1], branch_labels[2])
  annotation_colors=list("Cell Type"=branch_colors)
  names(annotation_colors$`Cell Type`) = c('Pre-branch', branch_labels)

  ph_res <- pheatmap::pheatmap(cds[, ], #ph$tree_row$order
                     useRaster = T,
                     cluster_cols=FALSE,
                     cluster_rows=FALSE,
                     show_rownames=show_rownames,
                     show_colnames=F,
                     #scale="row",
                     clustering_distance_rows=row_dist, #row_dist
                     clustering_method = hclust_method, #ward.D2
                     cutree_rows=num_clusters,
                     # cutree_cols = 2,
                     # annotation_row=annotation_row,
                     annotation_col=annotation_col,
                     annotation_colors=annotation_colors,
                     #gaps_col = col_gap_ind,
                     treeheight_row = 20,
                     breaks=bks,
                     fontsize = 6,
                     color=c(grDevices::colorRampPalette(colors = c("blue","white"))(length(bks)/2),grDevices::colorRampPalette(colors = c("white","red"))(length(bks)/2)),
                     border_color = NA,
                     silent=TRUE)
  grid::grid.rect(gp=grid::gpar("fill", col=NA))
  grid::grid.draw(ph_res$gtable)

}





