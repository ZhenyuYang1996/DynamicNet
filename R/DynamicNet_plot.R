
#' Title
#'
#' @param belong_data ???? row sample ???? matrix of features X
#' @param ori_data error
#' @param means_data Number of cpus for parallel computing, The default is 1
#' @param times Output file location and name prefix
#' @param J Output file location and name prefix
#' @param plot.control Whether to output running information
#'
#' @return NULL
#' @export
#'
#' @examples
#'
fclusty_plot <- function(belong_data, ori_data, means_data, times, J,
                         plot.control = list(y_cex = NULL, col_mean = "blue", main_label = NULL, no_out_pdf = TRUE)){
  #----------------------------------------------------------------------------
  y_cex <- plot.control$y_cex
  col_mean <- plot.control$col_mean
  model_text <- plot.control$model_text
  no_out_pdf <- plot.control$no_out_pdf
  #----------------------------------------------------------------------------
  cluster_real_index = which( belong_data[,2] == (J-1) )
  ori_data_1 = matrix(NA, nrow = length(times)*length(cluster_real_index), ncol = 3)
  for (i in 1:length(cluster_real_index))
  {
    k_1 = (i-1)*length(times) + 1
    ori_data_1_i = cbind(times, as.numeric(ori_data[cluster_real_index[i], 1:length(times)]), rep(i,length(times)))
    ori_data_1[k_1:(k_1+length(times)-1), ] = ori_data_1_i
  }
  ori_data_1 = as.data.frame(ori_data_1)
  colnames(ori_data_1) = c("time", "pheno", "index")
  #------------------------------------------------------------------------------------------
  if(length(cluster_real_index) > 100){
    cluster_i_index = sample(cluster_real_index, size = 100, replace = FALSE)
  }else{
    cluster_i_index = cluster_real_index
  }
  ori_data_p = matrix(NA, nrow = length(times)*length(cluster_i_index), ncol = 3)
  for (i in 1:length(cluster_i_index))
  {
    k_1 = (i-1)*length(times) + 1
    ori_data_1_i = cbind(times, as.numeric(ori_data[cluster_i_index[i], 1:length(times)]), rep(i,length(times)))
    ori_data_p[k_1:(k_1+length(times)-1), ] = ori_data_1_i
  }
  ori_data_p = as.data.frame(ori_data_p)
  colnames(ori_data_p) = c("time", "pheno", "index")
  #------------------------------------------------------------------------------------------
  #real_mean_data <- data.frame(x = 1:20, y = real_data[,J], index = rep("real",20))
  #------------------------------------------------------------------------------------------
  sim_mean_data <- data.frame(x = 1:20, y = means_data[,J], index = rep("sim",20))

  g <- ggplot() + theme_bw() + theme(panel.grid =element_blank()) + theme(text = element_text(family = windowsFonts()$serif))

  g <- g + geom_line(data = ori_data_p, aes(x = time, y = pheno, group = index), col = "gray", size = 0.5)
  #
  g <- g + geom_boxplot(data = ori_data_1, aes(x = time, y = pheno, group = time ), colour = "#000000")
  #g <- g + geom_line(data = real_mean_data , aes(x = x, y = y), col = "red", size = 1.5)
  g <- g + geom_line(data = sim_mean_data , aes(x = x, y = y), col = col_mean, size = 1.5)
  g <- g + scale_x_continuous(limits=c(1,length(times)),breaks = seq(1,length(times),2),labels = seq(1,length(times),2))
  g <- g + scale_y_continuous(limits=c(min(ori_data),max(ori_data)))
  g <- g + annotate("text", x = (times[1]+times[length(times)])/2, y = max(ori_data)-0.5, label = model_text, size = 5)
  g <- g + theme(plot.title = element_blank(),
                 axis.title.x = element_blank(),
                 axis.text.x= element_text(size = 14,vjust = 0.5,colour = "black"),
                 axis.ticks.length.x = unit(0.1,"cm"))
  if(is.null(NULL)) {
    g <- g + theme(axis.title.y = element_blank(),
                   axis.text.y = element_text(size = 14,vjust = 0.5,colour = "black"),
                   axis.ticks.length.y = unit(0.1,"cm"))
    g <-  g +   theme(plot.margin = unit(c(0.5,0,0.5,0.5),"lines"))
  }else{
    g <- g + theme( axis.title.y = element_blank(),
                    axis.ticks.y = element_blank(),
                    axis.ticks.length.y = unit(-0.1,"cm"),
                    axis.text.y = element_blank())
    g <-  g +   theme(plot.margin = unit(c(0.5,0.5,0.5,0),"lines"))
  }
  return(g)
}
#' Title
#'
#' @param ori_eff_data ???? row sample ???? matrix of features X
#' @param odee_matrix error
#' @param times Number of cpus for parallel computing, The default is 1
#' @param plot.control Whether to output running information
#'
#' @return NULL
#' @export
#'
#' @examples
#'
net_eff_plot <- function(ori_eff_data, odee_matrix, times,
                         plot.control = list(main_label = NULL)){
  cb_palette <- rainbow(length(ori_eff_data))

  main_label <- plot.control$main_label

  odee_col_names <- as.numeric(colnames(odee_matrix))

  ind_idx <- odee_col_names[1]
  dep_idx <- odee_col_names[-1]

  y_max <- max(max(odee_matrix), max(ori_eff_data))
  y_min <- min(min(odee_matrix), min(ori_eff_data))
  ori_eff_frame <- data.frame(times = times, eff = ori_eff_data, index = rep("G 0", length(times)))
  ind_eff_frame <- data.frame(times = times, eff = (odee_matrix[,1]+ori_eff_data[1]), index = rep(paste0("IND", ind_idx), length(times)))

  if(length(dep_idx)==0){
    legend_label <- c("ORI EFF", paste0("IND ",ind_idx), "NULL DEP")
    curve_col <- c("red", "blue", "grey")
    names(curve_col) <- c("G 0", paste0("IND",ind_idx), "NULL DEP")
    names(legend_label) <- c("G 0", paste0("IND",ind_idx), "NULL DEP")
    dep_eff_frame <- data.frame(times = times, eff = rep(0, length(ori_eff_data)), index = rep("NULL DEP", length(times)))
  }else{
    legend_label <- c("ORI EFF", paste0("IND ",ind_idx), paste0("DEP ",dep_idx))
    curve_col <- c("red", "blue", cb_palette[dep_idx])
    names(curve_col) <- c("G 0", paste0("IND",ind_idx), paste0("DEP",dep_idx))
    names(legend_label) <- c("G 0", paste0("IND",ind_idx), paste0("DEP",dep_idx))
    dep_eff_frame <- NULL
    for (i in 2:dim(odee_matrix)[2])
    {
      dep_eff_i <- data.frame(times = times, eff = (odee_matrix[,i]+ori_eff_data[1]), index = rep(paste0("DEP",dep_idx[i-1]), length(times)))
      dep_eff_frame <- rbind(dep_eff_frame, dep_eff_i)
    }
  }

  y_position = sort(odee_matrix[length(times),-1])

  dep_label <- names(y_position)

  factor_y <- 0.7 * (y_max - y_min ) / (length(dep_label)-1)
  adjust_x <- (y_max - y_min - factor_y*(length(dep_label)-1)) / 2
  annotate_y <- seq(y_min+adjust_x, y_max-adjust_x, length.out = length(dep_label))

  g <-  ggplot() + theme_bw() + theme(panel.grid =element_blank())
  g <- g + geom_line(data = dep_eff_frame , aes(x = times, y = eff, colour = index),
                     size = 1, alpha = 0.8)
  g <- g + geom_line(data = ind_eff_frame , aes(x = times, y = eff, colour = index),
                     size = 1.5)
  g <- g + geom_line(data = ori_eff_frame, aes(x = times, y = eff, colour = index),
                     size = 1.5)
  g <- g + scale_color_manual(name = "index",
                              values = curve_col,
                              label = legend_label)

  g <- g + scale_x_continuous(limits=c(1,(max(times+3.5))),breaks = seq(1,max(times),1),labels = seq(1,max(times),1))
  g <- g + labs(title = main_label)
  #=======================
  #g <- g + scale_y_continuous(limits=c(-1,1.2),breaks = seq(-1,1,0.2),labels = seq(-1,1,0.2))
  #=======================
  g <- g + annotate("segment", x = -Inf, xend =  Inf, y = 0, yend = 0, colour = "black", size = 0.5)

  g <- g + theme(plot.title = element_text(size = 22,hjust = 0.5,colour = "black"),
                 axis.title.x = element_blank(),
                 axis.text.x= element_text(size = 20,hjust = 0.5,colour = "black"),
                 axis.ticks.length.x = unit(0.1,"cm"))
  g <- g + theme(axis.title.y = element_blank(),
                 axis.text.y = element_text(size = 20,hjust = 0.5,colour = "black"),
                 axis.ticks.length.y = unit(0.1,"cm"))
  g <- g + theme(legend.title=element_blank(),
                 legend.position = c(0.98, 0.95),
                 legend.justification = c(0.98, 0.95),
                 legend.text = element_text(size = 15,hjust = 0.5,colour = "black"))
  return(g)
}
