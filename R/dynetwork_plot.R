#' Title
#'
#' @param Time ???? row sample ???? matrix of features X
#' @param ori_eff_mat error
#' @param times Number of cpus for parallel computing, The default is 1
#' @param connect Whether to output running information
#' @param net_odee Whether to output running information
#'
#' @return NULL
#' @export
#'
#' @examples
#'
dynetwork_plot <- function(Time, ori_eff_mat, times, connect, net_odee){


  t_idx <- which(times==Time)
  if(length(t_idx)==0){
    print("ERROR: Enter the parameter T corresponding to the actual time")
  }else{
    if(t_idx == 1){
      t <- 1:(1+2)
      res <- lapply(net_odee, t = t, function (net_odee,t){return(net_odee[t, ])})
      net_stat <- interType_2(con = connect, alle = res, sme = ori_eff_mat[,t])
    }else if(t_idx == length(times)){
      t <- (length(times)-2):length(times)
      res <- lapply(net_odee, t = t, function (net_odee,t){return(net_odee[t, ])})
      net_stat <- interType_2(con = connect, alle = res, sme = ori_eff_mat[,t])
    }else{
      t <- (t_idx-1):(t_idx+1)
      res <- lapply(net_odee, t = t, function (net_odee,t){return(net_odee[t, ])})
      net_stat <- interType_2(con = connect, alle = res, sme = ori_eff_mat[,t])
    }
    net <- as.matrix(connect)
    nn <- dim(net)[1]
    pt.cex_size <- seq(4, 2.5, length.out = 3)
    init_colors <- cm.colors(nn)
    init_colors <- colorRampPalette(c("#B94FFF", "#E8CCFF"))(nn)
    node_names <-  paste0("Node", 1:dim(ori_eff_mat)[1])
    rownames(ori_eff_mat) <- 1:dim(ori_eff_mat)[1]

    g1 <- graph.adjacency((net), mode="directed", weighted = TRUE)
    #-------------------------------------------------------
    OV_sorce_1 <- OV_sorce(net, net_stat$connfest)
    allsize <- OV_sorce_1/2
    #-------------------------------------------------------
    order_idx <- order(allsize,decreasing = TRUE)
    node_colors <- 1:nn
    node_colors[order_idx]<- init_colors
    legend_1_10_color <- init_colors[1:3]
    legend_1_10_names <- node_names[order_idx[1:3]]
    #------------------------------------------------------------------------------
    xx <- NULL
    for (i in 1:dim(net_stat$connfest)[1])
    {
      index_i <- which(net_stat$connfest[i,]!=0)
      xx_i <- cbind(rep(i, length(index_i)), index_i, net_stat$connfest[i,index_i])
      xx <- rbind(xx, xx_i)
    }
    Promote_index  <- which(xx[,3]>0)
    Inhibition_index <- which(xx[,3]<0)
    coll <- rep("#CCCCCC", length(E(g1)))
    coll[Promote_index] <- "blue"
    coll[Inhibition_index] <- "red"
    e.w <<- abs(xx[,3])*5
    nr <-  rep("triangle",length(e.w))
    nr[Inhibition_index] <- "T"
    rr1 <<- nr
    V(g1)$color <- node_colors
    V(g1)$label <- node_names

    x_label <- paste0("Time = ", Time)
    outfilename <- paste0("network_t",Time,".pdf")

    pdf(outfilename,height = 8,width=10)
    plot(g1,edge.arrow.size=1.5,
         layout=layout.circle,
         vertex.size=allsize,
         vertex.frame.color="transparent",
         vertex.label.color="black",
         vertex.label.cex=2,
         edge.curved=0,
         edge.color=coll)
    mtext(x_label,side = 1,line=1.5,cex=2,adj = 0.5)
    legend("right",    "center", title= "Rank", legend = legend_1_10_names, col = legend_1_10_color, pch=16, pt.cex = pt.cex_size, pt.bg = legend_1_10_color, cex = 1.5, bty = "n")
    dev.off()
  }
}




interType_2 <- function(con,alle,sme){

  diag(con) <- 0
  nn <- dim(con)[1]
  connfest <- matrix(0,nrow=nn,ncol=nn)
  indp <- c()
  inter <- list()
  for(i in 1:nn){
    al <- alle[[i]]
    index <- which(as.numeric(colnames(al))==i)
    if(is.matrix(al)){
      lindp <- al[,index]
      linter <- al[,-index]
      indp <- cbind(indp,lindp)
      inter[[i]] <- linter
      rcor <- cor(as.numeric(sme[i,]),linter)
    }else{
      indp <- cbind(indp,al)
      inter[[i]] <- 0
      rcor <- 0
    }
    connfest[i,which(con[i,]==1)] <- as.numeric(rcor)
  }
  return(list(connfest=connfest,connect=con,indp=indp,inter=inter))
}



OV_sorce <- function (net, connfest){
  n_node <- dim(net)[1]
  KS_NI_v <- data.frame(KS = rep(0, n_node), NI = rep(0, n_node))
  d_1 = 0
  iter = 1
  g <- graph.adjacency((net),mode="directed", weighted = TRUE)
  net_degree <- as.vector(igraph::degree(g))
  weightd <- rowSums(abs(connfest))
  net_n <- net
  rownames(net_n) <- 1:n_node
  colnames(net_n) <- 1:n_node
  while (!setequal(integer(0),which(KS_NI_v[, 1]==0))){
    d_0 <- d_1
    g <- graph.adjacency((net_n),mode="directed", weighted = TRUE)
    net_degree <- as.vector(igraph::degree(g))
    d_1 <- min(net_degree)
    min_degree_idx <- which(net_degree == d_1)
    idx_n <- as.numeric(rownames(net_n)[min_degree_idx])
    KS_NI_v[idx_n, 1] = d_1
    if(d_0 == d_1){
      iter <- iter + 1
    }else{
      iter <- 1
    }
    KS_NI_v[idx_n,2] <- iter
    net_n <- net_n[-min_degree_idx, -min_degree_idx]
  }
  tmp_1 <- rbind(degree = net_degree,KS = KS_NI_v[,1],NI = KS_NI_v[,2], weightd = weightd)
  tmp_2 <- tmp_1
  tmp_3 <- rep(NA, 4)
  for (i in 1:4){
    tmp_2[i,] <- scale(tmp_1[i,], center = FALSE, scale = TRUE)
    tmp_3[i] <- sum(tmp_2[i,]*log1p(tmp_2[i,]))/log(20)
  }
  tmp_4 <- (1-tmp_3)/(4-sum(tmp_3))
  tmp_5 <- colSums(tmp_1*tmp_4)
  ov <- (tmp_1[2,] + tmp_1[3,]/(tmp_1[2,]+tmp_1[3,]))*tmp_5 + tmp_1[4,]*tmp_5
}
#-----------------------------------------------------------------------------

