
#' Title
#'
#' @param X ???? row sample ???? matrix of features X
#' @param eps error
#' @param core_n Number of cpus for parallel computing, The default is 1
#' @param outfilename Output file location and name prefix
#' @param verbose Whether to output running information
#'
#' @return NULL
#' @export
#'
#' @examples
#'
Dnelastic <- function(X, eps, core_n = NULL, outfilename, verbose = FALSE){
  ori_data <- as.matrix(X)
  if(is.null(NULL)){
    core.number <- detectCores()
  }else{
    core.number <- 1
  }
  cl <- makeCluster(getOption("cl.cores", core.number))
  clusterExport(cl,c("cv.glmnet", "glmnet"))
  connect_matrix <- parSapply(cl,1:dim(ori_data)[1],function(col){  ## Linear regression
    #lasso_data <- scale(ori_data, center=T,scale=F)
    lasso_data <- ori_data
    res = rep(0,dim(lasso_data)[1])
    m <- as.numeric(lasso_data[col,])
    M <- lasso_data[-col,]
    n <- dim(M)[1]
    vec <- rep(NA,length(M[,1]))
    for (i in 1:length(M[,1]))
    {
      vec[i] <- cor(m,as.numeric(M[i,]))
    }
    #corre_index = which( vec %in% -sort(-vec)[1:(n/(1*log(n)))] )
    corre_index = which( vec %in% -sort(-vec)[1:n] )
    x_matrix <- t(M[corre_index,])
    y <- as.numeric(lasso_data[col,])
    cv.fit  <- cv.glmnet(x_matrix, y, family = "gaussian", type.measure = "mse")
    fit     <-    glmnet(x_matrix, y, family = "gaussian", lambda = cv.fit$lambda.min)

    coefficients <- coef(fit,s=cv.fit$lambda.min)
    Active.Index <- which(coefficients[-1] != 0)
    #rownames(coefficients[-1][which(coefficients[-1] != 0)])
    index_0 = corre_index[Active.Index]
    #rownames(M)[index_0]
    index_1 = index_0[which(index_0 >= col)]+1
    index_2 = index_0[which(index_0 < col)]
    res[c(index_1,index_2)] = 1
    #length(which(res!=0))
    return(res)
  })
  stopCluster(cl)
  write.table(t(connect_matrix), file = paste0(outfilename,"_connect.txt"), row.names = FALSE, col.names = FALSE)
}
