map_error_model <- function(input.map.list, error, tol, ncpu = NULL){
  
  if(inherits(input.map.list, "mappoly.map2")){
    if(is.null(ncpu))
      ncpu <- length(input.map.list$both$map.list)    
    if(parallel::detectCores() < ncpu) 
      ncpu <- parallel::detectCores()
    cl <- parallel::makeCluster(ncpu)
    parallel::clusterExport(cl, "est_full_hmm_with_global_error")
    parallel::clusterExport(cl, "dat")
    MAP.err <- parallel::parLapply(cl,input.map.list$both$map.list,
                                   est_full_hmm_with_global_error, 
                                   error = error, tol = tol)
    parallel::stopCluster(cl)
    input.map.list$both$map.list.err <- MAP.err
    return(input.map.list)    
  } else {
    if(is.null(ncpu))
      ncpu <- length(input.map.list)    
    if(parallel::detectCores() < ncpu) 
      ncpu <- parallel::detectCores()
    cl <- parallel::makeCluster(ncpu)
    parallel::clusterExport(cl, "est_full_hmm_with_global_error")
    parallel::clusterExport(cl, "dat")
    MAP.err <- parallel::parLapply(cl,input.map.list,
                                   est_full_hmm_with_global_error, 
                                   error = error, tol = tol)
    parallel::stopCluster(cl)
    return(MAP.err)    
  }
}

plot.mappoly.map2 <- function(x){
  {
    df1 <- t(sapply(x$both$map.list, function(x) summary(diff(extract_map(x)))))
    nit <- nrow(df1)
    df11 <- reshape2::melt(df1)
    df11$hmm <- "partial"
    if(!is.na(names(x$both)[4])){
      df2 <- t(sapply(x$both$map.list.err, function(x) summary(diff(extract_map(x)))))   
      df12 <- reshape2::melt(df2)
      df12$hmm <- "full"
      df <- rbind(df11, df12)
      p1<-ggplot(df, aes(x=as.factor(Var1), y=value, group=Var2)) +
        geom_line(aes(color=Var2)) +
        geom_point(aes(color=Var2)) + 
        ylab("nearest neighbor marker distance") +
        xlab("Iteration") + 
        #annotate(geom="label", x = 1:nit, y = df1[,6] + .5, 
        #         label=sapply(x$both$map.list, function(x) x$info$n.mrk),
        #         color="red") + 
        theme(legend.position = "none")  + facet_wrap(.~hmm)
    } else{
      df <- df11
      p1<-ggplot(df, aes(x=as.factor(Var1), y=value, group=Var2)) +
        geom_line(aes(color=Var2)) +
        geom_point(aes(color=Var2)) + 
        ylab("nearest neighbor marker distance") +
        xlab("Iteration") + 
        annotate(geom="label", x = 1:nit, y = df1[,6] + .5, 
                 label=sapply(x$both$map.list, function(x) x$info$n.mrk),
                 color="red") + 
        theme(legend.position = "none")
    }
  }  
  {
    df31 <- lapply(x$both$map.list, function(x) extract_map(x))
    df41 <- reshape2::melt(df31)
    df41$hmm <- "partial"
    colnames(df41)[1:2] <- c("position", "iteration")
    len.map1 <- sapply(df31, max)  
    if(!is.na(names(x$both)[4])){
      df32 <- lapply(x$both$map.list.err, function(x) extract_map(x))    
      df42 <- reshape2::melt(df32)
      df42$hmm <- "full"
      colnames(df42)[1:2] <- c("position", "iteration")
      len.map2 <- sapply(df32, max)  
      df4 <- rbind(df41, df42)
      #len.map <- c(len.map1, len.map2)
    } else{
      df4 <- df41
    }
    p2<-ggplot(df4, aes(y=-position, x = as.factor(iteration), color = hmm)) +
      geom_point(shape = 95, size = 5) +
      # annotate(geom="label", y = -len.map + 1, x = 1:length(len.map)  + 0.2, 
      #          label=round(len.map, 1),
      #          color="red", size = 3) + 
      xlab("iteration") + facet_wrap(.~hmm) + theme(legend.position = "none")+
      scale_color_manual(values=c('#E69F00', '#56B4E9'))
  }
  
  df5 <- lapply(x$single, function(x) extract_map(x))
  df6 <- reshape2::melt(df5)
  head(df6)
  colnames(df6) <- c("position", "parent")
  p3<-ggplot(df6, aes(x=position, y=parent, group = as.factor(parent))) +
    geom_point(ggplot2::aes(color = as.factor(parent)), shape = 108, size = 5, show.legend = FALSE) +
    scale_color_brewer(palette="Dark2")
  
  
  df7 <- reshape2::melt(lapply(x$both$calc.lod, na.omit))
  colnames(df7) <- c("LODScore", "iteration")
  p4<-ggplot(df7, aes(x=as.factor(iteration), y=LODScore, group = as.factor(iteration))) +
    geom_boxplot(fill='#A4A4A4', color="black") + xlab("iteration")
  
  f1 <- ggarrange(p1, p3 + font("x.text", size = 9),
                  ncol = 2, nrow = 1, widths = c(2,1))
  f2 <- ggarrange(p2, p4 + font("x.text", size = 9),
                  ncol = 2, nrow = 1, widths = c(2,1))
  ggarrange(f1, f2, ncol = 1, nrow = 2)
}








