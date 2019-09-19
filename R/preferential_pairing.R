#' Preferential pairing profiles
#'
#' Compute the probability profiles for all pairing configurations in both parents.
#'
#' @param input.genoprobs an object of class \code{"mappoly.genoprob"} 
#' @param x an object of class \code{"mappoly.prefpair.profiles"}
#' @param min.y.prof lower bound for y axis on the probability profile graphic.
#' @param max.y.prof upper bound for y axis on the probability profile graphic.
#' @param thresh line for chi-square test
#' @param ... unused arguments
#' 
#'@examples
#' \dontrun{
#'   ## hexaploid example
#'   w1 <- lapply(maps.hexafake, calc_genoprob)
#'   x1 <- calc_prefpair_profiles(w1)
#'   print(x1)
#'   plot(x1, min.y.prof = 0.05, max.y.prof = .15, thresh = 0.01)
#'   
#'   ## tetraploid example
#'   w2 <- lapply(solcap.err.map, calc_genoprob_error, error = 0.05)
#'   x2 <- calc_prefpair_profiles(w2)
#'   print(x2)
#'   plot(x2)
#'}
#'
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu}
#'
#' @references
#'     Mollinari, M. et al. (2019) Unraveling the hexaploid 
#'     sweetpotato inheritance using ultra-dense multilocus mapping, 
#'     _submited_. \url{https://doi.org/10.1101/689638}
#'     
#' @export
#' @importFrom ggplot2 ggplot geom_hline theme geom_smooth ggtitle facet_grid theme_minimal ylab xlab aes vars scale_color_manual
#' @importFrom gridExtra grid.arrange
#' @importFrom reshape melt
#' 
calc_prefpair_profiles<-function(input.genoprobs){
  if(class(input.genoprobs) == "mappoly.genoprob")
    input.genoprobs <- list(input.genoprobs)
  if(class(input.genoprobs) == "list"){
    if(!all(sapply(input.genoprobs, class) == "mappoly.genoprob"))
      stop(deparse(substitute(input.genoprobs)), " is not an object of class 'mappoly.sequence' neither a list containing 'mappoly.sequence' objects.")
  } 
  df.prefpair <- df.prefpair.pval <- NULL
  for(j in 1:length(input.genoprobs)){
    cat("\nLinkage group ", j, "...")
    Gnames<-dimnames(input.genoprobs[[j]]$probs)[[1]]
    x<-dim(input.genoprobs[[j]]$probs)
    nsta<-x[1]
    npos<-x[2]
    nind<-x[3]
    v<-c(2,4,6,8,10,12)
    names(v)<-choose(v,v/2)^2
    m<-v[as.character(nsta)]
    ## Possible pairing configurations (Psi's) for both parents
    Psi <- list(P = NULL, Q = NULL)
    aP<-combn(letters[1:m],m/2)
    aQ<-combn(letters[(1+m):(2*m)],m/2)
    for(i in 1:ncol(aP))
    {
      bP<-perm_tot(setdiff(letters[1:m], aP[,i]))
      bQ<-perm_tot(setdiff(letters[(1+m):(2*m)], aQ[,i]))
      for(k in 1:nrow(bP))
      {
        Psi$P<-rbind(Psi$P, sort(apply(cbind(aP[,i],bP[k,]), 1, function(x) paste(sort(x), collapse = ""))))
        Psi$Q<-rbind(Psi$Q, sort(apply(cbind(aQ[,i],bQ[k,]), 1, function(x) paste(sort(x), collapse = ""))))
      }
    }
    Psi$P<-unique(Psi$P)
    Psi$Q<-unique(Psi$Q)
    nbiv<-nrow(Psi$P)
    ## Probability of a pairing configuration (Psi) given a genitypic state (G) for both parents
    Psi_given_G<-list(P = matrix(0, length(Gnames), nrow(Psi$P), dimnames = list(Gnames, apply(Psi$P, 1, paste, collapse="/"))),
                      Q = matrix(0, length(Gnames), nrow(Psi$Q), dimnames = list(Gnames, apply(Psi$Q, 1, paste, collapse="/"))))
    ## Genotypic states
    G<-list(P = unlist(strsplit(substr(Gnames, 1, m/2), "")),
            Q = unlist(strsplit(substr(Gnames, 2+(m/2), m+1), "")))
    dim(G$Q) <- dim(G$P) <- c(m/2,length(Gnames))
    
    for(i in 1:nrow(Psi_given_G$P))
    {
      xP<-apply(Psi$P, 1, function(x,y) all(grepl(paste(y, collapse="|"), x)), y = G$P[,i])
      xQ<-apply(Psi$Q, 1, function(x,y) all(grepl(paste(y, collapse="|"), x)), y = G$Q[,i])
      #print(as.character(sum(xP)))
      #print(as.character(factorial(m/2)))
      #print("-----")
      Psi_given_G$Q[i,xQ] <- Psi_given_G$P[i,xP] <- 1/factorial(m/2)
    }
    ## Equation 3 in Mollinari et al. (2019)
    AQ <- AP <- array(NA, dim = c(nbiv,npos,nind))
    for(i in 1:nind){
      AP[,,i]<-crossprod(Psi_given_G$P,input.genoprobs[[j]]$probs[,,i])   
      AQ[,,i]<-crossprod(Psi_given_G$Q,input.genoprobs[[j]]$probs[,,i])   
    }
    ## 
    P<-apply(AP, MARGIN = c(1,2), mean)
    Q<-apply(AQ, MARGIN = c(1,2), mean)
    dimnames(P)<-list(colnames(Psi_given_G$P), names(input.genoprobs[[j]]$map))
    dimnames(Q)<-list(colnames(Psi_given_G$Q), names(input.genoprobs[[j]]$map))
    df.prefpair.temp<-rbind(data.frame(reshape::melt(P), parent = "P", lg = j),
                            data.frame(reshape::melt(Q), parent = "Q", lg = j))
    colnames(df.prefpair.temp) <- c("pair.conf", "marker", "probability", "parent", "LG")
    map<-data.frame(map.position = input.genoprobs[[j]]$map, marker = names(input.genoprobs[[j]]$map))
    df.prefpair.temp<-merge(df.prefpair.temp, map, sort = FALSE)
    df.prefpair <- rbind(df.prefpair, df.prefpair.temp)
    df.prefpair.pval.temp<-data.frame(p.val.P = apply(nind * P, 2, function(x) chisq.test(x)$p.value),
                                      p.val.Q = apply(nind * Q, 2, function(x) chisq.test(x)$p.value), 
                                      #p.val.P = apply((m/2) * nind * P, 2, function(x) chisq.test(x)$p.value),
                                      #p.val.Q = apply((m/2) * nind * Q, 2, function(x) chisq.test(x)$p.value), 
                                      LG = j,
                                      map.position = input.genoprobs[[j]]$map)
    df.prefpair.pval <- rbind(df.prefpair.pval, df.prefpair.pval.temp)
  }
  cat("\n")
  structure(list(info = list(m = m, nind = nind) , prefpair = df.prefpair, prefpair.pval = df.prefpair.pval), class = "mappoly.prefpair.profiles")
}
#' @rdname calc_prefpair_profiles
#' @keywords internal
#' @export
print.mappoly.prefpair.profiles <- function(x, ...){
  cat("  This is an object of class 'mappoly.prefpair.profiles'")
  cat("\n  -----------------------------------------------------\n")
  cat("  No. positions:                            ", nrow(x$prefpair.pval),"in", length(unique(x$prefpair$LG)),"LGs\n")
  cat("  No. individuals:                          ", x$info$nind, "\n")
  cat("  -----------------------------------------------------")
  ## printing summary
  parent <- NULL
  cat("\n  Pairing configurations:\n")
  a1<-as.character(unique(subset(x$prefpair, parent == "P")$pair.conf))
  a2<-as.character(unique(subset(x$prefpair, parent == "Q")$pair.conf))
  a1<-matrix(a1, nrow = 1)
  a2<-matrix(a2, nrow = 1)
  if(x$info$m == 6)
    dim(a2) <- dim(a1) <- c(3, 5)
  else if (x$info$m == 8)
    dim(a2) <- dim(a1) <- c(21, 5)
  else if (x$info$m == 10)
    dim(a2) <- dim(a1) <- c(189, 5)
  cat("     Parent 1:\n")
  for(i in 1:nrow(a1)){
    cat("       ", a1[i,], "\n", sep = "    ")
  }
  cat("     Parent 2:\n")
  for(i in 1:nrow(a2)){
    cat("       ", a2[i,], "\n", sep = "    ")
  }
  cat("  -----------------------------------------------------\n")
}
#' @rdname calc_prefpair_profiles
#' @keywords internal
#' @export
plot.mappoly.prefpair.profiles <- function(x, min.y.prof = 0, max.y.prof = 1, thresh = 0.01, ...){
  m<-x$info$m
  variable <- value <- map.position <- probability <- colour <- pair.conf <-NULL
  p1<-ggplot2::ggplot(x$prefpair) + 
    ggplot2::geom_smooth(aes(map.position, probability, colour = pair.conf), size = .5, se = FALSE) + 
    ggplot2::facet_grid(parent~LG, scales = "free_x", space = "free_x") +
    ggplot2::geom_hline(yintercept = 1/(prod(choose(seq(2,m,2),2))/factorial(m/2)), linetype="dashed") +
    ggplot2::ylim(min.y.prof,max.y.prof) + 
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) +
    ggplot2::ylab("Probability") + ggplot2::xlab("Distance (cM)") 
  DF<-reshape::melt(data = x$prefpair.pval, measure.vars = c("p.val.P", "p.val.Q"))
  p2<-ggplot2::ggplot(DF, aes(map.position, -log10(value), colour = variable)) +
    ggplot2::geom_point(alpha = .7, size = 1) +  
    ggplot2::facet_grid(.~LG, scales = "free_x", space = "free_x") +
    ggplot2::geom_hline(yintercept = -log10(thresh), linetype="dashed") +  
    ggplot2::scale_color_manual(values=c("#E69F00","#56B4E9"), name = "Parents") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) +
    ggplot2::ylab("-log_10(P)") + ggplot2::xlab("Distance (cM)")
  p<-gridExtra::grid.arrange(p1, p2, nrow = 2)
  invisible(p)
}

