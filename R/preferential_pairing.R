#' Preferential pairing profiles
#'
#' Compute the probability profiles for all pairing configurations in both parents.
#'
#' @param input.genoprobs an object of class \code{mappoly.genoprob}
#' 
#' @param type a character string indicating which type of graphic is ploted:
#'             \code{"pair.configs"} (default: plots the preferential pairing 
#'             profile for the pairing configurations) or \code{"hom.pairs"} (plots 
#'             the preferential pairing profile for the homolog pairs)
#'             
#' @param min.y.prof lower bound for y axis on the probability profile graphic (default = 0)
#' 
#' @param max.y.prof upper bound for y axis on the probability profile graphic (default = 1)
#' 
#' @param thresh threshold for chi-square test (default = 0.01)
#' 
#' @param x an object of class \code{mappoly.prefpair.profiles}
#' 
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
#'     _G3: Genes, Genomes, Genetics_. \url{http://dx.doi.org/10.1534/g3.119.400620}
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
  df.hompair.pval <- df.prefpair <- df.prefpair.pval <- NULL
  for(j in 1:length(input.genoprobs)){
    cat("\nLinkage group ", j, "...")
    # get names for all states
    Gnames<-dimnames(input.genoprobs[[j]]$probs)[[1]]
    x<-dim(input.genoprobs[[j]]$probs)
    # number of states, positions and individuals
    nsta<-x[1]
    npos<-x[2]
    nind<-x[3]
    # get ploidy level based on the number of states
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
    # number of possible bivalent configurations (size of Psi)
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
    ## Pr(psi_i | O_1, ..., Oz, lambda)
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
                                      LG = j,
                                      map.position = input.genoprobs[[j]]$map)
    ## Probability of bivalent pairs given the complete map
    h.h_prime.P <- apply(combn(letters[1:m], 2), 2, paste0, collapse = "")
    h.h_prime.Q <- apply(combn(letters[(m+1):(2*m)], 2), 2, paste0, collapse = "")
    p1<-get_w_m(m - 2)/get_w_m(m)
    p<-c(p1,1 - p1)
    df.hom.pair <- NULL
    for(k1 in h.h_prime.P){
      phh<-apply(P[grepl(k1, rownames(P)), , drop = FALSE], 2, sum)
      phh.P <- data.frame(marker = names(phh), phh = phh, hom.pair = k1, LG = j, parent = "P")
      x1 <- phh * nind
      x<-cbind(x1,nind - x1)
      p.val.P<-apply(x, 1, function(x) chisq.test(x, p = p)$p.value)
      p.val.P <- data.frame(marker = names(p.val.P), p.val = p.val.P)
      mtemp <- merge(phh.P, p.val.P)
      mtemp <- merge(mtemp, map)
      df.hom.pair<-rbind(df.hom.pair, mtemp)
    }
    for(k1 in h.h_prime.Q){
      phh <- apply(Q[grepl(k1, rownames(Q)), , drop = FALSE], 2, sum)
      phh.Q <- data.frame(marker = names(phh), phh = phh, hom.pair = k1, LG = j, parent = "Q")
      x1 <- phh * nind
      x<-cbind(x1, nind - x1)
      p.val.Q<-apply(x, 1, function(x) chisq.test(x, p = p)$p.value)
      p.val.Q <- data.frame(marker = names(p.val.Q), p.val = p.val.Q)
      mtemp <- merge(phh.Q, p.val.Q)
      mtemp <- merge(mtemp, map)
      df.hom.pair<-rbind(df.hom.pair, mtemp)
    }
    df.hompair.pval <- rbind(df.hompair.pval, df.hom.pair)
    df.prefpair.pval <- rbind(df.prefpair.pval, df.prefpair.pval.temp)
  }
  cat("\n")
  structure(list(info = list(m = m, nind = nind), 
                 prefpair.psi = df.prefpair, 
                 prefpair.psi.pval = df.prefpair.pval,
                 prefpair.homolog = df.hompair.pval),
            class = "mappoly.prefpair.profiles")
}
#' @rdname calc_prefpair_profiles
#' @keywords internal
#' @export
print.mappoly.prefpair.profiles <- function(x, ...){
  cat("  This is an object of class 'mappoly.prefpair.profiles'")
  cat("\n  -----------------------------------------------------\n")
  cat("  No. positions:                            ", nrow(x$prefpair.psi.pval),"in", length(unique(x$prefpair.psi$LG)),"LGs\n")
  cat("  No. individuals:                          ", x$info$nind, "\n")
  cat("  -----------------------------------------------------")
  ## printing summary
  parent <- NULL
  cat("\n  Pairing configurations:\n")
  a1<-as.character(unique(subset(x$prefpair.psi, parent == "P")$pair.conf))
  a1.h<-apply(combn(letters[1:x$info$m], 2), 2, paste0, collapse = "")
  a2<-as.character(unique(subset(x$prefpair.psi, parent == "Q")$pair.conf))
  a2.h<-apply(combn(letters[(x$info$m+1):(2*x$info$m)], 2), 2, paste0, collapse = "")
  a1<-matrix(a1, nrow = 1)
  a2<-matrix(a2, nrow = 1)
  a1.h<-matrix(a1.h, nrow = 1)
  a2.h<-matrix(a2.h, nrow = 1)
  if(x$info$m == 6){
    dim(a2) <- dim(a1) <- c(3, 5)    
    dim(a2.h) <- dim(a1.h) <- c(3, 5)    
  } else if (x$info$m == 8){
    dim(a2) <- dim(a1) <- c(21, 5)
    dim(a2.h) <- dim(a1.h) <- c(7, 4)
  } else if (x$info$m == 10){
    dim(a2) <- dim(a1) <- c(189, 5)
    dim(a2.h) <- dim(a1.h) <- c(9, 5)
  }
  cat("     Parent 1:\n")
  for(i in 1:nrow(a1)){
    cat("       ", a1[i,], "\n", sep = "    ")
  }
  cat("     Parent 2:\n")
  for(i in 1:nrow(a2)){
    cat("       ", a2[i,], "\n", sep = "    ")
  }
  cat("  -----------------------------------------------------\n")
  cat("  Homolog pairs:\n")
  cat("     Parent 1:\n")
  for(i in 1:nrow(a1.h)){
    cat("       ", a1.h[i,], "\n", sep = "    ")
  }
  cat("     Parent 2:\n")
  for(i in 1:nrow(a2.h)){
    cat("       ", a2.h[i,], "\n", sep = "    ")
  }
  cat("  -----------------------------------------------------\n")
}
#' @rdname calc_prefpair_profiles
#' @keywords internal
#' @export
plot.mappoly.prefpair.profiles <- function(x, type = c("pair.configs", "hom.pairs"), 
                                           min.y.prof = 0, max.y.prof = 1, 
                                           thresh = 0.01, ...){
  type <- match.arg(type)
  if(type == "pair.configs"){
    m<-x$info$m
    variable <- value <- map.position <- probability <- colour <- pair.conf <-NULL
    p1<-ggplot2::ggplot(x$prefpair.psi) + 
      ggplot2::geom_smooth(ggplot2::aes(map.position, probability, colour = pair.conf), size = 1, se = FALSE) + 
      ggplot2::facet_grid(parent~LG, scales = "free_x", space = "free_x") +
      ggplot2::geom_hline(yintercept = 1/(prod(choose(seq(2,m,2),2))/factorial(m/2)), linetype="dashed") +
      ggplot2::ylim(min.y.prof,max.y.prof) + 
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) +
      ggplot2::ylab("Probability") + ggplot2::xlab("Distance (cM)") 
    DF<-reshape::melt(data = x$prefpair.psi.pval, measure.vars = c("p.val.P", "p.val.Q"))
    p2<-ggplot2::ggplot(DF, ggplot2::aes(map.position, -log10(value), colour = variable)) +
      ggplot2::geom_point(alpha = .7, size = 1) +  
      ggplot2::facet_grid(.~LG, scales = "free_x", space = "free_x") +
      ggplot2::geom_hline(yintercept = -log10(thresh), linetype="dashed") +  
      ggplot2::scale_color_manual(values=c("#E69F00","#56B4E9"), name = "Parents") +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) +
      ggplot2::ylab("-log_10(P)") + ggplot2::xlab("Distance (cM)")
    p<-gridExtra::grid.arrange(p1, p2, nrow = 2)
    invisible(p)
  } else if(type == "hom.pairs"){
    m<-x$info$m
    hom.pair <- phh <- variable <- p.val <- map.position <- probability <- colour <- pair.conf <-NULL
    p1<-ggplot2::ggplot(x$prefpair.homolog) + 
      ggplot2::geom_smooth(ggplot2::aes(map.position, phh, colour = hom.pair), size = 1, se = FALSE) + 
      ggplot2::facet_grid(parent~LG, scales = "free_x", space = "free_x") +
      ggplot2::geom_hline(yintercept = get_w_m(m - 2)/get_w_m(m), linetype="dashed") +
      ggplot2::ylim(min.y.prof,max.y.prof) + 
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) +
      ggplot2::ylab("Probability") + ggplot2::xlab("Distance (cM)") 
    p2<-ggplot2::ggplot(x$prefpair.homolog) + 
      ggplot2::geom_point(ggplot2::aes(map.position, -log10(p.val), colour = hom.pair)) + 
      ggplot2::facet_grid(parent~LG, scales = "free_x", space = "free_x") +
      ggplot2::geom_hline(yintercept = -log10(thresh), linetype="dashed") +  
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) +
      ggplot2::ylab("-log10(P)") + ggplot2::xlab("Distance (cM)") 
    p<-gridExtra::grid.arrange(p1, p2, nrow = 2)
  }
  
}

