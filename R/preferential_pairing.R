#' Preferential pairing profiles
#'
#' Given the genotype conditional probabilities for a map, this function
#' computes the probability profiles for all possible homolog pairing 
#' configurations in both parents.
#'
#' @param input.genoprobs an object of class \code{mappoly.genoprob}
#'
#' @param verbose if \code{TRUE} (default), the current progress is shown; if
#'     \code{FALSE}, no output is produced
#' 
#'@examples
#' \donttest{
#'   ## tetraploid example
#'   w1 <- lapply(solcap.dose.map[1:12], calc_genoprob)
#'   x1 <- calc_prefpair_profiles(w1)
#'   print(x1)
#'   plot(x1)
#'}
#'
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu} and Guilherme Pereira, \email{g.pereira@cgiar.org}
#'
#' @references
#'     Mollinari M., Olukolu B. A.,  Pereira G. da S., 
#'     Khan A., Gemenet D., Yencho G. C., Zeng Z-B. (2020), 
#'     Unraveling the Hexaploid Sweetpotato Inheritance 
#'     Using Ultra-Dense Multilocus Mapping, 
#'     _G3: Genes, Genomes, Genetics_. 
#'     \url{https://doi.org/10.1534/g3.119.400620} 
#'     
#' @export
#' @importFrom ggplot2 ggplot geom_hline theme geom_smooth ggtitle facet_grid theme_minimal ylab xlab aes vars scale_color_manual
#' @importFrom reshape2 melt
#' @importFrom ggpubr ggarrange
#' @importFrom ggsci scale_color_d3 pal_d3
#' 
calc_prefpair_profiles<-function(input.genoprobs, verbose = TRUE){
  if(inherits(input.genoprobs, "mappoly.genoprob")) 
    input.genoprobs <- list(input.genoprobs)
  if(!inherits(input.genoprobs, "list"))
    stop(deparse(substitute(input.genoprobs)), 
         " is not an object of class 'mappoly.genoprob' neither a list containing 'mappoly.genoprob' objects.")
  if (any(!sapply(input.genoprobs, inherits, "mappoly.genoprob"))) 
    stop(deparse(substitute(input.genoprobs)), 
         " is not an object of class 'mappoly.genoprob' neither a list containing 'mappoly.genoprob' objects.")
  df.hompair.pval <- df.prefpair <- df.prefpair.pval <- NULL
  for(j in 1:length(input.genoprobs)){
    if (verbose) cat("\nLinkage group ", j, "...")
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
    df.prefpair.temp<-rbind(data.frame(reshape2::melt(P), parent = "P", lg = j),
                            data.frame(reshape2::melt(Q), parent = "Q", lg = j))
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
  if (verbose) cat("\n")
  structure(list(info = list(m = m, nind = nind), 
                 prefpair.psi = df.prefpair, 
                 prefpair.psi.pval = df.prefpair.pval,
                 prefpair.homolog = df.hompair.pval),
            class = "mappoly.prefpair.profiles")
}

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

#' Plots mappoly.prefpair.profiles
#' 
#' @param x an object of class \code{mappoly.prefpair.profiles}
#' 
#' @param type a character string indicating which type of graphic is plotted:
#'             \code{"pair.configs"} (default) plots the preferential pairing 
#'             profile for the pairing configurations or \code{"hom.pairs"} plots 
#'             the preferential pairing profile for the homolog pairs
#'             
#' @param min.y.prof lower bound for y axis on the probability profile graphic (default = 0)
#' 
#' @param max.y.prof upper bound for y axis on the probability profile graphic (default = 1)
#' 
#' @param thresh threshold for chi-square test (default = 0.01)
#' 
#' @param P a string containing the name of parent P
#' 
#' @param Q a string containing the name of parent Q
#' 
#' @param ... unused arguments
#' @export
plot.mappoly.prefpair.profiles <- function(x, type = c("pair.configs", "hom.pairs"), 
                                           min.y.prof = 0, max.y.prof = 1, 
                                           thresh = 0.01, P = "P", Q = "Q", ...){
  type <- match.arg(type)
  colnames(x$prefpair.psi.pval)[1:2] <- c(P, Q)
  if(type == "pair.configs"){
    m<-x$info$m
    variable <- value <- map.position <- probability <- colour <- pair.conf <-NULL
    p1<-ggplot2::ggplot(x$prefpair.psi) + 
      ggplot2::geom_smooth(ggplot2::aes(map.position, probability, colour = pair.conf), size = 1, se = FALSE) + 
      ggplot2::facet_grid(parent~LG, scales = "free_x", space = "free_x", labeller = ggplot2::labeller(parent=ggplot2::as_labeller(c(P=P, Q=Q)))) +
      ggplot2::geom_hline(yintercept = 1/(prod(choose(seq(2,m,2),2))/factorial(m/2)), linetype="dashed") +
      ggplot2::ylim(min.y.prof,max.y.prof) + 
      ggplot2::scale_x_continuous(breaks = seq(0, max(x$prefpair.psi$map.position), 100), expand = ggplot2::expansion(add = 15)) +
      ggplot2::labs(subtitle = "Linkage group", y = "Probability", x = ggplot2::element_blank(), col = "Pairing\nconfig.") + 
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1), panel.spacing = ggplot2::unit(0, "lines"), plot.subtitle = ggplot2::element_text(hjust = 0.5), axis.title.y = ggplot2::element_text(margin = ggplot2::margin(t = 0, r = 5, b = 0, l = 0))) 
    DF<-reshape2::melt(data = x$prefpair.psi.pval, measure.vars = c(P, Q))
    p2<-ggplot2::ggplot(DF, ggplot2::aes(map.position, -log10(value), colour = variable)) +
      ggplot2::geom_point(alpha = .7, size = 1) +  
      ggplot2::facet_grid(.~LG, scales = "free_x", space = "free_x") +
      ggplot2::geom_hline(yintercept = -log10(thresh), linetype="dashed") + 
      ggplot2::labs(y = expression(paste("-",log[10],"(",italic(P),")")), x = "Map Position (cM)") +
      ggplot2::scale_x_continuous(breaks = seq(0, max(x$prefpair.psi$map.position), 100), expand = ggplot2::expansion(add = 15)) +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1), panel.spacing = ggplot2::unit(0, "lines"), legend.box.margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 3)) 
    if(m <6 ){
      p1 <- p1 + ggsci::scale_color_d3()
      p2 <- p2 + ggplot2::scale_color_manual(values=c(ggsci::pal_d3(alpha = 0.8)(nlevels(x$prefpair.psi$pair.conf))[c(1,(nlevels(x$prefpair.psi$pair.conf)/2)+1)]), labels = c(P, Q), name = "Parents")
    }
    p <- ggpubr::ggarrange(p1, p2, ncol = 1, nrow = 2, labels = c("A", "B"), heights = c(1.5,1))
    print(p)
  } else if(type == "hom.pairs"){
    m<-x$info$m
    hom.pair <- phh <- variable <- p.val <- map.position <- probability <- colour <- pair.conf <-NULL
    p1<-ggplot2::ggplot(x$prefpair.homolog) + 
      ggplot2::geom_smooth(ggplot2::aes(map.position, phh, colour = hom.pair), size = 1, se = FALSE) + 
      ggplot2::facet_grid(parent~LG, scales = "free_x", space = "free_x", labeller = ggplot2::labeller(parent=ggplot2::as_labeller(c(P=P, Q=Q)))) +
      ggplot2::geom_hline(yintercept = get_w_m(m - 2)/get_w_m(m), linetype="dashed") +
      ggplot2::ylim(min.y.prof,max.y.prof) + 
      ggplot2::scale_x_continuous(breaks = seq(0, max(x$prefpair.psi$map.position), 100)) +
      ggplot2::labs(subtitle = "Linkage group", y = "Probability", x = ggplot2::element_blank(), col = "Homolog\npairs") +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1), panel.spacing = ggplot2::unit(0, "lines"), plot.subtitle = ggplot2::element_text(hjust = 0.5), axis.title.y = ggplot2::element_text(margin = ggplot2::margin(t = 0, r = 5, b = 0, l = 0))) 
    p2<-ggplot2::ggplot(x$prefpair.homolog) + 
      ggplot2::geom_point(ggplot2::aes(map.position, -log10(p.val), colour = hom.pair)) + 
      ggplot2::facet_grid(parent~LG, scales = "free_x", space = "free_x", labeller = ggplot2::labeller(parent=ggplot2::as_labeller(c(P=P, Q=Q)))) +
      ggplot2::geom_hline(yintercept = -log10(thresh), linetype="dashed") +  
      ggplot2::labs(y = expression(paste("-",log[10],"(",italic(P),")")), x = "Map Position (cM)") +
      ggplot2::scale_x_continuous(breaks = seq(0, max(x$prefpair.psi$map.position), 100)) +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1), panel.spacing = ggplot2::unit(0, "lines")) 
    # p<-gridExtra::grid.arrange(p1, p2, nrow = 2)
    p <- ggpubr::ggarrange(p1, p2, ncol = 1, nrow = 2, labels = c("A", "B"), common.legend = TRUE, legend = "right")
    print(p)
  }
}

