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
#'     \doi{10.1534/g3.119.400620} 
#'     
#' @export
#' @importFrom ggplot2 ggplot geom_hline theme geom_smooth ggtitle facet_grid theme_minimal ylab xlab aes vars scale_color_manual
#' @importFrom reshape2 melt
#' @importFrom ggpubr ggarrange
#' @importFrom ggsci scale_color_d3 pal_d3
#' 
calc_prefpair_profiles <- function(input.genoprobs, verbose = TRUE){
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
    Gnames <- dimnames(input.genoprobs[[j]]$probs)[[1]]
    x <- dim(input.genoprobs[[j]]$probs)
    # number of states, positions and individuals
    nsta <- x[1]
    npos <- x[2]
    n.ind <- x[3]
    # get ploidy level based on the number of states
    v <- c(2,4,6,8,10,12)
    names(v) <- choose(v,v/2)^2
    ploidy <- v[as.character(nsta)]
    ## Possible pairing configurations (Psi's) for both parents
    Psi <- list(P1 = NULL, P2 = NULL)
    aP1 <- combn(letters[1:ploidy],ploidy/2)
    aP2 <- combn(letters[(1+ploidy):(2*ploidy)],ploidy/2)
    for(i in 1:ncol(aP1))
    {
      bP1 <- perm_tot(setdiff(letters[1:ploidy], aP1[,i]))
      bP2 <- perm_tot(setdiff(letters[(1+ploidy):(2*ploidy)], aP2[,i]))
      for(k in 1:nrow(bP1))
      {
        Psi$P1 <- rbind(Psi$P1, sort(apply(cbind(aP1[,i],bP1[k,]), 1, function(x) paste(sort(x), collapse = ""))))
        Psi$P2 <- rbind(Psi$P2, sort(apply(cbind(aP2[,i],bP2[k,]), 1, function(x) paste(sort(x), collapse = ""))))
      }
    }
    Psi$P1 <- unique(Psi$P1)
    Psi$P2 <- unique(Psi$P2)
    # number of possible bivalent configurations (size of Psi)
    nbiv <- nrow(Psi$P1)
    ## Probability of a pairing configuration (Psi) given a genitypic state (G) for both parents
    Psi_given_G <- list(P1 = matrix(0, length(Gnames), nrow(Psi$P1), dimnames = list(Gnames, apply(Psi$P1, 1, paste, collapse = "/"))),
                      P2 = matrix(0, length(Gnames), nrow(Psi$P2), dimnames = list(Gnames, apply(Psi$P2, 1, paste, collapse = "/"))))
    ## Genotypic states
    G <- list(P1 = unlist(strsplit(substr(Gnames, 1, ploidy/2), "")),
            P2 = unlist(strsplit(substr(Gnames, 2+(ploidy/2), ploidy+1), "")))
    dim(G$P2) <- dim(G$P1) <- c(ploidy/2,length(Gnames))
    for(i in 1:nrow(Psi_given_G$P1))
    {
      xP1 <- apply(Psi$P1, 1, function(x,y) all(grepl(paste(y, collapse = "|"), x)), y = G$P1[,i])
      xP2 <- apply(Psi$P2, 1, function(x,y) all(grepl(paste(y, collapse = "|"), x)), y = G$P2[,i])
      #print(as.character(sum(xP)))
      #print(as.character(factorial(ploidy/2)))
      #print("-----")
      Psi_given_G$P2[i,xP2] <- Psi_given_G$P1[i,xP1] <- 1/factorial(ploidy/2)
    }
    ## Equation 3 in Mollinari et al. (2019)
    AP2 <- AP1 <- array(NA, dim = c(nbiv,npos,n.ind))
    for(i in 1:n.ind){
      AP1[,,i] <- crossprod(Psi_given_G$P1,input.genoprobs[[j]]$probs[,,i])   
      AP2[,,i] <- crossprod(Psi_given_G$P2,input.genoprobs[[j]]$probs[,,i])   
    }
    ## Pr(psi_i | O_1, ..., Oz, lambda)
    P1 <- apply(AP1, MARGIN = c(1,2), mean)
    P2 <- apply(AP2, MARGIN = c(1,2), mean)
    dimnames(P1) <- list(colnames(Psi_given_G$P1), names(input.genoprobs[[j]]$map))
    dimnames(P2) <- list(colnames(Psi_given_G$P2), names(input.genoprobs[[j]]$map))
    df.prefpair.temp <- rbind(data.frame(reshape2::melt(P1), parent = "P1", lg = j),
                            data.frame(reshape2::melt(P2), parent = "P2", lg = j))
    colnames(df.prefpair.temp) <- c("pair.conf", "marker", "probability", "parent", "LG")
    map <- data.frame(map.position = input.genoprobs[[j]]$map, marker = names(input.genoprobs[[j]]$map))
    df.prefpair.temp <- merge(df.prefpair.temp, map, sort = FALSE)
    df.prefpair <- rbind(df.prefpair, df.prefpair.temp)
    df.prefpair.pval.temp <- data.frame(p.val.P1 = apply(n.ind * P1, 2, function(x) chisq.test(x)$p.value),
                                      p.val.P2 = apply(n.ind * P2, 2, function(x) chisq.test(x)$p.value),
                                      LG = j,
                                      map.position = input.genoprobs[[j]]$map)
    ## Probability of bivalent pairs given the complete map
    h.h_prime.P1 <- apply(combn(letters[1:ploidy], 2), 2, paste0, collapse = "")
    h.h_prime.P2 <- apply(combn(letters[(ploidy+1):(2*ploidy)], 2), 2, paste0, collapse = "")
    p1 <- get_w_m(ploidy - 2)/get_w_m(ploidy)
    p <- c(p1,1 - p1)
    df.hom.pair <- NULL
    for(k1 in h.h_prime.P1){
      phh <- apply(P1[grepl(k1, rownames(P1)), , drop = FALSE], 2, sum)
      phh.P1 <- data.frame(marker = names(phh), phh = phh, hom.pair = k1, LG = j, parent = "P1")
      x1 <- phh * n.ind
      x <- cbind(x1,n.ind - x1)
      p.val.P1 <- apply(x, 1, function(x) chisq.test(x, p = p)$p.value)
      p.val.P1 <- data.frame(marker = names(p.val.P1), p.val = p.val.P1)
      mtemp <- merge(phh.P1, p.val.P1)
      mtemp <- merge(mtemp, map)
      df.hom.pair <- rbind(df.hom.pair, mtemp)
    }
    for(k1 in h.h_prime.P2){
      phh <- apply(P2[grepl(k1, rownames(P2)), , drop = FALSE], 2, sum)
      phh.P2 <- data.frame(marker = names(phh), phh = phh, hom.pair = k1, LG = j, parent = "P2")
      x1 <- phh * n.ind
      x <- cbind(x1, n.ind - x1)
      p.val.P2 <- apply(x, 1, function(x) chisq.test(x, p = p)$p.value)
      p.val.P2 <- data.frame(marker = names(p.val.P2), p.val = p.val.P2)
      mtemp <- merge(phh.P2, p.val.P2)
      mtemp <- merge(mtemp, map)
      df.hom.pair <- rbind(df.hom.pair, mtemp)
    }
    df.hompair.pval <- rbind(df.hompair.pval, df.hom.pair)
    df.prefpair.pval <- rbind(df.prefpair.pval, df.prefpair.pval.temp)
  }
  if (verbose) cat("\n")
  structure(list(info = list(ploidy = ploidy, n.ind = n.ind), 
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
  cat("  No. individuals:                          ", x$info$n.ind, "\n")
  cat("  -----------------------------------------------------")
  ## printing summary
  parent <- NULL
  cat("\n  Pairing configurations:\n")
  a1 <- as.character(unique(subset(x$prefpair.psi, parent  ==  "P1")$pair.conf))
  a1.h <- apply(combn(letters[1:x$info$ploidy], 2), 2, paste0, collapse = "")
  a2 <- as.character(unique(subset(x$prefpair.psi, parent  ==  "P2")$pair.conf))
  a2.h <- apply(combn(letters[(x$info$ploidy+1):(2*x$info$ploidy)], 2), 2, paste0, collapse = "")
  a1 <- matrix(a1, nrow = 1)
  a2 <- matrix(a2, nrow = 1)
  a1.h <- matrix(a1.h, nrow = 1)
  a2.h <- matrix(a2.h, nrow = 1)
  if(x$info$ploidy  ==  6){
    dim(a2) <- dim(a1) <- c(3, 5)    
    dim(a2.h) <- dim(a1.h) <- c(3, 5)    
  } else if (x$info$ploidy  ==  8){
    dim(a2) <- dim(a1) <- c(21, 5)
    dim(a2.h) <- dim(a1.h) <- c(7, 4)
  } else if (x$info$ploidy  ==  10){
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
#' @param P1 a string containing the name of parent P1
#' 
#' @param P2 a string containing the name of parent P2
#' 
#' @param ... unused arguments
#' @export
plot.mappoly.prefpair.profiles <- function(x, type = c("pair.configs", "hom.pairs"), 
                                           min.y.prof = 0, max.y.prof = 1, 
                                           thresh = 0.01, P1 = "P1", P2 = "P2", ...){
  type <- match.arg(type)
  colnames(x$prefpair.psi.pval)[1:2] <- c(P1, P2)
  if(type  ==  "pair.configs"){
    ploidy <- x$info$ploidy
    variable <- value <- map.position <- probability <- colour <- pair.conf  <- NULL
    p1 <- ggplot2::ggplot(x$prefpair.psi) + 
      ggplot2::geom_smooth(ggplot2::aes(map.position, probability, colour = pair.conf), size = 1, se = FALSE) + 
      ggplot2::facet_grid(parent~LG, scales = "free_x", space = "free_x", labeller = ggplot2::labeller(parent = ggplot2::as_labeller(c(P1 = P1, P2 = P2)))) +
      ggplot2::geom_hline(yintercept = 1/(prod(choose(seq(2,ploidy,2),2))/factorial(ploidy/2)), linetype = "dashed") +
      ggplot2::ylim(min.y.prof,max.y.prof) + 
      ggplot2::scale_x_continuous(breaks = seq(0, max(x$prefpair.psi$map.position), 100), expand = ggplot2::expansion(add = 15)) +
      ggplot2::labs(subtitle = "Linkage group", y = "Probability", x = ggplot2::element_blank(), col = "Pairing\nconfig.") + 
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1), panel.spacing = ggplot2::unit(0, "lines"), plot.subtitle = ggplot2::element_text(hjust = 0.5), axis.title.y = ggplot2::element_text(margin = ggplot2::margin(t = 0, r = 5, b = 0, l = 0))) 
    DF <- reshape2::melt(data = x$prefpair.psi.pval, measure.vars = c(P1, P2))
    p2 <- ggplot2::ggplot(DF, ggplot2::aes(map.position, -log10(value), colour = variable)) +
      ggplot2::geom_point(alpha = .7, size = 1) +  
      ggplot2::facet_grid(.~LG, scales = "free_x", space = "free_x") +
      ggplot2::geom_hline(yintercept = -log10(thresh), linetype = "dashed") + 
      ggplot2::labs(y = expression(paste("-",log[10],"(",italic(P1),")")), x = "Map Position (cM)") +
      ggplot2::scale_x_continuous(breaks = seq(0, max(x$prefpair.psi$map.position), 100), expand = ggplot2::expansion(add = 15)) +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1), panel.spacing = ggplot2::unit(0, "lines"), legend.box.margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 3)) 
    if(ploidy <6 ){
      p1 <- p1 + ggsci::scale_color_d3()
      p2 <- p2 + ggplot2::scale_color_manual(values = c(ggsci::pal_d3(alpha = 0.8)(nlevels(x$prefpair.psi$pair.conf))[c(1,(nlevels(x$prefpair.psi$pair.conf)/2)+1)]), labels = c(P1, P2), name = "Parents")
    }
    p <- ggpubr::ggarrange(p1, p2, ncol = 1, nrow = 2, labels = c("A", "B"), heights = c(1.5,1))
    print(p)
  } else if(type  ==  "hom.pairs"){
    ploidy <- x$info$ploidy
    hom.pair <- phh <- variable <- p.val <- map.position <- probability <- colour <- pair.conf  <- NULL
    p1 <- ggplot2::ggplot(x$prefpair.homolog) + 
      ggplot2::geom_smooth(ggplot2::aes(map.position, phh, colour = hom.pair), size = 1, se = FALSE) + 
      ggplot2::facet_grid(parent~LG, scales = "free_x", space = "free_x", labeller = ggplot2::labeller(parent = ggplot2::as_labeller(c(P1 = P1, P2 = P2)))) +
      ggplot2::geom_hline(yintercept = get_w_m(ploidy - 2)/get_w_m(ploidy), linetype = "dashed") +
      ggplot2::ylim(min.y.prof,max.y.prof) + 
      ggplot2::scale_x_continuous(breaks = seq(0, max(x$prefpair.psi$map.position), 100)) +
      ggplot2::labs(subtitle = "Linkage group", y = "Probability", x = ggplot2::element_blank(), col = "Homolog\npairs") +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1), panel.spacing = ggplot2::unit(0, "lines"), plot.subtitle = ggplot2::element_text(hjust = 0.5), axis.title.y = ggplot2::element_text(margin = ggplot2::margin(t = 0, r = 5, b = 0, l = 0))) 
    p2 <- ggplot2::ggplot(x$prefpair.homolog) + 
      ggplot2::geom_point(ggplot2::aes(map.position, -log10(p.val), colour = hom.pair)) + 
      ggplot2::facet_grid(parent~LG, scales = "free_x", space = "free_x", labeller = ggplot2::labeller(parent = ggplot2::as_labeller(c(P1 = P1, P2 = P2)))) +
      ggplot2::geom_hline(yintercept = -log10(thresh), linetype = "dashed") +  
      ggplot2::labs(y = expression(paste("-",log[10],"(",italic(P1),")")), x = "Map Position (cM)") +
      ggplot2::scale_x_continuous(breaks = seq(0, max(x$prefpair.psi$map.position), 100)) +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1), panel.spacing = ggplot2::unit(0, "lines")) 
    # p <- gridExtra::grid.arrange(p1, p2, nrow = 2)
    p <- ggpubr::ggarrange(p1, p2, ncol = 1, nrow = 2, labels = c("A", "B"), common.legend = TRUE, legend = "right")
    print(p)
  }
}

