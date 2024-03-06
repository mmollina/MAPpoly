#' Extract the LOD Scores in a \code{'mappoly.map'} object
#' @param x an object of class \code{mappoly.map}
#' @param sorted logical. if \code{TRUE}, the LOD Scores are displayed
#'     in a decreasing order
#' @return a numeric vector containing the LOD Scores
#' @keywords internal
#' @export
get_LOD <- function(x, sorted = TRUE) {
  w <- sapply(x$maps, function(x) x$loglike)
  LOD <- w - max(w)
  if (sorted)
    LOD <- sort(LOD, decreasing = TRUE)
  abs(LOD)
}

#' Get recombination fraction from a matrix
#' @keywords internal
get_rf_from_mat <- function(M){
  r <- numeric(nrow(M)-1)
  for(i in 1:(nrow(M)-1)){
    r[i] <- M[i, i+1]
  }
  r
}


#' Check if Object is a Probability Dataset in MAPpoly
#'
#' Determines whether the specified object is a probability dataset
#' by checking for the existence of the 'geno' component within a
#' `"mappoly.data"` object. 
#'
#' @param x An object of class `"mappoly.data"`
#'
#' @return A logical value: `TRUE` if the 'geno' component exists within `x`,
#' indicating it is a valid probability dataset for genetic analysis; `FALSE`
#' otherwise.
#'
#' @keywords internal
#' @export
is.prob.data <- function(x){
  # Verify if 'x' is indeed an object of class 'mappoly.data'
  if (!inherits(x, "mappoly.data")) {
    stop("Input is not an object of class 'mappoly.data'")
  }
  
  # Check for the existence of 'geno' within the 'mappoly.data' object
  exists('geno', where = x)
}


#' Get the number of bivalent configurations
#' @keywords internal
get_w_m <- function(ploidy){
  if(ploidy%%2 != 0) stop("ploidy level should be an even number") 
  if(ploidy <= 0) stop("ploidy level should be greater than zero")
  1/factorial((ploidy/2)) * prod(choose(seq(2, ploidy, 2),2))
}

#' Reverse map
#'
#' Provides the reverse of a given map.
#'
#' @param input.map an object of class \code{mappoly.map}
#' 
#' @examples
#' plot_genome_vs_map(solcap.mds.map[[1]])
#' plot_genome_vs_map(rev_map(solcap.mds.map[[1]]))
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu}
#'
#'@export
#'
rev_map <- function(input.map)
{
  output.map <- input.map
  output.map$info <- lapply(output.map$info, rev)
  for(i in 1:length(output.map$maps))
  {
    output.map$maps[[i]]$seq.num <- rev(input.map$maps[[i]]$seq.num)
    output.map$maps[[i]]$seq.rf <- rev(input.map$maps[[i]]$seq.rf)
    output.map$maps[[i]]$seq.ph$P <- rev(input.map$maps[[i]]$seq.ph$P)
    output.map$maps[[i]]$seq.ph$Q <- rev(input.map$maps[[i]]$seq.ph$Q)
  }
  return(output.map)
}

#' Returns the class with the highest probability in 
#' a genotype probability distribution
#'
#' @param geno the probabilistic genotypes contained in the object
#'     \code{'mappoly.data'}
#' @param prob.thres probability threshold to select the genotype. 
#'     Values below this genotype are assumed as missing data
#' @return a matrix containing the doses of each genotype and
#'     marker. Markers are disposed in rows and individuals are 
#'     disposed in columns. Missing data are represented by NAs
#' @keywords internal
#' @examples
#' \donttest{
#' geno.dose <- dist_prob_to_class(hexafake.geno.dist$geno)
#' geno.dose$geno.dose[1:10,1:10]
#'}   
#' @importFrom magrittr "%>%"
#' @importFrom reshape2 melt dcast
#' @importFrom dplyr group_by filter arrange
#' @export
dist_prob_to_class <- function(geno, prob.thres = 0.9) {
  a <- reshape2::melt(geno, id.vars = c("mrk", "ind"))
  mrk <- ind <- value <- variable <- NULL # Setting the variables to NULL first
  a$variable <- as.numeric(levels(a$variable))[a$variable]
  b <- a %>%
    dplyr::group_by(mrk, ind) %>%
    dplyr::filter(value > prob.thres) %>%
    dplyr::arrange(mrk, ind, variable)
  z <- reshape2::dcast(data = b[,1:3], formula = mrk ~ ind, value.var = "variable")
  rownames(z) <- z[,"mrk"]
  z <- data.matrix(frame = z[,-1])
  n <- setdiff(unique(geno$mrk), rownames(z))
  if(length(n) > 0)
  {
    ploidy <- matrix(NA, nrow = length(n), ncol = ncol(z), dimnames = list(n, colnames(z)))
    z <- rbind(z,ploidy)
  }
  rm.ind <- setdiff(unique(geno$ind), colnames(z))
  flag <- FALSE
  if(length(rm.ind) > 0){
    flag <- TRUE
    warning("Inividual(s) ", paste(rm.ind, collapse = " "), 
            "\n  did not meet the 'prob.thres' criteria for any of\n  the markers and was (were) removed.")
    geno <- geno %>% dplyr::filter(ind %in% colnames(z))
  }
  z <- z[as.character(unique(geno$mrk)), as.character(unique(geno$ind))]
  list(geno.dose = z, geno = geno, flag = flag)
}

#' Export data to \code{polymapR}
#' 
#' See examples at \url{https://rpubs.com/mmollin/tetra_mappoly_vignette}.
#' 
#' @param data.in an object of class \code{mappoly.data}
#' @return a dosage \code{matrix} 
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu}
#' @export export_data_to_polymapR
export_data_to_polymapR <- function(data.in)
{
  data.out <- as.matrix(data.frame(P1 = data.in$dosage.p1,     
                                   P2 = data.in$dosage.p2,
                                   data.in$geno.dose))
  data.out[data.out  ==  (data.in$ploidy + 1)] <- NA
  return(data.out)
}

#' Msg function
#' @keywords internal
#' @importFrom cli rule
msg <- function(text, line = 1){
  cli::rule(line = line, right = text) %>%
    text_col() %>%
    message()
}

#' @importFrom rstudioapi isAvailable hasFun getThemeInfo
#' @importFrom crayon white black
text_col <- function(x) {
  # If RStudio not available, messages already printed in black
  if (!rstudioapi::isAvailable()) {
    return(x)
  }
  if (!rstudioapi::hasFun("getThemeInfo")) {
    return(x)
  }
  theme <- rstudioapi::getThemeInfo()
  if (isTRUE(theme$dark)) crayon::white(x) else crayon::black(x)
}

#' Genetic Mapping Functions
#'
#' These functions facilitate the conversion between recombination fractions (r) and genetic distances (d) 
#' using various mapping models. The functions starting with `mf_` convert recombination fractions to genetic distances,
#' while those starting with `imf_` convert genetic distances back into recombination fractions.
#'
#' @name genetic-mapping-functions
#' @aliases mf_k mf_h mf_m imf_k imf_h imf_m
#' @usage mf_k(d)
#' @usage mf_h(d)
#' @usage mf_m(d)
#' @usage imf_k(r)
#' @usage imf_h(r)
#' @usage imf_m(r)
#' @param d Numeric or numeric vector, representing genetic distances in centiMorgans (cM) for direct functions (mf_k, mf_h, mf_m).
#' @param r Numeric or numeric vector, representing recombination fractions for inverse functions (imf_k, imf_h, imf_m).
#' @details
#' The `mf_` prefixed functions apply different models to convert recombination fractions into genetic distances:
#' \itemize{
#'   \item \code{mf_k}: Kosambi mapping function.
#'   \item \code{mf_h}: Haldane mapping function.
#'   \item \code{mf_m}: Morgan mapping function.
#'}
#' The `imf_` prefixed functions convert genetic distances back into recombination fractions:
#' \itemize{
#'   \item \code{imf_k}: Inverse Kosambi mapping function.
#'   \item \code{imf_h}: Inverse Haldane mapping function.
#'   \item \code{imf_m}: Inverse Morgan mapping function.
#'}
#' @references
#' Kosambi, D.D. (1944). The estimation of map distances from recombination values. Ann Eugen., 12, 172-175.
#' Haldane, J.B.S. (1919). The combination of linkage values, and the calculation of distances between the loci of linked factors. J Genet, 8, 299-309.
#' Morgan, T.H. (1911). Random segregation versus coupling in Mendelian inheritance. Science, 34(873), 384.
#' @keywords genetics
#' @export
mf_k <- function(d) 0.5 * tanh(d / 50)
mf_h <- function(d) 0.5 * (1 - exp(-d / 50))
mf_m <- function(d) sapply(d, function(a) min(a / 100, 0.5))
imf_k <- function(r) {
  r[r >= 0.5] <- 0.5 - 1e-14
  50 * atanh(2 * r)
}
imf_h <- function(r) {
  r[r >= 0.5] <- 0.5 - 1e-14
  -50 * log(1 - 2 * r)
}
imf_m <- function(r) sapply(r, function(a) min(a * 100, 50))

#' Compare two polyploid haplotypes stored in list format
#'
#' @param ploidy ploidy level
#' @param h1 homology group 1
#' @param h2 homology group 2
#' @return a numeric vector of size \code{ploidy} indicating which
#'     homolog in h2 represents the homolog in h1. If there is
#'     no correspondence, i.e. different homolog, it returns NA for
#'     that homolog.
#' @keywords internal
#' @export compare_haplotypes
#' 
compare_haplotypes <- function(ploidy, h1, h2) {
  I1 <- matrix(0, ploidy, length(h1))
  I2 <- matrix(0, ploidy, length(h2))
  for (i in 1:length(h1)) {
    I1[h1[[i]], i] <- 1
    I2[h2[[i]], i] <- 1
  }
  a <- apply(I1, 1, paste, collapse = "")
  b <- apply(I2, 1, paste, collapse = "")
  haplo.ord <- match(a, b)
  list(is.same.haplo = !any(is.na(haplo.ord)), haplo.ord = haplo.ord)
}

#' Genotypic information content 
#' 
#' This function plots the genotypic information content given 
#' an object of class \code{mappoly.homoprob}.
#' 
#' @param hprobs an object of class \code{mappoly.homoprob}
#' 
#' @param P a string containing the name of parent P
#' 
#' @param Q a string containing the name of parent Q
#' 

#' 
#' @examples
#' \donttest{
#'      w <- lapply(solcap.err.map[1:3], calc_genoprob)
#'      h.prob <- calc_homologprob(w)
#'      plot_GIC(h.prob)
#' }
#' 
#' @importFrom dplyr mutate 
#' @importFrom ggplot2 facet_wrap element_text ylim scale_color_discrete

#' @export plot_GIC 
plot_GIC <- function(hprobs, P = "P1", Q = "P2"){
  if(!inherits(hprobs, "mappoly.homoprob"))
    stop(" 'hprobs' should be of class 'mappoly.homoprob'")
  LG <- map.position <- GIC <- p1 <- m1 <- homolog <- marker <- probability <- NULL
  DF <- hprobs$homoprob %>% 
    dplyr::mutate(p1 = probability * (1-probability)) %>%
    group_by(marker, homolog, map.position, LG) %>%
    summarise(m1 = sum(p1)) %>%
    mutate(GIC = 1-(4/hprobs$info$n.ind)*m1 , parent = ifelse(homolog%in%letters[1:hprobs$info$ploidy], P, Q))
  head(as.data.frame(DF))
  
  print(ggplot(DF, aes(map.position, GIC, colour = homolog)) +
          geom_smooth(alpha = .8, se = FALSE) + facet_wrap(~LG, nrow = 3, ncol = 5) +
          facet_grid(parent~LG, scales = "free_x", space = "free_x") +
          geom_hline(yintercept = .8, linetype = "dashed") + ylim(0,1) +
          theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
          scale_color_discrete(name = "Homologs") +
          ylab("Genotypic Information Content") + xlab("Distance (cM)"))
  return(invisible(DF))
}

#' Plot Two Overlapped Haplotypes
#'
#' This function plots two sets of haplotypes for comparison, allowing for visual
#' inspection of homologous allele patterns across two groups or conditions.
#' It is designed to handle and display genetic data for organisms with varying ploidy levels.
#'
#' @param ploidy Integer, specifying the ploidy level of the organism being represented.
#' @param hom.allele.p1 A list where each element represents the alleles for a marker in the first haplotype group, for 'p' parent.
#' @param hom.allele.q1 A list where each element represents the alleles for a marker in the first haplotype group, for 'q' parent.
#' @param hom.allele.p2 Optionally, a list where each element represents the alleles for a marker in the second haplotype group, for 'p' parent.
#' @param hom.allele.q2 Optionally, a list where each element represents the alleles for a marker in the second haplotype group, for 'q' parent.
#'
#' @details The function creates a graphical representation of haplotypes, where each marker's alleles are plotted
#' along a line for each parent ('p' and 'q'). If provided, the second set of haplotypes (for comparison) are overlaid
#' on the same plot. This allows for direct visual comparison of allele presence or absence across the two sets.
#' Different colors are used to distinguish between the first and second sets of haplotypes.
#'
#' The function uses several internal helper functions (`ph_list_to_matrix` and `ph_matrix_to_list`) to manipulate
#' haplotype data. These functions should correctly handle the conversion between list and matrix representations
#' of haplotype data.
#'
#' @return Invisible. The function primarily generates a plot for visual analysis and does not return any data.
#' @export
#' @keywords internal
plot_compare_haplotypes <- function(ploidy, hom.allele.p1, hom.allele.q1, hom.allele.p2 = NULL, hom.allele.q2 = NULL) {
  nmmrk <- names(hom.allele.p1)
  oldpar <- par(mar = c(5.1, 4.1, 4.1, 2.1))
  on.exit(par(oldpar))
  o1 <- order(apply(ph_list_to_matrix(hom.allele.p1, ploidy), 2, paste, collapse = ""), decreasing = TRUE)
  hom.allele.p1 <- ph_matrix_to_list(ph_list_to_matrix(hom.allele.p1, ploidy)[, o1])
  o2 <- order(apply(ph_list_to_matrix(hom.allele.q1, ploidy), 2, paste, collapse = ""), decreasing = TRUE)
  hom.allele.q1 <- ph_matrix_to_list(ph_list_to_matrix(hom.allele.q1, ploidy)[, o2])
  
  if (!is.null(hom.allele.p2)) {
    o3 <- order(apply(ph_list_to_matrix(hom.allele.p2, ploidy), 2, paste, collapse = ""), decreasing = TRUE)
    hom.allele.p2 <- ph_matrix_to_list(ph_list_to_matrix(hom.allele.p2, ploidy)[, o3])
  }
  if (!is.null(hom.allele.q2)) {
    o4 <- order(apply(ph_list_to_matrix(hom.allele.q2, ploidy), 2, paste, collapse = ""), decreasing = TRUE)
    hom.allele.q2 <- ph_matrix_to_list(ph_list_to_matrix(hom.allele.q2, ploidy)[, o4])
  }
  col2 <- 0  #rgb(1,0,.5,0.5)
  col1 <- rgb(1, 0, 0, 0.5)
  col3 <- rgb(0, 0, 1, 0.5)
  n.mrk <- length(hom.allele.p1)
  plot(c(0, 22), c(0, -(ploidy + 15)), type = "n", axes = FALSE, xlab = "", main = "", ylab = "")
  for (i in -(1:ploidy)) {
    lines(c(0, 10), c(i, i), lwd = 1, col = "darkgray", lty = 2)
    lines(c(12, 22), c(i, i), lwd = 1, col = "darkgray", lty = 2)
  }
  pos.p <- cumsum(c(0, rep(1, n.mrk - 1)/sum(rep(1, n.mrk - 1))) * 10)
  for (i in 1:n.mrk) {
    points(x = rep(pos.p[i], ploidy), y = -c(1:ploidy), pch = 15, col = col2, cex = 2)
    if (any(hom.allele.p1[[i]] != 0))
      points(x = rep(pos.p[i], length(hom.allele.p1[[i]])), y = -hom.allele.p1[[i]], col = col1, pch = 15, cex = 2)
    if (any(hom.allele.p2[[i]] != 0))
      if (!is.null(hom.allele.p2))
        points(x = rep(pos.p[i], length(hom.allele.p2[[i]])), y = -hom.allele.p2[[i]], col = col3, pch = 15, cex = 2)
  }
  text(x = pos.p, y = rep(-(ploidy+1), length(pos.p)), labels = nmmrk, cex = .5)
  pos.q <- pos.p + 12
  for (i in 1:n.mrk) {
    points(x = rep(pos.q[i], ploidy), y = -c(1:ploidy), col = col2, pch = 15, cex = 2)
    if (any(hom.allele.q1[[i]] != 0))
      points(x = rep(pos.q[i], length(hom.allele.q1[[i]])), y = -hom.allele.q1[[i]], col = col1, pch = 15, cex = 2)
    if (any(hom.allele.q2[[i]] != 0))
      if (!is.null(hom.allele.q2))
        points(x = rep(pos.q[i], length(hom.allele.q2[[i]])), y = -hom.allele.q2[[i]], col = col3, pch = 15, cex = 2)
  }
  text(x = pos.q, y = rep(-(ploidy+1) , length(pos.q)), labels = nmmrk, cex = .5)
  text(x = 11, y = -(ploidy + 1)/2, labels = "X", cex = 1)
}

#' Check if it is possible to estimate the recombination
#' fraction between neighbor markers using two-point
#' estimation
#' @keywords internal
check_if_rf_is_possible <- function(input.seq){
  dp <- abs(abs(input.seq$seq.dose.p1-(input.seq$ploidy/2))-(input.seq$ploidy/2))
  dq <- abs(abs(input.seq$seq.dose.p2-(input.seq$ploidy/2))-(input.seq$ploidy/2))
  y <- numeric(length(input.seq$seq.num)-1)
  for(i in 1:length(y))
    y[i] <- (dp[i]  ==  0 && dq[i+1]  ==  0) ||
    (dp[i+1]  ==  0 && dq[i]  ==  0) ||
    (dp[i]  ==  0 && dq[i]  ==  0) ||
    (dp[i+1]  ==  0 && dq[i+1]  ==  0)
  y <- as.logical(y)
  !y
}

#' N! combination
#' @keywords internal
perm_tot <- function(v) {
  n <- length(v)
  result <- v
  if (n > 1) {
    M <- perm_tot(v[2:n])
    result <- cbind(v[1], M)
    if (n > 2) {
      for (i in 2:(n - 1)) {
        N <- cbind(M[, 1:(i - 1)], v[1], M[, i:(n - 1)])
        result <- rbind(result, N)
      }
    }
    N <- cbind(M, v[1])
    result <- rbind(result, N)
  }
  return(result)
}

#' N!/2 combination
#' @keywords internal
perm_pars <- function(v) {
  n <- length(v)
  result <- v
  if (n > 2) {
    Mt <- perm_tot(v[2:n])
    result <- cbind(v[1], Mt)
    f <- floor(n/2)
    c <- ceiling(n/2)
    if (n > 3) {
      for (i in 2:f) {
        N <- cbind(Mt[, 1:(i - 1)], v[1], Mt[, i:(n - 1)])
        result <- rbind(result, N)
      }
    }
    if (c > f) {
      Ms <- perm_pars(v[2:n])
      if (n > 3) {
        N <- cbind(Ms[, 1:f], v[1], Ms[, c:(n - 1)])
      } else {
        N <- cbind(Ms[1:f], v[1], Ms[c:(n - 1)])
      }
      result <- rbind(result, N)
    }
  }
  return(result)
}

#' Color pallet ggplot-like
#' @keywords internal
#' @importFrom grDevices hcl col2rgb hsv rgb2hsv
gg_color_hue <- function(n) {
  x <- rgb2hsv(col2rgb("steelblue"))[, 1]
  cols = seq(x[1], x[1] + 1, by = 1/n)
  cols = cols[1:n]
  cols[cols > 1] <- cols[cols > 1] - 1
  return(hsv(cols, x[2], x[3]))
}

#' MAPpoly Color Palettes
#'
#' Provides a set of color palettes designed for use with MAPpoly, 
#' a package for genetic mapping in polyploids. These palettes are 
#' intended to enhance the visual representation of genetic data.
#'
#' The available palettes are:
#' \describe{
#'   \item{\code{mp_pallet1}}{A palette with warm colors ranging from yellow to dark red and brown.}
#'   \item{\code{mp_pallet2}}{A palette with cool colors, including purples, blues, and green.}
#'   \item{\code{mp_pallet3}}{A comprehensive palette that combines colors from both \code{mp_pallet1} and \code{mp_pallet2}, offering a broad range of colors.}
#'}
#'
#' Each palette function returns a function that can generate color vectors of variable length, suitable for mapping or plotting functions in R.
#'
#' @name mappoly-color-palettes
#' @aliases mp_pallet1 mp_pallet2 mp_pallet3
#' @keywords internal
#' @export mp_pallet1 mp_pallet2 mp_pallet3
#' @examples
#' # Generate a palette of 5 colors from mp_pallet1
#' pal1 <- mp_pallet1(5)
#' plot(1:5, pch=19, col=pal1)
#'
#' # Generate a palette of 10 colors from mp_pallet2
#' pal2 <- mp_pallet2(10)
#' plot(1:10, pch=19, col=pal2)
#'
#' # Generate a palette of 15 colors from mp_pallet3
#' pal3 <- mp_pallet3(15)
#' plot(1:15, pch=19, col=pal3)

mp_pallet1 <- colorRampPalette(c("#ffe119", "#f58231","#e6194b","#808000","#9a6324", "#800000"))
mp_pallet2 <- colorRampPalette(c("#911eb4", "#000075","#4363d8","#42d4f4","#469990", "#3cb44b"))
mp_pallet3 <- colorRampPalette(c("#ffe119", "#f58231","#e6194b","#808000","#9a6324", "#800000","#911eb4", "#000075","#4363d8","#42d4f4","#469990", "#3cb44b"))

#' Update missing information
#' @keywords internal
update_missing <- function(input.data, prob.thres = 0.95){
  geno.dose <- dist_prob_to_class(geno = input.data$geno, prob.thres = prob.thres)
  if(geno.dose$flag)
  {
    geno <- geno.dose$geno
    geno.dose <- geno.dose$geno.dose
  } else {
    geno.dose <- geno.dose$geno.dose
  }
  geno.dose[is.na(geno.dose)] <- input.data$ploidy + 1
  input.data$geno.dose <- geno.dose
  input.data$prob.thres <- prob.thres
  return(input.data)
}

#' Chi-square test
#' @keywords internal
mrk_chisq_test <- function(x, ploidy){
  y <- x[-c(1:(ploidy+1))]
  y[y == ploidy+1] <- NA
  y <- table(y, useNA = "always")
  names(y) <- c(names(y)[-length(y)], "NA") 
  seg.exp <- x[0:(ploidy+1)]
  seg.exp <- seg.exp[seg.exp != 0]
  seg.obs <- seg.exp
  seg.obs[names(y)[-length(y)]] <- y[-length(y)]
  pval <- suppressWarnings(stats::chisq.test(x = seg.obs, p = seg.exp[names(seg.obs)])$p.value)
  pval
}

#' Get the genomic position of markers in a sequence
#'
#' This functions gets the genomic position of markers in a sequence and
#' return an ordered data frame with the name and position of each marker
#' @param input.seq a sequence object of class \code{mappoly.sequence}
#' @param verbose if \code{TRUE} (default), the current progress is shown; if
#'     \code{FALSE}, no output is produced
#' @param x 	an object of the class mappoly.geno.ord
#' @param ... 	currently ignored   
#'   
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu}
#' @examples
#' s1 <- make_seq_mappoly(tetra.solcap, "all")
#' o1 <- get_genomic_order(s1)
#' plot(o1)
#' s.geno.ord <- make_seq_mappoly(o1)
#' @export get_genomic_order
get_genomic_order <- function(input.seq, verbose = TRUE){
  if (!inherits(input.seq, "mappoly.sequence")) {
    stop(deparse(substitute(input.seq)), 
         " is not an object of class 'mappoly.sequence'")
  }
  if(all(is.na(input.seq$genome.pos))){
    if(all(is.na(input.seq$chrom))) 
      stop("No sequence or sequence position information found.")
    else{
      if (verbose) message("Ordering markers based on chromosome information")
      M <- data.frame(seq = input.seq$chrom, row.names = input.seq$seq.mrk.names)
      M.out <- M[order(M[,1]),]
    }
  } else if(all(is.na(input.seq$chrom))){
    if(all(is.na(input.seq$genome.pos))) 
      stop("No sequence or sequence position information found.")
    else{
      if (verbose) message("Ordering markers based on sequence position information")
      M <- data.frame(seq.pos = input.seq$genome.pos, row.names = input.seq$seq.mrk.names)
      M.out <- M[order(as.numeric(M[,1])),]
    }
  } else{
    M <- data.frame(seq = input.seq$chrom, 
                    seq.pos = input.seq$genome.pos, 
                    row.names = input.seq$seq.mrk.names)
    M.out <- M[order(as.numeric(M[,1]), as.numeric(M[,2])),]
  }
  structure(list(data.name = input.seq$data.name, ord = M.out), class = "mappoly.geno.ord")
}

#' @rdname get_genomic_order
#' @export
print.mappoly.geno.ord <- function(x, ...){
  print(head(x$ord))
}

#' @rdname get_genomic_order
#' @export
#' @importFrom ggplot2 ggplot aes geom_point xlab ylab theme
plot.mappoly.geno.ord <- function(x, ...){
  seq <- seq.pos <- NULL
  w <- x$ord
  w$seq <- as.factor(w$seq)
  w$seq = with(w, reorder(seq))
  ggplot2::ggplot(w, 
                  ggplot2::aes(x = seq.pos, y = seq, group = as.factor(seq))) +
    ggplot2::geom_point(ggplot2::aes(color = as.factor(seq)), shape = 108, size = 5, show.legend = FALSE) +
    ggplot2::xlab("Position") + 
    ggplot2::ylab("Chromosome") 
}

#' Data sanity check
#'
#' Checks the consistency of a dataset
#'
#' @param x an object of class \code{mappoly.data}
#' 
#' @return if consistent, returns 0. If not consistent, returns a 
#'         vector with a number of tests, where \code{TRUE} indicates
#'         a failed test.
#' 
#' @examples
#' check_data_sanity(tetra.solcap)
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu}
#'
#' @references
#'     Mollinari, M., and Garcia, A.  A. F. (2019) Linkage
#'     analysis and haplotype phasing in experimental autopolyploid
#'     populations with high ploidy level using hidden Markov
#'     models, _G3: Genes, Genomes, Genetics_. 
#'     \doi{10.1534/g3.119.400378}
#'     
#' @export check_data_sanity
check_data_sanity <- function(x){
  if(exists('geno', where = x)){
    check_data_dist_sanity(x)
  } else if (exists('geno.dose', where = x)){
    check_data_dose_sanity(x)
  } else
    stop("Inconsistent genotypic information.")
}

#' Checks the consistency of dataset (dosage)
#' @keywords internal
check_data_dose_sanity <- function(x){
  test <- logical(24L)
  names(test) <- 1:24
  
  # ploidy
  test[1] <- x$ploidy%%2 != 0 #is ploidy even?
  test[2] <- any(sapply(x$dosage.p1, function(y) max(y) > x$ploidy | min(y) < 0)) #are dosages in P higher than ploidy?
  test[3] <- any(sapply(x$dosage.p2, function(y) max(y) > x$ploidy | min(y) < 0)) #are dosages in Q higher than ploidy?
  test[4] <- max(x$geno.dose) > x$ploidy + 1 #is there any dose in offspring higher than ploidy?
  test[5] <- min(x$geno.dose) < 0 #is there any negative dose in offspring?
  
  # number of individuals
  test[6] <- x$n.ind < 0 #is the number of individuals greater than zero?
  test[7] <- length(x$ind.names) != x$n.ind #is the number of individual names equal to the number of individuals?
  test[8] <- ncol(x$geno.dose) != x$n.ind #is the number of columns in the dosage matrix equal to the number of individuals?
  
  # number of markers
  test[9] <- x$n.mrk < 0 #is the number of markers greater than zero?
  test[10] <- length(x$mrk.names) != x$n.mrk #is the number of marker names equal to the number of markers?
  test[11] <- length(x$dosage.p1) != x$n.mrk #is the number of marker dosages in P equal to the number of markers?
  test[12] <- length(x$dosage.p2) != x$n.mrk #is the number of marker dosages in Q equal to the number of markers?
  if(length(x$chrom) > 0)
    test[13] <- length(x$chrom) != x$n.mrk #is the number of sequences equal to the number of markers?
  if(length(x$genome.pos) > 0)
    test[14] <- length(x$genome.pos) != x$n.mrk #is the number of sequence positions equal to the number of markers?
  test[15] <- nrow(x$geno.dose) != x$n.mrk #is the number of rows in the dosage matrix equal to the number of markers?
  if(length(x$chisq.pval) > 0)
    test[16] <- length(x$chisq.pval) != x$n.mrk #is the number of chi-square tests equal to the number of markers?
  
  # individual names in the dosage dataset
  test[17] <- !is.character(x$ind.names) # are individual's names characters
  test[18] <- !identical(colnames(x$geno.dose), x$ind.names) # are column names in dosage matrix identical to individual names?
  
  # markers names in the dosage dataset
  test[19] <- !is.character(x$mrk.names)# are marker's names characters
  test[20] <- !identical(rownames(x$geno.dose), x$mrk.names)# are row names in dosage matrix identical to marker names?
  
  # dosage in both parents
  test[21] <- !is.integer(x$dosage.p1) # are dosages in P numeric
  test[22] <- !is.integer(x$dosage.p2) # are dosages in Q numeric 
  test[23] <- !identical(names(x$dosage.p1), x$mrk.names) # are names in P's dosage vector identical to marker names?
  test[24] <- !identical(names(x$dosage.p2), x$mrk.names) # are names in Q's dosage vector identical to marker names?
  
  if(any(test))
    return(test)
  else
    return(0)
}

#' Checks the consistency of dataset (probability distribution)
#' @keywords internal
check_data_dist_sanity <- function(x){
  test <- logical(29L)
  names(test) <- 1:29
  
  # ploidy
  test[1] <- x$ploidy%%2 != 0 #is ploidy even?
  test[2] <- any(sapply(x$dosage.p1, function(y) max(y) > x$ploidy | min(y) < 0)) #are dosages in P higher than ploidy?
  test[3] <- any(sapply(x$dosage.p2, function(y) max(y) > x$ploidy | min(y) < 0)) #are dosages in Q higher than ploidy?
  test[4] <- ncol(x$geno) > x$ploidy + 3 #is the number of columns in the probability data frame correct? (ploidy + 3)
  test[5] <- max(x$geno.dose) > x$ploidy + 1 #is there any dose in offspring higher than ploidy?
  test[6] <- min(x$geno.dose) < 0 #is there any negative dose in offspring?
  
  # number of individuals
  test[7] <- x$n.ind < 0 #is the number of individuals greater than zero?
  test[8] <- length(x$ind.names) != x$n.ind #is the number of individual names equal to the number of individuals?
  test[9] <- length(unique(x$geno$ind)) != x$n.ind #is the number of individuals in the probability data frame equal to the number of individuals?
  test[10] <- ncol(x$geno.dose) != x$n.ind #is the number of columns in the dosage matrix equal to the number of individuals?
  
  # number of markers
  test[11] <- x$n.mrk < 0 #is the number of markers greater than zero?
  test[12] <- length(x$mrk.names) != x$n.mrk #is the number of marker names equal to the number of markers?
  test[13] <- length(x$dosage.p1) != x$n.mrk #is the number of marker dosages in P equal to the number of markers?
  test[14] <- length(x$dosage.p2) != x$n.mrk #is the number of marker dosages in Q equal to the number of markers?
  if(length(x$chrom) > 0)
    test[15] <- length(x$chrom) != x$n.mrk #is the number of sequences equal to the number of markers?
  if(length(x$genome.pos) > 0)
    test[16] <- length(x$genome.pos) != x$n.mrk #is the number of sequence positions equal to the number of markers?
  test[17] <- nrow(x$geno)/x$n.ind != x$n.mrk#is the number of rows in the probability data frame divided by n.ind equal to the number of markers?
  test[18] <- nrow(x$geno.dose) != x$n.mrk#is the number of rows in the dosage matrix equal to the number of markers?
  if(length(x$chisq.pval) > 0)
    test[19] <- length(x$chisq.pval) != x$n.mrk#is the number of chi-square tests equal to the number of markers?
  
  # individual names in the dosage and probability dataset
  test[20] <- !is.character(x$ind.names) # are individual's names characters
  test[21] <- !identical(x$geno$ind, rep(x$ind.names, each = x$n.mrk)) #are individual's names in the probability data frame properly 
  #arranged and consistent with the informed individual's names 
  test[22] <- !identical(colnames(x$geno.dose), x$ind.names)# are column names in dosage matrix identical to individual names?
  
  # marker names in the dosage and probability dataset
  test[23] <- !is.character(x$mrk.names)# are marker's names characters
  test[24] <- !identical(x$geno$mrk, rep(x$mrk.names, x$n.ind))#are marker names in the probability data frame properly 
  #arranged and consistent with the informed marker names 
  test[25] <- !identical(rownames(x$geno.dose), x$mrk.names)# are row names in dosage matrix identical to marker names?
  
  # dosage in both parents
  test[26] <- !is.integer(x$dosage.p1)# are dosages in P numeric
  test[27] <- !is.integer(x$dosage.p2)# are dosages in Q numeric 
  test[28] <- !identical(names(x$dosage.p1), x$mrk.names)# are names in P's dosage vector identical to marker names?
  test[29] <- !identical(names(x$dosage.p2), x$mrk.names)# are names in Q's dosage vector identical to marker names?
  
  if(any(test))
    return(test)
  else
    return(0)
}

#' Merge datasets
#'
#' This function merges two datasets of class \code{mappoly.data}. This can be useful
#' when individuals of a population were genotyped using two or more techniques
#' and have datasets in different files or formats. Please notice that the datasets 
#' should contain the same number of individuals and they must be represented identically 
#' in both datasets  (e.g. \code{Ind_1} in both datasets, not \code{Ind_1}
#' in one dataset and \code{ind_1} or \code{Ind.1} in the other).
#'
#' @param dat.1 the first dataset of class \code{mappoly.data} to be merged
#'
#' @param dat.2 the second dataset of class \code{mappoly.data} to be merged (default = NULL);
#' if \code{dat.2 = NULL}, the function returns \code{dat.1} only
#'
#' @return An object of class \code{mappoly.data} which contains all markers
#' from both datasets. It will be a list with the following components:
#'     \item{ploidy}{ploidy level}
#'     \item{n.ind}{number individuals}
#'     \item{n.mrk}{total number of markers}
#'     \item{ind.names}{the names of the individuals}
#'     \item{mrk.names}{the names of the markers}
#'     \item{dosage.p1}{a vector containing the dosage in
#'       parent P for all \code{n.mrk} markers}
#'     \item{dosage.p2}{a vector containing the dosage in
#'       parent Q for all \code{n.mrk} markers}
#'     \item{chrom}{a vector indicating which sequence each marker
#'       belongs. Zero indicates that the marker was not assigned to any
#'       sequence}
#'     \item{genome.pos}{Physical position of the markers into the
#'       sequence}
#'     \item{seq.ref}{if one or both datasets originated from read_vcf, it keeps
#' reference alleles from sequencing platform, otherwise is NULL}
#'     \item{seq.alt}{if one or both datasets originated from read_vcf, it keeps
#' alternative alleles from sequencing platform, otherwise is NULL}
#'     \item{all.mrk.depth}{if one or both datasets originated from read_vcf, it keeps
#' marker read depths from sequencing, otherwise is NULL}
#'     \item{prob.thres}{(unused field)}
#'     \item{geno.dose}{a matrix containing the dosage for each markers (rows) 
#'       for each individual (columns). Missing data are represented by 
#'       \code{ploidy_level + 1}}
#'     \item{geno}{if both datasets contain genotype distribution information,
#' the final object will contain 'geno'. This is set to NULL otherwise}
#'     \item{nphen}{(0)}
#'     \item{phen}{(NULL)}
#'     \item{chisq.pval}{a vector containing p-values related to the chi-squared 
#'     test of Mendelian segregation performed for all markers in both datasets}
#'     \item{kept}{if elim.redundant = TRUE when reading any dataset, holds all non-redundant markers}
#'     \item{elim.correspondence}{if elim.redundant = TRUE when reading any dataset,
#' holds all non-redundant markers and its equivalence to the redundant ones}
#' 
#' @author Gabriel Gesteira, \email{gdesiqu@ncsu.edu}
#' @examples
#' \donttest{
#' ## Loading a subset of SNPs from chromosomes 3 and 12 of sweetpotato dataset 
#' ## (SNPs anchored to Ipomoea trifida genome)
#' dat <- NULL
#' for(i in c(3, 12)){
#'   cat("Loading chromosome", i, "...\n")
#'     tempfl <- tempfile(pattern = paste0("ch", i), fileext = ".vcf.gz")
#'     x <- "https://github.com/mmollina/MAPpoly_vignettes/raw/master/data/sweet_sample_ch"
#'     address <- paste0(x, i, ".vcf.gz")
#'     download.file(url = address, destfile = tempfl)
#'     dattemp <- read_vcf(file = tempfl, parent.1 = "PARENT1", parent.2 = "PARENT2",
#'                         ploidy = 6, verbose = FALSE)
#'     dat <- merge_datasets(dat, dattemp)
#'   cat("\n")
#' }
#' dat
#' plot(dat)
#'}
#'
#' @references
#'     Mollinari, M., and Garcia, A.  A. F. (2019) Linkage
#'     analysis and haplotype phasing in experimental autopolyploid
#'     populations with high ploidy level using hidden Markov
#'     models, _G3: Genes, Genomes, Genetics_. 
#'     \doi{10.1534/g3.119.400378} 
#'
#' @export merge_datasets
#' @importFrom dplyr bind_rows arrange
merge_datasets = function(dat.1 = NULL, dat.2 = NULL){
  ## Check objects class
  if (is.null(dat.1)){
    if (is.null(dat.2)) return(dat.1)
    if (!inherits(dat.2, "mappoly.data")) {
      stop(deparse(substitute(dat.2)), " is not an object of class 'mappoly.data'")
    } else {
      return(dat.2)
    }
  } 
  if (!inherits(dat.1, "mappoly.data")) {
    stop(deparse(substitute(dat.1)), " is not an object of class 'mappoly.data'")
  }
  if (is.null(dat.2)) return(dat.1)
  if (!inherits(dat.2, "mappoly.data")) {
    stop(deparse(substitute(dat.2)), " is not an object of class 'mappoly.data'")
  }
  ## Check ploidy
  if (dat.1$ploidy != dat.2$ploidy){
    stop("The ploidy levels in the datasets do not match. Please check your files and try again.")
  }
  ## Check individuals
  if (dat.1$n.ind != dat.2$n.ind){
    stop("The datasets have different number of individuals. Please check your files and try again.")
  }
  ## Check individual names
  if (!all(dat.1$ind.names %in% dat.2$ind.names)){
    nmi.1 = dat.1$ind.names[!(dat.1$ind.names %in% dat.2$ind.names)]
    nmi.2 = dat.2$ind.names[!(dat.2$ind.names %in% dat.1$ind.names)]
    cat("Some individuals' names don't match:\n")
    cat("Dataset 1\tDataset 2\n")
    for (i in 1:length(nmi.1)) cat(nmi.1[i], "\t\t", nmi.2[i], "\n")
    stop("Your datasets have the same number of individuals, but some of them are not the same. Please check the list above, fix your files and try again.")
  }
  ## Checking (and fixing) individuals order in geno.dose
  if (!all(colnames(dat.1$geno.dose)  ==  colnames(dat.2$geno.dose))){
    dat.2$geno.dose = dat.2$geno.dose[,colnames(dat.1$geno.dose)]
  }
  ## Merging all items
  dat.1$geno.dose = rbind(dat.1$geno.dose, dat.2$geno.dose)
  dat.1$dosage.p1 = c(dat.1$dosage.p1, dat.2$dosage.p1)
  dat.1$dosage.p2 = c(dat.1$dosage.p2, dat.2$dosage.p2)
  ##dat.1$mrk.names = c(dat.1$mrk.names, dat.2$mrk.names) ## Fixing possible name incompatibilities  
  dat.1$mrk.names = rownames(dat.1$geno.dose) ## Fixing possible name incompatibilities  
  dat.1$n.mrk = dat.1$n.mrk + dat.2$n.mrk
  dat.1$chrom = c(dat.1$chrom, dat.2$chrom)
  dat.1$genome.pos = c(dat.1$genome.pos, dat.2$genome.pos)
  # VCF info
  ## seq.ref and seq.alt
  if (!is.null(dat.1$seq.ref) && !is.null(dat.2$seq.ref)){
    dat.1$seq.ref = c(dat.1$seq.ref,dat.2$seq.ref)
    dat.1$seq.alt = c(dat.1$seq.alt,dat.2$seq.alt)  
  }  else if (is.null(dat.1$seq.ref) && is.null(dat.2$seq.ref)){
    dat.1$seq.ref = dat.1$seq.alt = NULL
  }  else if (!is.null(dat.1$seq.ref) && is.null(dat.2$seq.ref)){
    dat.1$seq.ref = c(dat.1$seq.ref,rep(NA,length(dat.2$n.mrk)))
    dat.1$seq.alt = c(dat.1$seq.alt,rep(NA,length(dat.2$n.mrk)))
  } else if (is.null(dat.1$seq.ref) && !is.null(dat.2$seq.ref)){
    dat.1$seq.ref = c(rep(NA,length(dat.1$n.mrk)),dat.2$seq.ref)
    dat.1$seq.alt = c(rep(NA,length(dat.1$n.mrk)),dat.2$seq.alt)
  }
  ## all.mrk.depth
  if (!is.null(dat.1$all.mrk.depth) && !is.null(dat.2$all.mrk.depth)){
    dat.1$all.mrk.depth = c(dat.1$all.mrk.depth,dat.2$all.mrk.depth)
  }  else if (is.null(dat.1$all.mrk.depth) && is.null(dat.2$all.mrk.depth)){
    dat.1$all.mrk.depth = NULL
  }  else if (!is.null(dat.1$all.mrk.depth) && is.null(dat.2$all.mrk.depth)){
    dat.1$all.mrk.depth = c(dat.1$all.mrk.depth,rep(NA,length(dat.2$n.mrk)))
  } else if (is.null(dat.1$all.mrk.depth) && !is.null(dat.2$all.mrk.depth)){
    dat.1$all.mrk.depth = c(rep(NA,length(dat.1$n.mrk)),dat.2$all.mrk.depth)
  }
  # Chisq info
  if (!is.null(dat.1$chisq.pval) && !is.null(dat.2$chisq.pval)){
    dat.1$chisq.pval = c(dat.1$chisq.pval,dat.2$chisq.pval)
  }  else if (is.null(dat.1$chisq.pval) && is.null(dat.2$chisq.pval)){
    dat.1$chisq.pval = NULL
  }  else if (!is.null(dat.1$chisq.pval) && is.null(dat.2$chisq.pval)){
    dat.1$chisq.pval = c(dat.1$chisq.pval,rep(NA,length(dat.2$n.mrk)))
  } else if (is.null(dat.1$chisq.pval) && !is.null(dat.2$chisq.pval)){
    dat.1$chisq.pval = c(rep(NA,length(dat.1$n.mrk)),dat.2$chisq.pval)
  }
  # Redundant info
  ## Kept info
  if (!is.null(dat.1$kept) && !is.null(dat.2$kept)){
    dat.1$kept = c(dat.1$kept,dat.2$kept)
  }  else if (is.null(dat.1$kept) && is.null(dat.2$kept)){
    dat.1$kept = NULL
  }  else if (!is.null(dat.1$kept) && is.null(dat.2$kept)){
    dat.1$kept = c(dat.1$kept, dat.2$mrk.names)
  } else if (is.null(dat.1$kept) && !is.null(dat.2$kept)){
    dat.1$kept = c(dat.1$mrk.names,dat.2$kept)
  }  
  ## elim.correspondence info
  if (!is.null(dat.1$elim.correspondence) && !is.null(dat.2$elim.correspondence)){
    dat.1$elim.correspondence = rbind(dat.1$elim.correspondence,dat.2$elim.correspondence)
  }  else if (is.null(dat.1$elim.correspondence) && is.null(dat.2$elim.correspondence)){
    dat.1$elim.correspondence = NULL
  }  else if (!is.null(dat.1$elim.correspondence) && is.null(dat.2$elim.correspondence)){
    dat.1$elim.correspondence = dat.1$elim.correspondence
  } else if (is.null(dat.1$elim.correspondence) && !is.null(dat.2$elim.correspondence)){
    dat.1$elim.correspondence = dat.2$elim.correspondence
  }  
  ## geno dist info (just keep if both datasets contain this information)
  if (exists('geno', where = dat.1) && exists('geno', where = dat.2)){
    #dat.1$geno = rbind(dat.1$geno,dat.2$geno)
    ind <- NULL
    dat.1$geno <- dplyr::bind_rows(dat.1$geno, dat.2$geno) %>% arrange(ind)
  } else dat.1$geno = NULL
  if(!is.null(dat.1$chisq.pval) | !any(is.na(dat.1$chisq.pval)))
    dat.1$chisq.pval <- dat.1$chisq.pval[dat.1$mrk.names]
  ## Fixing possible issues with marker names
  names(dat.1$dosage.p1) = names(dat.1$dosage.p2) = names(dat.1$genome.pos) = names(dat.1$chrom) = names(dat.1$chisq.pval) = dat.1$mrk.names
  return(dat.1)
}

#' Summary maps
#'
#' This function generates a brief summary table of a list of \code{mappoly.map} objects
#' @param map.list a list of objects of class \code{mappoly.map}
#' @param verbose if \code{TRUE} (default), the current progress is shown; if
#'     \code{FALSE}, no output is produced
#' @return a data frame containing a brief summary of all maps contained in \code{map.list}
#' @examples
#' tetra.sum <- summary_maps(solcap.err.map)
#' tetra.sum
#' @author Gabriel Gesteira, \email{gdesiqu@ncsu.edu}
#' @export summary_maps
summary_maps = function(map.list, verbose = TRUE){
  ## Check data
  if (any(!sapply(map.list, inherits, "mappoly.map"))) 
    stop(deparse(substitute(map.list)), 
         " is not a list containing 'mappoly.map' objects.")
  md <- unlist(lapply(map.list, function(x) sum(get_tab_mrks(x)) - sum(get_tab_mrks(x)[rbind(c(2,2),c(x$info$ploidy,x$info$ploidy))]) - sum(get_tab_mrks(x)[rbind(c(1,2),c(2,1),c(x$info$ploidy,(x$info$ploidy+1)),c((x$info$ploidy+1),x$info$ploidy))])))
  if(map.list[[1]]$info$ploidy == 2)
    md[] <- 0
  results = data.frame("LG" = as.character(seq(1,length(map.list),1)),
                       "Genomic sequence" = as.character(unlist(lapply(map.list, function(x) paste(unique(x$info$chrom), collapse = "-")))),
                       "Map length (cM)" = unlist(lapply(map.list, function(x) round(sum(c(0, imf_h(x$maps[[1]]$seq.rf))), 2))),
                       "Markers/cM" = round(unlist(lapply(map.list, function(x) x$info$n.mrk/(round(sum(c(0, imf_h(x$maps[[1]]$seq.rf))), 2)))),2),
                       "Simplex" = unlist(lapply(map.list, function(x) sum(get_tab_mrks(x)[rbind(c(1,2),c(2,1),c(x$info$ploidy,(x$info$ploidy+1)),c((x$info$ploidy+1),x$info$ploidy))]))),
                       "Double-simplex" = unlist(lapply(map.list, function(x) sum(get_tab_mrks(x)[rbind(c(2,2),c(x$info$ploidy,x$info$ploidy))]))),
                       "Multiplex" = md,
                        "Total" = unlist(lapply(map.list, function(x) x$info$n.mrk)),
                       "Max gap" = unlist(lapply(map.list, function(x) round(imf_h(max(x$maps[[1]]$seq.rf)),2))),
                       check.names = FALSE, stringsAsFactors = F)
  results = rbind(results, c('Total', NA, sum(as.numeric(results$`Map length (cM)`)), round(mean(as.numeric(results$`Markers/cM`)),2), sum(as.numeric(results$Simplex)), sum(as.numeric(results$`Double-simplex`)), sum(as.numeric(results$Multiplex)), sum(as.numeric(results$Total)), round(mean(as.numeric(results$`Max gap`)),2)))
  if (verbose){
    all.mrks = unlist(lapply(map.list, function(x) return(x$info$mrk.names)))
    if (!any(get(map.list[[1]]$info$data.name, pos = 1)$elim.correspondence$elim %in% all.mrks))
      message("\nYour dataset contains removed (redundant) markers. Once finished the maps, remember to add them back with the function 'update_map'.\n")
  }
  return(results)
}

#' Get table of dosage combinations
#' 
#' Internal function
#' @param x an object of class \code{mappoly.map}
#' @author Gabriel Gesteira, \email{gdesiqu@ncsu.edu}
get_tab_mrks = function(x){
  tab = table(get(x$info$data.name, pos = 1)$dosage.p1[which(get(x$info$data.name, pos = 1)$mrk.names %in% x$info$mrk.names)], get(x$info$data.name, pos = 1)$dosage.p2[which(get(x$info$data.name, pos = 1)$mrk.names %in% x$info$mrk.names)])
  doses = as.character(seq(0,x$info$ploidy,1))
  # checking dimensions
  if (!all(doses %in% colnames(tab))){
    newmat = matrix(NA, nrow(tab), length(doses))
    newmat[,which(doses %in% colnames(tab))] = tab[,which(colnames(tab) %in% doses)]
    newmat[,which(!doses %in% colnames(tab))] = 0
    rownames(newmat) = rownames(tab)
    colnames(newmat) = doses
    tab = newmat
  }
  if (!all(doses %in% rownames(tab))){
    newmat = matrix(NA, length(doses), ncol(tab))
    newmat[which(doses %in% rownames(tab)),] = tab[which(rownames(tab) %in% doses),]
    newmat[which(!doses %in% rownames(tab)),] = 0
    colnames(newmat) = colnames(tab)
    rownames(newmat) = doses
    tab = newmat
  }
  return(tab)
}

#' Update map
#' 
#' This function takes an object of class \code{mappoly.map} and checks for
#' removed redundant markers in the original dataset. Once redundant markers
#' are found, they are re-added to the map in their respective equivalent positions
#' and another HMM round is performed.
#' 
#' @param input.maps a single map or a list of maps of class \code{mappoly.map}
#' @param verbose if TRUE (default), shows information about each update process
#' @return an updated map (or list of maps) of class \code{mappoly.map}, containing the original map(s) plus redundant markers
#' @author Gabriel Gesteira, \email{gdesiqu@ncsu.edu}
#' @examples
#' orig.map <- solcap.err.map
#' up.map <- lapply(solcap.err.map, update_map)
#' summary_maps(orig.map)
#' summary_maps(up.map)
#' @export update_map
#' 
update_map = function(input.maps, verbose = TRUE){
  ## Checking object
  if (inherits(input.maps, "mappoly.map"))
    input.maps = list(input.maps)
  if(!inherits(input.maps, "list"))
    stop(deparse(substitute(input.maps)), 
         " is not an object of class 'mappoly.map' neither a list containing 'mappoly.map' objects.")
  if (any(!sapply(input.maps, inherits, "mappoly.map"))) 
    stop(deparse(substitute(input.maps)), 
         " is not an object of class 'mappoly.map' neither a list containing 'mappoly.map' objects.")
  ## Checking the existence of redundant markers
  if (is.null(get(input.maps[[1]]$info$data.name, pos = 1)$elim.correspondence))
    stop('Your dataset does not contain redundant markers. Please check it and try again.')
  ## Creating list to handle results
  results = list()
  for (i in 1:length(input.maps)){
    if (verbose) cat("Updating map", i , "\n")    
    input.map = input.maps[[i]]
    ## Checking if redundant markers belong to the informed map
    map.kept.mrks = which(as.character(get(input.map$info$data.name, pos = 1)$elim.correspondence$kept) %in% input.map$info$mrk.names)
    if (is.null(map.kept.mrks)){
      if (verbose) cat("There is no redundant marker on map ",i,". Skipping it.\n")
      next
    }
    corresp = get(input.map$info$data.name, pos = 1)$elim.correspondence[which(as.character(get(input.map$info$data.name, pos = 1)$elim.correspondence$kept) %in% input.map$info$mrk.names),]
    ## Check if redundant markers were not already added to the map
    if (any(corresp$elim %in% input.map$info$mrk.names)) {
      if (verbose) cat("Some redundant markers were already added to the map ", i,". These markers will be skipped.\n")
      corresp = corresp[!(corresp$elim %in% input.map$info$mrk.names),]
    }
    ## Updating number of markers
    input.map$info$n.mrk = input.map$info$n.mrk + nrow(corresp)
    ## Adding markers to the sequence
    while (nrow(corresp) > 0){
      pos.kep = match(as.character(corresp$kept), get(input.map$info$data.name, pos = 1)$mrk.names)
      input.map$info$seq.num = append(input.map$info$seq.num, NA, after = which(input.map$info$seq.num  ==  pos.kep[1]))
      input.map$info$chisq.pval = append(input.map$info$chisq.pval, input.map$info$chisq.pval[which(input.map$info$seq.num  ==  pos.kep[1])], after = which(input.map$info$seq.num  ==  pos.kep[1]))
      input.map$info$seq.dose.p1 = append(input.map$info$seq.dose.p1, input.map$info$seq.dose.p1[which(input.map$info$seq.num  ==  pos.kep[1])], after = which(input.map$info$seq.num  ==  pos.kep[1]))
      input.map$info$seq.dose.p2 = append(input.map$info$seq.dose.p2, input.map$info$seq.dose.p2[which(input.map$info$seq.num  ==  pos.kep[1])], after = which(input.map$info$seq.num  ==  pos.kep[1]))
      input.map$info$chrom = as.numeric(append(input.map$info$chrom, as.character(corresp$chrom[1]), after = which(input.map$info$seq.num  ==  pos.kep[1])))
      input.map$info$genome.pos = as.numeric(append(input.map$info$genome.pos, as.character(corresp$genome.pos[1]), after = which(input.map$info$seq.num  ==  pos.kep[1])))
      input.map$info$mrk.names = append(input.map$info$mrk.names, as.character(corresp$elim[1]), after = which(input.map$info$mrk.names  ==  as.character(corresp$kept[1])))
      if (!is.null(input.map$info$seq.ref) && !is.null(input.map$info$seq.alt)){
        input.map$info$seq.ref = append(input.map$info$seq.ref, as.character(corresp$seq.ref[1]), after = which(input.map$info$seq.num  ==  pos.kep[1]))
        input.map$info$seq.alt = append(input.map$info$seq.alt, as.character(corresp$seq.alt[1]), after = which(input.map$info$seq.num  ==  pos.kep[1]))
      }
      for (j in 1:length(input.map$maps)){
        input.map$maps[[j]]$seq.rf = append(input.map$maps[[j]]$seq.rf, 0.000001, after = which(input.map$info$seq.num  ==  pos.kep[1]))
        input.map$maps[[j]]$seq.ph$P = append(input.map$maps[[j]]$seq.ph$P, input.map$maps[[j]]$seq.ph$P[paste0(pos.kep[1])], after = which(input.map$info$seq.num  ==  pos.kep[1]))
        input.map$maps[[j]]$seq.ph$Q = append(input.map$maps[[j]]$seq.ph$Q, input.map$maps[[j]]$seq.ph$Q[paste0(pos.kep[1])], after = which(input.map$info$seq.num  ==  pos.kep[1]))
      }
      corresp = corresp[-1,]
    }
    ## Renaming objects and updating map information
    names(input.map$info$chrom) = names(input.map$info$genome.pos) = names(input.map$info$chisq.pval) = names(input.map$info$seq.dose.p1) = names(input.map$info$seq.dose.p2) = names(input.map$info$seq.num) = input.map$info$mrk.names
    for (j in 1:length(input.map$maps)){
      names(input.map$maps[[j]]$seq.ph$P) = names(input.map$maps[[j]]$seq.ph$Q) = as.character(input.map$info$seq.num)
      input.map$maps[[j]]$seq.num = input.map$info$seq.num
    }
    if (!is.null(input.map$info$seq.ref) && !is.null(input.map$info$seq.alt))
      names(input.map$info$seq.ref) = names(input.map$info$seq.alt) = input.map$info$mrk.names
    results[[i]] = input.map
  }
  if (length(results)  ==  1) return(results[[1]])
  else return(results)
}

#' Random sampling of dataset
#' @param input.data an object  of class \code{mappoly.data}
#' @param n number of individuals or markers to be sampled
#' @param percentage if \code{n == NULL}, the percentage of individuals or markers to be sampled
#' @param type should sample individuals or markers?
#' @param selected.ind a vector containing the name of the individuals to select. Only has effect 
#' if \code{type = "individual"}, \code{n = NULL} and \code{percentage = NULL}
#' @param selected.mrk a vector containing the name of the markers to select. Only has effect 
#' if \code{type = "marker"}, \code{n = NULL} and \code{percentage = NULL}
#' @return an object  of class \code{mappoly.data}
#' @keywords internal
#' @export
#' @importFrom magrittr "%>%"
#' @importFrom dplyr filter
sample_data <- function(input.data, n = NULL, 
                        percentage = NULL, 
                        type = c("individual", "marker"),
                        selected.ind = NULL,
                        selected.mrk = NULL){
  if(!is.null(selected.mrk))
    type <- "marker"
  else type <- match.arg(type)
  if(type  ==  "individual"){
    if(!is.null(n)){
      selected.ind.id <- sort(sample(input.data$n.ind, n))
    } else if(!is.null(percentage)){
      selected.ind.id <- sample(input.data$n.ind, ceiling(input.data$n.ind * percentage/100))
    } else if(!is.null(selected.ind)){
      selected.ind.id <- match(selected.ind, input.data$ind.names)
    } else {
      stop("Inform 'n', 'percentage' or selected.ind.")
    }
    ind <- mrk <- NULL
    selected.ind <- input.data$ind.names[selected.ind.id]
    if(length(selected.ind.id) >= input.data$n.ind) return(input.data)
    if(is.prob.data(input.data))
      input.data$geno <-  input.data$geno %>%
      dplyr::filter(ind%in%selected.ind)
    input.data$geno.dose <- input.data$geno.dose[,selected.ind.id]
    input.data$ind.names <- input.data$ind.names[selected.ind.id]
    input.data$n.ind <- ncol(input.data$geno.dose)
    ##Computing chi-square p.values
    if(!is.null(input.data$chisq.pval)){
      ploidy <- input.data$ploidy
      Ds <- array(NA, dim = c(ploidy+1, ploidy+1, ploidy+1))
      for(i in 0:ploidy)
        for(j in 0:ploidy)
          Ds[i+1,j+1,] <- segreg_poly(ploidy = ploidy, dP = i, dQ = j)
      Dpop <- cbind(input.data$dosage.p1, input.data$dosage.p2)
      M <- t(apply(Dpop, 1, function(x) Ds[x[1]+1, x[2]+1,]))
      dimnames(M) <- list(input.data$mrk.names, c(0:ploidy))
      M <- cbind(M, input.data$geno.dose)
      input.data$chisq.pval <- apply(M, 1, mrk_chisq_test, ploidy = ploidy)
    }
    return(input.data)
  } 
  else if(type  ==  "marker"){
    if(!is.null(n)){
      selected.mrk.id <- sort(sample(input.data$n.mrk, n))
    } else if(!is.null(percentage)){
      selected.mrk.id <- sort(sample(input.data$n.mrk, ceiling(input.data$n.mrk * percentage/100)))
    } else if(!is.null(selected.mrk)){
      selected.mrk.id <- match(selected.mrk, input.data$mrk.names)
    } else{
      stop("Inform 'n', 'percentage' or selected.mrk.")
    }
    selected.mrk <- input.data$mrk.names[selected.mrk.id]
    if(length(selected.mrk.id) >= input.data$n.mrk) return(input.data)
    if(is.prob.data(input.data))
      input.data$geno <-  input.data$geno %>%
      dplyr::filter(mrk%in%selected.mrk)
    input.data$geno.dose <- input.data$geno.dose[selected.mrk.id, , drop = FALSE]
    input.data$n.mrk <- nrow(input.data$geno.dose)
    input.data$mrk.names <- input.data$mrk.names[selected.mrk.id]
    input.data$dosage.p1 <- input.data$dosage.p1[selected.mrk.id]
    input.data$dosage.p2 <- input.data$dosage.p2[selected.mrk.id]
    input.data$chrom <- input.data$chrom[selected.mrk.id]
    input.data$genome.pos <- input.data$genome.pos[selected.mrk.id]
    input.data$seq.ref <- input.data$seq.ref[selected.mrk.id]
    input.data$seq.alt <- input.data$seq.alt[selected.mrk.id]
    input.data$all.mrk.depth <- input.data$all.mrk.depth[selected.mrk.id]
    input.data$kept <- intersect(input.data$mrk.names, input.data$kept)
    input.data$elim.correspondence <- input.data$elim.correspondence[input.data$elim.correspondence$kept%in%input.data$mrk.names, , drop = FALSE]
    input.data$chisq.pval <- input.data$chisq.pval[names(input.data$chisq.pval)%in%input.data$mrk.names]
    if(!is.null(input.data$chisq.pval)) 
      input.data$chisq.pval <- input.data$chisq.pval[names(input.data$chisq.pval)%in%input.data$mrk.names]
    return(input.data)  
  }
  else stop("Inform type")
}


#' Get weighted ordinary least squared map give a sequence and rf matrix
#' @keywords internal
#' @importFrom zoo na.approx
get_ols_map <- function(input.seq, input.mat, weight = TRUE){
  id <- input.seq$seq.mrk.names
  rf <- imf_h(get_rf_from_mat(input.mat$rec.mat[id,id]))
  pos <- c(0, cumsum(ifelse(is.na(rf), 0, rf)) + rf*0)
  names(pos) <- input.seq$seq.mrk.names
  pos <- rev(rev(pos)[!cumprod(is.na(rev(pos)))])
  id <- names(pos)
  z <- zoo::na.approx(pos)  
  names(z) <- id
  y <- as.numeric((imf_h(as.dist(input.mat$rec.mat[id,id]))))
  w <- as.numeric((imf_h(as.dist(input.mat$lod.mat[id,id]))))
  v <- t(combn(id,2))
  x <- numeric(nrow(v))
  names(x) <- names(y) <- apply(v, 1, paste0, collapse = "-")
  for(i in 1:nrow(v))
    x[i] <- z[v[i,2]]-z[v[i,1]]
  if(weight)
    model <- lm(y ~ x-1, weights = w)
  else
    model <- lm(y ~ x-1)
  new <- data.frame(x = z)
  u <- predict(model, new)
  d <- cumsum(imf_h(c(0, mf_h(diff(u)))))
  names(d) <- id
  d
}

#' Get Dosage Type in a Sequence
#'
#' Analyzes a genomic sequence object to categorize markers based on their dosage type.
#' The function calculates the dosage type by comparing the dosage of two parental sequences
#' (p1 and p2) against the ploidy level. It categorizes markers into simplex for parent 1 (simplex.p),
#' simplex for parent 2 (simplex.q), double simplex (ds), and multiplex based on the calculated dosages.
#'
#' @param input.seq An object of class \code{"mappoly.sequence"}:
#' 
#' @return A list with four components categorizing marker names into:
#' \describe{
#'   \item{simplex.p}{Markers with a simplex dosage from parent 1.}
#'   \item{simplex.q}{Markers with a simplex dosage from parent 2.}
#'   \item{double.simplex}{Markers with a double simplex dosage.}
#'   \item{multiplex}{Markers not fitting into the above categories, indicating a multiplex dosage.}
#' }
#' @keywords internal
#' @export
get_dosage_type <- function(input.seq){
  p <- abs(abs(input.seq$seq.dose.p1 - input.seq$ploidy/2) - input.seq$ploidy/2)
  q <- abs(abs(input.seq$seq.dose.p2 - input.seq$ploidy/2) - input.seq$ploidy/2)
  s.p <- p  ==  1 & q  ==  0
  s.q <- p  ==  0 & q  ==  1
  ds <- p  ==  1 & q  ==  1
  list(simplex.p = input.seq$seq.mrk.names[s.p],
       simplex.q = input.seq$seq.mrk.names[s.q], 
       double.simplex = input.seq$seq.mrk.names[ds],
       multiplex = input.seq$seq.mrk.names[!(s.p | s.q | ds)])
}

#' Aggregate matrix cells (lower the resolution by a factor)
#' @keywords internal
aggregate_matrix <- function(M, fact){
  id <- seq(1,ncol(M), by = fact)
  id <- cbind(id, c(id[-1]-1, ncol(M)))
  R <- matrix(NA, nrow(id), nrow(id))
  for(i in 1:(nrow(id)-1)){
    for(j in (i+1):nrow(id)){
      R[j,i] <-  R[i,j] <- mean(M[id[i,1]:id[i,2], id[j,1]:id[j,2]], na.rm = TRUE)
    }
  }
  R
}

#' Get states and emission in one informative parent
#' @keywords internal
get_states_and_emission_single_parent <- function(ploidy, ph, global.err, D, dose.notinf.P){
  n.mrk <- nrow(D)
  n.ind <- ncol(D)
  A <- matrix(0, nrow = choose(ploidy, ploidy/2), ncol = length(ph))
  for(i in 1:length(ph)){
    id1 <- numeric(ploidy)
    id1[ph[[i]]] <- 1
    A[,i] <- apply(combn(id1, ploidy/2), 2, sum) + dose.notinf.P[i]/2
  }
  if(round(global.err, 4) == 0.0){
    e <- h <- vector("list", n.mrk)
    for(j in 1:n.mrk){
      e.temp <- h.temp <- vector("list", n.ind)
      for(i in 1:n.ind){
        h.temp[[i]] <- which(A[,j] == D[j,i]) - 1
          if(length(h.temp[[i]]) == 0)
            h.temp[[i]] <- 1:nrow(A) - 1
        #e.temp[[i]] <- rep(1, length(h.temp[[i]]))/length(h.temp[[i]])
         e.temp[[i]] <- rep(1, length(h.temp[[i]]))
      }
      e[[j]] <- e.temp
      h[[j]] <- h.temp
    }
  } else if(round(global.err, 4) > 0.0){
    e <- h <- vector("list", n.mrk)
    for(j in 1:n.mrk){
      e.temp <- h.temp <- vector("list", n.ind)
      for(i in 1:n.ind){
        h.temp[[i]] <- 1:nrow(A) - 1
        s <- which(A[,j] == D[j,i])
        if(length(s) == 0)
          e.temp[[i]] <- rep(1/nrow(A), nrow(A))
        else{
          e.temp0 <- numeric(nrow(A))
          e.temp0[-s] <- global.err/(nrow(A) - length(s))
          e.temp0[s] <- (1- global.err)/length(s)
          e.temp[[i]] <- e.temp0
        }
      }
      e[[j]] <- e.temp
      h[[j]] <- h.temp
    }
  }
  list(states = h, emission = e)    
}


#' Conversion: vector to matrix
#' @keywords internal
v_2_m <- function(x, n){
  y <- base::matrix(NA, n, n) 
  y[base::lower.tri(y)] <- as.numeric(x)
  y[base::upper.tri(y)] <- t(y)[base::upper.tri(y)]
  y
}

#' Compare a list of maps
#' 
#' Compare lengths, density, maximum gaps and log likelihoods in a list of maps.
#' In order to make the maps comparable, the function uses the intersection of 
#' markers among maps. 
#' 
#' @param ... a list of objects of class \code{mappoly.map}
#' 
#' @return A data frame where the lines correspond to the maps in the order provided in input list list
#' 
#' @export
compare_maps <- function(...){
  l <- list(...)
  id <- lapply(l, function(x) x$info$mrk.names)
  id <- Reduce(intersect, id)
  if(length(id) < 2)
    stop("At least two markers must be present in all maps")
  map.out <- vector("list", length(l))
  for(i in 1:length(l)){
    id.temp <- which(l[[i]]$info$mrk.names%in%id)
    map.temp <- filter_map_at_hmm_thres(get_submap(l[[i]], mrk.pos = id.temp, 
                                                   reestimate.rf = FALSE, verbose = FALSE), 10e-5)
    map.out[[i]] <-loglike_hmm(map.temp)
  }
  a <- summary_maps(map.out, verbose = FALSE)
  a <- cbind(a[-nrow(a), -c(1,5,6,7,8)], sapply(map.out, function(x) x$maps[[1]]$loglike))
  a <- cbind(a, rep("-", nrow(a)))
  a[which.max(a[,5]),6] <- "*"
  colnames(a)[c(5,6)] <- c("loglike", "max_likelog")
  return(as.data.frame(a))
}

#' Detects which parent is informative
#'
#' @param x an object of class \code{mappoly.sequence} or \code{mappoly.map}
#' 
#' @export
detect_info_par<-function(x){
  ## checking for correct object
  if (inherits(x, "mappoly.sequence")) {
   if(all(x$seq.dose.p2 == 0 | x$seq.dose.p2 == x$ploidy))
     return("p1")
  else if(all(x$seq.dose.p1 == 0 | x$seq.dose.p1 == x$ploidy))
    return("p2")
  else
    return("both")
  } else if (inherits(x, "mappoly.map")) {
    if(all(x$info$seq.dose.p2 == 0 | x$info$seq.dose.p2 == x$info$ploidy))
      return("p1")
    else if(all(x$info$seq.dose.p1 == 0 | x$info$seq.dose.p1 == x$info$ploidy))
      return("p2")
    else
      return("both")
  } else{
    stop(deparse(substitute(x)), " is not an object of class 'mappoly.map' or 'mappoly.sequence'")    
  }
}

# Skeleton to test CPP functions
# @export
# test_CPP<-function(m, rf)
#   .Call("rec_number", as.integer(m), as.numeric(rf), PACKAGE = "mappoly")
#   


.mappoly_data_skeleton<-function(ploidy =NULL,
                                 n.ind =NULL,
                                 n.mrk =NULL,
                                 ind.names =NULL,
                                 mrk.names =NULL,
                                 dosage.p1 =NULL,
                                 dosage.p2 =NULL,
                                 chrom =NULL,
                                 genome.pos =NULL,
                                 seq.ref =NULL,
                                 seq.alt =NULL,
                                 all.mrk.depth =NULL,
                                 prob.thres =NULL,
                                 geno.dose =NULL,
                                 nphen =NULL,
                                 phen =NULL,
                                 kept =NULL,
                                 chisq.pval =NULL,
                                 elim.correspondence =NULL)
  structure(list(ploidy = ploidy,
                 n.ind = n.ind,
                 n.mrk = n.mrk,
                 ind.names = ind.names,
                 mrk.names = mrk.names,
                 dosage.p1 = dosage.p1,
                 dosage.p2 = dosage.p2,
                 chrom = chrom,
                 genome.pos = genome.pos,
                 seq.ref = seq.ref,
                 seq.alt = seq.alt,
                 all.mrk.depth = all.mrk.depth,
                 prob.thres = prob.thres,
                 geno.dose = geno.dose,
                 nphen = nphen,
                 phen = phen,
                 kept = kept,
                 chisq.pval = chisq.pval,
                 elim.correspondence = elim.correspondence),
            class = "mappoly.data")


.mappoly_map_skeleton<-function(ploidy=NULL,
                                n.mrk=NULL,
                                seq.num=NULL,
                                mrk.names=NULL,
                                seq.dose.p1=NULL,
                                seq.dose.p2=NULL,
                                chrom=NULL,
                                genome.pos=NULL,
                                seq.ref=NULL,
                                seq.alt=NULL,
                                chisq.pval=NULL,
                                data.name=NULL,
                                ph.thresh =NULL,
                                seq.rf =NULL,
                                seq.ph = list(P =NULL, Q =NULL),
                                loglike = 0)
  structure(list(info = list(ploidy=ploidy,
                             n.mrk=n.mrk,
                             seq.num=seq.num,
                             mrk.names=mrk.names,
                             seq.dose.p1=seq.dose.p1,
                             seq.dose.p2=seq.dose.p2,
                             chrom=chrom,
                             genome.pos=genome.pos,
                             seq.ref=seq.ref,
                             seq.alt=seq.alt,
                             chisq.pval=chisq.pval,
                             data.name=data.name,
                             ph.thresh = ph.thresh),
                 maps = list(list(seq.num = seq.num,
                                  seq.rf = seq.rf,
                                  seq.ph = seq.ph,
                                  loglike = loglike))),
            class = "mappoly.map")
  
  
#' Split map into sub maps given a gap threshold
#' @keywords internal
split_mappoly <- function(input.map,
                          gap.threshold = 5, 
                          size.rem.cluster = 1,
                          phase.config = "best",
                          tol.final = 10e-4,
                          verbose = TRUE){
  if (!inherits(input.map, "mappoly.map")) {
    stop(deparse(substitute(input.map)), " is not an object of class 'mappoly.map'")
  }
  ## choosing the linkage phase configuration
  LOD.conf <- get_LOD(input.map, sorted = FALSE)
  if(phase.config  ==  "best") {
    i.lpc <- which.min(LOD.conf)
  } else if (phase.config > length(LOD.conf)) {
    stop("invalid linkage phase configuration")
  } else i.lpc <- phase.config
  id <- which(imf_h(input.map$maps[[i.lpc]]$seq.rf) > gap.threshold)
  if(length(id) == 0){
    if(verbose) cat("no submaps found\n")
    return(input.map)
  } 
  id <- cbind(c(1, id+1), c(id, input.map$info$n.mrk))
  temp.map <- input.map
  ## Selecting map segments larger then the specified threshold
  segments <- id[apply(id, 1, diff) > size.rem.cluster - 1, , drop = FALSE]
  if(length(segments) == 0) stop("all markers were eliminated\n")
  ## Dividing map in sub-maps
  temp.maps <- vector("list", nrow(segments))
  if (verbose) {
    ns <- nrow(segments)
    if(ns == 1){
      cat("one submap found ...\n")
      map <- get_submap(input.map, c(segments[1, 1]:segments[1, 2]), tol.final = tol.final, verbose = FALSE)
      return(filter_map_at_hmm_thres(map, 10e-4))
    } 
    else cat(ns, "submaps found ...\n")
  }
  for(i in 1:length(temp.maps)){
    temp.id <- c(segments[i, 1]:segments[i, 2])
    if(length(temp.id) > 1)
      temp.maps[[i]] <- get_submap(input.map, temp.id, reestimate.rf = FALSE, verbose = FALSE)
    else
      temp.maps[[i]] <- input.map$info$mrk.names[temp.id]    
  }
  temp.maps
}
