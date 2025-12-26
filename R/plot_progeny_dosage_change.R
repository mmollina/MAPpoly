#' Plot progeny dosage changes after HMM-based correction with a global error model
#'
#' Computes genotype probabilities under a global genotyping error rate, derives
#' homolog probabilities, compares the most likely HMM-implied dosages to the
#' original dosage matrix, and plots which entries were unchanged, imputed
#' (originally missing), or changed.
#'
#' Missing dosages in \code{dat$geno.dose} are assumed to be encoded as
#' \code{ploidy + 1}. The function prints a heatmap-style \code{ggplot2} tile
#' plot and, optionally, returns a corrected dosage matrix with marker metadata.
#'
#' @param map_list A non-empty list of objects of class \code{mappoly.map}.
#' @param error Numeric scalar. Global genotyping error rate passed to
#'   \code{\link[mappoly]{calc_genoprob_error}}.
#' @param verbose Logical. If \code{TRUE} (default), progress is printed by
#'   underlying routines; if \code{FALSE}, suppresses progress messages where
#'   supported.
#' @param output_corrected Logical. If \code{FALSE} (default), the function only
#'   prints the plot and returns \code{invisible(NULL)}. If \code{TRUE}, returns
#'   a matrix containing marker metadata columns followed by the corrected
#'   dosages (individuals in columns).
#'
#' @return If \code{output_corrected = FALSE}, returns \code{invisible(NULL)} and
#'   prints a \code{ggplot}. If \code{output_corrected = TRUE}, returns a matrix
#'   with columns \code{P1}, \code{P2}, \code{sequence}, \code{sequence_position},
#'   followed by one column per individual containing corrected dosages.
#'
#' @details
#' The \code{mappoly.data} object is retrieved by name from
#' \code{map_list[[1]]$info$data.name}. The function:
#' \enumerate{
#'   \item Runs \code{calc_genoprob_error()} for each map in \code{map_list}.
#'   \item Computes homolog probabilities with \code{calc_homologprob()}.
#'   \item For each marker and individual, selects the \code{ploidy} most likely
#'         homologs and converts them to a dosage implied by the phased map.
#'   \item Compares implied dosages to \code{dat$geno.dose} and summarizes the
#'         fraction of entries that were imputed (original missing) vs changed.
#'   \item Produces a tile plot showing unchanged/imputed/changed cells.
#' }
#'
#' @examples
#' \dontrun{
#' x <- get_submap(solcap.err.map[[1]], 1:5, reestimate.rf = FALSE)
#'
#' corrected_matrix <- plot_progeny_dosage_change(list(x), error = 0.05,
#'                                                output_corrected = TRUE)
#' }
#'
#' @author Jeekin Lau, with optimization by Cristiane Taniguti
#'
#' @seealso \code{\link[mappoly]{calc_genoprob_error}},
#'   \code{\link[mappoly]{calc_homologprob}}
#'
#' @import ggplot2
#' @importFrom reshape2 melt
#' @export
plot_progeny_dosage_change <- function(map_list,
                                       error,
                                       verbose = TRUE,
                                       output_corrected = FALSE) {
  Var1 <- Var2 <- value <- NULL
  
  map <- map_list
  if (!is.list(map) || length(map) == 0L) {
    stop("'map_list' must be a non-empty list of 'mappoly.map' objects.")
  }
  
  # Retrieve mappoly.data by name (search calling env first, then global)
  data_name <- map[[1]]$info$data.name
  data_env <- parent.frame()
  if (!exists(data_name, envir = data_env, inherits = TRUE)) data_env <- .GlobalEnv
  if (!exists(data_name, envir = data_env, inherits = TRUE)) stop("mappoly.data object not here")
  dat <- get(data_name, envir = data_env, inherits = TRUE)
  
  ploidy <- as.integer(map[[1]]$info$ploidy)
  if (is.na(ploidy) || ploidy <= 0L) stop("Invalid ploidy in map object.")
  if (2L * ploidy > length(letters)) stop("Ploidy too large for letter-based labels.")
  
  print("calculating genoprob error")
  
  genoprob <- vector("list", length(map))
  for (i in seq_along(map)) {
    genoprob[[i]] <- calc_genoprob_error(input.map = map[[i]],
                                         error = error,
                                         verbose = verbose)
  }
  
  print("calculating homologprob")
  homoprobs <- calc_homologprob(genoprob, verbose = verbose)
  
  print("comparing to orginal")
  
  # P/Q lists for each marker
  P <- unlist(lapply(map, function(x) x$maps[[1]]$seq.ph$P), recursive = FALSE)
  Q <- unlist(lapply(map, function(x) x$maps[[1]]$seq.ph$Q), recursive = FALSE)
  
  n_mrk <- length(P)
  P_matrix <- matrix(0, nrow = n_mrk, ncol = ploidy)
  Q_matrix <- matrix(0, nrow = n_mrk, ncol = ploidy)
  
  for (i in seq_len(n_mrk)) {
    pi <- as.integer(P[[i]])
    qi <- as.integer(Q[[i]])
    pi <- pi[!is.na(pi)]
    qi <- qi[!is.na(qi)]
    if (length(pi)) P_matrix[i, pi] <- 1
    if (length(qi)) Q_matrix[i, qi] <- 1
  }
  
  mrks_mapped <- unlist(lapply(map, function(x) x$info$mrk.names), use.names = FALSE)
  
  PQ_matrix <- cbind(P_matrix, Q_matrix)
  rownames(PQ_matrix) <- mrks_mapped
  colnames(PQ_matrix) <- letters[seq_len(2L * ploidy)]
  
  homoprob <- homoprobs$homoprob
  if (!is.data.frame(homoprob)) homoprob <- as.data.frame(homoprob)
  if (ncol(homoprob) < 4L) stop("Unexpected 'homoprob' format.")
  names(homoprob)[1:4] <- c("marker", "homolog", "individual", "prob")
  
  inds <- sort(unique(homoprob$individual))
  mrks <- unique(homoprob$marker)
  
  # Build marker x individual matrix of concatenated top-ploidy homolog labels
  by_ind <- split(homoprob, homoprob$individual)
  test_list <- lapply(by_ind, function(df_i) {
    by_mrk <- split(df_i, df_i$marker)
    vapply(mrks, function(m) {
      df_m <- by_mrk[[m]]
      if (is.null(df_m)) return(NA_character_)
      df_m <- df_m[order(df_m$prob, decreasing = TRUE), , drop = FALSE]
      take <- seq_len(min(ploidy, nrow(df_m)))
      paste0(as.character(df_m$homolog[take]), collapse = "")
    }, FUN.VALUE = character(1))
  })
  
  test <- do.call(cbind, test_list)
  colnames(test) <- names(test_list)
  rownames(test) <- mrks
  test <- test[, inds, drop = FALSE]
  
  # Ensure markers exist in PQ_matrix
  if (!all(mrks %in% rownames(PQ_matrix))) {
    missing <- setdiff(mrks, rownames(PQ_matrix))
    stop("Markers missing from phase matrices: ",
         paste(missing[seq_len(min(10, length(missing)))], collapse = ", "),
         if (length(missing) > 10) " ...")
  }
  
  # Convert homolog strings -> dosage by summing membership in P/Q
  finished <- matrix(NA_real_,
                     nrow = length(mrks),
                     ncol = length(inds),
                     dimnames = list(mrks, inds))
  
  for (a in seq_along(mrks)) {
    m <- mrks[a]
    s <- test[m, ]
    
    # n_individuals x ploidy matrix with the k-th character of each string
    letter_mat <- vapply(seq_len(ploidy),
                         function(k) substring(s, k, k),
                         FUN.VALUE = character(length(s)))
    letter_mat[letter_mat == ""] <- NA_character_
    
    vals <- vapply(seq_len(ploidy),
                   function(k) PQ_matrix[m, letter_mat[, k]],
                   FUN.VALUE = numeric(length(s)))
    
    finished[a, ] <- rowSums(vals, na.rm = TRUE)
  }
  
  # Pull original dosages in exactly the same order
  if (!all(mrks %in% rownames(dat$geno.dose))) {
    missing <- setdiff(mrks, rownames(dat$geno.dose))
    stop("Markers missing from 'geno.dose': ",
         paste(missing[seq_len(min(10, length(missing)))], collapse = ", "),
         if (length(missing) > 10) " ...")
  }
  if (!all(inds %in% colnames(dat$geno.dose))) {
    missing <- setdiff(inds, colnames(dat$geno.dose))
    stop("Individuals missing from 'geno.dose': ",
         paste(missing[seq_len(min(10, length(missing)))], collapse = ", "),
         if (length(missing) > 10) " ...")
  }
  
  original_geno <- as.matrix(dat$geno.dose[mrks, inds, drop = FALSE])
  
  # Correct comparisons (the old code used !original_geno==finished, which is not "!=")
  same <- (original_geno == finished) & !is.na(original_geno) & !is.na(finished)
  imputed <- (!same) & (original_geno == (ploidy + 1L))
  changed <- (!same) & (original_geno != (ploidy + 1L))
  
  percent_imputed <- sum(imputed, na.rm = TRUE) / length(original_geno) * 100
  percent_changed <- sum(changed, na.rm = TRUE) / length(original_geno) * 100
  
  status <- matrix("unchanged",
                   nrow = nrow(original_geno),
                   ncol = ncol(original_geno),
                   dimnames = dimnames(original_geno))
  status[imputed] <- "imputed"
  status[changed] <- "changed"
  
  empty_matrix_melt <- reshape2::melt(status)
  
  cols <- c(changed = "red", imputed = "chartreuse", unchanged = "black")
  
  plot1 <- ggplot2::ggplot(
    empty_matrix_melt,
    ggplot2::aes(Var1, Var2,
                 fill = factor(value, levels = c("changed", "imputed", "unchanged")))
  ) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_manual(values = cols, drop = FALSE) +
    ggplot2::xlab("Markers") +
    ggplot2::ylab("Individuals") +
    ggplot2::ggtitle(paste0(
      "changed = ", round(percent_changed, digits = 3), "% ",
      "imputed = ", round(percent_imputed, digits = 3), "%"
    )) +
    ggplot2::guides(fill = ggplot2::guide_legend(title = "Dosage change"))
  
  print("done")
  print(plot1)
  
  if (isTRUE(output_corrected)) {
    mrk_names <- rownames(finished)
    
    # Keep order aligned with 'finished'
    mrk_index <- match(mrk_names, dat$mrk.names)
    if (anyNA(mrk_index)) stop("Some markers were not found in 'dat$mrk.names'.")
    
    P1 <- dat$dosage.p1[mrk_index]
    P2 <- dat$dosage.p2[mrk_index]
    sequence <- dat$chrom[mrk_index]
    sequence_position <- dat$genome.pos[mrk_index]
    
    hmm_imputed <- cbind(P1, P2, sequence, sequence_position, finished)
    return(hmm_imputed)
  }
  
  invisible(NULL)
}
