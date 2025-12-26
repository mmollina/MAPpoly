.safe_require_namespace <- function(pkg) {
  tryCatch(requireNamespace(pkg, quietly = TRUE), error = function(e) FALSE)
}

.can_use_plot3D <- function() {
  # On macOS without XQuartz, loading tcltk fails, which breaks misc3d/plot3D.
  # capabilities("X11") is usually FALSE in that situation.
  if (!isTRUE(capabilities("X11"))) return(FALSE)
  .safe_require_namespace("plot3D")
}

.plot_prob_2d_fallback <- function(prob_mat, prob_thres, pal, main = "Genotype probability (2D fallback)") {
  prob_mat <- as.matrix(prob_mat)
  storage.mode(prob_mat) <- "double"
  
  if (nrow(prob_mat) == 0L || ncol(prob_mat) == 0L) {
    graphics::plot.new()
    graphics::title(main = main)
    return(invisible(NULL))
  }
  
  ploidy <- ncol(prob_mat) - 1L
  nind <- nrow(prob_mat)
  
  # Flatten to match dose-major coloring (dose 0..ploidy, each across individuals)
  p <- as.vector(t(prob_mat))
  p[is.na(p)] <- 0
  p <- pmax(0, pmin(1, p))
  
  # Base colors per dose, repeated for each individual
  base_cols <- rep(pal, each = nind)
  base_cols <- base_cols[seq_along(p)]
  
  # Build RGBA colors without adjustcolor()
  rgbm <- grDevices::col2rgb(base_cols)
  alpha <- as.integer(round(p * 255))
  cols <- grDevices::rgb(rgbm[1, ], rgbm[2, ], rgbm[3, ], alpha = alpha, maxColorValue = 255)
  
  # Gray out below threshold
  cols[p < prob_thres] <- "#404040"
  
  x <- rep(seq_len(nind), times = ploidy + 1L)
  y <- rep(0:ploidy, each = nind)
  
  graphics::plot(
    NA,
    xlim = c(1, nind),
    ylim = c(-0.5, ploidy + 0.5),
    xlab = "Offspring (sorted)",
    ylab = "Dose",
    yaxt = "n"
  )
  graphics::axis(2, at = 0:ploidy, labels = 0:ploidy, las = 2)
  graphics::points(x, y, pch = 15, cex = 0.45, col = cols)
  graphics::mtext(main, side = 3, line = 0.5, cex = 0.85)
  
  invisible(NULL)
}


#' Plot marker information
#'
#' Plots summary statistics for a given marker.
#'
#' @param input.data an object of class \code{mappoly.data}
#' @param mrk marker name or position in the dataset
#'
#' @examples
#' plot_mrk_info(tetra.solcap.geno.dist, 2680)
#' plot_mrk_info(tetra.solcap.geno.dist, "solcap_snp_c2_23828")
#'
#' @export
#' @importFrom graphics barplot layout mtext legend par plot text axis points
#' @importFrom stats chisq.test
plot_mrk_info <- function(input.data, mrk) {
  if (!inherits(input.data, "mappoly.data")) {
    stop(deparse(substitute(input.data)), " is not an object of class 'mappoly.data'")
  }
  
  # Resolve marker name
  if (is.numeric(mrk)) {
    mrk <- mrk[1]
    if (mrk > input.data$n.mrk) stop(mrk, " exceeds the number of markers in the dataset")
    mrk <- input.data$mrk.names[mrk]
  }
  
  idx <- match(mrk, input.data$mrk.names)
  if (is.na(idx)) {
    stop(deparse(substitute(mrk)), " is not present in ", deparse(substitute(input.data)), " dataset")
  }
  
  dp <- input.data$dosage.p1[idx]
  dq <- input.data$dosage.p2[idx]
  
  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldpar), add = TRUE)
  on.exit(graphics::layout(1), add = TRUE)
  
  graphics::par(mar = c(2, 2, 5, 2), bg = "gray98")
  
  has_prob <- is.prob.data(input.data)
  
  suppressWarnings({
    if (!has_prob) {
      graphics::layout(matrix(c(1, 2), ncol = 2), widths = c(1, 2))
      
      xdose <- input.data$geno.dose[mrk, ]
      xdose[xdose == input.data$ploidy + 1] <- NA
      xtab <- table(as.numeric(xdose), useNA = "always")
      names(xtab) <- c(names(xtab)[-length(xtab)], "NA")
      
      graphics::plot(0:100, type = "n", axes = FALSE, xlab = "", ylab = "")
      graphics::mtext(side = 3, text = mrk, adj = 0, cex = 1.2, font = 3)
      graphics::text(x = 0, y = 90, labels = paste0("marker #: ", idx), adj = 0)
      graphics::text(x = 0, y = 80, labels = paste0("Dose in P1: ", dp), adj = 0)
      graphics::text(x = 0, y = 70, labels = paste0("Dose in P2: ", dq), adj = 0)
      graphics::text(x = 0, y = 60, labels = paste0("Missing: ", round(100 * tail(xtab, 1) / sum(xtab), 1), "%"), adj = 0)
      
      seg.exp <- segreg_poly(ploidy = input.data$ploidy, dP = dp, dQ = dq)
      names(seg.exp) <- 0:input.data$ploidy
      seg.exp <- seg.exp[seg.exp != 0]
      seg.obs <- seg.exp
      seg.obs[names(xtab)[-length(xtab)]] <- xtab[-length(xtab)]
      
      pval <- tryCatch(
        stats::chisq.test(x = seg.obs, p = seg.exp[names(seg.obs)])$p.value,
        error = function(e) NA_real_
      )
      graphics::text(
        x = 0, y = 50,
        labels = if (is.na(pval)) "p-value: NA" else paste0("p-value: ", formatC(pval, format = "e", digits = 2)),
        adj = 0
      )
      
      graphics::text(x = 0, y = 40, labels = paste0("sequence: ", input.data$chrom[idx]), adj = 0)
      graphics::text(x = 0, y = 30, labels = paste0("seq. position: ", input.data$genome.pos[idx]), adj = 0)
      
      graphics::barplot(
        xtab,
        col = c(mp_pallet2(input.data$ploidy + 1)[1:(length(xtab) - 1)], "#404040")
      )
      
    } else {
      graphics::layout(matrix(c(1, 2, 3, 3), ncol = 2, nrow = 2), widths = c(1, 2))
      
      xdose <- input.data$geno.dose[mrk, ]
      xdose[xdose == input.data$ploidy + 1] <- NA
      xtab <- table(xdose, useNA = "always")
      names(xtab) <- c(names(xtab)[-length(xtab)], "NA")
      
      graphics::plot(0:100, type = "n", axes = FALSE, xlab = "", ylab = "")
      graphics::mtext(side = 3, text = mrk, adj = 0, cex = 1.2, font = 3)
      graphics::text(x = 0, y = 90, labels = paste0("marker #: ", idx), adj = 0)
      graphics::text(x = 0, y = 80, labels = paste0("Dose in P1: ", dp), adj = 0)
      graphics::text(x = 0, y = 70, labels = paste0("Dose in P2: ", dq), adj = 0)
      graphics::text(x = 0, y = 60, labels = paste0("Missing: ", round(100 * tail(xtab, 1) / sum(xtab), 1), "%"), adj = 0)
      
      seg.exp <- segreg_poly(ploidy = input.data$ploidy, dP = dp, dQ = dq)
      names(seg.exp) <- 0:input.data$ploidy
      seg.exp <- seg.exp[seg.exp != 0]
      seg.obs <- seg.exp
      seg.obs[names(xtab)[-length(xtab)]] <- xtab[-length(xtab)]
      
      pval <- tryCatch(
        stats::chisq.test(x = seg.obs, p = seg.exp[names(seg.obs)])$p.value,
        error = function(e) NA_real_
      )
      graphics::text(
        x = 0, y = 50,
        labels = if (is.na(pval)) "p-value: NA" else paste0("p-value: ", formatC(pval, format = "e", digits = 2)),
        adj = 0
      )
      
      graphics::text(x = 0, y = 40, labels = paste0("sequence: ", input.data$chrom[idx]), adj = 0)
      graphics::text(x = 0, y = 30, labels = paste0("seq. position: ", input.data$genome.pos[idx]), adj = 0)
      graphics::text(x = 0, y = 20, labels = paste0("prob. threshold: ", input.data$prob.thres), adj = 0)
      
      pal <- mp_pallet2(input.data$ploidy + 1)
      names(pal) <- 0:input.data$ploidy
      op1 <- graphics::par(mar = c(5, 3, 0, 2), cex = 0.7)
      on.exit(graphics::par(op1), add = TRUE)
      
      graphics::barplot(xtab, col = c(na.omit(pal[names(xtab)]), "#404040"))
      
      # Build probability matrix for plotting
      df <- input.data$geno[input.data$geno$mrk == mrk, c(as.character(0:input.data$ploidy), "ind")]
      rownames(df) <- df$ind
      
      prob_mat <- as.matrix(df[, as.character(0:input.data$ploidy), drop = FALSE])
      
      # Sort offspring so high-probability individuals are grouped (similar intent as original)
      ord <- rev(do.call(order, as.data.frame(prob_mat)))
      prob_mat <- prob_mat[ord, , drop = FALSE]
      
      # Optional filter: remove rows that match exact expected segregation vector
      expv <- segreg_poly(input.data$ploidy, dp, dq)
      keep <- apply(prob_mat, 1, function(v) !all(round(v, 5) == round(expv, 5)))
      prob_mat <- prob_mat[keep, , drop = FALSE]
      
      # Plot 3D if possible, otherwise fall back to 2D
      dose_pal <- mp_pallet2(input.data$ploidy + 1)
      
      if (.can_use_plot3D()) {
        nind <- nrow(prob_mat)
        w <- expand.grid(seq_len(nind), 0:input.data$ploidy)
        xx <- w[, 1]
        yy <- w[, 2]
        pp <- as.vector(prob_mat[, as.character(0:input.data$ploidy), drop = FALSE])
        
        cols <- rep(dose_pal, each = nind)
        cols[is.na(pp) | pp < input.data$prob.thres] <- "#404040"
        
        # plot3D call guarded by tryCatch in case something still goes wrong
        tryCatch(
          plot3D::scatter3D(
            x = xx, y = yy, z = pp,
            theta = 30, phi = 30,
            bty = "g",
            type = "h", lwd = 0.3,
            ticktype = "detailed",
            pch = 19, cex = 0.5,
            colvar = NULL,
            col = cols,
            xlab = "Offspring (sorted)",
            ylab = "Dose",
            zlab = "Genotype probability",
            cex.axis = 0.7, cex.lab = 0.7,
            clab = ""
          ),
          error = function(e) {
            warning("3D probability plot disabled (plot3D/tcltk issue). Using 2D fallback.\n  ", conditionMessage(e))
            .plot_prob_2d_fallback(
              prob_mat = prob_mat,
              prob_thres = input.data$prob.thres,
              pal = dose_pal
            )
          }
        )
      } else {
        .plot_prob_2d_fallback(
          prob_mat = prob_mat,
          prob_thres = input.data$prob.thres,
          pal = dose_pal
        )
      }
    }
  })
  
  invisible(NULL)
}
