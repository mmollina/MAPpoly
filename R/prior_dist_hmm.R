#' Estimate genetic map using as input the probability distribution of
#' genotypes (wrapper function to C++)
#'
#' @param void internal function to be documented
#' @keywords internal
#'
poly_hmm_est <-
  function(m,
           n.mrk,
           n.ind,
           p,
           dp,
           q,
           dq,
           g,
           rf,
           verbose = TRUE,
           tol = 0.001) {
      ## Checking capabilities
      if (verbose && !capabilities("long.double")){
        cat("This function uses high precision calculations, but your system's architecture doesn't support long double allocation ('capabilities('long.double') = FALSE'). Running in low precision mode.\n")
      }
    res <-
      .Call(
        "poly_hmm_est_CPP",
        as.numeric(m),
        as.numeric(n.mrk),
        as.numeric(n.ind),
        as.numeric(p),
        as.numeric(dp),
        as.numeric(q),
        as.numeric(dq),
        as.double(g),
        as.double(rf),
        as.numeric(rep(0, choose(m, m / 2) ^ 2)),
        as.double(0),
        as.numeric(verbose),
        as.double(tol),
        PACKAGE = "mappoly"
      )
    names(res) <- c("loglike", "rf")
    return(res)
  }
