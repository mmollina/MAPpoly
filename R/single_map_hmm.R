#' Multipoint analysis using Hidden Markov Models (single phase)
#' @keywords internal
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu}
est_rf_hmm_single_phase <- function(input.seq,
                            input.ph.single,
                            rf.temp = NULL,
                            tol,
                            verbose = FALSE,
                            ret.map.no.rf.estimation = FALSE,
                            high.prec = TRUE,
                            max.rf.to.break.EM = 0.5)
{
  input_classes <- c("mappoly.sequence")
  if (!inherits(input.seq, input_classes[1])) {
    stop(deparse(substitute(input.seq)), " is not an object of class 'mappoly.sequence'")
  }
  if(length(input.seq$seq.num)  ==  1)
    stop("Input sequence contains only one marker.", call. = FALSE)
  if (verbose && !capabilities("long.double") && high.prec){
    cat("You've requested high precision calculations ('high.prec = TRUE'), but your system's architecture doesn't support long double allocation ('capabilities('long.double') = FALSE'). Running in low precision mode.\n")
  }
  if(is.null(rf.temp))
    rf.temp <- rep(0.001, length(input.seq$seq.num)-1)
  if(!ret.map.no.rf.estimation)
  {
    D <- get(input.seq$data.name, pos = 1)$geno.dose[input.seq$seq.num,]
    dp <- get(input.seq$data.name)$dosage.p1[input.seq$seq.num]
    dq <- get(input.seq$data.name)$dosage.p2[input.seq$seq.num]
    for (j in 1:nrow(D))
      D[j, D[j, ]  ==  input.seq$ploidy + 1] <- dp[j] + dq[j] + 1 + as.numeric(dp[j] == 0 || dq[j] == 0)

    if(high.prec)
    {
      res.temp <- .Call("est_map_hmm_highprec",
                      input.seq$ploidy,
                      t(D),
                      lapply(input.ph.single$P, function(x) x-1),
                      lapply(input.ph.single$Q, function(x) x-1),
                      rf.temp,
                      verbose = verbose,
                      max.rf.to.break.EM,
                      tol,
                      PACKAGE = "mappoly")
    } else{
      res.temp <- .Call("est_map_hmm",
                      input.seq$ploidy,
                      t(D),
                      lapply(input.ph.single$P, function(x) x-1),
                      lapply(input.ph.single$Q, function(x) x-1),
                      rf.temp,
                      verbose = verbose,
                      max.rf.to.break.EM,
                      tol,
                      PACKAGE = "mappoly")
    }
  } else res.temp <- list(1, rf.temp)
  map <- list(seq.num = input.seq$seq.num,
            seq.rf = res.temp[[2]],
            seq.ph = input.ph.single,
            loglike = res.temp[[1]])
  map
}
