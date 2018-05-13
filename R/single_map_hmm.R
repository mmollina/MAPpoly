#' Multipoint analysis using Hidden Markov Models (single phase)
#'
#' @param void interfunction to be documented
#' @keywords internal
#' @export est_rf_hmm_single
est_rf_hmm_single<-function(input.seq,
                            input.ph.single,
                            rf.temp = NULL,
                            tol,
                            verbose = FALSE,
                            ret.map.no.rf.estimation = FALSE,
                            high.prec = TRUE)
{
  if(is.null(rf.temp))
    rf.temp<-rep(0.001, length(input.seq$seq.num)-1)
  if(!ret.map.no.rf.estimation)
  {
    D<-get(input.seq$data.name, pos=1)$geno.dose[input.seq$seq.num,]
    dp <- get(input.seq$data.name)$dosage.p[input.seq$seq.num]
    dq <- get(input.seq$data.name)$dosage.q[input.seq$seq.num]
    for (j in 1:nrow(D))
      D[j, D[j, ] == input.seq$m + 1] <- dp[j] + dq[j] + 1 + as.numeric(dp[j]==0 || dq[j]==0)

    if(high.prec)
    {
      res.temp<-.Call("est_map_hmm_highprec",
                      input.seq$m,
                      t(D),
                      lapply(input.ph.single$P, function(x) x-1),
                      lapply(input.ph.single$Q, function(x) x-1),
                      rf.temp,
                      verbose=verbose,
                      0.5,
                      tol,
                      PACKAGE = "mappoly")
    } else{
      res.temp<-.Call("est_map_hmm",
                      input.seq$m,
                      t(D),
                      lapply(input.ph.single$P, function(x) x-1),
                      lapply(input.ph.single$Q, function(x) x-1),
                      rf.temp,
                      verbose=verbose,
                      0.5,
                      tol,
                      PACKAGE = "mappoly")
    }
  } else res.temp<-list(1, rf.temp)
  map<-list(seq.num=input.seq$seq.num,
            seq.rf=res.temp[[2]],
            seq.ph=input.ph.single,
            loglike=res.temp[[1]])
  map
}
