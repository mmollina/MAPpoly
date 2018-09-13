#' Compute the conditional QTL probabilities using genotype distribution
#'
#' @param void interfunction to be documented
#' @keywords internal
#' @export
#'
calc_genoprob_dist<-function(input.map, dat.dist,  phase.config = "best", verbose = TRUE)
{
  if (!inherits(input.map, "mappoly.map")) {
    stop(deparse(substitute(input.map)), " is not an object of class 'mappoly.map'")
  }
  ## choosing the linkage phase configuration
  LOD.conf <- get_LOD(input.map, sorted = FALSE)
  if(phase.config == "best") {
    i.lpc <- which.min(LOD.conf)
  } else if (phase.config > length(LOD.conf)) {
    stop("invalid linkage phase configuration")
  } else i.lpc <- phase.config
  mrk<-NULL
  original.map.mrk<-get(input.map$info$data.name, pos=1)$mrk.names[input.map$maps[[i.lpc]]$seq.num]
  dat.dist.pos<-match(original.map.mrk, dat.dist$mrk.names)
  which.is.na<-which(is.na(dat.dist.pos))
  if(length(which.is.na) > 0)
    stop("Markers", original.map.mrk[which.is.na], "are not present in the 'dat.dist' object")
  temp.map<-input.map
  temp.map$info$data.name<-as.character(sys.call())[3]
  temp.map$maps[[i.lpc]]$seq.num<-dat.dist.pos
  names(temp.map$maps[[i.lpc]]$seq.ph$P)<-names(temp.map$maps[[i.lpc]]$seq.ph$Q)<-dat.dist.pos
  if(!all(sort(get(temp.map$info$data.name, pos = 1)$ind.names) %in% sort(get(input.map$info$data.name, pos = 1)$ind.names)))
    stop("The individuals in the new data set are not contained in the original data set")
  geno<-subset(get(temp.map$info$data.name, pos = 1)$geno, mrk%in%original.map.mrk)
  geno.new<-NULL
  for(i in unique(geno$ind))
    geno.new<-rbind(geno.new, geno[geno[,"ind"] == i, ][match(original.map.mrk, geno[,"mrk"]),])
  g <- as.double(t(geno.new[, -c(1:2)]))
  m = as.numeric(temp.map$info$m)
  n.mrk = as.numeric(temp.map$info$n.mrk)
  n.ind = get(temp.map$info$data.name, pos = 1)$n.ind
  p = as.numeric(unlist(temp.map$maps[[1]]$seq.ph$P))
  dp = as.numeric(cumsum(c(0, sapply(temp.map$maps[[1]]$seq.ph$P, function(x) sum(length(x))))))
  q = as.numeric(unlist(temp.map$maps[[1]]$seq.ph$Q))
  dq = as.numeric(cumsum(c(0, sapply(temp.map$maps[[1]]$seq.ph$Q, function(x) sum(length(x))))))
  rf = temp.map$maps[[1]]$seq.rf
  indnames<-get(input.map$info$data.name, pos=1)$ind.names
  res.temp <-
    .Call(
      "calc_genoprob_prior",
      as.numeric(m),
      as.numeric(n.mrk),
      as.numeric(n.ind),
      as.numeric(p),
      as.numeric(dp),
      as.numeric(q),
      as.numeric(dq),
      as.double(g),
      as.double(rf),
      as.numeric(rep(0, choose(m, m/2)^2 * n.mrk * n.ind)),
      as.double(0),
      as.numeric(verbose),
      PACKAGE = "mappoly"
    )
  dim(res.temp[[1]])<-c(choose(m,m/2)^2,n.mrk,n.ind)
  dimnames(res.temp[[1]])<-list(kronecker(apply(combn(letters[1:m],m/2),2, paste, collapse=""),
                                          apply(combn(letters[(m+1):(2*m)],m/2),2, paste, collapse=""), paste, sep=":"),
                                original.map.mrk, indnames)
  structure(list(probs = res.temp[[1]], map = create_map(input.map)), class="mappoly.genoprob")
}
