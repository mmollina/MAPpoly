#' Reestimate map using genotype distribution
#'
#' @param void interfunction to be documented
#' @keywords internal
#' @export
#'
reest_map_geno_dist<-function(input.map, dat.dist,  phase.config = "best", verbose = TRUE, tol = 10e-3)
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
  temp.map$info$data.name<-deparse(substitute(dat.dist))
  #temp.map$info$data.name<-as.character(sys.call())[3]
  temp.map$maps[[i.lpc]]$seq.num<-dat.dist.pos
  names(temp.map$maps[[i.lpc]]$seq.ph$P)<-names(temp.map$maps[[i.lpc]]$seq.ph$Q)<-dat.dist.pos
  if(!all(sort(get(temp.map$info$data.name, pos = 1)$ind.names) %in% sort(get(input.map$info$data.name, pos = 1)$ind.names)))
    stop("The individuals in the new data set are not contained in the original data set")
  geno<-subset(get(temp.map$info$data.name, pos = 1)$geno, mrk%in%original.map.mrk)
  geno.new<-NULL
  for(i in unique(geno$ind))
    geno.new<-rbind(geno.new, geno[geno[,"ind"] == i, ][match(original.map.mrk, geno[,"mrk"]),])
  g <- as.double(t(geno.new[, -c(1:2)]))
  map.res<-poly_hmm_est(m = as.numeric(temp.map$info$m),
                        n.mrk = as.numeric(temp.map$info$n.mrk),
                        n.ind = dat.dist$n.ind,
                        p = as.numeric(unlist(temp.map$maps[[1]]$seq.ph$P)),
                        dp = as.numeric(cumsum(c(0, sapply(temp.map$maps[[1]]$seq.ph$P, function(x) sum(length(x)))))),
                        q = as.numeric(unlist(temp.map$maps[[1]]$seq.ph$Q)),
                        dq = as.numeric(cumsum(c(0, sapply(temp.map$maps[[1]]$seq.ph$Q, function(x) sum(length(x)))))),
                        g = g,
                        rf = temp.map$maps[[1]]$seq.rf,
                        verbose = verbose,
                        tol = tol)
  temp.map$maps[[1]]$seq.rf<-map.res$rf
  temp.map$maps[[1]]$loglike<-map.res$loglike
  return(list(original.map = input.map, reestimated.map = temp.map))
}
