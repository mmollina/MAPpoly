#' Generates the states the model should visit including global error rate
#'
#' @param void internal function to be documented
#' @examples
#'
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu}
#' @export
#'
states_to_visit <- function(input.data){
  ####States to visit####
  pd <- input.data$pedigree
  pd$id <- apply(pd[,1:2], 1, paste0, collapse = "x")
  upd <- unique(pd)
  unique.cross <- apply(upd[,1:2], 1, paste0, collapse = "x")
  unique.cross.id <- seq_along(unique.cross)
  names(unique.cross.id) <- unique.cross
  pd$id <- unique.cross.id[pd$id]
  upd <- cbind(upd, unique.cross.id)
  mrk.names <- rownames(input.data$offspring[[1]])
  ind.names <- rownames(input.data$pedigree)
  h <- vector("list", length(mrk.names))
  names(h) <- mrk.names
  for(i in names(h)){
    h[[i]] <- vector("list", length(ind.names))
    names(h[[i]]) <- ind.names
  }
  emit <- h
  for(k in unique.cross.id){
    ngam1 <- choose(upd$pl1[k], upd$pl1[k]/2)
    ngam2 <- choose(upd$pl2[k], upd$pl2[k]/2)
    S <- as.matrix(expand.grid(0:(ngam2-1), 0:(ngam1-1))[,2:1])
    P1 <- upd[k,"Par1"]
    P2 <- upd[k,"Par2"]
    pl1 <- upd$pl1[k]
    pl2 <- upd$pl2[k]
    cur.ind.names <- rownames(pd)[pd$id == k]
    for(j in mrk.names){
      a1 <- input.data$phases[[P1]][j,]
      a2 <- input.data$phases[[P2]][j,]
      if(any(is.na(a1)) | any(is.na(a2))){
        for(i in cur.ind.names){
          h[[j]][[i]] <- S
          #emit[[j]][[i]] <- matrix(rep(1/nrow(h[[j]][[i]]), nrow(h[[j]][[i]])), ncol = 1)
          emit[[j]][[i]] <- matrix(rep(1, nrow(h[[j]][[i]])), ncol = 1)
        }
      } else {
        A1 <- combn(a1, pl1/2) ## Ordered vector ---> phased
        A2 <- combn(a2, pl2/2) ## Ordered vector ---> phased
        w <- kronecker(apply(A1, 2, paste0, collapse = "_"),
                       apply(A2, 2, paste0, collapse = "_"), paste, sep = "_")
        u <- str_split_fixed(w,"_", (pl1 + pl2)/2)
        storage.mode(u) <- "integer"
        v <- apply(u,1,sort)
        for(i in cur.ind.names){
          id <- which(apply(v, 2, function(x) all(x == sort(input.data$offspring[[i]][j,]))))
          h[[j]][[i]] <- S[id, , drop = FALSE]
          #emit[[j]][[i]] <- matrix(rep(1/nrow(h[[j]][[i]]), nrow(h[[j]][[i]])), ncol = 1)
          emit[[j]][[i]] <- matrix(rep(1, nrow(h[[j]][[i]])), ncol = 1)
        }
      }
    }
  }
  ploidy <-input.data$pedigree[, c("pl1", "pl2")]
  list( n.mrk = length(h), n.ind = length(h[[1]]), states = h,
        emit = emit, ploidy = ploidy)
}

#' Generates the states the model should visit for mappoly legacy
#'
#' @param void internal function to be documented
#' @examples
#' require(mappoly)
#' n.mrk <- 2
#' h.temp <- sim_homologous(ploidy = 4,
#'                          n.mrk = n.mrk,
#'                          max.d = 2,
#'                          max.ph = 0,
#'                          seed = 1)
#'
#' for(i in 1:n.mrk){
#'   h.temp$hom.allele.p[[i]] <- 1:2
#'   h.temp$hom.allele.q[[i]] <- 1
#'   h.temp$p[i] <- 2
#'   h.temp$q[i] <- 1
#' }
#' dat <- poly_cross_simulate(ploidy = 4,
#'                            rf.vec = .01,
#'                            n.mrk = n.mrk,
#'                            n.ind = 10,
#'                            h.temp,
#'                            seed = 8532,
#'                            draw = TRUE)
#' sim.map<-cumsum(c(0,rep(imf_h(.01), (n.mrk - 1))))
#' plot(dat)
#' s <- make_seq_mappoly(dat, "all")
#' tpt <- est_pairwise_rf(s)
#' tpt$pairwise$`1-2`
#' map <- est_rf_hmm_sequential(input.seq = s, twopt = tpt)
#' plot(map)
#' mp <- round(cumsum(mappoly::imf_h(c(0, map$maps[[1]]$seq.rf))),2)
#' ph <- list(ph_list_to_matrix(h.temp$hom.allele.p, dat$ploidy),
#'            ph_list_to_matrix(h.temp$hom.allele.q, dat$ploidy))
#' states.hmm <- states_to_visit_mp1(dat, ph, is.log = TRUE)
#' x1 <- est_map_R(states.hmm,tol = 1e-3, verbose = FALSE)
#' mp2 <- round(cumsum(mappoly::imf_h(c(0, x1[[2]]))), 2)
#' Y1 <- rbind(sim.map,mp,mp2)
#' rownames(Y1) <- c("simulation", "mappoly", "mappoly2")
#' Y1
#'
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu}
#' @export
#'
states_to_visit_mp1 <- function(dat, ph, is.log = TRUE){
  ngam <- choose(dat$ploidy, dat$ploidy/2)
  S <- data.frame(expand.grid(0:(ngam-1),
                              0:(ngam-1))[,2:1],
                  as.numeric(1/(ngam * ngam)))
  Y <- vector("list", length(dat$geno.dose))
  cte<-1
  for(i in 1:nrow(dat$geno.dose))
  {
    a1<-apply(combn(ph[[1]][i,], ncol(ph[[1]])/2), 2, sum)
    a2<-apply(combn(ph[[2]][i,], ncol(ph[[2]])/2), 2, sum)
    st<-apply(expand.grid(a2, a1)[2:1], 1, sum)
    for(j in 1:ncol(dat$geno.dose))
    {
      Y[[cte]] <- data.frame(ind = j,
                             mrk = i,
                             pl.1 = dat$ploidy,
                             pl.2 = dat$ploidy,
                             pl.id = 1,
                             st.p1 = S[st == dat$geno.dose[i,j],1],
                             st.p2 = S[st == dat$geno.dose[i,j],2],
                             st.all = which(st == dat$geno.dose[i,j])-1,
                             emit = S[st == dat$geno.dose[i,j],3],
                             row.names = NULL)
      if(is.log) Y[[cte]]$emit <- log(Y[[cte]]$emit)
      cte <- cte + 1
    }
  }
  Y <- Y %>% bind_rows %>% arrange(ind)
  Y$ind <- as_factor(Y$ind)
  Y$mrk <- as_factor(Y$mrk)
  pls <- matrix(c(dat$ploidy, dat$ploidy), nrow = 1, dimnames = list(paste(dat$ploidy, dat$ploidy, sep = "x"), NULL))
  list(hmm.info = Y, err = 0.0, is.log = TRUE, ploidy.cross.id = pls)
}





