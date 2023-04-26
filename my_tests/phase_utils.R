perm_rank <- function(mat.hom.k,mat.hom.k1,m){
  #### Stationary homologs ####
  hom.k <- apply(mat.hom.k, 2, paste, collapse = "")
  dimnames(mat.hom.k) = list(rownames(mat.hom.k), hom.k)
  #### Permutation homologs ####
  hom.k1 <- apply(mat.hom.k1, 2, paste, collapse = "")
  dimnames(mat.hom.k1) = list(rownames(mat.hom.k1), hom.k1)
  #### Unique permutations ####
  i.k1 <- perm_tot(seq_along(hom.k1))
  id.k1 <- !duplicated(perm_tot(hom.k1))
  i.k1 <- i.k1[id.k1,]
  hom.k1 <- perm_tot(hom.k1)
  hom.k1 <- hom.k1[id.k1,]
  rownames(i.k1) <- rownames(hom.k1) <- 1:nrow(hom.k1)
  #### Unique homologs ####
  a <- unique(t(apply(apply(hom.k1, 1, function(x) paste(hom.k, x, sep = "-")), 2, sort, decreasing = TRUE)))
  hom.k1 <- hom.k1[rownames(a),]
  i.k1 <- i.k1[rownames(a),]
  #### Merging each permuted homologs with stationary
  #x <- numeric(nrow(hom.k1))
  x <- matrix(NA, nrow(hom.k1), 1 + ncol(mat.hom.k)/2)
  colnames(x) <- c(0:(ncol(mat.hom.k)/2))
  for(i in 1:nrow(hom.k1)){
    #### Shared homologs per markers between blocks
    sh <- tcrossprod(mat.hom.k[, drop = FALSE], 
                     mat.hom.k1[,hom.k1[i,], drop = FALSE])
    
    sh.obs <- m$ShP[rownames(sh), colnames(sh)]
    y <- table(abs(sh - sh.obs), useNA = "no")
    x[i,names(y)] <- y
    #u <- expand.grid(match(rownames(sh), dat$mrk.names),
    #                 match(colnames(sh), dat$mrk.names))
    #id <- apply(u, 1, function(x) paste(sort(x), collapse = "-"))
    #z <- tpt$pairwise[id]
    #w <- as.numeric(sh)
    #for(j in 1:length(z))
    #  w[j] <- prob_connections(v = z[[j]], n = w[j])
    #x[i] <- mean(w)
  }
  id <- do.call(order, c(decreasing = TRUE, data.frame(x)))
  list(homolog.permutations = i.k1[id,], matches = x[id,])
}
expected_max_connections <- function(v){
  x <- exp(-v[2,]/log10(exp(1)))
  x <- x/sum(x)
  w <- as.numeric(str_split_fixed(string = names(x), pattern = "-", 2)[,1])
  return(as.numeric(x[which.max(w)] * w[which.max(w)]))
}
prob_connections <- function(v, n){
  x <- exp(-v[2,]/log10(exp(1)))
  x <- x/sum(x)
  w <- as.numeric(str_split_fixed(string = names(x), pattern = "-", 2)[,1])
  return(x[match(n, w)])
}
matrix_to_map <- function(ph.mat = NULL, cur.ph.mat, 
                          ploidy, seq.temp, m, dat, 
                          na.perc = 0.2, top.eval = 10){
  if(is.null(ph.mat))
    return(cur.ph.mat)
  perm.ph <- perm_rank(mat.hom.k = ph.mat, 
                       mat.hom.k1 = cur.ph.mat, 
                       m = m)
  v <- is.na(perm.ph$matches)
  if(sum(v)/length(v) < na.perc | nrow(v) < top.eval){
    which.test <- 1:nrow(v)          
  } else {
    which.test <- 1:top.eval      
  } 
  ct<-1
  if(length(which.test) > 1){
    x <- numeric(length(which.test))
    for(l in which.test){
      PH1 <- rbind(ph.mat, cur.ph.mat[,perm.ph$homolog.permutations[l,], drop = FALSE])
      PH1 <- PH1[intersect(seq.temp$seq.mrk.names, rownames(PH1)),]
      s1 <- make_seq_mappoly(dat, rownames(PH1), data.name = seq.temp$data.name)
      PH2 <- PH1 <- ph_matrix_to_list(PH1)
      PH2[] <- 0
      ph <- list(PH1, PH2)
      map1<-est_rf_hmm_single_phase_single_parent(s1, ph, tol = 10e-2)
      x[ct] <- map1$maps[[1]]$loglike
      ct <- 1 + ct
    }
    plot(x, type = "b", col = 2)
    ph.mat <- rbind(ph.mat, cur.ph.mat[,perm.ph$homolog.permutations[which.test[which.max(x)],], drop = FALSE])
  } else {
    ph.mat <- rbind(ph.mat, cur.ph.mat[,perm.ph$homolog.permutations[which.test,], drop = FALSE])
  }
  ph.mat1 <- ph.mat[intersect(seq.temp$seq.mrk.names, rownames(ph.mat)),]
  s1 <- make_seq_mappoly(dat, rownames(ph.mat1), data.name = seq.temp$data.name)
  ph.mat2 <- ph.mat1 <- ph_matrix_to_list(ph.mat1)
  ph.mat2[] <- 0
  ph <- list(ph.mat1, ph.mat2)
  map1<-est_rf_hmm_single_phase_single_parent(s1, ph, tol = 10e-4)
  plot(map1, mrk.names = T)
  return(ph.mat[s1$seq.mrk.names,])
}
list_to_map <- function(L, ploidy, seq.temp, m, 
                        dat, na.perc = 0.2, 
                        top.eval = 10){
  ph.mat <- NULL
  for(i in 1:length(L)){
    cat("\nSuperMarker", i, "\n")
    for(j in 1:length(L[[i]])){
      cat(".")
      cur.ph.mat <- matrix(0, 
                           length(L[[i]][[j]]), 
                           ploidy, 
                           dimnames = list(L[[i]][[j]], NULL))
      z <- unique(seq.temp$seq.dose.p1[match(rownames(cur.ph.mat), seq.temp$seq.mrk.names)])
      cur.ph.mat[,1:z] <- 1 
      if(is.null(ph.mat)){
        ph.mat <- cur.ph.mat
        next()
      }
      perm.ph <- perm_rank(mat.hom.k = ph.mat, 
                           mat.hom.k1 = cur.ph.mat, 
                           m = m)
      v <- is.na(perm.ph$matches)
      if(sum(v)/length(v) < na.perc | nrow(v) < top.eval){
        which.test <- 1:nrow(v)          
      } else {
        which.test <- 1:top.eval      
      } 
      ct<-1
      if(length(which.test) > 1){
        x <- numeric(length(which.test))
        for(l in which.test){
          PH1 <- rbind(ph.mat, cur.ph.mat[,perm.ph$homolog.permutations[l,], drop = FALSE])
          PH1 <- PH1[intersect(seq.temp$seq.mrk.names, rownames(PH1)),]
          s1 <- make_seq_mappoly(dat, rownames(PH1), data.name = seq.temp$data.name)
          PH2 <- PH1 <- ph_matrix_to_list(PH1)
          PH2[] <- 0
          ph <- list(PH1, PH2)
          map1<-est_rf_hmm_single_phase_single_parent(s1, ph, tol = 10e-2)
          x[ct] <- map1$maps[[1]]$loglike
          ct <- 1 + ct
        }
        #plot(x, type = "b", col = 2)
        ph.mat <- rbind(ph.mat, cur.ph.mat[,perm.ph$homolog.permutations[which.test[which.max(x)],], drop = FALSE])
      } else {
        ph.mat <- rbind(ph.mat, cur.ph.mat[,perm.ph$homolog.permutations[which.test,], drop = FALSE])
      }
      ph.mat1 <- ph.mat[intersect(seq.temp$seq.mrk.names, rownames(ph.mat)),]
      s1 <- make_seq_mappoly(dat, rownames(ph.mat1), data.name = seq.temp$data.name)
      ph.mat2 <- ph.mat1 <- ph_matrix_to_list(ph.mat1)
      ph.mat2[] <- 0
      ph <- list(ph.mat1, ph.mat2)
      map1<-est_rf_hmm_single_phase_single_parent(s1, ph, tol = 10e-4)
      plot(map1, mrk.names = T)
    }
  }
  return(ph.mat[s1$seq.mrk.names,])
}
make_mat_prob_connections <- function(tpt.ll, dat){
  M <- matrix(NA, tpt.ll$n.mrk, tpt.ll$n.mrk, 
              dimnames =  list(dat$mrk.names[tpt.ll$seq.num],
                               dat$mrk.names[tpt.ll$seq.num]))
  A<-sapply(tpt.ll$pairwise, expected_max_connections)
  M[lower.tri(M)] <- as.numeric(A)
  M[upper.tri(M)] <- t(M)[upper.tri(M)]
  M
}
split_in_phase_blocks <- function(input.seq, M.ll){
  L <- vector("list", length(unique(input.seq$seq.dose.p1)))
  for(i in sort(unique(input.seq$seq.dose.p1))){
    cat("\n", i)
    id <- input.seq$seq.mrk.names[which(input.seq$seq.dose.p1 == i)]
    if(length(id) == 1){
      L[[i]][[1]] <- id
      next()
    }
    Mnew <- round(M.ll[id, id],2)
    if(any(Mnew != 0, na.rm = T))
      Mnew <- Mnew/max(Mnew, na.rm = T)
    Mnew <- 1 - Mnew
    Mnew[Mnew >= 0.5] <- 1
    Mnew[Mnew < 0.5] <- 0
    diag(Mnew) <- 0
    hc <- hclust(dist(Mnew), "complete")
    #plot(hc)
    #abline(h = 0.05)
    ct <- cutree(hc, h = 0.05)
    L[[i]]<-vector("list", length(unique(ct)))
    names(L[[i]]) <- paste0("SupMrk",i,"_", unique(ct))
    id2 <- match(names(sort(ct)), map.orig$info$mrk.names)
    map <- get_submap(map.orig, mrk.pos = id2, reestimate.rf = F, verbose = FALSE)
    plot(map, mrk.names = TRUE)
    for(j in unique(ct))
      L[[i]][[j]] <- names(which(ct == j))
  }
  L
}