get_codom_level<-function(map){
  A<-cbind(ph_list_to_matrix(L = map$maps[[1]]$seq.ph$P, m = map$info$m),
           ph_list_to_matrix(L = map$maps[[1]]$seq.ph$Q, m = map$info$m))
  codom<-length(unique(apply(A, 2, paste, collapse  = "")))
  return(codom)
}


## hap1: firts haplotype block
## hap2: second haplotype block

## Function
my_map_func_test<-function(hap1, hap2,
                           probs.hap1,
                           thresh.cut.path = 10e-3, 
                           thresh.twopt = 5, 
                           tol = 10e-4,
                           rf.matrix){
  ## Hash table: homolog combination --> states to visit in both parents
  A<-as.matrix(expand.grid(0:(sqrt(dim(probs.hap1$probs)[1])-1), 
                           0:(sqrt(dim(probs.hap1$probs)[1])-1))[,2:1])
  rownames(A) <- dimnames(probs.hap1$probs)[[1]]
  ## h: states to visit in both parents
  ## e: probability distribution 
  e.fix <- h.fix <- vector("list", dim(probs.hap1$probs)[3])
  for(i in 1:dim(probs.hap1$probs)[3]){
    a <- probs.hap1$probs[,dim(probs.hap1$probs)[2],i]  
    e.fix[[i]] <- a[a > thresh.cut.path]
    h.fix[[i]] <- A[names(e.fix[[i]]), , drop = FALSE]
  }
  ## Second block
  w<-generate_all_link_phases_elim_equivalent_haplo(block1 = hap1$maps[[1]], 
                                                    block2 = hap2, 
                                                    rf.matrix = rf.matrix, 
                                                    m = hap1$info$m)
  
  ## get testmaps and conditional probabilities
  testmaps<-p<-vector("list", length(w))
  flag <- FALSE
  if(is.numeric(hap2)){
    for(i in 1:length(w))
    {
      suppressMessages(hap.temp <- get_submap(hap1, 
                                              c(hap1$info$n.mrk,hap1$info$n.mrk), 
                                              reestimate.rf = FALSE))
      hap.temp$maps[[1]]$seq.num<-rep(hap2, 2)
      hap.temp$maps[[1]]$seq.ph <- list(P = c(w[[i]]$P, w[[i]]$P),
                                        Q = c(w[[i]]$Q, w[[i]]$Q))
      hap.temp$maps[[1]]$seq.rf <- 10e-6
      testmaps[[i]] <- hap.temp
      
      
      
      ##############################################
      ##############################################
      ##############################################
      ## Implement usage of index function CPP to R
      p[[i]] <- calc_genoprob(testmaps[[i]], verbose = FALSE)
    }
    flag <- TRUE
  } else { 
    for(i in 1:length(w))
    {
      hap2$maps[[1]]$seq.ph <- w[[i]]
      testmaps[[i]] <- hap2
      p[[i]] <- calc_genoprob(testmaps[[i]], verbose = FALSE)
    }
  }
  ## h: states to visit in both parents
  ## e: probability distribution 
  h<-e<-vector("list", length(w))
  for(j in 1:length(w)){
    etemp<-htemp<-vector("list", dim(probs.hap1$probs)[3])
    for(i in 1:dim(probs.hap1$probs)[3]){
      a <- p[[j]]$probs[,1,i]  
      etemp[[i]] <- a[a > thresh.cut.path]
      htemp[[i]] <- A[names(etemp[[i]]), , drop = FALSE]
    }
    h[[j]] <- htemp
    e[[j]] <- etemp
  }
  configs<-vector("list", length(testmaps))
  names(configs) <- paste0("Phase_config.", 1:length(testmaps))
  res<-matrix(NA, nrow = length(h), ncol = 2, dimnames = list(names(configs), c("log_like", "rf")))
  probs.temp<-NULL
  ## HMM
  for(i in 1:length(h)){
    #cat("testing", i, "of", length(h), "\n")
    cat(".")
    h.test<-c(list(h.fix), h[i])
    names(h.test) <- c("hap1", "hap2")
    e.test<-c(list(e.fix), e[i])
    restemp<-est_haplo_hmm(m = hap1$info$m, 
                           n.mrk = length(h.test), 
                           n.ind = dim(probs.hap1$probs)[3], 
                           haplo = h.test, 
                           #emit = e.test, 
                           rf_vec = rep(0.01, length(h.test)-1), 
                           verbose = FALSE, 
                           use_H0 = FALSE, 
                           tol = tol) 
    res[i,]<-unlist(restemp)
    probs.temp[[i]]<-calc_genoprob_haplo(m = hap1$info$m, n.mrk = length(h.test), 
                                     n.ind = dim(probs.hap1$probs)[3], 
                                     haplo = h.test,
                                     indnames = BT.trifida.triloba$ind.names,
                                     rf_vec = res[i, 2], verbose = FALSE) 
  }
  for(i in 1:length(testmaps)){
    if(flag){
      P<-c(hap1$maps[[1]]$seq.ph$P, testmaps[[i]]$maps[[1]]$seq.ph$P[1])
      Q<-c(hap1$maps[[1]]$seq.ph$Q, testmaps[[i]]$maps[[1]]$seq.ph$Q[1])
      names(P)<-names(Q)<-c(hap1$maps[[1]]$seq.num, testmaps[[i]]$maps[[1]]$seq.num[1])
      configs[[i]]<-list(P = P, Q = Q)
    } else {
      P<-c(hap1$maps[[1]]$seq.ph$P, testmaps[[i]]$maps[[1]]$seq.ph$P)
      Q<-c(hap1$maps[[1]]$seq.ph$Q, testmaps[[i]]$maps[[1]]$seq.ph$Q)
      names(P)<-names(Q)<-c(hap1$maps[[1]]$seq.num, testmaps[[i]]$maps[[1]]$seq.num)
      configs[[i]]<-list(P = P, Q = Q)
    }
  }
  cat("\n")
  res<-res[order(res[,"log_like"], decreasing = TRUE),,drop = FALSE]
  return(list(stats = res, phase.config  = configs, probs = probs.temp))  
}
