#' Look at genotypes that were imputed or changed by the HMM chain given a level of global genotypic error
#'
#' Outputs a graphical representation ggplot with the percent of data changed
#'
#' @param map_list a list of multiple \code{mappoly.map.list}
#'
#' @param error error rate used in global error in the `calc_genoprob_error()`
#' @param verbose T or F for `calc_genoprob_error()` and `calc_homologprob()`
#' 
#'
#' @return A ggplot of the changed and imputed genotypic dosages
#'
#' @examples
#'     plot_progeny_dosage_change(map_list=MAPs, error=0.05)
#'
#' @author Jeekin Lau, \email{jzl0026@tamu.edu}, with optimization by Cristiane Taniguti, \email{chtaniguti@tamu.edu}
#'
#' @import ggplot2 
#' @import reshape2
#' @export plot_progeny_dosage_change

plot_progeny_dosage_change <- function(map_list, error, verbose=T){
  map=map_list
  if(!exists(map[[1]]$info$data.name)) stop("mappoly.data object not here")
  
  dat <- get(map[[1]]$info$data.name)
  ploidy <- map[[1]]$info$ploidy
  
  print("calculating genoprob error")
  
  genoprob  <- vector("list", 7)
  for(i in 1:7){
    genoprob[[i]] <- calc_genoprob_error(input.map = map[[i]], error = error, verbose = verbose)
  }
  
  print("calculating homologprob")
  
  homoprobs = calc_homologprob(genoprob, verbose=verbose)
  
  print("comparing to orginal")
  
  P=unlist(lapply(map, function(x) x$maps[[1]]$seq.ph$P), recursive=F)
  Q=unlist(lapply(map, function(x) x$maps[[1]]$seq.ph$Q), recursive=F)
  
  
  
  P_matrix=matrix(0,nrow=length(P),ncol=ploidy)
  Q_matrix=matrix(0,nrow=length(Q),ncol=ploidy)
  #colnames(P_matrix)=c("a","b","c","d")
  #colnames(Q_matrix)=c("e","f","g","h")
  
  
  for(i in 1:nrow(P_matrix)){
    P_matrix[i,unlist(P[i])]=1
  }
  
  for(i in 1:nrow(Q_matrix)){
    Q_matrix[i,unlist(Q[i])]=1
  }
  
  mrks_mapped=unlist(lapply(map, function(x) x$info$mrk.names) )
  
  PQ_matrix=cbind(P_matrix,Q_matrix)
  rownames(PQ_matrix)=mrks_mapped
  colnames(PQ_matrix)=c(letters[1:(ploidy*2)])
  
  
  homoprob=homoprobs$homoprob
  inds = unique(homoprob$individual)
  mrks = unique(homoprob$marker)
  
  
  
  temp = matrix(NA,length(mrks),length(inds))
  rownames(temp)=mrks
  colnames(temp)=inds
  
  
  homoprob=as.matrix(homoprob)
  
  
  
  homoprob_ind <- split.data.frame(homoprob, homoprob[,3])
  
  test <- sapply(homoprob_ind, function(x) {
    by_marker <- split.data.frame(x, x[,1]) 
    final_matrix <- t(sapply(by_marker, function(y) y[order(y[,4], decreasing = T),][1:4,2]))
    final_vector <- apply(final_matrix, 1, function(w) paste0(w, collapse = ""))
    return(final_vector)
  })
  
  
  
  test=test[,order(colnames(test))]
  test=test[match(mrks,rownames(test)),]
  
  finished = test

  
  
  for(a in 1:length(mrks)){
    for(b in 1:length(inds)){
      finished[a,b]=sum(PQ_matrix[a,substr(test[a,b], 1,1)],
                        PQ_matrix[a,substr(test[a,b], 2,2)],
                        PQ_matrix[a,substr(test[a,b], 3,3)],
                        PQ_matrix[a,substr(test[a,b], 4,4)])
      
    }
    #print(a)
  }
  
  
  
  original_geno=as.matrix(dat$geno.dose[which(rownames(dat$geno.dose)%in%mrks),])
  original_geno=original_geno[,order(colnames(original_geno))]
  identical(colnames(finished), colnames(original_geno))
  
  identical(original_geno, finished)
  which(!original_geno==finished)
  original_geno[which(!original_geno==finished)]
  
  
  percent_imputed=length(which(!original_geno==finished&original_geno==5))/length(original_geno)*100
  percent_changed=length(which(!original_geno==finished&!original_geno==5))/length(original_geno)*100
  
  
  
  
  colors = c("red","chartreuse","black")
  empty_matrix=matrix(0,nrow(original_geno),ncol(original_geno))
  empty_matrix[which(original_geno==finished)]="unchanged"
  empty_matrix[which(!original_geno==finished&original_geno==5)]="imputed"
  empty_matrix[which(!original_geno==finished&!original_geno==5)]="changed"
  empty_matrix_melt=melt(empty_matrix)
  plot1<-ggplot(empty_matrix_melt, aes(Var1, Var2, fill= factor(value))) + 
    geom_tile()+scale_fill_manual(values=colors)+
    xlab("Markers")+
    ylab("Individuals")+
    ggtitle(paste0("changed = ",round(percent_changed, digits=3),"% ","imputed = ", round(percent_imputed, digits=3),"%"))
  
  print("done")
  
  return(plot1) 
  
}


