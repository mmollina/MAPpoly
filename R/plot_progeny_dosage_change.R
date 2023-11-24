#' Display genotypes imputed or changed by the HMM chain given a global genotypic error
#'
#' Outputs a graphical representation ggplot with the percent of data changed. 
#'  
#' @param map_list a list of multiple \code{mappoly.map.list}
#'
#' @param error error rate used in global error in the `calc_genoprob_error()`
#' @param verbose if TRUE (default), current progress is shown; if FALSE, 
#' no output is produced
#' @param output_corrected logical. if FALSE only the ggplot of the changed 
#' dosage is printed, if TRUE then a new corrected dosage matrix is output. 
#'
#' @return A ggplot of the changed and imputed genotypic dosages
#'
#' @examples
#'       x <- get_submap(solcap.err.map[[1]], 1:30, reestimate.rf = FALSE)   
#'       plot_progeny_dosage_change(list(x), error=0.05, output_corrected=FALSE) 
#'       corrected_matrix <- plot_progeny_dosage_change(list(x), error=0.05, 
#'       output_corrected=FALSE) #output corrected
#'
#' @author Jeekin Lau, \email{jzl0026@tamu.edu}, with optimization by Cristiane 
#' Taniguti, \email{chtaniguti@tamu.edu}
#'
#' @import ggplot2 
#' @import reshape2
#' @export plot_progeny_dosage_change

plot_progeny_dosage_change <- function(map_list, 
                                       error, 
                                       verbose = TRUE, 
                                       output_corrected = FALSE){
  Var1 <- Var2 <- value <- NULL
  map=map_list
  if(!exists(map[[1]]$info$data.name)) stop("mappoly.data object not here")
  
  dat <- get(map[[1]]$info$data.name)
  ploidy <- map[[1]]$info$ploidy
  
  print("calculating genoprob error")
  
  genoprob  <- vector("list", length(map))
  for(i in 1:length(map)){
    genoprob[[i]] <- calc_genoprob_error(input.map = map[[i]], error = error, 
                                         verbose = verbose)
  }
  
  print("calculating homologprob")
  
  homoprobs = calc_homologprob(genoprob, verbose=verbose)
  
  print("comparing to orginal")
  
  P=unlist(lapply(map, function(x) x$maps[[1]]$seq.ph$P), recursive=FALSE)
  Q=unlist(lapply(map, function(x) x$maps[[1]]$seq.ph$Q), recursive=FALSE)
  
  
  
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
    final_matrix <- t(sapply(by_marker, function(y) y[order(y[,4], decreasing = TRUE),][1:ploidy,2]))
    final_vector <- apply(final_matrix, 1, function(w) paste0(w, collapse = ""))
    return(final_vector)
  })
  
  
  
  test=test[,order(colnames(test))]
  test=test[match(mrks,rownames(test)),]
  
  finished = test
  
  if (ploidy==2){
    for(a in 1:length(mrks)){
      for(b in 1:length(inds)){
        finished[a,b]=sum(PQ_matrix[a,substr(test[a,b], 1,1)],
                          PQ_matrix[a,substr(test[a,b], 2,2)])}}
  }  
  
  if (ploidy==4){
    for(a in 1:length(mrks)){
      for(b in 1:length(inds)){
        finished[a,b]=sum(PQ_matrix[a,substr(test[a,b], 1,1)],
                          PQ_matrix[a,substr(test[a,b], 2,2)],
                          PQ_matrix[a,substr(test[a,b], 3,3)],
                          PQ_matrix[a,substr(test[a,b], 4,4)])}}
  }  
  
  if (ploidy==6){
    for(a in 1:length(mrks)){
      for(b in 1:length(inds)){
        finished[a,b]=sum(PQ_matrix[a,substr(test[a,b], 1,1)],
                          PQ_matrix[a,substr(test[a,b], 2,2)],
                          PQ_matrix[a,substr(test[a,b], 3,3)],
                          PQ_matrix[a,substr(test[a,b], 4,4)], 
                          PQ_matrix[a,substr(test[a,b], 5,5)],
                          PQ_matrix[a,substr(test[a,b], 6,6)])}}
  }
  
  if (ploidy==8){
    for(a in 1:length(mrks)){
      for(b in 1:length(inds)){
        finished[a,b]=sum(PQ_matrix[a,substr(test[a,b], 1,1)],
                          PQ_matrix[a,substr(test[a,b], 2,2)],
                          PQ_matrix[a,substr(test[a,b], 3,3)],
                          PQ_matrix[a,substr(test[a,b], 4,4)], 
                          PQ_matrix[a,substr(test[a,b], 5,5)],
                          PQ_matrix[a,substr(test[a,b], 6,6)],
                          PQ_matrix[a,substr(test[a,b], 7,7)],
                          PQ_matrix[a,substr(test[a,b], 8,8)])}}
  }
  
  
  
  
  original_geno=as.matrix(dat$geno.dose[which(rownames(dat$geno.dose)%in%mrks),])
  original_geno=original_geno[,order(colnames(original_geno))]
  identical(colnames(finished), colnames(original_geno))
  
  identical(original_geno, finished)
  which(!original_geno==finished)
  original_geno[which(!original_geno==finished)]
  
  
  percent_imputed=length(which(!original_geno==finished&original_geno==ploidy+1))/length(original_geno)*100
  percent_changed=length(which(!original_geno==finished&!original_geno==ploidy+1))/length(original_geno)*100
  
  
  
  
  colors = c("red","chartreuse","black")
  empty_matrix=matrix(0,nrow(original_geno),ncol(original_geno))
  empty_matrix[which(original_geno==finished)]="unchanged"
  empty_matrix[which(!original_geno==finished&original_geno==ploidy+1)]="imputed"
  empty_matrix[which(!original_geno==finished&!original_geno==ploidy+1)]="changed"
  empty_matrix_melt=melt(empty_matrix)
  plot1<-ggplot(empty_matrix_melt, aes(Var1, Var2, fill= factor(value, levels=c("changed","imputed","unchanged")))) + 
    geom_tile()+scale_fill_manual(values=colors, drop=FALSE)+
    xlab("Markers")+
    ylab("Individuals")+
    ggtitle(paste0("changed = ",round(percent_changed, digits=3),"% ","imputed = ", round(percent_imputed, digits=3),"%"))+
    guides(fill=guide_legend(title="Dosage change"))
  print("done")
  
  print(plot1)
  if (output_corrected ==TRUE){
  mrk_names=rownames(finished)
  mrk_index=which(dat$mrk.names%in%mrk_names)
  P1 = dat$dosage.p1[mrk_index]
  P2 = dat$dosage.p2[mrk_index]
  sequence = dat$chrom[mrk_index]
  sequence_position = dat$genome.pos[mrk_index]
  
  hmm_imputed = cbind(P1,P2,sequence,sequence_position, finished)
  return(hmm_imputed)}
}


