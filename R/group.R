#' Assign markers to linkage groups
#'
#' Identifies linkage groups of markers using the results of two-point
#' (pairwise) analysis.
#'
#' @param input.mat an object of class \code{mappoly.rf.matrix}
#'
#' @param expected.groups when available, inform the number of expected 
#' linkage groups (i.e. chromosomes) for the species
#'
#' @param inter if \code{TRUE} (default), plots a dendrogram highlighting the
#'    expected groups before continue
#'
#' @param comp.mat if \code{TRUE}, shows a comparison between the reference
#'     based and the linkage based grouping, if the sequence information is
#'     available (default = FALSE)
#'
#' @param verbose logical. If \code{TRUE} (default), current progress is shown;
#'     if \code{FALSE}, no output is produced
#'
#' @return Returns an object of class \code{mappoly.group}, which is a list
#'     containing the following components:
#'     \item{data.name}{the referred dataset name}
#'     \item{hc.snp}{a list containing information related to 
#'     the UPGMA grouping method}
#'     \item{expected.groups}{the number of expected linkage groups}
#'     \item{groups.snp}{the groups to which each of the markers belong}
#'     \item{seq.vs.grouped.snp}{comparison between the genomic group information
#'     (when available) and the groups provided by \code{group_mappoly}}
#'     \item{chisq.pval.thres}{the threshold used on the segregation test when reading the dataset}
#'     \item{chisq.pval}{the p-values associated with the segregation test for all markers in the sequence}
#'
#' @examples
#'     ## Getting first 20 markers from two linkage groups
#'     all.mrk <- make_seq_mappoly(hexafake, c(1:20,601:620))
#'     red.mrk <- elim_redundant(all.mrk)
#'     unique.mrks <- make_seq_mappoly(red.mrk)
#'     counts <- cache_counts_twopt(unique.mrks, cached = TRUE)
#'     all.pairs <- est_pairwise_rf(input.seq = unique.mrks,
#'                                  count.cache = counts,
#'                                  ncpus = 1,
#'                                  verbose=TRUE)
#'
#'     ## Full recombination fraction matrix
#'     mat.full<-rf_list_to_matrix(input.twopt=all.pairs)
#'     plot(mat.full, index = FALSE)
#'
#'     lgs <- group_mappoly(input.mat = mat.full,
#'                          expected.groups = 2,
#'                          inter = TRUE,
#'                          comp.mat = TRUE, #this data has physical information
#'                          verbose = TRUE)
#'     lgs
#'     plot(lgs)
#'    
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu}
#'
#' @references
#'     Mollinari, M., and Garcia, A.  A. F. (2019) Linkage
#'     analysis and haplotype phasing in experimental autopolyploid
#'     populations with high ploidy level using hidden Markov
#'     models, _G3: Genes, Genomes, Genetics_. 
#'     \url{https://doi.org/10.1534/g3.119.400378} 
#'
#' @importFrom graphics abline pie
#' @importFrom stats as.dendrogram as.dist cutree hclust lm predict quantile rect.hclust
#' @importFrom dendextend color_branches
#' @export group_mappoly

group_mappoly <- function(input.mat, expected.groups = NULL,
                          inter = TRUE, comp.mat = FALSE, verbose = TRUE)
  {
    ## checking for correct object
    input_classes <- c("mappoly.rf.matrix")
    if (!inherits(input.mat, input_classes)) {
      stop(deparse(substitute(input.mat)), " is not an object of class 'mappoly.rf.matrix'")
    }
    # if(is.null(input.seq))
      input.seq<-make_seq_mappoly(input.obj = get(input.mat$data.name, pos=1), 
                                  arg = rownames(input.mat$rec.mat), 
                                  data.name = input.mat$data.name)
    if(!setequal(intersect(input.seq$seq.mrk.names,colnames(input.mat$rec.mat)), input.seq$seq.mrk.names)){
      stop(deparse(substitute(input.mat)), " does not contain all markers present in", deparse(substitute(input.seq)))
    }
    MSNP <- input.mat$rec.mat[input.seq$seq.mrk.names, input.seq$seq.mrk.names]
    mn<-input.seq$sequence
    mn[is.na(mn)]<-"NH"
    dimnames(MSNP)<-list(mn, mn)
    diag(MSNP)<-0
    MSNP[is.na(MSNP)]<-.5
    hc.snp<-hclust(as.dist(MSNP), method = "average")
    ANSWER <- "flag"
    if(interactive() && inter)
    {
      dend.snp <- as.dendrogram(hc.snp)
      while(substr(ANSWER, 1, 1) != "y" && substr(ANSWER, 1, 1) != "yes" && substr(ANSWER, 1, 1) != "Y" && ANSWER !="")
      {
        dend1 <- dendextend::color_branches(dend.snp, k = expected.groups)
        plot(dend1, leaflab = "none")
        if(is.null(expected.groups))
          expected.group <- as.numeric(readline("Enter the number of expected groups: "))
        z<-rect.hclust(hc.snp, k = expected.groups, border = "red")
        groups.snp  <- cutree(tree = hc.snp, k = expected.groups)
        xy<-sapply(z, length)
        xt<-as.numeric(cumsum(xy)-ceiling(xy/2))
        yt<-.1
        points(x = xt, y = rep(yt, length(xt)), cex = 6, pch = 20, col = "lightgray")
        text(x = xt, y = yt, labels = pmatch(xy, table(groups.snp, useNA = "ifany")), adj = .5)
        ANSWER <- readline("Enter 'Y/n' to proceed or update the number of expected groups: ")
        if(substr(ANSWER, 1, 1) == "n" | substr(ANSWER, 1, 1) == "no" | substr(ANSWER, 1, 1) == "N")
            stop("You decided to stop the function.")
        if(substr(ANSWER, 1, 1) != "y" && substr(ANSWER, 1, 1) != "yes" && substr(ANSWER, 1, 1) != "Y" && ANSWER !="")
          expected.groups <- as.numeric(ANSWER)
      }
    }
    if(is.null(expected.groups))
      stop("Inform the 'expected.groups' or use 'inter = TRUE'")

    # Distribution of SNPs into linkage groups
    seq.vs.grouped.snp <- NULL
    if(all(unique(mn) == "NH") & comp.mat)
    {
      comp.mat <- FALSE
      seq.vs.grouped.snp <- NULL
      warning("There is no physical reference to generate a comparison matrix")
    }
    groups.snp  <- cutree(tree = hc.snp, k = expected.groups)
    if(comp.mat){
      seq.vs.grouped.snp<-matrix(0, expected.groups, length(na.omit(unique(input.seq$sequence)))+1,
                                 dimnames = list(1:expected.groups, c(na.omit(unique(input.seq$sequence)),"NH")))
      for(i in 1:expected.groups)
      {
        x<-table(names(which(groups.snp==i)))
        seq.vs.grouped.snp[i,names(x)]<-x
      }
      idtemp2 <- apply(seq.vs.grouped.snp, 1, which.max)
      seq.vs.grouped.snp<-cbind(seq.vs.grouped.snp[,unique(idtemp2)], seq.vs.grouped.snp[,"NH"])
      cnm<-colnames(seq.vs.grouped.snp)
      cnm[cnm==""]<-"NoChr"
      colnames(seq.vs.grouped.snp)<-cnm
    } else {
      seq.vs.grouped.snp <- NULL
    }
    names(groups.snp)<-input.seq$seq.num
    structure(list(data.name = input.mat$data.name, 
                   hc.snp = hc.snp, 
                   expected.groups = expected.groups,
                   groups.snp = groups.snp, 
                   seq.vs.grouped.snp = seq.vs.grouped.snp, 
                   chisq.pval.thres = input.seq$chisq.pval.thres, 
                   chisq.pval = input.seq$chisq.pval),
                   class = "mappoly.group")
 }

#' @export
print.mappoly.group <- function(x, detailed = TRUE, ...) {
    cat("  This is an object of class 'mappoly.group'\n  ------------------------------------------\n")
    ## criteria
    cat("  Criteria used to assign markers to groups:\n\n")
    cat("    - Number of markers:         ", length(x$groups.snp), "\n")
    cat("    - Number of linkage groups:  ", length(unique(x$groups.snp)), "\n")
    cat("    - Number of markers per linkage groups: \n")
    w<-data.frame(table(x$groups.snp, useNA = "ifany"))
    colnames(w) = c("   group", "n.mrk")
    print (w, row.names = FALSE)
    cat("  ------------------------------------------\n")
    
    ## printing summary
    if(!is.null(x$seq.vs.grouped.snp)){
      print(x$seq.vs.grouped.snp)
      cat("  ------------------------------------------\n")
    }
}

#' @export
plot.mappoly.group <- function(x, ...) {
  dend <- as.dendrogram(x$hc.snp)
  dend1 <- dendextend::color_branches(dend, k = x$expected.groups)
  plot(dend1, leaflab = "none")
  z<-rect.hclust(x$hc.snp, k = x$expected.groups, border = "red")
  xy<-sapply(z, length)
  xt<-as.numeric(cumsum(xy)-ceiling(xy/2))
  yt<-.1
  points(x = xt, y = rep(yt, length(xt)), cex = 6, pch = 20, col = "lightgray")
  text(x = xt, y = yt, labels = pmatch(xy, table(x$groups.snp, useNA = "ifany")), adj = .5)
}
