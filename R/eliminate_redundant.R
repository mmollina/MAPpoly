#' Eliminate redundant markers
#'
#' Eliminate markers with identical dosage infomation throughout individuals.
#'
#' @param input.seq an object of class \code{mappoly.sequence}
#' 
#' @param ... currently ignored
#'
#' @return An object of class \code{mappoly.unique.seq} which
#'     is a list  containing the following components:
#'     \item{unique.seq}{an object of class \code{mappoly.sequence}
#'           with the redundant markers removed}
#'     \item{kept}{a vector containing the name of the informative markers}
#'     \item{eliminated}{a vector containing the name of the non-informative (eliminated) markers}
#'
#' @examples
#'   \dontrun{
#'     data(hexafake)
#'     all.mrk<-make_seq_mappoly(hexafake, 'all')
#'     red.mrk<-elim_redundant(all.mrk)
#'     plot(red.mrk)
#'     unique.mrks<-make_seq_mappoly(red.mrk)
#'    }
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
#' @export elim_redundant
#'
elim_redundant<-function(input.seq)
{
  x<-get(input.seq$data.name, pos = 1)$geno.dose[input.seq$seq.num, ]
    dat_temp <- unique(x, dimnames = TRUE)
  output.seq <- make_seq_mappoly(get(input.seq$data.name, pos = 1), rownames(dat_temp), data.name = input.seq$data.name)
  output.seq$chisq.pval.thres <-input.seq$chisq.pval.thres
  output.seq$chisq.pval <-input.seq$chisq.pval
  mrknames <- get(input.seq$data.name, pos = 1)$mrk.names
  elim<-setdiff(input.seq$seq.num,output.seq$seq.num)
  n1<-apply(dat_temp, 1, paste, collapse="")
  n2<-apply(x[mrknames[elim],], 1, paste, collapse="")
  elim.out <- data.frame(kept = rownames(dat_temp)[match(n2,n1)], elim = mrknames[elim])
  structure(list(unique.seq = output.seq, kept = mrknames[output.seq$seq.num],
                 elim.correspondence =  elim.out),
            class = "mappoly.unique.seq")
}

#' @export
plot.mappoly.unique.seq<-function(x, ...)
{
  slc <- c(nrow(x$elim.correspondence), length(x$kept))
  lbls <- c("eliminated", "kept")
  pct <- round(slc/sum(slc)*100)
  pct <- paste0(slc, " (", pct, "%)")
  lbls <- paste(lbls, pct, sep = "\n")
  pie(slc, labels = lbls)
}

#' @export
print.mappoly.unique.seq<-function(x, ...)
{
  print(x$unique.seq)
  cat("------------\n")
  cat("Eliminated markers: ", nrow(x$elim.correspondence), "\n")
}

