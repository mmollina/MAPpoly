% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/make_seq.R
\name{make_seq_mappoly}
\alias{make_seq_mappoly}
\alias{print.mappoly.sequence}
\alias{plot.mappoly.sequence}
\title{Create a Sequence of Markers}
\usage{
make_seq_mappoly(
  input.obj,
  arg = NULL,
  data.name = NULL,
  info.parent = c("all", "p1", "p2"),
  genomic.info = NULL
)

\method{print}{mappoly.sequence}(x, ...)

\method{plot}{mappoly.sequence}(x, ...)
}
\arguments{
\item{input.obj}{An object belonging to one of the specified classes: \code{mappoly.data},
\code{mappoly.map}, \code{mappoly.sequence}, \code{mappoly.group}, \code{mappoly.unique.seq},
\code{mappoly.pcmap}, \code{mappoly.pcmap3d}, \code{mappoly.geno.ord}, or \code{mappoly.edit.order}.}

\item{arg}{Specifies the markers to include in the sequence, accepting several formats: a string 'all' for all
markers; a string or vector of strings 'seqx' where x is the sequence number (0 for unassigned markers); a
vector of integers indicating specific markers; or a vector of integers representing linkage group numbers if
\code{input.obj} is of class \code{mappoly.group}. For certain classes (\code{mappoly.pcmap}, \code{mappoly.pcmap3d},
\code{mappoly.unique.seq}, or \code{mappoly.geno.ord}), \code{arg} can be \code{NULL}.}

\item{data.name}{Name of the \code{mappoly.data} class object.}

\item{info.parent}{Selection criteria based on parental information: \code{'all'} for all dosage combinations,
\code{'P1'} for markers informative in parent 1, or \code{'P2'} for markers informative in parent 2. Default
is \code{'all'}.}

\item{genomic.info}{Optional and applicable only to \code{mappoly.group} objects. Specifies the use of genomic
information in sequence creation. With \code{NULL} (default), all markers defined by the grouping function are
included. Numeric values indicate the use of specific sequences from genomic information, aiming to match the
maximum number of markers with the group. Supports single values or vectors for multiple sequence consideration.}

\item{x}{An object of class \code{mappoly.sequence}.}

\item{...}{Currently ignored.}
}
\value{
Returns an object of class `mappoly.sequence`, comprising:
  \item{"seq.num"}{Ordered vector of marker indices according to the input.}
  \item{"seq.phases"}{List of linkage phases between markers; -1 for undefined phases.}
  \item{"seq.rf"}{Vector of recombination frequencies; -1 for not estimated frequencies.}
  \item{"loglike"}{Log-likelihood of the linkage map.}
  \item{"data.name"}{Name of the `mappoly.data` object with raw data.}
  \item{"twopt"}{Name of the `mappoly.twopt` object with 2-point analyses; -1 if not computed.}
}
\description{
Constructs a sequence of markers based on an object belonging to various specified classes. This
function is versatile, supporting multiple input types and configurations for generating marker sequences.
}
\examples{
all.mrk <- make_seq_mappoly(hexafake, 'all')
seq1.mrk <- make_seq_mappoly(hexafake, 'seq1')
plot(seq1.mrk)
some.mrk.pos <- c(1,4,28,32,45)
some.mrk.1 <- make_seq_mappoly(hexafake, some.mrk.pos)
plot(some.mrk.1)

}
\references{
Mollinari, M., and Garcia, A. A. F. (2019). Linkage analysis and haplotype phasing in experimental
autopolyploid populations with high ploidy level using hidden Markov models. _G3: Genes|Genomes|Genetics_,
\doi{10.1534/g3.119.400378}.
}
\author{
Marcelo Mollinari \email{mmollin@ncsu.edu}, with modifications by Gabriel Gesteira 
\email{gdesiqu@ncsu.edu}
}
