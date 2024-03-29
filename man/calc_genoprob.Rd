% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calc_genoprob.R
\name{calc_genoprob}
\alias{calc_genoprob}
\title{Compute conditional probabilities of the genotypes}
\usage{
calc_genoprob(input.map, step = 0, phase.config = "best", verbose = TRUE)
}
\arguments{
\item{input.map}{An object of class \code{mappoly.map}}

\item{step}{Maximum distance (in cM) between positions at which
the genotype probabilities are calculated, though for
step = 0, probabilities are calculated only at the
marker locations.}

\item{phase.config}{which phase configuration should be used. "best" (default)
will choose the phase configuration associated with the
maximum likelihood}

\item{verbose}{if \code{TRUE} (default), current progress is shown; if
\code{FALSE}, no output is produced}
}
\value{
An object of class 'mappoly.genoprob' which has two elements: a tridimensional
array containing the probabilities of all possible genotypes for each individual
in each marker position; and the marker sequence with it's recombination frequencies
}
\description{
Conditional genotype probabilities are calculated for each marker
position and each individual given a map.
}
\examples{
 ## tetraploid example
 probs.t <- calc_genoprob(input.map = solcap.dose.map[[1]],
                        verbose = TRUE)
 probs.t
 ## displaying individual 1, 36 genotypic states
 ## (rows) across linkage group 1 (columns)                          
 image(t(probs.t$probs[,,1]))

}
\references{
Mollinari, M., and Garcia, A.  A. F. (2019) Linkage
    analysis and haplotype phasing in experimental autopolyploid
    populations with high ploidy level using hidden Markov
    models, _G3: Genes, Genomes, Genetics_. 
    \doi{10.1534/g3.119.400378}
}
\author{
Marcelo Mollinari, \email{mmollin@ncsu.edu}
}
