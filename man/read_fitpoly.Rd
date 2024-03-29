% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_fitpoly.R
\name{read_fitpoly}
\alias{read_fitpoly}
\title{Data Input in fitPoly format}
\usage{
read_fitpoly(
  file.in,
  ploidy,
  parent1,
  parent2,
  offspring = NULL,
  filter.non.conforming = TRUE,
  elim.redundant = TRUE,
  parent.geno = c("joint", "max"),
  thresh.parent.geno = 0.95,
  prob.thres = 0.95,
  file.type = c("table", "csv"),
  verbose = TRUE
)
}
\arguments{
\item{file.in}{a character string with the name of (or full path to) the input file}

\item{ploidy}{the ploidy level}

\item{parent1}{a character string containing the name (or pattern of genotype IDs) of parent 1}

\item{parent2}{a character string containing the name (or pattern of genotype IDs) of parent 2}

\item{offspring}{a character string containing the name (or pattern of genotype IDs) of the offspring 
individuals. If \code{NULL} (default) it considers all individuals as offsprings, except 
\code{parent1} and \code{parent2}.}

\item{filter.non.conforming}{if \code{TRUE} (default) converts data points with unexpected 
genotypes (i.e. no double reduction) to 'NA'. See function \code{\link[mappoly]{segreg_poly}} 
for information on expected classes and their respective frequencies.}

\item{elim.redundant}{logical. If \code{TRUE} (default), removes redundant markers
during map construction, keeping them annotated to in order to include them in the final map.}

\item{parent.geno}{indicates whether to use the joint probability \code{'joint'} (default) or the 
maximum probability of multiple replicates (if available) to assign dosage to parents. 
If there is one observation per parent, both options will yield the same results.}

\item{thresh.parent.geno}{threshold probability to assign a dosage to parents. If the probability 
is smaller than \code{thresh.parent.geno}, the marker is discarded.}

\item{prob.thres}{threshold probability to assign a dosage to offspring. If the probability 
is smaller than \code{prob.thres}, the data point is converted to 'NA'.}

\item{file.type}{indicates whether the characters in the input file are separated by 
'white spaces' ("table") or by commas ("csv").}

\item{verbose}{if \code{TRUE} (default), the current progress is shown; if
\code{FALSE}, no output is produced}
}
\value{
An object of class \code{mappoly.data} which contains a
    list with the following components:
    \item{ploidy}{ploidy level}
    \item{n.ind}{number individuals}
    \item{n.mrk}{total number of markers}
    \item{ind.names}{the names of the individuals}
    \item{mrk.names}{the names of the markers}
    \item{dosage.p1}{a vector containing the dosage in
      parent P for all \code{n.mrk} markers}
    \item{dosage.p2}{a vector containing the dosage in
      parent Q for all \code{n.mrk} markers}
    \item{chrom}{a vector indicating which sequence each marker
      belongs. Zero indicates that the marker was not assigned to any
      sequence}
    \item{genome.pos}{Physical position of the markers into the
      sequence}
    \item{seq.ref}{NULL (unused in this type of data)}
    \item{seq.alt}{NULL (unused in this type of data)}
    \item{all.mrk.depth}{NULL (unused in this type of data)}
    \item{geno.dose}{a matrix containing the dosage for each markers (rows) 
      for each individual (columns). Missing data are represented by 
      \code{ploidy_level + 1}}
    \item{n.phen}{number of phenotypic traits}
    \item{phen}{a matrix containing the phenotypic data. The rows
                correspond to the traits and the columns correspond
                to the individuals}
    \item{kept}{if elim.redundant = TRUE, holds all non-redundant markers}
    \item{elim.correspondence}{if elim.redundant = TRUE, holds all non-redundant markers and
its equivalence to the redundant ones}
}
\description{
Reads an external data file generated as output of \code{\link[fitPoly]{saveMarkerModels}}. 
This function creates an object of class \code{mappoly.data}.
}
\examples{
\donttest{
#### Tetraploid Example
ft <- "https://raw.githubusercontent.com/mmollina/MAPpoly_vignettes/master/data/fitpoly.dat"
tempfl <- tempfile()
download.file(ft, destfile = tempfl)
fitpoly.dat <- read_fitpoly(file.in = tempfl, ploidy = 4, 
                            parent1 = "P1", parent2 = "P2", 
                            verbose = TRUE)
print(fitpoly.dat, detailed = TRUE)
plot(fitpoly.dat)
plot_mrk_info(fitpoly.dat, 37)
}

}
\references{
Voorrips, R.E., Gort, G. & Vosman, B. (2011) Genotype calling 
    in tetraploid species from bi-allelic marker data using mixture 
    models. _BMC Bioinformatics_.
    \doi{10.1186/1471-2105-12-172}
}
\author{
Marcelo Mollinari, \email{mmollin@ncsu.edu}
}
