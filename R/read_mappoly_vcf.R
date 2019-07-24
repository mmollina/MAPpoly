#' Data Input VCF
#'
#' Reads an external VCF file. This function accepts version 4.0 and higher, and creates an object of class \code{mappoly.data}
#' 
#' This function can handle .vcf files of version 4.0 and higher. The ploidy can be automatically detected, but you
#' should inform it to check mismatches. All individual and marker names will be kept as they are in the .vcf file.
#'
#' @param file.in name of input file which contains the data to
#'     be read.
#'
#' @param ploidy the species ploidy (optional, it will be automatically detected)
#'
#' @param filter.non.conforming if \code{TRUE} (default) exclude samples with non 
#'     expected genotypes under random chromosome pairing and no double reduction 
#'     
#' @param thresh.line threshold used for p-values on segregation test
#' 
#' @param parent.1 name of parent 1
#' 
#' @param parent.2 name of parent 2
#' 
#' @param update.prob Logical. Should MAPpoly update genotype probabilities based on HMM? (default = FALSE)
#' 
#' @param output if \code{update.prob = TRUE}, defines the output of the updated .vcf file
#' 
#' @param ... currently ignored
#'
#' @return An object of class \code{mappoly.data} which contains a
#'     list with the following components:
#'     \item{m}{ploidy level}
#'     \item{n.ind}{number individuals}
#'     \item{n.mrk}{total number of markers}
#'     \item{ind.names}{the names of the individuals}
#'     \item{mrk.names}{the names of the markers}
#'     \item{dosage.p}{a vector containing the dosage in
#'       parent P for all \code{n.mrk} markers}
#'     \item{dosage.q}{a vector containing the dosage in
#'       parent Q for all \code{n.mrk} markers}
#'     \item{sequence}{a vector indicating which sequence each marker
#'       belongs. Zero indicates that the marker was not assigned to any
#'       sequence.}
#'     \item{sequence.pos}{Physical position of the markers into the
#'       sequence}
#'     \item{geno.dose}{a matrix containing the dosage for each markers (rows) 
#'       for each individual (columns). Missing data are represented by 
#'       \code{ploidy_level + 1}.}
#'     \item{input.file}{Full path to input file, used when \code{update.prob = TRUE}}
#' @examples
#' \dontrun{
#'     mydata <- read_vcf(hexasubset, parent.1 = "P1", parent.2 = "P2")
#'     print(mydata, detailed = TRUE)
#'}
#' @author Gabriel Gesteira, \email{gabrielgesteira@usp.br}
#'
#' @references
#'     Mollinari, M., and Garcia, A.  A. F. (2019) Linkage
#'     analysis and haplotype phasing in experimental autopolyploid
#'     populations with high ploidy level using hidden Markov
#'     models, _submited_. \url{https://doi.org/10.1101/415232}
#'
#' @export read_vcf

read_vcf <- function(file.in, filter.non.conforming = TRUE, parent.1, parent.2, ploidy = NA,
                     thresh.line = 0.05, update.prob = FALSE, output = NULL) {
  # Checking even ploidy
  if(!is.na(ploidy) && (ploidy %% 2) != 0){
    stop("MAPpoly only supports even ploidy levels. Please check the 'ploidy' parameter and try again.")
  }
  input.file = normalizePath(file.in) # Getting full input file path
  input.size = file.size(file.in) / 1024000 # Getting input file size in MB
  if (input.size > 3000){
    warning("Your VCF file is greater than 3 GB. Check for available RAM memory.")
  }
  cat("Reading data...\n")
  input.data = vcfR::read.vcfR(file.in) # Reading file
  ind.names = colnames(input.data@gt)[-1]
  n.mrk = dim(input.data@gt)[1] # Getting number of markers
  n.ind = length(ind.names) - 2 # Number of individuals excepting two parents
  sequence = input.data@fix[,1] # Getting chromosome information
  sequence.pos = input.data@fix[,2] # Getting positions
  mrk.names = input.data@fix[,3] # Getting marker names
  cat("Processing genotypes... (this may take some minutes)\n")
  cname = which(unlist(strsplit(unique(input.data@gt[,1]), ":")) == "GT") # Defining GT position
  geno.dose = matrix(unlist(lapply(strsplit(input.data@gt[,-1], ":"), "[", cname)), nrow = n.mrk, byrow = F) # Selecting genotypes
  file.ploidy = length(unlist(strsplit(unique(geno.dose[,1])[1], "/"))) # Checking ploidy
  if (sum(grepl("\\|", geno.dose)) > 0){ # Checking phased data
    input.phased = TRUE # Treat this in a different way when reading file
    warning("Phased genotypes detected. Should MAPpoly consider this information?")
  } else {input.phased = FALSE}
  if ((file.ploidy %% 2) != 0){ # Checking odd ploidy level
    stop("Your VCF file shows an odd ploidy level, but MAPpoly only supports even ploidy levels. Please check your VCF file and try again.")
  }
  if (!is.na(ploidy) && (file.ploidy != ploidy)){ # Checking informed and file ploidy
    warning("Informed ploidy doesn't match the detected ploidy. Using detected ploidy instead.")
  }
  if (!is.na(ploidy) && (ploidy > 2) && (file.ploidy == 2)){ # Checking absence of dosages
    stop("Informed ploidy is ",ploidy, ", but detected ploidy is ", file.ploidy, ".\nYou should provide allelic dosages for all individuals. You can estimate allelic dosages using packages such as 'SuperMASSA', 'updog', 'fitTetra' and others. We are working for a integrated function to estimate dosages, which will be available soon.")
  }
  m = file.ploidy # Setting ploidy level
  if (!(parent.1 %in% ind.names) | !(parent.2 %in% ind.names)){
    stop("Provided parents were not found in VCF file. Please check it and try again.")
  }
  lost.data = paste(c(rep("./", (m-1)), "."), collapse = "")
      for (i in 1:length(geno.dose)){
          if (geno.dose[i] == lost.data){
            geno.dose[i] = NA
        } else {geno.dose[i] = length(which(unlist(strsplit(geno.dose[i], "/")) == 0))}
       }
  geno.dose = matrix(as.numeric(geno.dose), nrow = n.mrk, byrow = F)
  colnames(geno.dose) = ind.names
  rownames(geno.dose) = mrk.names
  dosage.p = geno.dose[,which(colnames(geno.dose) == parent.1)] # Selecting dosages for parent 1
  dosage.q = geno.dose[,which(colnames(geno.dose) == parent.2)] # Selecting dosages for parent 2
  geno.dose = geno.dose[, -c(which(colnames(geno.dose) %in% c(parent.1, parent.2)))] # Updating geno.dose matrix
  ind.names = ind.names[-c(which(ind.names %in% c(parent.1, parent.2)))] # Updating individual names
  geno.dose = data.frame(geno.dose)
  
  ## monomorphic markers
  dp = abs(abs(dosage.p-(m/2))-(m/2))
  dq = abs(abs(dosage.q-(m/2))-(m/2))
  id = dp+dq!=0
  
  cat("Done!\n")
  cat("Read the following data:")
  cat("\n    Ploidy level:", m)
  cat("\n    No. individuals: ", n.ind)
  cat("\n    No. markers: ", n.mrk) 
  cat("\n    No. informative markers:  ", sum(id), " (", round(100*sum(id)/n.mrk,1), "%)", sep = "")
  # if (all(unique(nphen) != 0))
  #   cat("\n    This dataset contains phenotypic information.")
  
  if (length(sequence) > 1)
    cat("\n    This dataset contains sequence information.")
  cat("\n    ...")
  ## get genotypic info --------------------
  if(nrow(geno.dose)!=length(mrk.names))
    stop("\n\t\t-------------------------------------
         Number of marker names is different from
         the number of markeres in the dataset.
         Please, check data.
         ------------------------------------------\n")
  if(ncol(geno.dose)!=length(ind.names))
    stop("\n\t\t-------------------------------------
         Number of individual names is different from
         the number of individuals in the dataset.
         Please, check data.
         ------------------------------------------\n")
  dimnames(geno.dose)<-list(mrk.names, ind.names)
  geno.dose[is.na(geno.dose)] <- m + 1
  ## returning the 'mappoly.data' object
  cat("\n    Done with reading.\n")
  
  
  res = structure(list(m = m,
                        n.ind = n.ind,
                        n.mrk = sum(id),
                        ind.names = ind.names,
                        mrk.names = mrk.names[id],
                        dosage.p = dosage.p[id],
                        dosage.q = dosage.q[id],
                        sequence = sequence[id],
                        sequence.pos = sequence.pos[id],
                        prob.thres = NULL,
                        geno.dose = geno.dose,
                        nphen = 0,
                        phen = NULL,
                        input.file = input.file,
                        input.phased = input.phased
                        ),
                   class = "mappoly.data")
  
# Not working: fix it  
  # if(filter.non.conforming){
  #   cat("    Filtering non-conforming markers.\n    ...")
  #   res<-filter_non_conforming_classes(res)
  #   cat("\n    Performing chi-square test.\n    ...")
  #   ##Computing chi-square p.values
  #   Ds <- array(NA, dim = c(m+1, m+1, m+1))
  #   for(i in 0:m)
  #     for(j in 0:m)
  #       Ds[i+1,j+1,] <- segreg_poly(m = m, dP = i, dQ = j)
  #   Dpop<-cbind(res$dosage.p, res$dosage.q)
  #   M<-t(apply(Dpop, 1, function(x) Ds[x[1]+1, x[2]+1,]))
  #   dimnames(M)<-list(res$mrk.names, c(0:m))
  #   M<-cbind(M, res$geno.dose)
  #   res$chisq.pval<-apply(M, 1, mrk_chisq_test, m = m)
  #   cat("\n    Done.\n")
  #   return(res)
  # }
  return(res)
}
