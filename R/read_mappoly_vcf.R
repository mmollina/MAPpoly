#' Data Input VCF
#'
#' Reads an external VCF file and creates an object of class \code{mappoly.data}
#' 
#' This function can handle .vcf files versions 4.0 or higher. The ploidy 
#' can be automatically detected, but it is highly recommended that you 
#' inform it to check for mismatches. All individual and marker names 
#' will be kept as they are in the .vcf file.
#'
#' @param file.in a character string with the name of (or full path to) the 
#'   input file which contains the data (VCF format)
#' 
#' @param parent.1 a character string containing the name of parent 1
#' 
#' @param parent.2 a character string containing the name of parent 2
#'
#' @param ploidy the species ploidy (optional, it will be automatically detected)
#' 
#' @param filter.non.conforming if \code{TRUE} (default) converts data points with unexpected 
#'        genotypes (i.e. no double reduction) to 'NA'. See function \code{\link[mappoly]{segreg_poly}} 
#'        for information on expected classes and their respective frequencies.  
#'     
#' @param thresh.line threshold used for p-values on segregation test (default = 0.05)     
#' 
#' @param min.gt.depth minimum genotype depth to keep information. 
#'  If the genotype depth is below \code{min.gt.depth},
#'  it will be replaced with NA (default = 0)
#'
#' @param min.av.depth minimum average depth to keep markers (default = 0)
#'
#' @param max.missing maximum proportion of missing data to keep markers (range = 0-1; default = 1)
#'
#' @param elim.redundant logical. If \code{TRUE} (default), removes redundant markers
#' during map construction, keeping them annotated to export to the final map.
#'
#' @param verbose if \code{TRUE} (default), the current progress is shown; if
#'     \code{FALSE}, no output is produced
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
#'       sequence}
#'     \item{sequence.pos}{Physical position of the markers into the
#'       sequence}
#'     \item{seq.ref}{Reference base used for each marker (i.e. A, T, C, G)}
#'     \item{seq.alt}{Alternative base used for each marker (i.e. A, T, C, G)}
#'     \item{prob.thres}{(unused field)}
#'     \item{geno.dose}{a matrix containing the dosage for each markers (rows) 
#'       for each individual (columns). Missing data are represented by 
#'       \code{ploidy_level + 1}}
#'     \item{nphen}{(unused field)}
#'     \item{phen}{(unused field)}
#'     \item{all.mrk.depth}{DP information for all markers on VCF file}
#'     \item{chisq.pval}{a vector containing p-values related to the chi-squared 
#'     test of Mendelian segregation performed for all markers}
#'     \item{kept}{if elim.redundant=TRUE, holds all non-redundant markers}
#'     \item{elim.correspondence}{if elim.redundant=TRUE, holds all non-redundant markers and
#' its equivalence to the redundant ones}
#' 
#' @examples
#' \donttest{
#' ## Hexaploid sweetpotato: Subset of chromosome 3
#' fl = "https://github.com/mmollina/MAPpoly_vignettes/raw/master/data/sweet_sample_ch3.vcf.gz"
#' tempfl <- tempfile(pattern = 'chr3_', fileext = '.vcf.gz')
#' download.file(fl, destfile = tempfl)
#' dat.dose.vcf = read_vcf(file = tempfl, parent.1 = "PARENT1", parent.2 = "PARENT2")
#' print(dat.dose.vcf)
#' plot(dat.dose.vcf)
#'}
#' 
#' @author Gabriel Gesteira, \email{gabrielgesteira@usp.br}
#' 
#' @references
#'     Mollinari M., Olukolu B. A.,  Pereira G. da S., 
#'     Khan A., Gemenet D., Yencho G. C., Zeng Z-B. (2020), 
#'     Unraveling the Hexaploid Sweetpotato Inheritance 
#'     Using Ultra-Dense Multilocus Mapping, 
#'     _G3: Genes, Genomes, Genetics_. 
#'     \url{https://doi.org/10.1534/g3.119.400620} 
#'     
#'     Mollinari, M., and Garcia, A.  A. F. (2019) Linkage
#'     analysis and haplotype phasing in experimental autopolyploid
#'     populations with high ploidy level using hidden Markov
#'     models, _G3: Genes, Genomes, Genetics_. 
#'     \url{https://doi.org/10.1534/g3.119.400378} 
#'
#' @export read_vcf
#' 
read_vcf = function(file.in, parent.1, parent.2, ploidy = NA, 
                    filter.non.conforming = TRUE, thresh.line = 0.05, 
                    min.gt.depth = 0, min.av.depth = 0, max.missing = 1, 
                    elim.redundant = TRUE, verbose = TRUE) {
  # Checking even ploidy
  if(!is.na(ploidy) && (ploidy %% 2) != 0){
    stop("MAPpoly only supports even ploidy levels. Please check the 'ploidy' parameter and try again.")
  }
  input.file = normalizePath(file.in) # Getting full input file path
  input.size = file.size(file.in) / 1024000 # Getting input file size in MB
  if (input.size > 3000){
    warning("Your VCF file is greater than 3 GB. Check for available RAM memory.")
  }
  if (verbose) cat("Reading data...\n")
  input.data = vcfR::read.vcfR(file.in, verbose = verbose) # Reading file
  ind.names = colnames(input.data@gt)[-1]
  n.mrk = dim(input.data@gt)[1] # Getting number of markers
  n.ind = length(ind.names) - 2 # Number of individuals excepting two parents
  mrk.names = mrk.names.all = input.data@fix[,3] # Getting marker names
  sequence = input.data@fix[,1] # Getting chromosome information
  sequence.pos = as.numeric(input.data@fix[,2]) # Getting positions
  if (any(is.na(unique(mrk.names)))){
      if (verbose) cat("No named markers. Using genome information instead.\n")
      no_name = sum(is.na(mrk.names))
      ##mrk.names[which(is.na(mrk.names))] = paste0("no_name_", seq(1, no_name, 1))
      mrk.names[which(is.na(mrk.names))] = paste0(sequence[which(is.na(mrk.names))],"_", sequence.pos[which(is.na(mrk.names))])
  }
  names(sequence)  = mrk.names
  names(sequence.pos)  = mrk.names
  seq.ref = input.data@fix[,4] # Getting reference alleles
  names(seq.ref)  = mrk.names
  seq.alt = input.data@fix[,5] # Getting alternative alleles
  names(seq.alt)  = mrk.names
  if (verbose) cat("Processing genotypes...")
  cname = which(unlist(strsplit(unique(input.data@gt[,1]), ":")) == "GT")[1] # Defining GT position
  dname = which(unlist(strsplit(unique(input.data@gt[,1]), ":")) == "DP")[1] # Defining DP position
  ## file.ploidy = length(unlist(strsplit(unique(input.data@gt[,2])[1], "/"))) # Checking ploidy (old)
  geno.ploidy = .vcf_get_ploidy(input.data@gt[,-1], cname) # Getting all ploidy levels
  geno.depth = .vcf_get_depth(input.data@gt[,-1], dname) # Getting all depths
  file.ploidy = unique(c(geno.ploidy)) # Getting different ploidy levels
  geno.dose = .vcf_transform_dosage(input.data@gt[,-1], cname) # Accounting for allele dosages
  geno.dose[which(geno.dose == -1)] = NA # Filling NA values
  if (verbose) cat("Done!\n")
  # geno.dose = matrix(unlist(lapply(strsplit(input.data@gt[,-1], ":"), "[", cname)), nrow = n.mrk, byrow = F) # Selecting genotypes (time consuming step)
  # if (sum(grepl("\\|", geno.dose)) > 0){ # Checking phased data
  #   input.phased = TRUE # Treat this in a different way when reading file
  #   warning("Phased genotypes detected. Should MAPpoly consider this information?")
  # } else {input.phased = FALSE}
  if (any((file.ploidy[-which(file.ploidy == -1)] %% 2) != 0)){ # Checking odd ploidy level
    stop("Your VCF file shows an odd ploidy level, but MAPpoly only supports even ploidy levels. Please check your VCF file and try again.")
  }
  if (!is.na(ploidy) && !(ploidy %in% file.ploidy)){ # Checking informed and file ploidy
    stop("Informed ploidy doesn't match any detected ploidy level. Detected ploidy level(s): ", paste0(file.ploidy, ' '))
  }
  if (!is.na(ploidy) && (ploidy > 2) && all(file.ploidy == 2)){ # Checking absence of dosages
    warning("Informed ploidy is ",ploidy, ", but detected ploidy is ", file.ploidy, ".\nIf your species is polyploid, you should provide allelic dosages for all individuals. You can estimate allelic dosages using packages such as 'SuperMASSA', 'updog', 'fitTetra', 'polyRAD' and others. We are working on an integrated function to estimate dosages, which will be available soon.\nUsing ploidy = 2 instead.")
  } # Allow option for building genetic maps for diploid species
    if (!is.na(ploidy)){ # If ploidy is informed and passed previous checks, then use it
        m = ploidy
    } else { # Else, use the first ploidy level detected on file
        m = file.ploidy[1]
    }
    if (verbose) cat("Selected ploidy:", m, "\n")
  if (!(parent.1 %in% ind.names) | !(parent.2 %in% ind.names)){
    stop("Provided parents were not found in VCF file. Please check it and try again.")
  }
  # lost.data = paste(c(rep("./", (m-1)), "."), collapse = "")
  #     for (i in 1:length(geno.dose)){
  #         if (geno.dose[i] == lost.data){
  #           geno.dose[i] = NA
  #       } else {geno.dose[i] = length(which(unlist(strsplit(geno.dose[i], "/")) == 0))} # Transforming dosages (time consuming step)
  #      }
  # geno.dose = matrix(as.numeric(geno.dose), nrow = n.mrk, byrow = F)

    ## Updating some info
    dif_ploidy = which(rowSums(geno.ploidy == m) == (n.ind+2)) # Markers with different ploidy levels
    all_mrk_depth = rowMeans(geno.depth)
    av_depth = which(all_mrk_depth >= min.av.depth) # Markers with average depths below threshold
    max_miss = which(rowSums(is.na(geno.dose))/dim(geno.dose)[2] <= max.missing) # Markers with missing data above the threshold
    selected_markers = intersect(intersect(dif_ploidy,av_depth),max_miss) # Selecting markers that passed all thresholds
    geno.dose = geno.dose[selected_markers,] # Selecting markers
    geno.dose[which(geno.dose < min.gt.depth)] = NA # removing genotypes with depths below the threshold
    all_mrk_depth = all_mrk_depth[selected_markers]
    n.mrk = nrow(geno.dose)
    sequence = sequence[selected_markers]
    sequence.pos = sequence.pos[selected_markers]
    seq.ref = seq.ref[selected_markers]
    seq.alt = seq.alt[selected_markers]
    mrk.names = mrk.names[selected_markers]
    ## sequence = sequence[which(rowSums(geno.ploidy == m) == (n.ind+2))]
    ## sequence.pos = sequence.pos[which(rowSums(geno.ploidy == m) == (n.ind+2))]
    ## seq.ref = seq.ref[which(rowSums(geno.ploidy == m) == (n.ind+2))]
    ## seq.alt = seq.alt[which(rowSums(geno.ploidy == m) == (n.ind+2))]
    ## mrk.names = mrk.names[which(rowSums(geno.ploidy == m) == (n.ind+2))]
    colnames(geno.dose) = ind.names
    rownames(geno.dose) = mrk.names
    dosage.p = as.integer(geno.dose[,which(colnames(geno.dose) == parent.1)]) # Selecting dosages for parent 1
    dosage.q = as.integer(geno.dose[,which(colnames(geno.dose) == parent.2)]) # Selecting dosages for parent 2
    names(dosage.p) <- names(dosage.q) <- mrk.names
    geno.dose = geno.dose[, -c(which(colnames(geno.dose) %in% c(parent.1, parent.2)))] # Updating geno.dose matrix
    ind.names = ind.names[-c(which(ind.names %in% c(parent.1, parent.2)))] # Updating individual names
    geno.dose = data.frame(geno.dose)
  
  ## monomorphic markers
  dp = abs(abs(dosage.p-(m/2))-(m/2))
  dq = abs(abs(dosage.q-(m/2))-(m/2))
  id = dp+dq!=0
  id[which(is.na(id))] = FALSE

  if (verbose){
      cat("Done!\n")
      cat("Read the following data:")
      cat("\n    Ploidy level:", m)
      cat("\n    No. individuals: ", n.ind)
      cat("\n    No. markers: ", n.mrk) 
      cat("\n    No. informative markers:  ", sum(id), " (", round(100*sum(id)/n.mrk,1), "%)", sep = "")
      ## if (all(unique(nphen) != 0))
      ##   cat("\n    This dataset contains phenotypic information.")        
      if (length(sequence) > 1)
          cat("\n    This dataset contains sequence information.")
      cat("\n    ...")
  }

  ## get genotypic info --------------------
  if(nrow(geno.dose)!=length(mrk.names))
    stop("\n\t\t-------------------------------------
         Number of marker names is different from
         the number of markers in the dataset.
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
  if (verbose) cat("\n    Done with reading.\n")
    
    res = structure(list(m = m,
                         n.ind = n.ind,
                         n.mrk = sum(id),
                         ind.names = ind.names,
                         mrk.names = mrk.names[id],
                         dosage.p = dosage.p[id],
                         dosage.q = dosage.q[id],
                         sequence = sequence[id],
                         sequence.pos = sequence.pos[id],
                         seq.ref = seq.ref[id],
                         seq.alt = seq.alt[id],
                         prob.thres = NULL,
                         geno.dose = geno.dose[id,],
                         nphen = 0,
                         phen = NULL,
                         all.mrk.depth = all_mrk_depth[id],
                         kept = NULL,
                         elim.correspondence = NULL
                         ),
                    class = "mappoly.data")
  
  if(filter.non.conforming){
    if (verbose) cat("    Filtering non-conforming markers.\n    ...")
    res<-filter_non_conforming_classes(res)
    if (verbose) cat("\n    Performing chi-square test.\n    ...")
    ##Computing chi-square p.values
    Ds <- array(NA, dim = c(m+1, m+1, m+1))
    for(i in 0:m)
      for(j in 0:m)
        Ds[i+1,j+1,] <- segreg_poly(m = m, dP = i, dQ = j)
    Dpop<-cbind(res$dosage.p, res$dosage.q)
    M<-t(apply(Dpop, 1, function(x) Ds[x[1]+1, x[2]+1,]))
    dimnames(M)<-list(res$mrk.names, c(0:m))
    M<-cbind(M, res$geno.dose)
    res$chisq.pval<-apply(M, 1, mrk_chisq_test, m = m)
    if (verbose) cat("\n    Done.\n")
  }
    if (elim.redundant){
      seqred = make_seq_mappoly(res, arg = 'all', data.name = res)
      redun = elim_redundant(seqred, data = res)
      if (nrow(redun$elim.correspondence) < 1) return(res)
      res$kept = redun$kept
      res$elim.correspondence = redun$elim.correspondence
      mrks.rem = match(res$elim.correspondence$elim, res$mrk.names)
      res$elim.correspondence$sequence = res$sequence[c(mrks.rem)]
      res$elim.correspondence$sequence.pos = res$sequence.pos[c(mrks.rem)]
      res$elim.correspondence$seq.ref = res$seq.ref[c(mrks.rem)]
      res$elim.correspondence$seq.alt = res$seq.alt[c(mrks.rem)]
      res$elim.correspondence$all.mrk.depth = res$all.mrk.depth[c(mrks.rem)]
      res$n.mrk = length(res$kept)
      res$mrk.names = res$mrk.names[-c(mrks.rem)]
      res$geno.dose = res$geno.dose[-c(mrks.rem),]
      res$dosage.p = res$dosage.p[-c(mrks.rem)]
      res$dosage.q = res$dosage.q[-c(mrks.rem)]
      res$sequence = res$sequence[-c(mrks.rem)]
      res$sequence.pos = res$sequence.pos[-c(mrks.rem)]
      res$seq.ref = res$seq.ref[-c(mrks.rem)]
      res$seq.alt = res$seq.alt[-c(mrks.rem)]
      res$all.mrk.depth = res$all.mrk.depth[-c(mrks.rem)]
      res$chisq.pval = res$chisq.pval[-c(mrks.rem)]
    }
  return(res)
}
