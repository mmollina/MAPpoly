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
#' @param read.geno.prob If genotypic probabilities are available (PL field),
#' generates a probability-based dataframe (default = \code{FALSE}).
#'
#' @param prob.thres probability threshold to associate a marker call to a 
#'     dosage. Markers with maximum genotype probability smaller than \code{prob.thres} 
#'     are considered as missing data for the dosage calling purposes (default = 0.95)
#' 
#' @return An object of class \code{mappoly.data} which contains a
#'     list with the following components:
#'     \item{ploidy}{ploidy level}
#'     \item{n.ind}{number individuals}
#'     \item{n.mrk}{total number of markers}
#'     \item{ind.names}{the names of the individuals}
#'     \item{mrk.names}{the names of the markers}
#'     \item{dosage.p1}{a vector containing the dosage in
#'       parent P for all \code{n.mrk} markers}
#'     \item{dosage.p2}{a vector containing the dosage in
#'       parent Q for all \code{n.mrk} markers}
#'     \item{chrom}{a vector indicating which sequence each marker
#'       belongs. Zero indicates that the marker was not assigned to any
#'       sequence}
#'     \item{genome.pos}{Physical position of the markers into the
#'       sequence}
#'     \item{seq.ref}{Reference base used for each marker (i.e. A, T, C, G)}
#'     \item{seq.alt}{Alternative base used for each marker (i.e. A, T, C, G)}
#'     \item{prob.thres}{(unused field)}
#'     \item{geno.dose}{a matrix containing the dosage for each markers (rows) 
#'       for each individual (columns). Missing data are represented by 
#'       \code{ploidy_level + 1}}
#'     \item{geno}{a dataframe containing all genotypic probabilities columns for each
#'       marker and individual combination (rows). Missing data are represented by 
#'       \code{ploidy_level + 1}}
#'     \item{nphen}{(unused field)}
#'     \item{phen}{(unused field)}
#'     \item{all.mrk.depth}{DP information for all markers on VCF file}
#'     \item{chisq.pval}{a vector containing p-values related to the chi-squared 
#'     test of Mendelian segregation performed for all markers}
#'     \item{kept}{if elim.redundant = TRUE, holds all non-redundant markers}
#'     \item{elim.correspondence}{if elim.redundant = TRUE, holds all non-redundant markers and
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
#' @author Gabriel Gesteira, \email{gdesiqu@ncsu.edu}
#' 
#' @references
#'     Mollinari M., Olukolu B. A.,  Pereira G. da S., 
#'     Khan A., Gemenet D., Yencho G. C., Zeng Z-B. (2020), 
#'     Unraveling the Hexaploid Sweetpotato Inheritance 
#'     Using Ultra-Dense Multilocus Mapping, 
#'     _G3: Genes, Genomes, Genetics_. 
#'     \doi{10.1534/g3.119.400620} 
#'     
#'     Mollinari, M., and Garcia, A.  A. F. (2019) Linkage
#'     analysis and haplotype phasing in experimental autopolyploid
#'     populations with high ploidy level using hidden Markov
#'     models, _G3: Genes, Genomes, Genetics_. 
#'     \doi{10.1534/g3.119.400378} 
#'
#' @export read_vcf
#' 
read_vcf = function(file.in, parent.1, parent.2, ploidy = NA, 
                    filter.non.conforming = TRUE, thresh.line = 0.05, 
                    min.gt.depth = 0, min.av.depth = 0, max.missing = 1, 
                    elim.redundant = TRUE, verbose = TRUE, read.geno.prob = FALSE, prob.thres = 0.95) {
  # Checking even ploidy
  if(!is.na(ploidy) && (ploidy %% 2) != 0){
    stop("MAPpoly only supports even ploidy levels. Please check the 'ploidy' parameter and try again.")
  }
  input.file = normalizePath(file.in) # Getting full input file path
  input.size = file.size(file.in) / 1024000 # Getting input file size in MB
  if (verbose && input.size > 3000){
    warning("Your VCF file is greater than 3 GB. Check for available RAM memory.")
  }
  if (verbose) cat("Reading data...\n")
  input.data = vcfR::read.vcfR(file.in, verbose = verbose) # Reading file
  ind.names = colnames(input.data@gt)[-1]
  n.mrk = dim(input.data@gt)[1] # Getting number of markers
  n.ind = length(ind.names) - 2 # Number of individuals excepting two parents
  mrk.names = mrk.names.all = input.data@fix[,3] # Getting marker names
  chrom = input.data@fix[,1] # Getting chromosome information
  genome.pos = as.numeric(input.data@fix[,2]) # Getting positions
  if (any(is.na(unique(mrk.names)))){
    if (verbose) cat("No named markers. Using genome information instead.\n")
    no_name = sum(is.na(mrk.names))
    mrk.names[which(is.na(mrk.names))] = paste0(chrom[which(is.na(mrk.names))],"_", genome.pos[which(is.na(mrk.names))])
  }
  names(chrom)  = mrk.names
  names(genome.pos)  = mrk.names
  seq.ref = input.data@fix[,4] # Getting reference alleles
  names(seq.ref)  = mrk.names
  seq.alt = input.data@fix[,5] # Getting alternative alleles
  names(seq.alt)  = mrk.names
  if (verbose) cat("Processing genotypes...")
  cname = which(unlist(strsplit(unique(input.data@gt[,1]), ":"))  ==  "GT")[1] # Defining GT position
  dname = which(unlist(strsplit(unique(input.data@gt[,1]), ":"))  ==  "DP")[1] # Defining DP position
  geno.ploidy = .vcf_get_ploidy(input.data@gt[,-1], cname) # Getting all ploidy levels
  geno.depth = .vcf_get_depth(input.data@gt[,-1], dname) # Getting all depths
  file.ploidy = unique(c(geno.ploidy)) # Getting different ploidy levels
  geno.dose = .vcf_transform_dosage(input.data@gt[,-1], cname) # Accounting for allele dosages
  ## Checking for PL information
  if (read.geno.prob){
    pname = which(unlist(strsplit(unique(input.data@gt[,1]), ":"))  ==  "PL")[1] # Defining PL position  
    if (isTRUE(pname > 0)){
      if (verbose) cat('\nPL information available. Getting probabilities...')
      geno = .vcf_get_probabilities(input.data@gt[,-1], pname)
      geno = lapply(geno, function(x){rownames(x) = ind.names;return(x)})
      names(geno) = mrk.names
    }
  }  
  geno.dose[which(geno.dose  ==  -1)] = NA # Filling NA values
  if (verbose) cat("Done!\n")
  # geno.dose = matrix(unlist(lapply(strsplit(input.data@gt[,-1], ":"), "[", cname)), nrow = n.mrk, byrow = F) # Selecting genotypes (time consuming step)
  # if (sum(grepl("\\|", geno.dose)) > 0){ # Checking phased data
  #   input.phased = TRUE # Treat this in a different way when reading file
  #   warning("Phased genotypes detected. Should MAPpoly consider this information?")
  # } else {input.phased = FALSE}
  if (any((file.ploidy[-which(file.ploidy  ==  -1)] %% 2) != 0)){ # Checking odd ploidy level
    stop("Your VCF file shows an odd ploidy level, but MAPpoly only supports even ploidy levels. Please check your VCF file and try again.")
  }
  if (!is.na(ploidy) && !(ploidy %in% file.ploidy)){ # Checking informed and file ploidy
    stop("Informed ploidy doesn't match any detected ploidy level. Detected ploidy level(s): ", paste0(file.ploidy, ' '))
  }
  if (!is.na(ploidy) && (ploidy > 2) && all(file.ploidy  ==  2)){ # Checking absence of dosages
    warning("Informed ploidy is ",ploidy, ", but detected ploidy is ", file.ploidy, ".\nIf your species is polyploid, you should provide allelic dosages for all individuals. You can estimate allelic dosages using packages such as 'SuperMASSA', 'updog', 'fitTetra', 'polyRAD' and others. We are working on an integrated function to estimate dosages, which will be available soon.\nUsing ploidy = 2 instead.")
  } # Allow option for building genetic maps for diploid species
  if (!is.na(ploidy)){ # If ploidy is informed and passed previous checks, then use it
    ploidy = ploidy
  } else { # Else, use the first ploidy level detected on file
    ploidy = file.ploidy[1]
  }
  if (verbose) cat("Selected ploidy:", ploidy, "\n")
  if (!(parent.1 %in% ind.names) | !(parent.2 %in% ind.names)){
    stop("Provided parents were not found in VCF file. Please check it and try again.")
  }

  ## Updating some info
  dif_ploidy = which(rowSums(geno.ploidy  ==  ploidy)  ==  (n.ind+2)) # Markers with different ploidy levels
  all_mrk_depth = rowMeans(geno.depth)
  av_depth = which(all_mrk_depth >= min.av.depth) # Markers with average depths below threshold
  max_miss = which(rowSums(is.na(geno.dose))/dim(geno.dose)[2] <= max.missing) # Markers with missing data above the threshold
  ## Filtering non-biallelic markers
  if (exists('geno')) {
    biallelic = unlist(lapply(geno, ncol))
    biallelic = which(biallelic  ==  (ploidy+1))
    selected_markers = intersect(intersect(intersect(dif_ploidy,av_depth),max_miss),biallelic) # Selecting markers that passed all thresholds
  } else {
    selected_markers = intersect(intersect(dif_ploidy,av_depth),max_miss) # Selecting markers that passed all thresholds
  }
  geno.dose = geno.dose[selected_markers,] # Selecting markers
  geno.dose[which(geno.dose < min.gt.depth)] = NA # removing genotypes with depths below the threshold
  all_mrk_depth = all_mrk_depth[selected_markers]
  n.mrk = nrow(geno.dose)
  chrom = chrom[selected_markers]
  genome.pos = genome.pos[selected_markers]
  seq.ref = seq.ref[selected_markers]
  seq.alt = seq.alt[selected_markers]
  mrk.names = mrk.names[selected_markers]
  ## chrom = chrom[which(rowSums(geno.ploidy  ==  ploidy)  ==  (n.ind+2))]
  ## genome.pos = genome.pos[which(rowSums(geno.ploidy  ==  ploidy)  ==  (n.ind+2))]
  ## seq.ref = seq.ref[which(rowSums(geno.ploidy  ==  ploidy)  ==  (n.ind+2))]
  ## seq.alt = seq.alt[which(rowSums(geno.ploidy  ==  ploidy)  ==  (n.ind+2))]
  ## mrk.names = mrk.names[which(rowSums(geno.ploidy  ==  ploidy)  ==  (n.ind+2))]
  colnames(geno.dose) = ind.names
  rownames(geno.dose) = mrk.names
  dosage.p1 = as.integer(geno.dose[,which(colnames(geno.dose)  ==  parent.1)]) # Selecting dosages for parent 1
  dosage.p2 = as.integer(geno.dose[,which(colnames(geno.dose)  ==  parent.2)]) # Selecting dosages for parent 2
  names(dosage.p1) <- names(dosage.p2) <- mrk.names
  geno.dose = geno.dose[, -c(which(colnames(geno.dose) %in% c(parent.1, parent.2)))] # Updating geno.dose matrix
  if (exists('geno')){
    geno2 = geno[c(selected_markers)] # Selecting markers
    geno2 = lapply(geno2, as.data.frame)
    for (i in 1:length(geno2)){
      geno2[[i]] = geno2[[i]][-c(which(rownames(geno2[[i]]) %in% c(parent.1,parent.2))),]
      geno2[[i]] = cbind(names(geno2)[i], rownames(geno2[[i]]), geno2[[i]])
    }
    geno2 = do.call(rbind.data.frame, geno2)
    colnames(geno2) = c('mrk','ind', as.character(seq(0,(ploidy),1)))
    geno2 = geno2[order(geno2$ind),]
    rownames(geno2) = seq(1,nrow(geno2),1)
    geno = geno2
  }
  ind.names = ind.names[-c(which(ind.names %in% c(parent.1, parent.2)))] # Updating individual names
  geno.dose = data.frame(geno.dose)
  
  ## monomorphic markers
  dp = abs(abs(dosage.p1-(ploidy/2))-(ploidy/2))
  dq = abs(abs(dosage.p2-(ploidy/2))-(ploidy/2))
  id = dp+dq != 0
  id[which(is.na(id))] = FALSE

  if (verbose){
    cat("Done!\n")
    cat("Read the following data:")
    cat("\n    Ploidy level:", ploidy)
    cat("\n    No. individuals: ", n.ind)
    cat("\n    No. markers: ", n.mrk) 
    cat("\n    No. informative markers:  ", sum(id), " (", round(100*sum(id)/n.mrk,1), "%)", sep = "")
    ## if (all(unique(nphen) != 0))
    ##   cat("\n    This dataset contains phenotypic information.")        
    if (length(chrom) > 1)
      cat("\n    This dataset contains chromosome information.")
    cat("\n    ...")
  }

  ## get genotypic info --------------------
  if(nrow(geno.dose) != length(mrk.names))
    stop("\n\t\t-------------------------------------
         Number of marker names is different from
         the number of markers in the dataset.
         Please, check data.
         ------------------------------------------\n")
  if(ncol(geno.dose) != length(ind.names))
    stop("\n\t\t-------------------------------------
         Number of individual names is different from
         the number of individuals in the dataset.
         Please, check data.
         ------------------------------------------\n")
  dimnames(geno.dose) <- list(mrk.names, ind.names)
  geno.dose[is.na(geno.dose)] <- ploidy + 1
  ## returning the 'mappoly.data' object
  if (verbose) cat("\n    Done with reading.\n")

  if (exists('geno')){
    res = structure(list(ploidy = ploidy,
                       n.ind = n.ind,
                       n.mrk = sum(id),
                       ind.names = ind.names,
                       mrk.names = mrk.names[id],
                       dosage.p1 = dosage.p1[id],
                       dosage.p2 = dosage.p2[id],
                       chrom = chrom[id],
                       genome.pos = genome.pos[id],
                       seq.ref = seq.ref[id],
                       seq.alt = seq.alt[id],
                       prob.thres = prob.thres,
                       geno = subset(geno, geno$mrk%in%mrk.names[id]),
                       geno.dose = geno.dose[id,],
                       nphen = 0,
                       phen = NULL,
                       all.mrk.depth = all_mrk_depth[id],
                       chisq.pval = NULL,
                       kept = NULL,
                       elim.correspondence = NULL
                       ),
                  class = "mappoly.data")
  } else {
    res = structure(list(ploidy = ploidy,
                       n.ind = n.ind,
                       n.mrk = sum(id),
                       ind.names = ind.names,
                       mrk.names = mrk.names[id],
                       dosage.p1 = dosage.p1[id],
                       dosage.p2 = dosage.p2[id],
                       chrom = chrom[id],
                       genome.pos = genome.pos[id],
                       seq.ref = seq.ref[id],
                       seq.alt = seq.alt[id],
                       prob.thres = NULL,
                       geno.dose = geno.dose[id,],
                       nphen = 0,
                       phen = NULL,
                       all.mrk.depth = all_mrk_depth[id],
                       chisq.pval = NULL,
                       kept = NULL,
                       elim.correspondence = NULL
                       ),
                  class = "mappoly.data")
  }
    
  if(filter.non.conforming){
    if (verbose) cat("    Filtering non-conforming markers.\n    ...")
    res <- filter_non_conforming_classes(res)
    if (verbose) cat("\n    Performing chi-square test.\n    ...")
    ##Computing chi-square p.values
    Ds <- array(NA, dim = c(ploidy+1, ploidy+1, ploidy+1))
    for(i in 0:ploidy)
      for(j in 0:ploidy)
        Ds[i+1,j+1,] <- segreg_poly(ploidy = ploidy, dP = i, dQ = j)
    Dpop <- cbind(res$dosage.p1, res$dosage.p2)
    M <- t(apply(Dpop, 1, function(x) Ds[x[1]+1, x[2]+1,]))
    dimnames(M) <- list(res$mrk.names, c(0:ploidy))
    M <- cbind(M, res$geno.dose)
    res$chisq.pval <- apply(M, 1, mrk_chisq_test, ploidy = ploidy)
    if (verbose) cat("\n    Done.\n")
  }
  if (elim.redundant){
    seqred = make_seq_mappoly(res, arg = 'all', data.name = res)
    redun = elim_redundant(seqred, data = res)
    if (nrow(redun$elim.correspondence) < 1) return(res)
    res$kept = redun$kept
    res$elim.correspondence = redun$elim.correspondence
    mrks.rem = match(res$elim.correspondence$elim, res$mrk.names)
    res$elim.correspondence$chrom = res$chrom[c(mrks.rem)]
    res$elim.correspondence$genome.pos = res$genome.pos[c(mrks.rem)]
    res$elim.correspondence$seq.ref = res$seq.ref[c(mrks.rem)]
    res$elim.correspondence$seq.alt = res$seq.alt[c(mrks.rem)]
    res$elim.correspondence$all.mrk.depth = res$all.mrk.depth[c(mrks.rem)]
    res$n.mrk = length(res$kept)
    res$mrk.names = res$mrk.names[-c(mrks.rem)]
    if (exists('geno')) res$geno = subset(res$geno, res$geno$mrk%in%res$mrk.names) ## res$geno = res$geno[-c(mrks.rem)]
    res$geno.dose = res$geno.dose[-c(mrks.rem),]
    res$dosage.p1 = res$dosage.p1[-c(mrks.rem)]
    res$dosage.p2 = res$dosage.p2[-c(mrks.rem)]
    res$chrom = res$chrom[-c(mrks.rem)]
    res$genome.pos = res$genome.pos[-c(mrks.rem)]
    res$seq.ref = res$seq.ref[-c(mrks.rem)]
    res$seq.alt = res$seq.alt[-c(mrks.rem)]
    res$all.mrk.depth = res$all.mrk.depth[-c(mrks.rem)]
    res$chisq.pval = res$chisq.pval[-c(mrks.rem)]
  }
  return(res)
}
