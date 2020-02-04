#' Plot marker information
#'
#' Plots summary statistics for a given marker
#' 
#' @param input.data an object of class \code{mappoly.data}
#' 
#' @param mrk marker name or position in the dataset
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
#' @export
#' @importFrom graphics barplot layout mtext legend 
#' @importFrom stats chisq.test
plot_mrk_info<-function(input.data, mrk)
  {
  input_classes <- c("mappoly.data")
  if (!inherits(input.data, input_classes)) {
    stop(deparse(substitute(input.data)), " is not an object of class 'mappoly.data'")
  }
  if(is.numeric(mrk))
  {
    mrk <- mrk[1]
    if(mrk > input.data$n.mrk)
      stop(mrk, " exceeds the number of markers in the dataset")
    mrk<-input.data$mrk.names[mrk]
  }
  
  if(!mrk%in%input.data$mrk.names)
    stop(deparse(substitute(mrk)), " is not presnet in ", deparse(substitute(input.data)), " dataset")
  ## Parents dosage
  dp<-input.data$dosage.p[input.data$mrk.names==mrk]
  dq<-input.data$dosage.q[input.data$mrk.names==mrk]
  suppressWarnings({
    ## if no probabilities
    if(nrow(input.data$geno)==input.data$n.mrk){
      layout(matrix(c(1,2), ncol = 2), widths = c(1, 2))
      ## Genotype frequencies
      x<-input.data$geno.dose[mrk, ]
      x[x==input.data$m+1]<-NA
      x<-table(as.numeric(x), useNA = "always")
      names(x)<-c(names(x)[-length(x)], "NA") 
      ## empty plot area
      op<-par(mar = c(2,2,5,2))
      plot(0:100, type = "n", axes = FALSE, xlab = "", ylab = "")
      mtext(side = 3, text = mrk, adj = 0, cex = 1.2, font = 3)
      text(x = 0, y = 90 , labels = paste0("marker #: ", which(input.data$mrk.names==mrk)), adj = 0)
      ## parents dosage
      text(x = 0, y = 80, paste0("Dose in P1: ", dp), adj = 0)
      text(x = 0, y = 70, paste0("Dose in P2: ", dq), adj = 0)
      ## missing data 
      text(x = 0, y = 60 , labels = paste0("Missing: ", round(100*tail(x,1)/sum(x),1), "%"), adj = 0)
      ## p.value for chi-square test under Mendelian segregation
      seg.exp <- segreg_poly(m = input.data$m, 
                             dP = dp, 
                             dQ = dq)
      names(seg.exp) <- 0:input.data$m
      seg.exp <- seg.exp[seg.exp!=0]
      seg.obs <- seg.exp
      seg.obs[names(x)[-length(x)]]<-x[-length(x)]
      pval <- chisq.test(x = seg.obs, p = seg.exp[names(seg.obs)])$p.value
      text(x = 0, y = 50 , labels = paste0("p-value: ", formatC(pval, format = "e", digits = 2)), adj = 0)
      ##genomic position and sequence
      text(x = 0, y = 40 , labels = paste0("sequence: ", input.data$sequence[input.data$mrk.names==mrk]), adj = 0)
      text(x = 0, y = 30 , labels = paste0("seq. position: ", input.data$sequence.pos[input.data$mrk.names==mrk]), adj = 0)
      par(op)
      barplot(x, col = c(gg_color_hue(input.data$m + 1)[1:(length(x)-1)], "#404040"))
    } else if(nrow(input.data$geno)!=input.data$n.mrk){ ## if genotype probabilities are available
      layout(matrix(c(1,2,3,3), ncol = 2, nrow = 2), widths = c(1, 2))
      ## Genotype frequencies
      x<-input.data$geno.dose[mrk, ]
      x[x==input.data$m+1]<-NA
      x<-table(x, useNA = "always")
      names(x)<-c(names(x)[-length(x)], "NA") 
      ## empty plot area
      op<-par(mar = c(2,2,5,2))
      plot(0:100, type = "n", axes = FALSE, xlab = "", ylab = "")
      mtext(side = 3, text = mrk, adj = 0, cex = 1.2, font = 3)
      text(x = 0, y = 90 , labels = paste0("marker #: ", which(input.data$mrk.names==mrk)), adj = 0)
      ## parents dosage
      text(x = 0, y = 80, paste0("Dose in P1: ", 
                                 dp), adj = 0)
      text(x = 0, y = 70, paste0("Dose in P2: ",
                                 dq), adj = 0)
      ## missing data 
      text(x = 0, y = 60 , labels = paste0("Missing: ", round(100*tail(x,1)/sum(x),1), "%"), adj = 0)
      ## p.value for chi-square test under Mendelian segregation
      seg.exp <- segreg_poly(m = input.data$m, 
                             dP = dp, 
                             dQ = dq)
      names(seg.exp) <- 0:input.data$m
      seg.exp <- seg.exp[seg.exp!=0]
      seg.obs <- seg.exp
      seg.obs[names(x)[-length(x)]]<-x[-length(x)]
      pval <- chisq.test(x = seg.obs, p = seg.exp[names(seg.obs)])$p.value
      text(x = 0, y = 50 , labels = paste0("p-value: ", formatC(pval, format = "e", digits = 2)), adj = 0)
      
      ##genomic position and sequence
      text(x = 0, y = 40 , labels = paste0("sequence: ", input.data$sequence[input.data$mrk.names==mrk]), adj = 0)
      text(x = 0, y = 30 , labels = paste0("seq. position: ", input.data$sequence.pos[input.data$mrk.names==mrk]), adj = 0)
      text(x = 0, y = 20 , labels = paste0("prob. threshold: ", input.data$prob.thres), adj = 0)
      par(op)
      pal<-gg_color_hue(input.data$m + 1)
      names(pal)<-0:input.data$m 
      op<-par(mar = c(5,3,0,2), cex = .7)
      barplot(x, col = c(na.omit(pal[names(x)]), "#404040"))
      par(op)
      ## probability distribution of the genotypes
      mrk.name<-mrk
      df<-subset(input.data$geno, 
                 mrk == mrk.name, 
                 select = c(as.character(0:input.data$m), "ind"))
      rownames(df)<-df$ind
      z <- df[rev(do.call(order, as.data.frame(df[,1:(input.data$m+1)]))),]
      z<-z[apply(z[,1:(input.data$m + 1)], 1, function(x) !all(round(x,5)==round(segreg_poly(input.data$m, dp, dq),5))),]
      w<-expand.grid(1:nrow(z), 0:(ncol(z)-2))
      x<-w[,1]
      y<-w[,2]
      pal<-rep(gg_color_hue(input.data$m + 1), each = nrow(z))
      pal[as.numeric(as.matrix(z[,1:(input.data$m + 1)])) < input.data$prob.thres] <- "#404040"
      op<-par(mar = c(2,2,2,2))
      plot3D::scatter3D(x,y, as.matrix(z[,1:(input.data$m + 1)]),  theta = 30, phi = 30, bty = "g",  type = "h", lwd = .3 ,
                        ticktype = "detailed", pch = 19, cex = 0.5, 
                        colvar = NULL, 
                        col = pal,
                        xlab = "Offspring",
                        ylab ="Dose", 
                        zlab = "Genotype probability", 
                        cex.axis = .7, cex.lab = .7, clab = )
      par(op)
      par(mfrow=c(1,1))
    } else stop("Shouldn't get here.")
  })
}