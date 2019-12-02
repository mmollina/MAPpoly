#' Generates all possible linkage phases between two blocks of markers,
#' eliminating equivalent configurations, i.e., configurations with the
#' same likelihood
#'
#' @param block1 submap with markers of the first block
#' @param block2 submap with markers of the second block, or just a single marker identified by its position in the \code{mappoly.data} object
#' @param rf.matrix matrix obtained with the function \code{rf_list_to_matrix}
#' using the parameter \code{shared.alleles = TRUE}
#' @param m ploidy level (i.e. 4, 6 and so on)
#' @param max.inc maximum number of allowed inconsistencies
#' @keywords internal
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu} and Gabriel Gesteira, \email{gabrielgesteira@usp.br}
#' @export generate_all_link_phases_elim_equivalent_haplo
generate_all_link_phases_elim_equivalent_haplo <-
    function(block1, block2, rf.matrix, m, max.inc = NULL) {
        ## Check block2 class (block or single marker)
        if (is.integer(block2)){
            seq.ph = list(P = get(rf.matrix$data.name)$dosage.p[block2],
                          Q = get(rf.matrix$data.name)$dosage.q[block2])
            block2 = list(seq.num = block2, seq.ph = seq.ph)
        }
        
        ## Getting M matrix
        M = list(P = rf.matrix$ShP[as.character(block1$seq.num),as.character(block2$seq.num)], Q = rf.matrix$ShQ[as.character(block1$seq.num),as.character(block2$seq.num)])
        
        ## Parent P: all permutations between blocks
        hP1 <- ph_list_to_matrix(L = block1$seq.ph$P, m = m)
        p1 <- apply(hP1, 2, paste, collapse = "")
        hP2 <- ph_list_to_matrix(L = block2$seq.ph$P, m = m)
        p2 <- apply(hP2, 2, paste, collapse = "")
        dimnames(hP1) <- list(block1$seq.num, p1)
        dimnames(hP2) <- list(block2$seq.num, p2)
        p2 <- perm_tot(p2)
        
        ## Parent Q: all permutations between blocks
        hQ1 <- ph_list_to_matrix(L = block1$seq.ph$Q, m = m)
        q1 <- apply(hQ1, 2, paste, collapse = "")
        hQ2 <- ph_list_to_matrix(L = block2$seq.ph$Q, m = m)
        q2 <- apply(hQ2, 2, paste, collapse = "")
        dimnames(hQ1) <- list(block1$seq.num, q1)
        dimnames(hQ2) <- list(block2$seq.num, q2)
        q2 <- perm_tot(q2)

        ## WP: removing redundancy and accounting for shared elleles
        wp <- NULL
        for (i in 1:nrow(p2))
            wp <- rbind(wp, paste(p1, p2[i, ], sep = "-"))
        wp.n <- unique(t(apply(unique(wp), 1, sort)))
        yp <- apply(t(apply(wp, 1, function(x) sort(x))), 1, paste, collapse = "|")
        yp.n <- apply(wp.n, 1, function(x) paste(sort(x), collapse = "|"))
        wp <- wp[match(yp.n, yp), , drop = FALSE]
        ct <- numeric(nrow(wp))
        for (i in 1:nrow(wp)){
            a = matrix(unlist(strsplit(wp[i, ], "-")), ncol = 2, byrow = TRUE)
            sharedP = tcrossprod(hP1[,a[,1]], t(hP2[,a[,2]]))            
            ct[i] = sum((M$P != sharedP), na.rm = T)
        }
        ## Checking inconsistency
        if(is.null(max.inc)){
            id <- rep(TRUE, length(ct))      
        } else
            id <- ct <= max.inc
        ## Maximum inconsistency
        if (sum(id) == 0)
            id <- which.min(ct)
        wp <- matrix(wp[id, ], ncol = m)

        ## WQ: removing redundancy and accounting for shared elleles
        wq <- NULL
        for (i in 1:nrow(q2))
            wq <- rbind(wq, paste(q1, q2[i, ], sep = "-"))
        wq.n <- unique(t(apply(unique(wq), 1, sort)))
        yq <- apply(t(apply(wq, 1, function(x) sort(x))), 1, paste, collapse = "|")
        yq.n <- apply(wq.n, 1, function(x) paste(sort(x), collapse = "|"))
        wq <- wq[match(yq.n, yq), , drop = FALSE]
        ct <- numeric(nrow(wq))
        for (i in 1:nrow(wq)){
            a = matrix(unlist(strsplit(wq[i, ], "-")), ncol = 2, byrow = TRUE)
            sharedQ = tcrossprod(hQ1[,a[,1]], t(hQ2[,a[,2]]))
            ct[i] = sum((M$Q != sharedQ), na.rm = T)
        }
        ## Checking inconsistency
        if(is.null(max.inc)){
            id <- rep(TRUE, length(ct))      
        } else
            id <- ct <= max.inc
        ## Maximum inconsistency
        if (sum(id) == 0)
            id <- which.min(ct)
        wq <- matrix(wq[id, ], ncol = m)
        
        ## Re-arranging phases
        phase.to.test <- vector("list", nrow(wp) * nrow(wq))
        cte <- 1
        for(i in 1:nrow(wp)){
            for(j in 1:nrow(wq)){
                phase.to.test[[cte]] <- list(P = ph_matrix_to_list(hP2[,sapply(strsplit(wp[i,], "-"), function(x) x[2])]),
                                             Q = ph_matrix_to_list(hQ2[,sapply(strsplit(wq[j,], "-"), function(x) x[2])]))
                cte <- cte + 1
            }
        }
        phase.to.test <- unique(phase.to.test)
        names(phase.to.test) <- paste0("config.", 1:length(phase.to.test))
        return(phase.to.test)
    }
