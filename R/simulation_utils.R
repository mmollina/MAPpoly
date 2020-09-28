#' Simulate mapping population (one parent)
#'
#' This function simulates a polyploid mapping population
#' under random chromosome segregation
#' with one informative parent. This function is not to be
#' directly called by the user
#'
#' @param void internal function to be documented
#' @keywords internal
sim_cross_one_informative_parent<-function(m,
                                           n.mrk,
                                           rf.vec,
                                           hom.allele,
                                           n.ind,
                                           seed = NULL,
                                           prob = NULL){
    if(!is.null(seed)) set.seed(seed)
    if(length(rf.vec)==1) rf.vec<-rep(rf.vec, n.mrk-1)
    res<-matrix(NA,n.ind,n.mrk)
    rf.res <- numeric(n.mrk-1)

    ## Listing all possible bivalent configurations
    a<-perm_tot(1:m)
    bv.conf<-vector("list", nrow(a))
    for(i in 1:nrow(a))
    {
      temp <- apply(matrix(a[i,], 2, m/2), 2, sort)
      bv.conf[[i]] <- temp[,order(temp[1,])]
    }
    bv.conf <- unique(bv.conf)
    names(bv.conf) <- sapply(bv.conf, function(x) paste(apply(x, 2,
                                                  function(x)
                                                    paste0("[", paste0(x, collapse = ""), "]", collapse="")),
                                            collapse = ""))
    if(is.null(prob))
      prob<-rep(1/length(bv.conf), length(bv.conf))

    for(k in 1:n.ind){              #for each individual
        gen.1<-matrix(1:m,m,n.mrk)  #simulates the chromosomes multiallelic markers in 'n.mrk' positions
        id <- sample(x = 1:length(bv.conf), size = 1, prob = prob) #sampling one bivalent configuration based on given probabilities
        choosed_biv <- bv.conf[[id]]
        choosed_biv <- choosed_biv[,sample(1:(m/2))]
        for(i in 1:ncol(choosed_biv))
        {
          choosed_biv[,i]<-sample(choosed_biv[,i])
        }
        pole.1<-choosed_biv[1,]
        pole.2<-choosed_biv[2,]
        set.2<-gen.1[pole.1,]      #allocating the chromosomes on the variables set.1 and set.2, thus (set.1[i], set.2[i]) represents a bivalent
        set.1<-gen.1[pole.2,]
        for(i in 1:(m/2)){         #for each one of the m/2 chromosome pair (bivalents)
            a<-set.1[i,]
            b<-set.2[i,]
            for(j in 1:(n.mrk-1)){             #for each adjacent interval between.mrkkers
                if(runif(1)  < rf.vec[j]){       #if a random number drawn from the interval [0,1] (according a uniform distribution)
                                        #is less than the recombination fraction for that interval
                    which.swap<-c((j+1):n.mrk)     #the alleles for that interval and bivalent are swapped
                    temp<-a[which.swap]
                    a[which.swap]<-b[which.swap]
                    b[which.swap]<-temp
                }
            }             #this completes the whole bivalent
            set.1[i,]<-a  #attributing the resulting vector to the initial variables
            set.2[i,]<-b
        }               #for all bivalents
        if(sample(0:1,1)) gam <- set.1 #sample one of the meiotic products
        else gam<-set.2
        for(i in 1:(m/2)){ #counting the recombinant chromosomes in their multiallelic form
            for(j in 2:ncol(gam)){
                if(!gam[i,j]==gam[i,j-1])
                    rf.res[j-1]<-rf.res[j-1]+1
            }
        }
        for(i in 1:n.mrk)
            gam[,i]<-as.numeric(!is.na(match(gam[,i], hom.allele[[i]])))
        res[k,]<-apply(gam, 2, sum)
    }
    rf.calc<-rf.res/(n.ind*m/2)  #computing the recombination fraction
    dimnames(res)<-list(paste("Ind",1:n.ind, sep="_"), paste("M",1:n.mrk, sep="_"))
    list(data.sim.one.parent=res, cross.count.one.parent=rf.res, rf.calc.one.parent=rf.calc)
}

#' Simulate mapping population (tow parents)
#'
#' @param void internal function to be documented
#' @keywords internal
sim_cross_two_informative_parents<-function(m,
                                           n.mrk,
                                           rf.vec,
                                           n.ind,
                                           hom.allele.p,
                                           hom.allele.q,
                                           prob.P = NULL,
                                           prob.Q = NULL,
                                           seed = NULL){
    if(!is.null(seed)) set.seed(seed)
    dose.p<-unlist(lapply(hom.allele.p, function(x) sum(as.logical(x))))
    dose.q<-unlist(lapply(hom.allele.q, function(x) sum(as.logical(x))))
    if(any(apply(rbind(dose.p,dose.q),2,sum)==0)) stop("Found zero doses in both parents at the same marker")
    if(!is.null(seed)) set.seed(seed)
    ##Parent 1
    data.P<-sim_cross_one_informative_parent(m=m, n.mrk=n.mrk, rf.vec=rf.vec,
                           hom.allele=hom.allele.p, n.ind=n.ind, prob = prob.P)
    ##Parent 2
    data.Q<-sim_cross_one_informative_parent(m=m, n.mrk=n.mrk, rf.vec=rf.vec,
                           hom.allele=hom.allele.q, n.ind=n.ind, prob = prob.Q)
    rf.calc<-(data.P$cross.count.one.parent + data.Q$cross.count.one.parent)/(n.ind*m)
    data.sim.two.parents<-data.P$data.sim.one.parent + data.Q$data.sim.one.parent
    list(data.sim.two.parents=data.sim.two.parents, rf.calc=rf.calc)
}

#' Draw simple parental linkage phase configurations
#'
#' This function draws the parental map (including the linkage
#' phase configuration) in a pdf output. This function is not to
#' be directly called by the user
#'
#' @param void internal function to be documented
#' @importFrom grDevices pdf dev.off
#' @keywords internal
draw_cross<-function(m,rf.vec=NULL,hom.allele.p,hom.allele.q, file=NULL, width=12, height=6){
    if(!is.null(file))
        pdf(file, width=width, height=height)
    oldpar <- par(xaxt="n")
    on.exit(par(oldpar))
    plot(c(0,22), c(0,-(m+10)), type="n", axes=FALSE, xlab="Partental homology groups", main=paste("Ploidy: ", m), ylab="")
    for(i in -(1:m)){
        lines(c(0,10), c(i,i))
        lines(c(12,22), c(i,i))
    }
    pos.p<-cumsum(c(0,rf.vec/sum(rf.vec)))*10
    for(i in 1:length(hom.allele.p)){
        abline(v=pos.p[i], lty=2, lwd=.5)
        text(pos.p[i], 0, i, cex=.7)
        if(any(hom.allele.p[[i]]!=0))
            points(x=rep(pos.p[i],length(hom.allele.p[[i]])), y=-hom.allele.p[[i]], col=2, pch=20, cex=2)
        points(pos.p[i] , -(m+10), pch="|")
        text(pos.p[i]+diff(pos.p)[i]/2, -(m+10)+.8, labels=rf.vec[i], srt=90)
    }
    pos.q<-pos.p+12
    for(i in 1:length(hom.allele.q)){
        abline(v=pos.q[i], lty=2, lwd=.5)
        text(pos.q[i], 0, i, cex=.7)
        if(any(hom.allele.q[[i]]!=0))
            points(x=rep(pos.q[i],length(hom.allele.q[[i]])), y=-hom.allele.q[[i]], col=2, pch=20, cex=2)
        points(pos.q[i] , -(m+10), pch="|")
        text(pos.q[i]+diff(pos.q)[i]/2, -(m+10)+.8, labels=rf.vec[i], srt=90)
    }
    text(x=11,y=-(m+1)/2,labels="X", cex=2)
    lines(c(0,10), c(-(m+10),-(m+10)))
    lines(c(12,22), c(-(m+10),-(m+10)))
    if(!is.null(file)) dev.off()
}
