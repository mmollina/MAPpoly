m<-4
ch<-1
probs<-genoprob.err[[1]]




pref_pair<-function(m, ch, probs)
{
  Gnames<-dimnames(probs$probs)[[1]]
  x<-dim(probs$probs)
  nsta<-x[1]
  npos<-x[2]
  nind<-x[3]
  ##Parent P
  PsiP<-NULL
  a<-combn(letters[1:m],m/2)
  for(i in 1:ncol(a))
  {
    b<-perm_tot(setdiff(letters[1:m], a[,i]))
    for(j in 1:nrow(b))
    {
      PsiP<-rbind(PsiP, sort(apply(cbind(a[,i],b[j,]), 1, function(x) paste(sort(x), collapse = ""))))
    }
  }
  PsiP<-unique(PsiP)
  nbiv<-nrow(PsiP)
  PsiP_given_GP<-matrix(0, length(Gnames), nrow(PsiP), dimnames = list(Gnames, apply(PsiP, 1, paste, collapse="/")))
  GP<-unlist(strsplit(substr(Gnames, 1, m/2), ""))
  dim(GP)<-c(m/2,length(Gnames))
  for(i in 1:nrow(PsiP_given_GP))
  {
    x<-apply(PsiP, 1, function(x,y) all(grepl(paste(y, collapse="|"), x)), y = GP[,i])
    PsiP_given_GP[i,x]<-1/sum(x)
  }
  A<-array(NA, dim = c(nbiv,npos,nind))
  for(i in 1:nind)
    A[,,i]<-crossprod(PsiP_given_GP,probs$probs[,,i])
  resfinalP<-apply(A, MARGIN = c(1,2), mean)
  ##Parent Q
  PsiQ<-NULL
  a<-combn(letters[(1+m):(2*m)],m/2)
  for(i in 1:ncol(a))
  {
    b<-perm_tot(setdiff(letters[(1+m):(2*m)], a[,i]))
    for(j in 1:nrow(b))
    {
      PsiQ<-rbind(PsiQ, sort(apply(cbind(a[,i],b[j,]), 1, function(x) paste(sort(x), collapse = ""))))
    }
  }
  PsiQ<-unique(PsiQ)
  PsiQ_given_GQ<-matrix(0, length(Gnames), nrow(PsiQ), dimnames = list(Gnames, apply(PsiQ, 1, paste, collapse="/")))
  GQ<-unlist(strsplit(substr(Gnames, 2+(m/2), m+1), ""))
  dim(GQ)<-c(m/2,length(Gnames))
  for(i in 1:nrow(PsiQ_given_GQ))
  {
    x<-apply(PsiQ, 1, function(x,y) all(grepl(paste(y, collapse="|"), x)), y = GQ[,i])
    PsiQ_given_GQ[i,x]<-1/sum(x)
  }
  A<-array(NA, dim = c(nbiv,npos,nind))
  for(i in 1:nind)
    A[,,i]<-crossprod(PsiQ_given_GQ,probs$probs[,,i])
  resfinalQ<-apply(A, MARGIN = c(1,2), mean)
  dimnames(resfinalP)<-list(colnames(PsiP_given_GP), names(probs$map))
  dimnames(resfinalQ)<-list(colnames(PsiQ_given_GQ), names(probs$map))
  list(P = resfinalP, Q = resfinalQ)
}

DF2<-DF<-NULL
for(ch in 1:12)
{
  cat("\n~~~~~~~~~ ch: ", ch, "~~~~~~~~~~\n")
  pp<-pref_pair(m = 4, ch = ch, probs = genoprob.err[[ch]])
  names(pp)<-c("P1", "P2")
  tt<-as.data.frame(sapply(pp, function(y) apply(2*154*y, 2, function(x) chisq.test(x)$p.value)))
  tt$ch<-ch
  tt$map<-genoprob.err[[ch]]$map
  temp<-reshape2::melt(pp)
  head(temp)
  temp$map<-rep(genoprob.err[[ch]]$map, each = 3)
  temp$ch <- ch
  DF<-rbind(DF,temp)
  DF2<-rbind(DF2,tt)
}

DF3<-data.table::melt(data = DF2, measure.vars = c("P1", "P2"))
head(DF); head(DF2); head(DF3)

#DF<-subset(DF, ch == 2)
require(ggplot2)
ggplot(DF) + 
  geom_smooth(aes(map, value, colour = Var1), size = .5, se = FALSE) + 
  facet_grid(L1~ch, scales = "free_x", space = "free_x") +
  geom_hline(yintercept = 1/3, linetype="dashed") +
  ylim(0,1) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab("Probability") + xlab("Distance (cM)") 

ggplot(DF3, aes(map, -log10(value), colour = variable)) +
  geom_point(alpha = .7, size = 1) +  facet_grid(.~ch, scales = "free_x", space = "free_x") + 
  geom_hline(yintercept = -log10(0.01), linetype="dashed") +  scale_color_manual(values=c("#E69F00","#56B4E9"), name = "Parents") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="none") +
  ylab("-log_10(P)") + xlab("Distance (cM)") 


ggplotly(p)

