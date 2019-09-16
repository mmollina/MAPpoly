require(mappoly)
require(ggplot2)
require(plotly)

## hexaploid example
#input.probs<-calc_genoprob(hexafake[[1]])

## solcap example
load("~/repos/MAPpoly_vignettes/vignette_1/maps.rda")

input.probs<-genoprob.err[[1]]

homolog_prob<-function(input.probs){
  stt.names<-dimnames(input.probs$probs)[[1]] ## state names
  mrk.names<-dimnames(input.probs$probs)[[2]] ## mrk names
  ind.names<-dimnames(input.probs$probs)[[3]] ## individual names
  v<-c(2,4,6,8,10,12)
  names(v)<-choose(v,v/2)^2
  m<-v[as.character(length(stt.names))]
  hom.prob<-array(NA, dim = c(m*2, length(mrk.names), length(ind.names)))
  dimnames(hom.prob)<-list(letters[1:(2*m)], mrk.names, ind.names)
  for(i in letters[1:(2*m)])
    hom.prob[i,,] <- apply(input.probs$probs[grep(stt.names, pattern = i),,], c(2,3), function(x) round(sum(x, na.rm = TRUE),4))
  df.hom<-reshape2::melt(hom.prob)
  map<-data.frame(map.position = input.probs$map, marker = names(input.probs$map))
  colnames(df.hom)<-c("homolog", "marker", "individual", "probability")
  df.hom<-merge(df.hom, map, sort = FALSE)
  #class(df.hom) <- c("data.frame", "homolog.probability.mappoly")
  df.hom
}
print.homolog.probability.mappoly<-function(x, ...){
  print.data.frame(x, max = 50)
}
plot.homolog.probability.mappoly<-function(x, ind = NULL, use.plotly = TRUE){
  all.ind<-unique(df.pr$individual)
  if(is.null(ind)){
    ind<-as.character(all.ind[1])
    df.pr1<-subset(df.pr, individual == ind)  
  } else if(is.numeric(ind)) {
    if(ind > length(all.ind))
      stop("Please chose an individual number between 1 and ", length(all.ind))
    ind<-as.character(all.ind[ind])
    df.pr1<-subset(df.pr, individual == ind)  
  } else if (is.character(ind)){
    if(!ind%in%all.ind)
      stop("Invalid individual name")
  } else stop("Invalid individual name")
  p <- ggplot(df.pr1, aes(x = map.position, y = probability, fill = homolog, color  = homolog)) +
    geom_density(stat = "identity", alpha = 0.7) + 
    #scale_fill_manual(values = pal) + 
    #scale_color_manual(values = pal) +
    ggtitle(ind) + 
    facet_grid(rows = vars(homolog)) + 
    theme_minimal() + 
    ylab(label = "Homologs probabilty") +
    xlab(label = "Map position")
  if(use.plotly)
    p <- ggplotly(p)
  return(p)
}
df.pr<-homolog_prob(input.probs)

print.homolog.probability.mappoly(df.pr)
plot.homolog.probability.mappoly(df.pr, ind = 12)



