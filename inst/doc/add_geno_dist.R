myfunc<-function(x){
  if(is.na(x))
    return(rep(NA,7))
  v<-rep(0,7)
  v[x+1]<-1
  return(v)
}

x<-as.data.frame(as.table(as.matrix(hexafake$geno.dose)))
x$Freq[sample(1:nrow(x), 2*nrow(x)/100)]<-NA
y<-t(sapply(x$Freq, myfunc))
head(y)

z<-cbind(x[,1:2], y)
head(z)
write.table(x = z, file = "~/repos/bla.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
