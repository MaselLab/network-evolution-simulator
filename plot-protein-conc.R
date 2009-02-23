library(lattice)

df = NULL

for (i in seq(0,10)) {
  fs=paste("protein",i,".dat",sep="")
  ##print(fs)
  ##print(i)
  for (f in fs) {
    if (!file.exists(f)) {
      next
    }
    series=read.table(f)
    ## names(series)=c("time (min)", "protein molecules/cell")
    names(series)=c("t", "protein")

    df <- rbind(df, data.frame(gene=i, t=series$t, protein=series$protein))
  }
}

xyplot(protein~t, groups=gene, data=df,
       auto.key=list(space="right", type='l'), type="l", 
       xlab="time (min)", ylab="protein molecules/cell")
