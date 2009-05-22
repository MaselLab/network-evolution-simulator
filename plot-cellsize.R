library(lattice)

df = NULL

for (i in seq(0,499)) {
  fs=sprintf("cellsize-%03d.dat",i)
  print(fs)
  ##print(i)
  for (f in fs) {
    if (!file.exists(f)) {
      next
    }
    series=read.table(f)
    ## names(series)=c("time (min)", "protein molecules/cell")
    names(series)=c("t", "cellsize")

    df <- rbind(df, data.frame(gene=i, t=series$t, cellsize=series$cellsize))
  }
}

xyplot(cellsize~t, groups=gene, data=df,
       auto.key=list(space="right", type='l'), type="l", 
       xlab="time (min)", ylab="cellsize")
