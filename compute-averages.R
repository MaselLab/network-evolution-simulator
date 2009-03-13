library(lattice)

if (TRUE) {

  for (i in seq(00, 99)) {
    cmd=sprintf("time ../../netsim -g -r %i  -d output%02d", i, i)
    print(cmd)
    system(cmd)
    cmd2=sprintf("mv netsimerrors.txt output%02d", i, i)
    print(cmd2)
    system(cmd2)
  }
}

df = NULL

for (i in seq(0,9)) {
  fs=paste(dir(pattern="output[0-9]+",full.names=T),"/protein",i,".dat",sep="")
  ##print(fs)
  ##print(i)
  for (f in fs) {
    if (!file.exists(f)) {
      next
    }
    series=read.table(f)
    names(series)=c("time", "conc")

    ## get the average of the length of timesteps
    time.av = mean(diff(series$time))

    ## get the mean protein concentration over the time series
    conc.av = mean(series$conc)

    ## get the last protein concentration
    conc.last = tail(series,n=1)$conc
    ##print(time.av)
    ##print(conc.av)
    df <- rbind(df, data.frame(gene=i, time.av=time.av, conc.av=conc.av, conc.last=conc.last))
  }
}
print(df)

bwplot(gene~conc.av, data=df, scales=list(x=list(log=TRUE)))
