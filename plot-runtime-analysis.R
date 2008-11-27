library(lattice)
##runtimeall=read.table("/tmp/new-runtime/runtime-all.txt", head=T)
##runtimeall=read.table("/tmp/new-runtime/runtime-profile-double-gentime.txt", head=T)
runtimeall=read.table("/tmp/new-runtime/runtime-profile-small-burnin-new.txt", head=T)
outall=aggregate(runtimeall, by=list(runtimeall$kon), mean)
outall=transform(subset(outall, kon>1e-10), gen=77040000/(500*real*33))

##png("generations-vs-kon.png")
##postscript("generations-vs-kon.ps")
##postscript("generations-vs-kon-double.ps")
postscript("generations-vs-kon-new.ps")

print(outall)

xat=10^seq(-10,10)
xlabels=seq(-10,10)

yat=10^seq(-10,10,0.5)
ylabels=seq(-10,10,0.5)

xyplot(gen~kon,
       data=outall,
       scales=list(log=TRUE, cex=2.0,
         x=list(at=xat,labels=xlabels),
         y=list(at=yat,labels=ylabels)),
       xlab=list(expression(log~~italic(k[on])), cex=2.0),
       ylab=list(expression(log~~italic(generations)), cex=2.0),
       panel=function(x, y,...) {
         panel.xyplot(x, y, ...)
         sbs=subset(outall, kon>1e-8)
         regr=lm(log10(sbs$gen)~log10(sbs$kon))
         slope=print(regr[[1]][[2]])
         intercept=print(regr[[1]][[1]])
         cat(sprintf("gen for k=0.2225 = %g\n", 10^(intercept)*(0.2225^slope)))
         panel.abline(regr)
       })
dev.off()
