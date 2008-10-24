library(lattice)
runtimeall=read.table("/tmp/new-runtime/runtime-all.txt", head=T)
outall=aggregate(runtimeall, by=list(runtimeall$kon), mean)
outall=transform(subset(outall, kon>1e-10), gen=77040000/(500*real*33))

#png("generations-vs-kon.png")
postscript("generations-vs-kon.ps")

myat=10^seq(-10,10)
mylabels=seq(-10,10)

xyplot(gen~kon,
       data=outall,
       scales=list(log=TRUE, cex=2.0,
         x=list(at=myat,labels=mylabels),
         y=list(at=myat,labels=mylabels)),
       xlab=list(expression(log~~italic(k[on])), cex=2.0),
       ylab=list(expression(log~~italic(generations)), cex=2.0),
       panel=function(x, y,...) {
         panel.xyplot(x, y, ...)
         sbs=subset(outall, kon>1e-8)
         regr=lm(log10(sbs$gen)~log10(sbs$kon)); print(regr[[1]][[2]]); print(regr[[1]][[1]]); 
         panel.abline(regr)
       })
dev.off()
