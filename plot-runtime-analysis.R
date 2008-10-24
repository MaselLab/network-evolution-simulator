runtimeall=read.table("/tmp/new-runtime/runtime-all.txt", head=T)
outall=aggregate(runtimeall, by=list(runtimeall$kon), mean)
xyplot(gen~kon,
       scales=list(y=list(log=T),x=list(log=T)),
       data=outall,
       panel=function(x, y,...) {
         panel.xyplot(x, y, ...)
         sbs=subset(outall, kon>1e-8)
         regr=lm(log10(sbs$gen)~log10(sbs$kon)); print(regr[[1]][[2]]); print(regr[[1]][[2]]); 
         panel.abline(regr)
       })
