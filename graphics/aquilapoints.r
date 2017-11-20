#!/usr/bin/Rscript
#plot core positions

sample <- 'starless' #'all', 'prestellar', or 'starless'
var <- 'mass'
#var: (corei, corenum, corename, RA, dec,) rcore, robs, mcore, (merr,) Tdust, (Terr,)
#     ncolPeak, ncolObs, ncolCore, nvolPeak, nvolObs, nvolCore, mBE, (coretype)
n <- 20 # 10, 20, 50

masterdir <- '/media/claire/Elements/Work/aquila/data'
setwd(masterdir)

objfn <- file.path('..',paste0('cores_',sample),var,paste0('n',n,'_objpositions.dat'))
objdat <- read.table(objfn,col.names=c("x","y","z","val"))

starless <- 'starless.dat'
prestellar <- 'prestellar.dat'
protostellar <- 'protostellar.dat'

starlessdat <- read.table(starless)#,col.names=c("corei","corenum","corename",...))
prestellardat <- read.table(prestellar)
protostellardat <- read.table(protostellar)

RAstarless <- starlessdat$V4
decstarless <- starlessdat$V5
RApre <- prestellardat$V4
decpre <- prestellardat$V5
RAproto <- protostellardat$V4
decproto <- protostellardat$V5
xmax <- max(RAstarless,RApre,RAproto)
xmin <- min(RAstarless,RApre,RAproto)
ymax <- max(decstarless,decpre,decproto)
ymin <- min(decstarless,decpre,decproto)



#h <- 800
#w <- 16*h / ((ymax-ymin)/(xmax-xmin))
#png(filename='testplot.png',width = w, height = h)

plot.new()
plot(RAstarless, decstarless, pch=20,col='gray40',#axes = FALSE,#ann = FALSE,
     xlim=c(xmin,xmax), ylim=c(ymin,ymax),
     #xlim=c(18.35,18.6), ylim=c(-5.5,0),
     xlab='RA (deg)', ylab='dec (deg)', cex.lab=1)
title(main='Positions of object cores')
mtext(paste0(var,', n=',n))

points(RApre,decpre,pch=3,col='sienna4')
points(RAproto,decproto,pch=4,col='black')
points(objdat$x,objdat$y,pch=19,col='deeppink')

#dev.off() #close plot