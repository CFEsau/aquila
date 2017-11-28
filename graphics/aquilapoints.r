#!/usr/bin/Rscript
#plot core positions

plotdir <- getwd() #save current working directory

samples <- c('all','prestellar','starless') #'all': starless, prestellar & protostellar
                       #'prestellar': prestellar only
                       #'starless': starless & prestellar
varlist<- c('mass','mBE','T','rcore','ncolCore','nvolCore','invT')
#var: (corei, corenum, corename, RA, dec,) rcore, robs, mcore, (merr,) T, (Terr,)
#     ncolPeak, ncolObs, ncolCore, nvolPeak, nvolObs, nvolCore, mBE, (coretype)
nmst <- c(10, 20, 50) #list of numbers of MST object stars

masterdir <- '~/Documents/Work/aquila/data'
setwd(masterdir)

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

for (sample in samples){
  for (var in varlist){
    for (n in nmst){
      objfn <- file.path('..',paste0('cores_',sample),var,paste0('n',n,'_objpositions.dat'))
      objdat <- read.table(objfn,col.names=c("x","y","val"))
      
      h <- 800
      w <- 16*h / ((ymax-ymin)/(xmax-xmin))
      png(filename=file.path(plotdir,paste0(var,'_',sample,'_n',n,'.png')),
          width = w, height = h)
      
      #plot.new()
      #always plot prestellar cores:
      plot(RApre, decpre, pch=20,col='gray40',#axes = FALSE,#ann = FALSE,
          xlim=c(xmin,xmax), ylim=c(ymin,ymax),
          #xlim=c(18.35,18.6), ylim=c(-5.5,0),
          xlab='RA (deg)', ylab='dec (deg)', cex.lab=1)
      title(main='Positions of object cores')
      mtext(paste0(var,', n=',n))
      
      if (sample=='all'){ #also plot starless and protostellar cores:
        points(RAstarless,decstarless,pch=3,col='sienna4')
        points(RAproto,decproto,pch=4,col='black')
      } else if (sample=='starless'){ #plot starless cores:
        points(RAstarless,decstarless,pch=3,col='sienna4')
      }
      #(if sample is prestellar, only plot prestellar cores)
      
      #plot object cores:
      points(objdat$x,objdat$y,pch=17,col='red')
      
      dev.off() #close plotting device
    }
  }
}

setwd(plotdir) #go back to original directory