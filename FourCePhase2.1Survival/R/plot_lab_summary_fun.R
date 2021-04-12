plot_lab_summary_fun=function(lab.res, siteid, nm.lab, nm.setting="mean", ylim.range){
  mycol=rainbow(10)
  dat.junk=lab.res
      tmp=data.frame(cbind(dat.junk[[1]][,nm.lab],
                           dat.junk[[2]][,nm.lab],
                           dat.junk[[3]][,nm.lab]))
      
      nn=dat.junk[[1]][,"n"]
      
      rownames(tmp)=c("March-April", "May-June", "July-August", "September-October", "November-December")[1:dim(tmp)[1]]
      colnames(tmp)=c("day0", "day0-1", "day0-7")
      tmp=apply(tmp,2,as.numeric)
      ylim=ylim.range
      if(nm.setting=="obs"){ylim=c(0,1)}else{
      ylim[1]=0
      ylim[2]=ylim[2]*1.5
      }
      barplot(tmp, beside = TRUE, ylim=ylim,col=mycol[1:dim(tmp)[1]], main=ifelse(nm.setting=="obs",paste0("observation Rate for ", nm.lab, " [", siteid, "]"),paste0("mean measure for ", nm.lab, " [", siteid, "]")))
      legend("topleft", paste0(c("March-April", "May-June", "July-August", "September-October", "November-December")[1:dim(tmp)[1]], " [n=", nn,"]"), col=mycol[1:dim(tmp)[1]], pch=15, ncol=3)

  }


