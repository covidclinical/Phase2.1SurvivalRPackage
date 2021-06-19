obfuscation.fun=function(summary.report, KM, survfit.coxnet,lab.dist.original, lab.dist.log, lab.summary, lab.recover, lab.recover.rmDead,cls.obs.summary, cls.early, cls.late, obfuscation.level){
  obfuscation.level=as.numeric(obfuscation.level)
  summary.report[summary.report<obfuscation.level]=-99
  tryCatch({
  nm.check= c("n.risk", "n.event", "n.censor")
  for(nm in nm.check){
    tmp=as.numeric(as.character(KM[[nm]]))
    tmp[tmp<obfuscation.level]=-99
    KM[[nm]]=tmp
  }},error=function(e) NA)
  
  for(aa in ls(survfit.coxnet)){
    for(bb in ls(survfit.coxnet[[aa]])){
      for(cc in ls(survfit.coxnet[[aa]][[bb]])){
        for(dd in ls(survfit.coxnet[[aa]][[bb]][[cc]])){
          for(ee in ls(survfit.coxnet[[aa]][[bb]][[cc]][[dd]])){
            tmp=survfit.coxnet[[aa]][[bb]][[cc]][[dd]][[ee]]
            nm.check=names(tmp)[substr(names(tmp),1,5)=="score"]
            for(nm in nm.check){
            tmp[[nm]][['28']][tmp[[nm]][['28']]<obfuscation.level]=-99
            }
            survfit.coxnet[[aa]][[bb]][[cc]][[dd]][[ee]]=tmp
        }
        
      }
      
    }
    }
  }

  
  for(aa in ls(lab.dist.original$res.all)){
    tmp.aa=lab.dist.original$res.all[[aa]]
    for(bb in ls(tmp.aa)){
      tmp.bb=tmp.aa[[bb]]
      for(cc in ls(tmp.bb)){
        tmp.cc=tmp.bb[[cc]]
        for(dd in ls(tmp.cc)){
          tmp.dd=tmp.cc[[dd]]
          if(is.na(tmp.dd[1])!=1){
            tmp.dd$counts[tmp.dd$counts<obfuscation.level]=-99
            lab.dist.original$res.all[[aa]][[bb]][[cc]][[dd]]=tmp.dd}
        }
      }
    }
  }
  
  for(aa in ls(lab.dist.log$res.all)){
    tmp.aa=lab.dist.log$res.all[[aa]]
    for(bb in ls(tmp.aa)){
      tmp.bb=tmp.aa[[bb]]
      for(cc in ls(tmp.bb)){
        tmp.cc=tmp.bb[[cc]]
        for(dd in ls(tmp.cc)){
          tmp.dd=tmp.cc[[dd]]
          if(is.na(tmp.dd[1])!=1){
            tmp.dd$counts[tmp.dd$counts<obfuscation.level]=-99
            lab.dist.log$res.all[[aa]][[bb]][[cc]][[dd]]=tmp.dd}
        }
      }
    }
  }
  
  for(aa in ls(lab.dist.original$obs.all)){
    tmp.aa=lab.dist.original$obs.all[[aa]]
    for(bb in ls(tmp.aa)){
      tmp.bb=tmp.aa[[bb]]
      for(cc in ls(tmp.bb)){
        tmp.cc=tmp.bb[[cc]]
        for(dd in ls(tmp.cc)){
          tmp.dd=tmp.cc[[dd]]
          if(is.na(tmp.dd[1])!=1){
            tmp.dd[tmp.dd<obfuscation.level]=-99
            lab.dist.original$obs.all[[aa]][[bb]][[cc]][[dd]]=tmp.dd}
        }
      }
    }
  }
  
  for(aa in ls(lab.dist.log$obs.all)){
    tmp.aa=lab.dist.log$obs.all[[aa]]
    for(bb in ls(tmp.aa)){
      tmp.bb=tmp.aa[[bb]]
      for(cc in ls(tmp.bb)){
        tmp.cc=tmp.bb[[cc]]
        for(dd in ls(tmp.cc)){
          tmp.dd=tmp.cc[[dd]]
          if(is.na(tmp.dd[1])!=1){
            tmp.dd[tmp.dd<obfuscation.level]=-99
            lab.dist.log$obs.all[[aa]][[bb]][[cc]][[dd]]=tmp.dd}
        }
      }
    }
  }
  cls.obs.summary$obs_sum[cls.obs.summary$obs_sum<obfuscation.level]=-99
  cls.early$counts[cls.early$counts<obfuscation.level]=-99
  cls.late$counts[cls.late$counts<obfuscation.level]=-99
  
  for(aa in ls(lab.summary)){
    tmp.aa=lab.summary[[aa]]
    for(bb in ls(tmp.aa)){
      tmp.bb=tmp.aa[[bb]]
      for(cc in 1:length(tmp.bb)){
        tmp.cc=tmp.bb[[cc]]
        tmp.cc$n[tmp.cc$n<=obfuscation.level]=-99
        lab.summary[[aa]][[bb]][[cc]]=tmp.cc
      }
    }
  }
  
  for(aa in ls(lab.recover)){
    tmp.aa=lab.recover[[aa]]
    for(bb in ls(tmp.aa)){
      tmp.bb=tmp.aa[[bb]]
      for(cc in 1:length(tmp.bb)){
        tmp.cc=tmp.bb[[cc]]
        for(dd in ls(tmp.cc)){
          tmp.dd=tmp.cc[[dd]]
          for(ee in ls(tmp.dd)){
            if(ee!="resN"){
            tmp.ee=tmp.dd[[ee]]
            for(ff in ls(tmp.ee))
              tmp.ff=tmp.ee[[ff]]
            if(tmp.ff$n.early<=obfuscation.level){tmp.ff$n.early=-99}
            if(tmp.ff$n.late<=obfuscation.level){tmp.ff$n.late=-99}
            if(tmp.ff$n.lab.early<=obfuscation.level){tmp.ff$n.lab.early=-99}
            if(tmp.ff$n.lab.late<=obfuscation.level){tmp.ff$n.lab.late=-99}
            }
          }
          
        }
      }
    }
  }
  
  lab.recover$max_day$`0`$'Inf'$`1-14`$resN[lab.recover$max_day$`0`$'Inf'$`1-14`$resN<=obfuscation.level]=-99
  
  
  
  for(aa in ls(lab.recover.rmDead)){
    tmp.aa=lab.recover[[aa]]
    for(bb in ls(tmp.aa)){
      tmp.bb=tmp.aa[[bb]]
      for(cc in 1:length(tmp.bb)){
        tmp.cc=tmp.bb[[cc]]
        for(dd in ls(tmp.cc)){
          tmp.dd=tmp.cc[[dd]]
          for(ee in ls(tmp.dd)){
            if(ee!="resN"){
              tmp.ee=tmp.dd[[ee]]
              for(ff in ls(tmp.ee))
                tmp.ff=tmp.ee[[ff]]
              if(tmp.ff$n.early<=obfuscation.level){tmp.ff$n.early=-99}
              if(tmp.ff$n.late<=obfuscation.level){tmp.ff$n.late=-99}
              if(tmp.ff$n.lab.early<=obfuscation.level){tmp.ff$n.lab.early=-99}
              if(tmp.ff$n.lab.late<=obfuscation.level){tmp.ff$n.lab.late=-99}
            }
          }
          
        }
      }
    }
  }
  
  lab.recover.rmDead$max_day$`0`$'Inf'$`1-14`$resN[lab.recover$max_day$`0`$'Inf'$`1-14`$resN<=obfuscation.level]=-99
  
  return(list(summary.report=summary.report, 
              KM=KM,
              survfit.coxnet=survfit.coxnet,
              lab.dist.original=lab.dist.original, 
              lab.dist.log=lab.dist.log, 
              lab.summary=lab.summary, 
              lab.recover=lab.recover,
              lab.recover.rmDead=lab.recover.rmDead,
              cls.obs.summary=cls.obs.summary,
              cls.early=cls.early,
              cls.late=cls.late
              ))
  
}