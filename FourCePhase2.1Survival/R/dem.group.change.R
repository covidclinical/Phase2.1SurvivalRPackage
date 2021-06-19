dem.group.change.fun=function(dat.survival, nm.var, old.group.nm, new.group.nm){
  for(dat.nm in c("dat.analysis.deceased", "dat.analysis.severe", "dat.analysis.severedeceased")){
  tmp=dat.survival[[dat.nm]][,nm.var]
  levels(tmp)[match(old.group.nm,levels(tmp))] <- new.group.nm
  dat.survival[[dat.nm]][,nm.var]=tmp
  }
  dat.survival
}
