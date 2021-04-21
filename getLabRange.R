load(paste0("../Phase2.1SurvivalRSummariesPublic/BIDMC_Result.Rdata"))

nm.lab.all=unique(c(ls(lab.dist.original$res.all$all$'1'$res.early),
                 ls(lab.dist.original$res.all$all$'2'$res.early),
                 ls(lab.dist.original$res.all$all$'3'$res.early),
                 ls(lab.dist.original$res.all$all$'4'$res.early),
                 ls(lab.dist.original$res.all$all$'1'$res.late),
                 ls(lab.dist.original$res.all$all$'2'$res.late),
                 ls(lab.dist.original$res.all$all$'3'$res.late),
                 ls(lab.dist.original$res.all$all$'4'$res.late)))

lab.breaks.original=NULL
for(nm.lab in nm.lab.all){
  breaks.all=NULL
  for(period in c("early", "laste")){
    for(mylist in c(1:4)){
    junk=tryCatch(lab.dist.original$res.all$all[[mylist]][[paste0("res.", period)]][[nm.lab]]$breaks,error=function(e) NA)
    breaks.all=c(breaks.all, junk)
   }
  }
  res=data.frame(nm.lab, min(breaks.all, na.rm=T), max(breaks.all, na.rm=T))
  lab.breaks.original=rbind(lab.breaks.original, res)
}
colnames(lab.breaks.original)=c("labname", "break.lb", "break.ub")
lab.breaks.original$break.ub[lab.breaks.original$labname=="aspartate_aminotransferase_AST"]=2000
lab.breaks.original$break.ub[lab.breaks.original$labname=="Ferritin"]=10000
lab.breaks.original$break.ub[lab.breaks.original$labname=="lactate_dehydrogenase_LDH"]=2000
lab.breaks.original$break.lb=0


###
lab.breaks.log=NULL
for(nm.lab in nm.lab.all){
  breaks.all=NULL
  for(period in c("early", "laste")){
    for(mylist in c(1:4)){
      junk=tryCatch(lab.dist.log$res.all$all[[mylist]][[paste0("res.", period)]][[nm.lab]]$breaks,error=function(e) NA)
      breaks.all=c(breaks.all, junk)
    }
  }
  res=data.frame(nm.lab, min(breaks.all, na.rm=T), max(breaks.all, na.rm=T))
  lab.breaks.log=rbind(lab.breaks.log, res)
}
colnames(lab.breaks.log)=c("labname", "break.lb", "break.ub")
lab.breaks.original=rbind(lab.breaks.original, data.frame(labname="cardiac_troponin_High_Sensitivity", break.lb=0, break.ub=20))
lab.breaks.log=rbind(lab.breaks.log, data.frame(labname="cardiac_troponin_High_Sensitivity", break.lb=-10, break.ub=4))


save(lab.breaks.original, file="FourCePhase2.1Survival/data/lab.breaks.original.rda")
save(lab.breaks.log, file="FourCePhase2.1Survival/data/lab.breaks.log.rda")

