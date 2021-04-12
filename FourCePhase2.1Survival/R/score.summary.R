score.summary=function(dat.label, nm.event, score.cv, roc.cv, is.combine=F, is.sep=F){
t0.all=names(roc.cv)
score.summary=vector("list", length = length(t0.all))
names(score.summary)=t0.all


for(tt in t0.all){
if(is.sep==F){
score=score.cv[,c("patient_num",tt)]}else{
score=score.cv
}
colnames(score)[2]="score"
dat.score=left_join(dat.label, score, by="patient_num")
junk=matrix(roc.cv[[tt]]$roc, ncol=6)
colnames(junk)=c("cut", "p.pos", "fpr","tpr","ppv","npv")

cut1=junk[which.min(abs(junk[,"tpr"]-0.85)),"cut"]
cut2=junk[which.min(abs(junk[,"fpr"]-0.15)),"cut"]
score.cat=rep("M", length(dat.score$score))
score.cat[dat.score$score<cut1]="L"
score.cat[dat.score$score>cut2]="H"
score.cat=factor(score.cat, levels=c("L", "M", "H"))
score.cat=data.frame(y=dat.score[,nm.event],cat=score.cat, calendar_month=substr(dat.score$calendar_date,1,7))
score.cat$calendar_month=as.character(score.cat$calendar_month)
if(is.combine==T){score.cat[as.character(score.cat$calendar_month)>="2020-07","calendar_month"]="Since July"}
score.cat.all=score.cat
score.cat=table(score.cat)
score.summary[[tt]]=score.cat
}
score.summary
}