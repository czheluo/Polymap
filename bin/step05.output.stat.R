#!/usr/bin/env Rscript
library('getopt');
options(bitmapType='cairo')
spec = matrix(c(
	'out','o',1,'character',
	'help','h',0,'logical'
	), byrow=TRUE, ncol=4)
opt = getopt(spec);
print_usage <- function(spec=NULL){
	cat(getopt(spec, usage=TRUE));
	cat("Usage example: \n")
	cat("
Usage example:
	Rscript G.R --input --hb --lb --output --all --opt
Usage:
	--out output dir
	--help		usage
\n")
	q(status=1);
}

if (is.null(opt$out)){ print_usage(spec) }

times<-Sys.time()

library(polymapR)

setwd(paste(opu$out,"01.vcf-convert/",sep=""))
ALL_dos<-read.csv("dosage.matrix",header=T)
map<-read.csv("dosage.matrix.map",header=T)
load(paste(opu$out,"step4.qtlanalysis/","step4_qtlanalysis.RData",sep=""))
dir.create(file.path(opt$out, "step5.result.stat"),showWarnings = FALSE)
setwd(file.path(opt$out, "step5.result.stat"))

####output genetype

integrated.maplist_P1<-integrated.maplist_P1_chr1new

###FATHER MAP
MarN<-NULL;MapD<-NULL;AvD<-NULL;MaxG<-NULL;#Gap5cM<-NULL
Position_P1<-list();Genetype_P1<-list()
for(i in 1:7){
  MarN[i]<-dim(integrated.maplist_P1[[i]])[1]
  MapD[i]<-max(integrated.maplist_P1[[i]][,2])
  AvD[i]<-max(integrated.maplist_P1[[i]][,2])/dim(integrated.maplist_P1[[i]])[1]
  MaxG[i]<-max(integrated.maplist_P1[[i]]$position[c(-1)]-integrated.maplist_P1[[i]]$position[c(-length(integrated.maplist_P1[[i]]$position))])
  #Gap5cM
  Position_P1[[i]]<-map[which(ALL_dos$SNP %in% integrated.maplist_P1[[i]][,1]),]
  Genetype_P1[[i]]<-ALL_dos[which(ALL_dos$SNP %in% integrated.maplist_P1[[i]][,1]),]
  }

phymap<-do.call(rbind,Position_P1)
write.table(phymap,file="phymap_P1.txt",
           quote = F,row.names = F,col.names = T)

Genetype<-do.call(rbind,Genetype_P1)

write.table(Genetype,file="Genetype_P1.txt",
            quote = F,row.names = F,col.names = T)


all_stats<-as.matrix(cbind(MarN,MapD,AvD,MaxG))

colnames(all_stats)<-c("Mar Num","Map Distance","Aver Distance","Max Gap")

for(i in 1:7){
  rownames(all_stats)[i]<-paste("LG",i,sep = "")
}
write.csv(all_stats,file = "stats_P1.csv",quote=F,row.names=T)

integrated.maplist_P2<-integrated.maplist_P2_chr1new

####MOTHER MAP
MarN<-NULL;MapD<-NULL;AvD<-NULL;MaxG<-NULL;#Gap5cM<-NULL
Position_P2<-list();Genetype_P2<-list()
for(i in 1:7){
  MarN[i]<-dim(integrated.maplist_P2[[i]])[1]
  MapD[i]<-max(integrated.maplist_P2[[i]][,2])
  AvD[i]<-max(integrated.maplist_P2[[i]][,2])/dim(integrated.maplist_P2[[i]])[1]
  MaxG[i]<-max(integrated.maplist_P2[[i]]$position[c(-1)]-integrated.maplist_P2[[i]]$position[c(-length(integrated.maplist_P2[[i]]$position))])
  #Gap5cM
  Position_P2[[i]]<-map[which(ALL_dos$SNP %in% integrated.maplist_P2[[i]][,1]),]
  Genetype_P2[[i]]<-ALL_dos[which(ALL_dos$SNP %in% integrated.maplist_P2[[i]][,1]),]
}

phymap<-do.call(rbind,Position_P2)
write.table(phymap,file="phymap_P2.txt",
            quote = F,row.names = F,col.names = T)

Genetype<-do.call(rbind,Genetype_P2)

write.table(Genetype,file="Genetype_P2.txt",
            quote = F,row.names = F,col.names = T)


all_stats<-as.matrix(cbind(MarN,MapD,AvD,MaxG))

colnames(all_stats)<-c("Mar Num","Map Distance","Aver Distance","Max Gap")

for(i in 1:7){
  rownames(all_stats)[i]<-paste("LG",i,sep = "")
}
write.csv(all_stats,file = "stats_P2.csv",quote=F,row.names=T)

integrated.maplist<-integrated.maplist_chr1
###INTEGRATE MAP
MarN<-NULL;MapD<-NULL;AvD<-NULL;MaxG<-NULL;#Gap5cM<-NULL
Position_PP<-list();Genetype_PP<-list()
for(i in 1:7){
  MarN[i]<-dim(integrated.maplist[[i]])[1]
  MapD[i]<-max(integrated.maplist[[i]][,2])
  AvD[i]<-max(integrated.maplist[[i]][,2])/dim(integrated.maplist[[i]])[1]
  MaxG[i]<-max(integrated.maplist[[i]]$position[c(-1)]-integrated.maplist[[i]]$position[c(-length(integrated.maplist[[i]]$position))])
  #Gap5cM
  Position_PP[[i]]<-map[which(ALL_dos$SNP %in% integrated.maplist[[i]][,1]),]
  Genetype_PP[[i]]<-ALL_dos[which(ALL_dos$SNP %in% integrated.maplist[[i]][,1]),]
}

phymap<-do.call(rbind,Position_PP)
write.table(phymap,file="phymap_PP.txt",
            quote = F,row.names = F,col.names = T)

Genetype<-do.call(rbind,Genetype_PP)

write.table(Genetype,file="Genetype_PP.txt",
            quote = F,row.names = F,col.names = T)


all_stats<-as.matrix(cbind(MarN,MapD,AvD,MaxG))

colnames(all_stats)<-c("Mar Num","Map Distance","Aver Distance","Max Gap")

for(i in 1:7){
  rownames(all_stats)[i]<-paste("LG",i,sep = "")
}
write.csv(all_stats,file = "stats_PP.csv",quote=F,row.names=T)


escaptime=Sys.time()-times;
print("Done!");
print(escaptime)



