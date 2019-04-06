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
	--output	output dir
	--help		usage
\n")
	q(status=1);
}

if (is.null(opt$out)){ print_usage(spec) }

times<-Sys.time()

library(polymapR)
###load data
load(paste(opu$out,"step2.linkage/","step2_linkage.RData",sep=""))
dir.create(file.path(opt$out, "step3.evaluating.map"),showWarnings = FALSE)
setwd(file.path(opt$out, "step3.evaluating.map"))
source("/mnt/ilustre/users/meng.luo/Pipeline/10.Ploymap-Pipeline/bin/check-map.R")
for (i in 1:7){
	map<-integrated.maplist[[i]][c(1,2)]
	colnames(map)<-c("group",paste(i))
	write.table(map,file=paste("LG",i,".sexAver.map",sep=""),sep="\t",quote=F,row.names=F,col.names=T)
}

for (i in 1:7){
	map<-integrated.maplist_P1_chr1new[[i]][c(1,2)]
	colnames(map)<-c("group",paste(i))
	write.table(map,file=paste("LG",i,".male.map",sep=""),sep="\t",quote=F,row.names=F,col.names=T)
}

for (i in 1:7){
	map<-integrated.maplist_P2_chr1new[[i]][c(1,2)]
	colnames(map)<-c("group",paste(i))
	write.table(map,file=paste("LG",i,".female.map",sep=""),sep="\t",quote=F,row.names=F,col.names=T)
}

for (i in 1:7){
	maplist<-list(integrated.maplist[[i]])
	linkage_list<-list(linkages_chr1[[i]])
	LG<-i
	png(paste("checkmap","LG",LG,"_PP.png",sep=""),width=1000, height=800)
	checkmap(linkage_list, maplist, mapfn = "haldane", lod.thresh = 3,LG)
	dev.off()
	pdf(paste("checkmap","LG",LG,"_PP.pdf",sep=""),width=16, height=9)
	checkmap(linkage_list, maplist, mapfn = "haldane", lod.thresh = 3,LG)
	dev.off()

}

for (i in 1:7){
	maplist<-list(integrated.maplist_P1_chr1new[[i]])
	linkage_list<-list(all_linkages_list_P1_chr1[[i]][[1]])
	LG<-i
	png(paste("checkmap","LG",LG,"_P1.png",sep=""),width=1000, height=800)
	checkmap(linkage_list, maplist, mapfn = "haldane", lod.thresh = 5,LG)
	dev.off()
	pdf(paste("checkmap","LG",LG,"_P1.pdf",sep=""),width=16, height=9)
	checkmap(linkage_list, maplist, mapfn = "haldane", lod.thresh = 5,LG)
	dev.off()
}

for (i in 1:7){
	maplist<-list(integrated.maplist_P2_chr1new[[i]])
	linkage_list<-list(all_linkages_list_P2_chr1[[i]][[1]])
	LG<-i
	png(paste("checkmap","LG",LG,"_P2.png",sep=""),width=1000, height=800)
	checkmap(linkage_list, maplist, mapfn = "haldane", lod.thresh = 5,LG)
	dev.off()
	pdf(paste("checkmap","LG",LG,"_P2.pdf",sep=""),width=16, height=9)
	checkmap(linkage_list, maplist, mapfn = "haldane", lod.thresh = 5,LG)
	dev.off()

}

save.image("step3_evaluating_map.RData")

escaptime=Sys.time()-times;
print("Done!");
print(escaptime)
