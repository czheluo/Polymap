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
	--out	output dir
	--help		usage
\n")
	q(status=1);
}

if (is.null(opt$out)){ print_usage(spec) }

times<-Sys.time()

library(polymapR)

load(paste(opu$out,"step3.evaluating.map/","step3_evaluating_map.RData",sep=""))
dir.create(file.path(opt$out, "step4.qtlanalysis"),showWarnings = FALSE)
setwd(file.path(opt$out, "step4.qtlanalysis"))
source("/mnt/ilustre/users/meng.luo/Pipeline/10.Ploymap-Pipeline/bin/create_phased_maplist.R")


#phased.maplist

dir.create(file.path(getwd(), "TetraploidSNPMap_QTLfiles_pp"),showWarnings = FALSE)

setwd(file.path(getwd(), "TetraploidSNPMap_QTLfiles_pp"))

integrated.maplist<-integrated.maplist_chr1

marker_assignments<-marker_assignments_chr1

filtered_dosage<-do.call(rbind,unname(filtered_data))
phased.maplist<-list()
for(i in 1:7){
	phased.maplist [[i]]<- create_phased_maplist(maplist = integrated.maplist[i],
                                        dosage_matrix.conv = filtered_dosage,
					N_linkages = 7,
                                        ploidy = 4,
                                        marker_assignment.1 = marker_assignments[[i]]$P1,
                                        marker_assignment.2 = marker_assignments[[i]]$P2,
					)
}
phased.maplist_new<-list()

for (i in 1:7){

	phased.maplist_new[i]<-phased.maplist[[i]]

	}
names(phased.maplist_new)<-c("LG1","LG2","LG3","LG4","LG5","LG6","LG7")
write.TSNPM(phased.maplist = phased.maplist_new,ploidy=4)


setwd(file.path(opt$out, "step4.qtlanalysis"))

dir.create(file.path(getwd(), "TetraploidSNPMap_QTLfiles_p1"),showWarnings = FALSE)

setwd(file.path(getwd(), "TetraploidSNPMap_QTLfiles_p1"))

integrated.maplist_P1<-integrated.maplist_P1_chr1new

phased.maplist_P1<-list()

for(i in 1:7){
	phased.maplist_P1 [[i]]<- create_phased_maplist(maplist = integrated.maplist_P1[i],
                                        dosage_matrix.conv = filtered_dosage,
					N_linkages = 7,
                                        ploidy = 4,
                                        marker_assignment.1 = marker_assignments[[i]]$P1,
                                        marker_assignment.2 = marker_assignments[[i]]$P2,
					)
}

phased.maplist_P1NEW<-list()

for (i in 1:7){

	phased.maplist_P1NEW[i]<-phased.maplist[[i]]

	}

names(phased.maplist_P1NEW)<-c("LG1","LG2","LG3","LG4","LG5","LG6","LG7")
write.TSNPM(phased.maplist = phased.maplist_P1NEW,ploidy=4)



setwd(file.path(opt$out, "step4.qtlanalysis"))

dir.create(file.path(getwd(), "TetraploidSNPMap_QTLfiles_p2"),showWarnings = FALSE)

setwd(file.path(getwd(), "TetraploidSNPMap_QTLfiles_p2"))

integrated.maplist_P2<-integrated.maplist_P2_chr1new
phased.maplist_P2<-list()

for(i in 1:7){

	phased.maplist_P2 [[i]]<- create_phased_maplist(maplist = integrated.maplist_P2[i],
                                        dosage_matrix.conv = filtered_dosage,
					N_linkages = 7,
                                        ploidy = 4,
                                        marker_assignment.1 = marker_assignments[[i]]$P1,
                                        marker_assignment.2 = marker_assignments[[i]]$P2,
					)
}

phased.maplist_P2NEW<-list()

for (i in 1:7){

	phased.maplist_P2NEW[i]<-phased.maplist_P2[[i]]

	}
names(phased.maplist_P2NEW)<-c("LG1","LG2","LG3","LG4","LG5","LG6","LG7")

write.TSNPM(phased.maplist = phased.maplist_P2NEW,ploidy=4)


save.image("step4_qtlanalysis.RData")

escaptime=Sys.time()-times;

print("Done!");
print(escaptime)




