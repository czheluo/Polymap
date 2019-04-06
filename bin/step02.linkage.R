#!/usr/bin/env Rscript
library('getopt');
options(bitmapType='cairo')
spec = matrix(c(
	'out','o',1,'character',
	'pid','p',1,'character',
	'mid','m',1,'character',
	'lod','l',1,'character',
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
	--pid	<str>	paternal id
	--mid	<str>	maternal id
	--out	output dir
	--lod the threathold value for the select marker to calculate (default was 2)
	--help		usage
\n")
	q(status=1);
}
if (is.null(opt$pid)) { print_usage(spec)}
if (is.null(opt$mid)){ print_usage(spec) }
if (is.null(opt$lod)) { opt$lod=2;print(opt$lod))}else{print(opt$lod)}

times<-Sys.time()

library(polymapR)

load(paste(opu$out,"01.vcf-convert/","step1_filter.RData",sep=""))
dir.create(file.path(opt$out, "step2.linkage"),showWarnings = FALSE)
setwd(file.path(opt$out, "step2.linkage"))

##lod threshold for filter

lod<-opt$lod

SN_SS_P1_chr1<-list();P1_SxS_Assigned_chr1<-list();
SN_DN_P1_chr1<-list();P1_DxN_Assigned_chr1<-list();
SN_SS_P2_chr1<-list();P2_SxS_Assigned_chr1<-list();
SN_DN_P2_chr1<-list();P2_DxN_Assigned_chr1<-list();
marker_assignments_P1_chr1<-list();marker_assignments_P2_chr1<-list();
marker_assignments_chr1<-list();
integrated.maplist_P1_chr1<-list();integrated.maplist_P2_chr1<-list();
integrated.maplist_chr1<-list();
linkages_chr1 <- list();
all_linkages_list_P1_chr1<-list();all_linkages_list_P2_chr1<-list();
filtered_data<-all_dosages_list

for (i in 1:7){

	SN_SS_P1_chr1[[i]] <- linkage(dosage_matrix = filtered_data[[i]],
                    markertype1 = c(1,0),
                    markertype2 = c(1,1),
                    target_parent = P1,
                    other_parent = P2,
                    ploidy = 4,
                    pairing = "random")

	P1_SxS_Assigned_chr1[[i]] <- assign_linkage_group(linkage_df = SN_SS_P1_chr1[[i]],
                                        LG_hom_stack = LGHomDf_P1[[i]],
                                        SN_colname = "marker_a",
                                        unassigned_marker_name = "marker_b",
                                        phase_considered = "coupling",
                                        LG_number = 1,
                                        LOD_threshold = lod,
                                        ploidy = 4)

	pq_before_convert <- parental_quantities(dosage_matrix = filtered_data[[i]],parent1 = P1,parent2 = P2, las=2)
	labels = round(pq_before_convert,0)
	if(labels["2x0"]>0){
		SN_DN_P1_chr1[[i]] <- linkage(dosage_matrix = filtered_data[[i]],
                    markertype1 = c(1,0),
                    markertype2 = c(2,0),
                    target_parent = P1,
                    other_parent = P2,
                    ploidy = 4,
                    pairing = "random")
		}else{
		SN_DN_P1_chr1[[i]] <- linkage(dosage_matrix = filtered_data[[i]],
                    markertype1 = c(1,0),
                    #markertype2 = c(2,0),
                    target_parent = P1,
                    other_parent = P2,
                    ploidy = 4,
                    pairing = "random")
		}

	P1_DxN_Assigned_chr1[[i]]<- assign_linkage_group(linkage_df = SN_DN_P1_chr1[[i]],
                                        LG_hom_stack = LGHomDf_P1[[i]],
                                        SN_colname = "marker_a",
                                        unassigned_marker_name = "marker_b",
                                        phase_considered = "coupling",
                                        LG_number = 1,
                                        LOD_threshold = lod,
                                        ploidy = 4)

	SN_SS_P2_chr1[[i]] <- linkage(dosage_matrix = filtered_data[[i]],
                    markertype1 = c(1,0),
                    markertype2 = c(1,1),
                    target_parent = P2,
                    other_parent = P1,
                    ploidy = 4,
                    pairing = "random")

	P2_SxS_Assigned_chr1[[i]] <- assign_linkage_group(linkage_df = SN_SS_P2_chr1[[i]],
                                        LG_hom_stack = LGHomDf_P2[[i]],
                                        SN_colname = "marker_a",
                                        unassigned_marker_name = "marker_b",
                                        phase_considered = "coupling",
                                        LG_number = 1,
                                        LOD_threshold = lod,
                                        ploidy = 4)

	pq_before_convert <- parental_quantities(dosage_matrix = filtered_data[[i]],parent1 = P1,parent2 = P2, las=2)
	labels = round(pq_before_convert,0)
	if(labels["2x0"]>0){
		SN_DN_P2_chr1[[i]] <- linkage(dosage_matrix = filtered_data[[i]],
                    markertype1 = c(1,0),
                    markertype2 = c(2,0),
                    target_parent = P2,
                    other_parent = P1,
                    ploidy = 4,
                    pairing = "random")
		}else{
		SN_DN_P2_chr1[[i]] <- linkage(dosage_matrix = filtered_data[[i]],
                    markertype1 = c(1,0),
                    #markertype2 = c(2,0),
                    target_parent = P2,
                    other_parent = P1,
                    ploidy = 4,
                    pairing = "random")
		}

	P2_DxN_Assigned_chr1[[i]] <- assign_linkage_group(linkage_df = SN_DN_P2_chr1[[i]],
                                        LG_hom_stack =  LGHomDf_P2[[i]],
                                        SN_colname = "marker_a",
                                        unassigned_marker_name = "marker_b",
                                        phase_considered = "coupling",
                                        LG_number = 1,
                                        LOD_threshold =lod,
                                        ploidy = 4)

	marker_assignments_P1_chr1[[i]] <- homologue_lg_assignment(dosage_matrix = filtered_data[[i]],
                                                 assigned_list = list(P1_SxS_Assigned_chr1[[i]],
                                                                      P1_DxN_Assigned_chr1[[i]]),
                                                 assigned_markertypes = list(c(1,1),c(2,0)),
                                                 LG_hom_stack = LGHomDf_P1[[i]],
                                                 target_parent = P1,
                                                 other_parent = P2,
                                                 ploidy = 4,
                                                 pairing = "random",
                                                 convert_palindrome_markers = FALSE,
                                                 LG_number = 1,
                                                 LOD_threshold = lod,
                                                 write_intermediate_files = FALSE
	)

	marker_assignments_P2_chr1[[i]] <-homologue_lg_assignment(dosage_matrix = filtered_data[[i]],
                           assigned_list = list(P2_SxS_Assigned_chr1[[i]],P2_DxN_Assigned_chr1[[i]]),
                           assigned_markertypes = list(c(1,1), c(2,0)),
                           LG_hom_stack = LGHomDf_P2[[i]],
                           target_parent = P2,
                           other_parent = P1,
                           ploidy = 4,
                           pairing = "random",
                           convert_palindrome_markers = TRUE,
                           LG_number =1,
                           LOD_threshold =lod,
                           write_intermediate_files = FALSE
	)

	marker_assignments_chr1[[i]] <- check_marker_assignment(marker_assignments_P1_chr1[[i]],marker_assignments_P2_chr1[[i]])

	all_linkages_list_P1_chr1[[i]] <- finish_linkage_analysis(marker_assignment = marker_assignments_chr1[[i]]$P1,
                                                dosage_matrix = filtered_data[[i]],
                                                target_parent = P1,
                                                other_parent = P2,
                                                convert_palindrome_markers = FALSE,
                                                ploidy = 4,
                                                pairing = "random",
                                                LG_number = 1)

	all_linkages_list_P2_chr1[[i]] <- finish_linkage_analysis(marker_assignment = marker_assignments_chr1[[i]]$P2,
                                                dosage_matrix = filtered_data[[i]],
                                                target_parent = P2,
                                                other_parent = P1,
                                                convert_palindrome_markers = TRUE, # convert 3.1 markers
                                                ploidy = 4,
                                                pairing = "random",
                                                LG_number = 1)

	integrated.maplist_P1_chr1[[i]] <- MDSMap_from_list(all_linkages_list_P1_chr1[[i]],write_to_file = FALSE)

	integrated.maplist_P2_chr1[[i]] <- MDSMap_from_list(all_linkages_list_P2_chr1[[i]],write_to_file = FALSE)
}

linkages_chr1 <- list()
for(i in 1:7){
  linkages_chr1[[i]] <- rbind(all_linkages_list_P1_chr1[[i]][1]$LG1,all_linkages_list_P2_chr1[[i]][1]$LG1)
}
names(linkages_chr1)<-c("LG1","LG2","LG3","LG4","LG5","LG6","LG7")

integrated.maplist_chr1 <- MDSMap_from_list(linkages_chr1,write_to_file = FALSE)

integrated.maplist_P1_chr1new<-list()
for (i in 1:7){

	integrated.maplist_P1_chr1new[i]<-integrated.maplist_P1_chr1[[i]]

	}
names(integrated.maplist_P1_chr1new)<-c("LG1","LG2","LG3","LG4","LG5","LG6","LG7")

integrated.maplist_P2_chr1new<-list()
for (i in 1:7){

	integrated.maplist_P2_chr1new[i]<-integrated.maplist_P2_chr1[[i]]

	}
names(integrated.maplist_P2_chr1new)<-c("LG1","LG2","LG3","LG4","LG5","LG6","LG7")

##Plotting a map

png(paste("integrated.maplist",".png",sep=""),width=1600, height=900)

plot_map(maplist=integrated.maplist_chr1)

dev.off()

png(paste("integrated.maplist_P2",".png",sep=""),width=1600, height=900)



plot_map(maplist=integrated.maplist_P2_chr1new)

dev.off()

png(paste("integrated.maplist_P1",".png",sep=""),width=1600, height=900)

plot_map(maplist=integrated.maplist_P1_chr1new)

dev.off()


pdf(paste("integrated.maplist",".pdf",sep=""),width=16, height=9)


plot_map(maplist=integrated.maplist_chr1)

dev.off()

pdf(paste("integrated.maplist_P2",".pdf",sep=""),width=16, height=9)

plot_map(maplist=integrated.maplist_P2_chr1new)

dev.off()

pdf(paste("integrated.maplist_P1",".pdf",sep=""),width=16, height=9)

plot_map(maplist=integrated.maplist_P1_chr1new)

dev.off()

##KEEP THE RESULT

mappp<-do.call(rbind,unname(integrated.maplist_chr1))
write.csv(mappp,file="integrated.maplist.csv",quote=F,row.names=F)

mappp<-do.call(rbind,unname(integrated.maplist_P1_chr1new))
write.csv(mappp,file="integrated.maplist_P1.csv",quote=F,row.names=F)

mappp<-do.call(rbind,unname(integrated.maplist_P2_chr1new))
write.csv(mappp,file="integrated.maplist_P2.csv",quote=F,row.names=F)

save.image("step2_linkage.RData")


################################################################
#####################bin linkage analysis#######################
################################################################

dir.create(file.path(getwd(), "bin"),showWarnings = FALSE)
setwd(file.path(getwd(), "bin"))

all_linkages_list_P1_chr_bin<-list()
for (i in 1:7){
	all_linkages_list_P1_chr_bin[[i]]<-all_linkages_list_P1_chr1[[i]][1]$LG1
}
all_linkages_list_P2_chr_bin<-list()
for (i in 1:7){
	all_linkages_list_P2_chr_bin[[i]]<-all_linkages_list_P2_chr1[[i]][1]$LG1
}

genotype<-do.call(rbind,unname(filtered_data))

all_linkages_list_P1_bin<-list()

for (i in 1:7){

	bin<-marker_binning(dosage_matrix=genotype, linkage_df=all_linkages_list_P1_chr_bin[[i]],target_parent = P1, other_parent = P2)
	all_linkages_list_P1_bin[i]<-bin[1]
	rm(bin)
}
names(all_linkages_list_P1_bin)<-c("LG1","LG2","LG3","LG4","LG5","LG6","LG7")

#can use the lod_thresh and r_thresh to chage the number of bin marker

all_linkages_list_P2_bin<-list()

for (i in 1:7){

	bin<-marker_binning(dosage_matrix=genotype, linkage_df=all_linkages_list_P2_chr_bin[[i]],
		target_parent = P2, other_parent = P1)
	all_linkages_list_P2_bin[i]<-bin[1]
	rm(bin)
}

names(all_linkages_list_P2_bin)<-c("LG1","LG2","LG3","LG4","LG5","LG6","LG7")

linkages_bin <- list()

for(lg in names(all_linkages_list_P1_bin)){
  linkages_bin[[lg]] <- rbind(all_linkages_list_P1_bin[[lg]],all_linkages_list_P2_bin[[lg]])
}

integrated.maplist_P1_bin <- MDSMap_from_list(all_linkages_list_P1_bin,write_to_file = TRUE)

integrated.maplist_P2_bin <- MDSMap_from_list(all_linkages_list_P2_bin,write_to_file = TRUE)

integrated.maplist_bin <- MDSMap_from_list(linkages_bin,write_to_file = TRUE)

mapp1_bin<-do.call(rbind,unname(integrated.maplist_P1_bin))
write.csv(mapp1_bin,file="integrated.maplist_P1_bin.csv",quote=F,row.names=F)


mapp2_bin<-do.call(rbind,unname(integrated.maplist_P2_bin))
write.csv(mapp2_bin,file="integrated.maplist_P2_bin.csv",quote=F,row.names=F)


mappp_bin<-do.call(rbind,unname(integrated.maplist_bin))
write.csv(mappp_bin,file="integrated.maplist_bin.csv",quote=F,row.names=F)

save.image("step2_linkage_bin.RData")

escaptime=Sys.time()-times;
print("Done!");
print(escaptime)



