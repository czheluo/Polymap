
t<-t(dosage)
tr<-tr[names(sort(tr[,1],decreasing=F)),]
all_dosage<-data.frame(dosage[,c(1,189,190],tr)
all_dosage<-data.frame(dosage[,c(1,189,190],tr)
P1<-"F"
P2<-"M"
ALL_dosages<-read.csv("dosage_10x_snp.csv",header=T)
map<-read.csv("dosage_10x_pos.csv",header=T)
rownames(ALL_dosages)<-ALL_dosages[,1]
ALL_dosages <- as.matrix(ALL_dosages[,-1])

all_dosages_list<-list()
nmap<-unique(map$Chr)
for (i in 1:length(nmap)) {
  all_dosages_list[[i]]<-ALL_dosages[which(map$Chr %in% nmap[i]),]
}
#ALL_dosages <- as.matrix(ALL_dosages)
#plot four different genotype

pq_before_convert <- parental_quantities(dosage_matrix = ALL_dosages,parent1 = P1,parent2 = P2, las=2)


time<-Sys.time()
for (i in 1: 7){
	LGHomDf_P1[[i]]$LG<-c(matrix(1,length(LGHomDf_P1[[i]]$LG)))

}
for (i in 1: 7){
	LGHomDf_P2[[i]]$LG<-c(matrix(1,length(LGHomDf_P2[[i]]$LG)))

}


times<-Sys.time()

##lod threshold for filter

lod<-1

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

SN_DN_P1_chr1[[i]] <- linkage(dosage_matrix = filtered_data[[i]],
                    markertype1 = c(1,0),
                    #markertype2 = c(2,0),
                    target_parent = P1,
                    other_parent = P2,
                    ploidy = 4,
                    pairing = "random")

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


SN_DN_P2_chr1[[i]] <- linkage(dosage_matrix = filtered_data[[i]],
                    markertype1 = c(1,0),
                    #markertype2 = c(2,0),
                    target_parent = P2,
                    other_parent = P1,
                    ploidy = 4,
                    pairing = "random")


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

#times<-Sys.time()

integrated.maplist_chr1 <- MDSMap_from_list(linkages_chr1,write_to_file = FALSE)

cal_linkageall<-Sys.time()-times


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

linkages_chr1 <- list()
for(i in 1:7){
  linkages_chr1[[i]] <- rbind(all_linkages_list_P1_chr1[[i]][1]$LG1,all_linkages_list_P2_chr1[[i]][1]$LG1)
}

names(linkages_chr1)<-c("LG1","LG2","LG3","LG4","LG5","LG6","LG7")

times<-Sys.time()

integrated.maplist_chr1 <- MDSMap_from_list(linkages_chr1,write_to_file = FALSE)

cal_linkageall<-Sys.time()-times

integrated.maplist_chr1$LG4<-integrated.maplist_chr1$LG4[,2]-integrated.maplist_chr1$LG4[2,2]

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

####keep the result


integrated.maplist_P1_split<-list()
for (i in 1:7){
integrated.maplist_P1_split[[i]]<-integrated.maplist_P1_chr1[[i]]$LG1
}

names(integrated.maplist_P1_split)<-c("LG1","LG2","LG3","LG4","LG5","LG6","LG7")
mapp1<-do.call(rbind,unname(integrated.maplist_chr1))
write.csv(mapp1,file="integrated.maplist_P1_split.csv",quote=F,row.names=F)


integrated.maplist_P2_split<-list()
for (i in 1:7){
integrated.maplist_P2_split[[i]]<-integrated.maplist_P2_chr1[[i]]$LG1
}

names(integrated.maplist_P2_split)<-c("LG1","LG2","LG3","LG4","LG5","LG6","LG7")
mapp2<-do.call(rbind,unname(integrated.maplist_P2_split))
write.csv(mapp2,file="integrated.maplist_P2_split.csv",quote=F,row.names=F)


names(integrated.maplist_P2_split)<-c("LG1","LG2","LG3","LG4","LG5","LG6","LG7")
mappp<-do.call(rbind,unname(integrated.maplist_chr1))
write.csv(mappp,file="integrated.maplist_chr1.csv",quote=F,row.names=F)


##KEEP THE RESULT

mappp<-do.call(rbind,unname(integrated.maplist_chr1))
write.csv(mappp,file="integrated.maplist_chr1.csv",quote=F,row.names=F)

mappp<-do.call(rbind,unname(integrated.maplist_P1_chr1new))
write.csv(mappp,file="integrated.maplist_P1_chr1.csv",quote=F,row.names=F)

mappp<-do.call(rbind,unname(integrated.maplist_P2_chr1new))
write.csv(mappp,file="integrated.maplist_P2_chr1.csv",quote=F,row.names=F)

###bin
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

times<-Sys.time()

integrated.maplist_P1_bin <- MDSMap_from_list(all_linkages_list_P1_bin,write_to_file = TRUE)

cal_linkage<-Sys.time()-times
times<-Sys.time()
integrated.maplist_P2_bin <- MDSMap_from_list(all_linkages_list_P2_bin,write_to_file = TRUE)
cal_linkage<-Sys.time()-times
times<-Sys.time()

integrated.maplist_bin <- MDSMap_from_list(linkages_bin,write_to_file = TRUE)

cal_linkage<-Sys.time()-times

mapp1_bin<-do.call(rbind,unname(integrated.maplist_P1_bin))
write.csv(mapp1_bin,file="integrated.maplist_P1_bin.csv",quote=F,row.names=F)
mapp2_bin<-do.call(rbind,unname(integrated.maplist_P2_bin))
write.csv(mapp2_bin,file="integrated.maplist_P2_bin.csv",quote=F,row.names=F)
mappp_bin<-do.call(rbind,unname(integrated.maplist_bin))
write.csv(mappp_bin,file="integrated.maplist_bin.csv",quote=F,row.names=F)

