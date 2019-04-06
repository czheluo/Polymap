
library(polymapR)

map<-read.csv("poly10Xpos.csv",header=T)
ALL_dosages<-read.csv("poly10X.csv",header=T)

P1<-"F"
P2<-"M"

all_dosages_list<-list()

nmap<-unique(map$Chr)

for (i in 1:length(nmap)) {
  all_dosages_list[[i]]<-ALL_dosages[which(map$Chr %in% nmap[i]),]
}

##remove
for (i in 1:length(all_dosages_list)) {

  all_dosages_list[[i]]<-all_dosages_list[[i]][which(!is.na(all_dosages_list[[i]][,1])),]
  all_dosages_list[[i]]<-all_dosages_list[[i]][which(!is.na(all_dosages_list[[i]][,2])),]
  }

#str(all_dosages_list)

nsnp<-NULL
for (i in 1:length(all_dosages_list)) {
  nsnp[i]<-dim(all_dosages_list[[i]])[1]
}
print(paste("Total Number of SNP Markers:",sum(nsnp)))


##second step checkF1
##remove q mult=0

check<-list();F1checked<-list();

times<-Sys.time()
for (i in 1:length(all_dosages_list)) {

F1checked[[i]] <- checkF1(dosage_matrix = all_dosages_list[[i]],parent1 = P1,parent2 = P2,
                     F1 = colnames(all_dosages_list[[i]])[3:ncol(all_dosages_list[[i]])],
                     polysomic = TRUE, disomic = FALSE, mixed = FALSE, ploidy = 4)

check[[i]]<-all_dosages_list[[i]][-which(F1checked[[i]]$qall_mult==0),]
}

calF1check<-Sys.time()-times


nsnp<-NULL
for (i in 1:length(check)) {
  nsnp[i]<-dim(check[[i]])[1]
}
print(paste("Total Number of SNP Markers:",sum(nsnp)))

all_dosages_list_orig<-all_dosages_list
all_dosages_list<-check
#check[[7]]<-all_dosages_list[[7]]
####REMOVE THE MISSING AND DUPLICATE MARKERS

segregating_data<-list();screened_data2<-list()
screened_data3<-list();screened_data4<-list()
ALL_dosages<-list()
for (i in 1:length(all_dosages_list)) {
  #pq_before_convert <- parental_quantities(
   # dosage_matrix = ALL_dosages,parent1 = P1,parent2 = P2,
   # las = 2)

  segregating_data[[i]] <- convert_marker_dosages(
    dosage_matrix = all_dosages_list[[i]],parent1 = P1,parent2 = P2)


  #pq_after_convert <- parental_quantities(
   # dosage_matrix = segregating_data,parent1 = P1,parent2 = P2,las=2)

  ##here may didnot work just jump to the next
  ##remove the makers
  #screened_data <- screen_for_NA_values(dosage_matrix = segregating_data,
                                       # margin = 1, # margin 1 means markers
                                       # cutoff =  0.05,
                                       # print.removed = FALSE)
  screened_data2[[i]] <- screen_for_NA_values(dosage_matrix = segregating_data[[i]],
                                         cutoff = 0.1,
                                         margin = 2, #margin = 2 means columns
                                         print.removed = FALSE)
  screened_data3[[i]] <- screen_for_duplicate_individuals(dosage_matrix = screened_data2[[i]],
                                                     cutoff = 0.95,
                                                     plot_cor = T)
  screened_data4[[i]] <- screen_for_duplicate_markers(dosage_matrix = screened_data3[[i]])
  ALL_dosages[[i]] <- screened_data4[[i]]$filtered_dosage_matrix
}

all_dosages_list<-ALL_dosages
#rm(segregating_data,screened_data3,screened_data4,screened_data2,ALL_dosages)

nsnp<-NULL
for (i in 1:length(all_dosages_list)) {
  nsnp[i]<-dim(all_dosages_list[[i]])[1]
}
print(paste("Total Number of SNP Markers:",sum(nsnp)))

nsnp<-NULL
for (i in 1:length(filtered_data)) {
  nsnp[i]<-dim(filtered_data[[i]])[1]
}
print(paste("Total Number of SNP Markers:",sum(nsnp)))



### Simplex x nulliplex markers - defining chromosomes and homologues
##few "outliers" in the coupling pairwise data (identified by the stars)
##from the expected relationship between r and LOD


SN_SS_P1<-list();SN_SS_P2<-list()
SN_SS_P1_1<-list();SN_SS_P2_1<-list()
filtered_data<-list();P1deviations<-list()
P2deviations<-list()
LGHomDf_P1<-list();LGHomDf_P2<-list()

times<-Sys.time()
for (i in 1:length(all_dosages_list)) {

  filtered_data[[i]]<-all_dosages_list[[i]]

  SN_SS_P1[[i]] <- linkage(dosage_matrix = filtered_data[[i]],
                    markertype1 = c(1,0),
                    markertype2 = c(1,1),
                    target_parent = P1,
                    other_parent = P2,
                    ploidy = 4,
                    pairing = "random")

  SxS_Marker<-unique(SN_SS_P1[[i]]$marker_a)

  LG<-rep(i,length(unique(SN_SS_P1[[i]]$marker_a)))

  homologue<-NULL
  for(j in 1:length(unique(SN_SS_P1[[i]]$marker_a))){
    homologue[j]<-sample(1:4,1)
  }

  #homologue<-rep(sample(1:4),length(unique(SN_SN_P1[[i]]$marker_a)))[1:length(unique(SN_SN_P1[[i]]$marker_a))]

  LGHomDf_P1[[i]]<-data.frame(SxS_Marker,LG,homologue)

  SN_SS_P2[[i]] <- linkage(dosage_matrix = filtered_data[[i]],
                           markertype1 = c(1,0),
                           markertype2 = c(1,1),
                           target_parent = P2,
                           other_parent = P1,
                           ploidy = 4,###change the ploidy type
                           pairing = "random")

  SxS_Marker<-unique(SN_SS_P2[[i]]$marker_a)

  LG<-rep(i,length(unique(SN_SS_P2[[i]]$marker_a)))

  homologue<-NULL
  for(j in 1:length(unique(SN_SS_P2[[i]]$marker_a))){
    homologue[j]<-sample(1:4,1)
  }

  #[1:length(unique(SN_SN_P1[[i]]$marker_a))]

  LGHomDf_P2[[i]]<-data.frame(SxS_Marker,LG,homologue)

  ###REMOVE A few "outliers"
  #PL1<-plot_linkage_df(SN_SN_P1[[i]], r_max = 0.5)
  #PL2<-plot_linkage_df(SN_SN_P2[[i]], r_max = 0.5)
  #png(PL1,paste("SN_SN_P1",i,".png",sep = ''))
  #dev.off()
  #pdf(PL1,paste("SN_SN_P1",i,".pdf",sep = ''))
  #dev.off()

  #png(paste("SN_SN_P2",i,".png",sep = ''))
  #dev.off()
 # pdf(PL2,paste("SN_SN_P2",i,".pdf",sep = ''))
  #dev.off()

  #png(paste("SN_SN_LOD_DE_P1",i,".png",sep = ''))
  #P1deviations[[i]] <- SNSN_LOD_deviations(linkage_df = SN_SN_P1[[i]],
                                     # ploidy = 4,
                                    #  N = ncol(filtered_data[[i]]) - 2, #The F1 population size
                                    #  alpha = c(0.05,0.2),
                                    #  plot_expected = TRUE,
                                    #  phase="coupling")

  #dev.off()

  #png(paste("SN_SN_LOD_DE_P2",i,".png",sep = ''))
  #P2deviations[[i]] <- SNSN_LOD_deviations(linkage_df = SN_SN_P2[[i]],
                                      #     ploidy = 4,
                                      #     N = ncol(filtered_data[[i]]) - 2, #The F1 population size
                                      #     alpha = c(0.05,0.2),
                                      #     plot_expected = TRUE,
                                     #      phase="coupling")
  #dev.off()
  ###KEEP THEM TO THE OTHER RESULT AND SAVE
  #SN_SN_P1_1[[i]] <- SN_SN_P1[[i]][SN_SN_P1$phase == "coupling",][-which(P1deviations[[i]] > 0.2),]
  #SN_SN_P2_1[[i]] <- SN_SN_P2[[i]][SN_SN_P2$phase == "coupling",][-which(P2deviations[[i]] > 0.2),]

}

calculate_time<-Sys.time()-times


SN_SN_P1<-list()
SN_SN_P2<-list()
SN_SN_P1_1<-list()
SN_SN_P2_1<-list()
filtered_data<-list()
P1deviations<-list()
P2deviations<-list()
LGHomDf_P1<-list()
LGHomDf_P2<-list()

times<-Sys.time()
for (i in 1:length(all_dosages_list)) {

  filtered_data[[i]]<-all_dosages_list[[i]]

  SN_SN_P1[[i]] <- linkage(dosage_matrix = filtered_data[[i]],
                           markertype1 = c(1,0),
                           #markertype2 = c(2,0),
                           target_parent = P1,
                           other_parent = P2,
                           ploidy = 4,
                           pairing = "random")

  SxN_Marker<-unique(SN_SN_P1[[i]]$marker_a)

  LG<-rep(i,length(unique(SN_SN_P1[[i]]$marker_a)))

  homologue<-NULL
  for(j in 1:length(unique(SN_SN_P1[[i]]$marker_a))){
    homologue[j]<-sample(1:4,1)
  }

  #homologue<-rep(sample(1:4),length(unique(SN_SN_P1[[i]]$marker_a)))[1:length(unique(SN_SN_P1[[i]]$marker_a))]

  LGHomDf_P1[[i]]<-data.frame(SxN_Marker,LG,homologue)

  SN_SN_P2[[i]]<-linkage(dosage_matrix = filtered_data[[i]],
                           markertype1 = c(1,0),
                           #markertype2 = c(2,0),
                           target_parent = P2,
                           other_parent = P1,
                           ploidy = 4,###change the ploidy type
                           pairing = "random")

  SxN_Marker<-unique(SN_SN_P2[[i]]$marker_a)

  LG<-rep(i,length(unique(SN_SN_P2[[i]]$marker_a)))

  homologue<-NULL
  for(j in 1:length(unique(SN_SN_P2[[i]]$marker_a))){
    homologue[j]<-sample(1:4,1)
  }

  #[1:length(unique(SN_SN_P1[[i]]$marker_a))]

  LGHomDf_P2[[i]]<-data.frame(SxN_Marker,LG,homologue)

  ###REMOVE A few "outliers"
  #PL1<-plot_linkage_df(SN_SN_P1[[i]], r_max = 0.5)
  #PL2<-plot_linkage_df(SN_SN_P2[[i]], r_max = 0.5)
  #png(PL1,paste("SN_SN_P1",i,".png",sep = ''))
  #dev.off()
  #pdf(PL1,paste("SN_SN_P1",i,".pdf",sep = ''))
  #dev.off()

  #png(paste("SN_SN_P2",i,".png",sep = ''))
  #dev.off()
  # pdf(PL2,paste("SN_SN_P2",i,".pdf",sep = ''))
  #dev.off()

  #png(paste("SN_SN_LOD_DE_P1",i,".png",sep = ''))
  #P1deviations[[i]] <- SNSN_LOD_deviations(linkage_df = SN_SN_P1[[i]],
  # ploidy = 4,
  #  N = ncol(filtered_data[[i]]) - 2, #The F1 population size
  #  alpha = c(0.05,0.2),
  #  plot_expected = TRUE,
  #  phase="coupling")

  #dev.off()

  #png(paste("SN_SN_LOD_DE_P2",i,".png",sep = ''))
  #P2deviations[[i]] <- SNSN_LOD_deviations(linkage_df = SN_SN_P2[[i]],
  #     ploidy = 4,
  #     N = ncol(filtered_data[[i]]) - 2, #The F1 population size
  #     alpha = c(0.05,0.2),
  #     plot_expected = TRUE,
  #      phase="coupling")
  #dev.off()
  ###KEEP THEM TO THE OTHER RESULT AND SAVE
  #SN_SN_P1_1[[i]] <- SN_SN_P1[[i]][SN_SN_P1$phase == "coupling",][-which(P1deviations[[i]] > 0.2),]
  #SN_SN_P2_1[[i]] <- SN_SN_P2[[i]][SN_SN_P2$phase == "coupling",][-which(P2deviations[[i]] > 0.2),]
}


calculate_time<-Sys.time()-times

#SN_SN_P1<-SN_SN_P1_1;SN_SN_P2<-SN_SN_P2_1
#rm(SN_SN_P1_1,SN_SN_P2_1)

#names(SN_SN_P1) <- c("LG1", "LG2", "LG3","LG4","LG5","LG6","LG7")
#names(SN_SN_P2) <- c("LG1", "LG2", "LG3","LG4","LG5","LG6","LG7")



time<-Sys.time()
for (i in 1: 7){
	LGHomDf_P1[[i]]$LG<-c(matrix(1,length(LGHomDf_P1[[i]]$LG)))

}
for (i in 1: 7){
	LGHomDf_P2[[i]]$LG<-c(matrix(1,length(LGHomDf_P2[[i]]$LG)))

}

save.image("step1_filter.RData")




time<-Sys.time()

SN_SS_P1 <- linkage(dosage_matrix = ALL_dosages_a,
                    markertype1 = c(1,0),
                    markertype2 = c(1,1),
                    target_parent = P1,
                    other_parent = P2,
                    ploidy = 4,
                    pairing = "random")

caltim<-Sys.time()-times
####13 mintues

P1_SxS_Assigned <- assign_linkage_group(linkage_df = SN_SS_P1,
                                        LG_hom_stack = SN_LGHomDf_P1_2,
                                        SN_colname = "marker_a",
                                        unassigned_marker_name = "marker_b",
                                        phase_considered = "coupling",
                                        LG_number = length(nmap),
                                        LOD_threshold = 3,
                                        ploidy = 4)



###

SN_DN_P1 <- linkage(dosage_matrix = ALL_dosages_a,
                    markertype1 = c(1,0),
                    markertype2 = c(2,0),
                    target_parent = P1,
                    other_parent = P2,
                    ploidy = 4,
                    pairing = "random")


P1_DxN_Assigned <- assign_linkage_group(linkage_df = SN_DN_P1,
                                        LG_hom_stack = SN_LGHomDf_P1_2,
                                        SN_colname = "marker_a",
                                        unassigned_marker_name = "marker_b",
                                        phase_considered = "coupling",
                                        LG_number = length(nmap),
                                        LOD_threshold = 3,
                                        ploidy = 4)




SN_SS_P2 <- linkage(dosage_matrix = ALL_dosages_a,
                    markertype1 = c(1,0),
                    markertype2 = c(1,1),
                    target_parent = P2,
                    other_parent = P1,
                    ploidy = 4,
                    pairing = "random")


P2_SxS_Assigned <- assign_linkage_group(linkage_df = SN_SS_P2,
                                        LG_hom_stack = SN_LGHomDf_P2_2,
                                        SN_colname = "marker_a",
                                        unassigned_marker_name = "marker_b",
                                        phase_considered = "coupling",
                                        LG_number = length(nmap),
                                        LOD_threshold = 3,
                                        ploidy = 4)


SN_DN_P2 <- linkage(dosage_matrix = ALL_dosages_a,
                    markertype1 = c(1,0),
                    markertype2 = c(2,0),
                    target_parent = P2,
                    other_parent = P1,
                    ploidy = 4,
                    pairing = "random")


P2_DxN_Assigned <- assign_linkage_group(linkage_df = SN_DN_P2,
                                        LG_hom_stack = SN_LGHomDf_P2_2,
                                        SN_colname = "marker_a",
                                        unassigned_marker_name = "marker_b",
                                        phase_considered = "coupling",
                                        LG_number = length(nmap),
                                        LOD_threshold = 3,
                                        ploidy = 4)


times<-Sys.time()
marker_assignments_P1 <- homologue_lg_assignment(dosage_matrix = ALL_dosages_a,
                                                 assigned_list = list(P1_SxS_Assigned,
                                                                      P1_DxN_Assigned),
                                                 assigned_markertypes = list(c(1,1), c(2,0)),
                                                 LG_hom_stack = SN_LGHomDf_P1_2,
                                                 target_parent = P1,
                                                 other_parent = P2,
                                                 ploidy = 4,
                                                 pairing = "random",
                                                 convert_palindrome_markers = FALSE,
                                                 LG_number = length(nmap),
                                                 LOD_threshold = 3,
                                                 write_intermediate_files = FALSE
)

markercalcu<-Sys.time()-times


times<-Sys.time()

marker_assignments_P2 <-homologue_lg_assignment(dosage_matrix = ALL_dosages_a,
                           assigned_list = list(P2_SxS_Assigned,P2_DxN_Assigned),
                           assigned_markertypes = list(c(1,1), c(2,0)),
                           LG_hom_stack = SN_LGHomDf_P2_2,
                           target_parent = P2,
                           other_parent = P1,
                           ploidy = 4,
                           pairing = "random",
                           convert_palindrome_markers = TRUE,
                           LG_number = length(nmap),
                           LOD_threshold = 3,
                           write_intermediate_files = FALSE
)

markercalcu<-Sys.time()-times


marker_assignments <- check_marker_assignment(marker_assignments_P1,marker_assignments_P2)


##finish the linkage analysis

times<-Sys.time()
all_linkages_list_P1 <- finish_linkage_analysis(marker_assignment = marker_assignments$P1,
                                                dosage_matrix = ALL_dosages_a,
                                                target_parent = P1,
                                                other_parent = P2,
                                                convert_palindrome_markers = FALSE,
                                                ploidy = 4,
                                                pairing = "random",
                                                LG_number = 7)

cal_linkage<-Sys.time()-times


times<-Sys.time()

all_linkages_list_P2 <- finish_linkage_analysis(marker_assignment = marker_assignments$P2,
                                                dosage_matrix = ALL_dosages_a,
                                                target_parent = P2,
                                                other_parent = P1,
                                                convert_palindrome_markers = TRUE, # convert 3.1 markers
                                                ploidy = 4,
                                                pairing = "random",
                                                LG_number = 7)


cal_linkage<-Sys.time()-times



integrated.maplist <- MDSMap_from_list(linkages)



pq_before_convert <- parental_quantities(dosage_matrix = ALL_dosages_a,
                                         parent1 = P1,parent2 = P2,
                                         las = 2)



####split all homologes

all_linkages_list_P1_split <- split_linkage_info(all_linkages = all_linkages_list_P1,
                                                 marker_assignment = marker_assignments_P1,
                                                 ploidy = 4)

all_linkages_list_P2_split <- split_linkage_info(all_linkages = all_linkages_list_P2,
                                                 marker_assignment = marker_assignments_P2,
                                                 ploidy = 4)


##make the bin
all_linkages_list_P1.binned <- marker_binning_list(dosage_matrix = ALL_dosages_a,
                                            linkage_list = all_linkages_list_P1_split,
                                            target_parent = "P1",
                                            other_parent = "P2",
                                            return_removed_marker_info = TRUE)


all_linkages_list_P2.binned <- marker_binning_list(dosage_matrix = ALL_dosages_a,
                                            linkage_list = all_linkages_list_P2_split,
                                            target_parent = "P1",
                                            other_parent = "P2",
                                            return_removed_marker_info = TRUE)



##maker order(more time )

#maplist_P1_LG1 <- MDSMap_from_list(all_linkages_list_P1_split$LG1,write_to_file = FALSE)

#maplist_P2_LG1 <- MDSMap_from_list(all_linkages_list_P2_split$LG1,write_to_file = FALSE)




##Optional: Independent homologue map integration
library("LPmerge")

integrated_map_LG2 <- orient_and_merge_maps(maplist_P1_LG1,
                                            maplist_P2_LG1,
                                            connection_threshold = 3,
                                            plot_graph=TRUE)






#This results in a nested list, which you can write to files:
str(all_linkages_list_P1_split)
write_nested_list(nested_list = all_linkages_list_P1, directory = "linkages_P1")


#names(SN_SS_P1) <- c("LG1", "LG2", "LG3","LG4","LG5","LG6","LG7")
#names(SN_SS_P2) <- c("LG1", "LG2", "LG3","LG4","LG5","LG6","LG7")




###use the bin to calculate

##one for the homelogues

all_linkages_list_P1.bin<-list()
for (i in 1:7){
	all_linkages_list_P1.bin[i]<-all_linkages_list_P1.binned[[1]][i]
	}
names(all_linkages_list_P1.bin)<-c("LG1","LG2","LG3","LG4","LG5","LG6","LG7")



##two for the linkages map

all_linkages_list_P1_bin<-list()

for (i in 1:7){

	bin<-marker_binning(dosage_matrix=ALL_dosages_a, linkage_df=all_linkages_list_P1[[i]],
		target_parent = "P1", other_parent = "P2")
	all_linkages_list_P1_bin[i]<-bin[1]
	rm(bin)


}

names(all_linkages_list_P1_bin)<-c("LG1","LG2","LG3","LG4","LG5","LG6","LG7")

all_linkages_list_P2_bin<-list()

for (i in 1:7){

	bin<-marker_binning(dosage_matrix=ALL_dosages_a, linkage_df=all_linkages_list_P2[[i]],
		target_parent = "P2", other_parent = "P1")
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




#phased.maplist bin

phased.maplist_bin <- create_phased_maplist(maplist = integrated.maplist_bin,
                                        dosage_matrix.conv = ALL_dosages_a,
                                        N_linkages = 7,
                                        ploidy = 4,
                                        marker_assignment.1 = marker_assignments$P1,
                                        marker_assignment.2 = marker_assignments$P2)




phased.maplist_P1_bin <- create_phased_maplist(maplist = integrated.maplist_P1_bin,
                                        dosage_matrix.conv = ALL_dosages_a,
                                        N_linkages = 7,
                                        ploidy = 4,
                                        marker_assignment.1 = marker_assignments$P1,
                                        marker_assignment.2 = marker_assignments$P2)



phased.maplist_P2_bin <- create_phased_maplist(maplist = integrated.maplist_P2,
                                        dosage_matrix.conv = ALL_dosages_a,
                                        N_linkages = 7,
                                        ploidy = 4,
                                        marker_assignment.1 = marker_assignments$P1,
                                        marker_assignment.2 = marker_assignments$P2)





#Creating an integrated chromosomal linkage map
#more time

times<-Sys.time()

integrated.maplist_P1 <- MDSMap_from_list(all_linkages_list_P1,write_to_file = TRUE)

cal_linkage<-Sys.time()-times


times<-Sys.time()

integrated.maplist_P2 <- MDSMap_from_list(all_linkages_list_P2,write_to_file = TRUE)
`
cal_linkage<-Sys.time()-times

times<-Sys.time()-times


linkages <- list()
for(lg in names(all_linkages_list_P1)){
  linkages[[lg]] <- rbind(all_linkages_list_P1[[lg]], all_linkages_list_P2[[lg]])
}



times<-Sys.time()-times

integrated.maplist <- MDSMap_from_list(linkages,write_to_file = FALSE)

cal_linkageall<-Sys.time()-times

#phased.maplist

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


##output TetraploidSNPMap

#load("polymap_tetra.RData")


SNP_SSR<-list()

for (i in 1:length(phased.maplist)){

	SNP_SSR[[i]]<-ALL_dosages_a[pmatch(phased.maplist[[i]][,1],rownames(ALL_dosages_a),dup = FALSE),]

}


names(SNP_SSR)<-paste("LG",c(1:7),sep="")

snpnameall<-phased.maplist

for (i in 1: length(snpnameall)){
	for (j in 1:length(snpnameall[[i]][,1])) {
		if (is.na(pmatch(snpnameall[[i]][j,1],rownames(polymap_tetra[[1]]),dup = FALSE))) {
			next
	} else {
	new<-polymap_tetra[[2]][pmatch(snpnameall[[i]][j,1],rownames(polymap_tetra[[1]]),dup = FALSE),]
			rownames(SNP_SSR[[i]])[j]<-paste("c",new[,1],"_",new[,2],sep="")
			snpnameall[[i]][j,1]<-paste("c",new[,1],"_",new[,2],sep="")


	}

	}

	SNP_SSR[[i]]<-

	snpnameall[[i]][,1]<-rownames(SNP_SSR[[i]])

}



SNP_SSR<-list()

for (i in 1:length(phased.maplist_P1)){
	SNP_SSR[[i]]<-ALL_dosages_a[which(phased.maplist_P1[[i]][,1] %in% rownames(ALL_dosages_a)),]

}


for (i in 1: length(SNP_SSR)){
	rownames(SNP_SSR[[i]])<-paste("LG",i,"_",paste("SNP",c(1:length(rownames(SNP_SSR[[i]]))),sep=""),sep="")
}


for (i in 1: length(SNP_SSR)){
	SNP_SSR[[i]][is.na(SNP_SSR[[i]])]<-9
	write.table(SNP_SSR[[i]],file=paste("LGGEN",i,".txt",sep=""),quote=F,row.names=T,col.names=T)
	}

snpmap<-phased.maplist_P1
for (i in 1: length(SNP_SSR)){
	snpmap[[i]]$marker<-rownames(SNP_SSR[[i]])
}


write.TSNPM(phased.maplist = snpmap,ploidy=4)



write.TSNPM(phased.maplist = phased.maplist_P1,ploidy=4)

write.TSNPM(phased.maplist = phased.maplist_P2,ploidy=4)





#plot phased the map plot
plot_phased_maplist(phased.maplist = phased.maplist[1], #Can plot full list also, remove "[1]"
                    ploidy = 4,
                    cols = c("black","grey50","grey50"))




##Plotting a map

png(paste("integrated.maplist_P2",".png",sep=""),width=1600, height=900)

plot_map(maplist=integrated.maplist_P2)

dev.off()

pdf(paste("integrated.maplist_P2",".pdf",sep=""),width=16, height=9)

plot_map(maplist=integrated.maplist_P2)

dev.off()



png(paste("integrated.maplist_P1",".png",sep=""),width=1600, height=900)

plot_map(maplist=integrated.maplist_P1)

dev.off()

pdf(paste("integrated.maplist_P1",".pdf",sep=""),width=16, height=9)

plot_map(maplist=integrated.maplist_P1)

dev.off()




png(paste("integrated.maplist_all_new",".png",sep=""),width=1600, height=900)

plot_map(maplist=integrated.maplist)

dev.off()

pdf(paste("integrated.maplist_all_new",".pdf",sep=""),width=16, height=9)

plot_map(maplist=integrated.maplist)

dev.off()


###plot.phased.maplist


plot_phased_maplist(phased.maplist = phased.maplist[1], #Can plot full list also, remove "[1]"
                    ploidy = 4,
                    cols = c("black","grey50","grey50"))





#remove the unnormal data
#integrated.maplist_P2[[2]][,2]<-integrated.maplist_P2[[2]][,2]-integrated.maplist_P2[[2]][1,2]
#integrated.maplist[[2]][,2]<-integrated.maplist[[2]][,2]-integrated.maplist[[2]][2,2]


##Evaluating map quality

for (i in 1:length(SN_SN_P)) {

  check_map(linkage_list = list(SN_SN_P[[i]]), maplist = list(integrated.maplist[[i]]))

}







###make the collinerity

snpnameall<-maplist_bin

for (i in 1: length(snpnameall)){
	for (j in 1:length(snpnameall[[i]][,1])) {
		if (is.na(pmatch(snpnameall[[i]][j,1],rownames(polymap_tetra[[1]]),dup = FALSE))) {
			next
	} else {
	new<-polymap_tetra[[2]][pmatch(snpnameall[[i]][j,1],rownames(polymap_tetra[[1]]),dup = FALSE),]
			#rownames(SNP_SSR[[i]])[j]<-paste("c",new[,1],"_",new[,2],sep="")
			snpnameall[[i]][j,1]<-paste("chr",new[,1],"_",new[,2],sep="")


	}

	}

	#SNP_SSR[[i]]<-

	#snpnameall[[i]][,1]<-rownames(SNP_SSR[[i]])

}








setwd("H:\\PROJECT\\wanhuihua\\RESULT")
load("example_maplist.RDdata")

linkage_list<-list(linkages$LG1);maplist<-lis(integrated.maplist$LG1)


load("all_example_maplist.RDdata")
mapfn <- match.arg(mapfn, choices = c("haldane", "kosambi"))
rev.haldane <- function(d) (1 - exp(-d/50))/2
rev.kosambi <- function(d) ((exp(d/25) - 1)/(exp(d/25) +
                                               1))/2
orig.mar <- c(5.1, 4.1, 4.1, 2.1)
colbar.mar <- c(5.1, 2, 4.1, 0.5)
posmat <- matrix(c(maplist[[l]][match(linkage_list[[l]]$marker_a,
                                      maplist[[l]]$marker), ]$position,
                   maplist[[l]][match(linkage_list[[l]]$marker_b,
                                      maplist[[l]]$marker), ]$position, linkage_list[[l]]$r,
                   linkage_list[[l]]$LOD), ncol = 4)



if (mapfn == "haldane") {
  expected.recom <- rev.haldane(abs(posmat[, 1] -
                                      posmat[, 2]))
}else if (mapfn == "kosambi") {
  expected.recom <- rev.kosambi(abs(posmat[, 1] -
                                      posmat[, 2]))
}


dev <- abs(linkage_list[[l]]$r - expected.recom)
wRMSE <- sqrt(mean((dev * linkage_list[[l]]$LOD)^2))
layout(matrix(c(1, 1, 1, 1, 2, 3, 4, 5), ncol = 4, byrow = TRUE),
       widths = c(1, 0.2, 1, 0.2))



pdf(paste("ALL_sample_LG1",".pdf",sep=""),width=16, height=9)

png(paste("ALL_SAMPLE_lg1",".png",sep=""),width=1600, height=900)

par(oma = c(0, 0, 3, 0))

png(paste("ALL_SAMPLE_LG1",".png",sep=""),width=1600, height=900)

plot(dev, linkage_list[[l]]$LOD, ylab = "LOD", xlab = expression(delta(r)),
     cex.lab = 1.25, main = expression("|r"["pairwise"] ~
                                         "- r"["map"] ~ "|"))
dev.off()

legend("topright", legend = paste("Weighted RMSE =",
                                  round(wRMSE, 3)), bty = "n")
colours <- colorRampPalette(c("green", "yellow", "orange",
                              "red"))(100)



lod.thresh<-5
##remove all the NAN data

posmat<-posmat[-which(is.na(posmat[,3])),]
rcolmin <- min(posmat[, 3])
rcolmax <- max(posmat[, 3])
LODcolmin <- lod.thresh
LODcolmax <- max(posmat[, 4])
rcolbreaks <- seq(rcolmin, rcolmax, (rcolmax - rcolmin)/100)

rcols <- colours[findInterval(x = posmat[, 3], vec = rcolbreaks)]
rcols[is.na(rcols)] <- colours[100]
LODcolbreaks <- seq(LODcolmin, LODcolmax, (LODcolmax -
                                             LODcolmin)/100)
LODcols <- colours[findInterval(x = posmat[, 4], vec = LODcolbreaks)]
LODcols[is.na(LODcols)] <- colours[100]


png(paste("ALL_SAMPLE_LG_LOD",".png",sep=""),width=1600, height=900)


plot(posmat[, 1], posmat[, 2], pch = ".", col = rcols,
     main = "r", cex = 3, xlab = "cM", ylab = "cM")



dev.off()


par(mar = colbar.mar)


colour.bar(col.data = colours, min = rcolmin, max = rcolmax,
           nticks = 8, ticks = round(seq(rcolmin, rcolmax,
          len = 8), 2), cex.ticks = 0.8)


par(mar = orig.mar)

png(paste("ALL_SAMPLE_LG_LOD",".png",sep=""),width=1600, height=900)

plot(posmat[posmat[, 4] > lod.thresh, 1], posmat[posmat[,
                                                        4] > lod.thresh, 2], pch = ".", col = LODcols, main = "LOD",
     cex = 3, xlab = "cM", ylab = "cM")

dev.off()

par(mar = colbar.mar)
colour.bar(col.data = colours, min = LODcolmin, max = LODcolmax,
           nticks = 8, ticks = round(seq(LODcolmin, LODcolmax,
                                         len = 8), 2), cex.ticks = 0.8)
mtext(text = paste("LG", l, "map diagnostics"), side = 3,
      outer = TRUE, cex = 2)
par(mfrow = c(1, 1), oma = c(0, 0, 0, 0), mar = orig.mar)







