#!/usr/bin/env Rscript
library(<p></p>'getopt');
options(bitmapType='cairo')
spec = matrix(c(
	'out','o',1,'character',
	'pid','p',1,'character',
	'mid','m',1,'character',
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
	--help		usage
\n")
	q(status=1);
}
if (is.null(opt$pid)) { print_usage(spec)}
if (is.null(opt$mid)){ print_usage(spec) }
times<-Sys.time()

library(polymapR)

setwd(paste(opu$out,"01.vcf-convert/",sep=""))

ALL_dosages<-read.csv("dosage.matrix",header=T)

map<-read.csv("dosage.matrix.map",header=T)

P1<-opt$pid
P2<-opt$mid

all_dosages_list<-list()

nmap<-unique(map$Chr)
for (i in 1:length(nmap)) {
  all_dosages_list[[i]]<-ALL_dosages[which(map$Chr %in% nmap[i]),]
}

##remove missing value from parents

for (i in 1:length(all_dosages_list)) {

  all_dosages_list[[i]]<-all_dosages_list[[i]][which(!is.na(all_dosages_list[[i]][,P1])),]
  all_dosages_list[[i]]<-all_dosages_list[[i]][which(!is.na(all_dosages_list[[i]][,P2])),]
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

for (i in 1:length(all_dosages_list)) {

F1checked[[i]] <- checkF1(dosage_matrix = all_dosages_list[[i]],parent1 = P1,parent2 = P2,
                     F1 = colnames(all_dosages_list[[i]])[3:ncol(all_dosages_list[[i]])],
                     polysomic = TRUE, disomic = FALSE, mixed = FALSE, ploidy = 4)
check[[i]]<-all_dosages_list[[i]][-which(F1checked[[i]]$qall_mult==0),]
}

nsnp<-NULL
for (i in 1:length(check)) {
  nsnp[i]<-dim(check[[i]])[1]
}
print(paste("Total Number of SNP Markers:",sum(nsnp)))
#keep the original data
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

##plot the segregation summary barplot for the type of dosage matrix
library(RColorBrewer)
color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
rcolor<-color[sample(1:length(color),length(color))]

ALL_dosages<-do.call(rbind,unname(filtered_data))
pq_before_convert <- parental_quantities(dosage_matrix = ALL_dosages,parent1 = P1,parent2 = P2, las=2)

png(paste("ALL segregation summary",".png",sep=""),width=1600, height=900)
par(mfrow=c(1,1))
col<-length(pq_before_convert)
cor<-rcolor[1:col]
bp<-barplot(pq_before_convert,#,space = 0,
        col=cor,border = NULL,ylab = "Nr.Markers",
        main = "Marker segregation summary",cex.names =1.2,cex.axis = 1, cex.lab=1.2)
mtext(paste("Total number markers:", sum(pq_before_convert)),3, font = 3)
text(x=bp,y=pq_before_convert,labels = round(pq_before_convert,0),pos=3,xpd=NA)

dev.off()

pdf(paste("ALL segregation summary",".pdf",sep=""),width=16, height=9)

par(mfrow=c(1,1))
col<-length(pq_before_convert)
cor<-rcolor[1:col]
bp<-barplot(pq_before_convert,#,space = 0,
        col=cor,border = NULL,ylab = "Nr.Markers",
        main = "Marker segregation summary",cex.names =1.2,cex.axis = 1, cex.lab=1.2)
mtext(paste("Total number markers:", sum(pq_before_convert)),3, font = 3)

text(x=bp,y=pq_before_convert,labels = round(pq_before_convert,0),pos=3,xpd=NA)

dev.off()

filter<-do.call(rbind,unname(filtered_data))
pq_before_convert <- parental_quantities(dosage_matrix = filter,parent1 = P1,parent2 = P2, las=2)

png(paste("filtered segregation summary",".png",sep=""),width=1600, height=900)

par(mfrow=c(1,1))
col<-length(pq_before_convert)
cor<-rcolor[1:col]
bp<-barplot(pq_before_convert,#,space = 0,
        col=cor,border = NULL,ylab = "Nr.Markers",
        main = "Marker segregation summary",cex.names =1.2,cex.axis = 1, cex.lab=1.2)
mtext(paste("Total number markers:", sum(pq_before_convert)),3, font = 3)
text(x=bp,y=pq_before_convert,labels = round(pq_before_convert,0),pos=3,xpd=NA)

dev.off()

pdf(paste("filtered segregation summary",".pdf",sep=""),width=16, height=9)

par(mfrow=c(1,1))
col<-length(pq_before_convert)
cor<-rcolor[1:col]
bp<-barplot(pq_before_convert,#,space = 0,
        col=cor,border = NULL,ylab = "Nr.Markers",
        main = "Marker segregation summary",cex.names =1.2,cex.axis = 1, cex.lab=1.2)
mtext(paste("Total number markers:", sum(pq_before_convert)),3, font = 3)

text(x=bp,y=pq_before_convert,labels = round(pq_before_convert,0),pos=3,xpd=NA)

dev.off()



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

save.image("RData/step1_filter.RData")

escaptime=Sys.time()-times;
print("Done!");
print(escaptime)

