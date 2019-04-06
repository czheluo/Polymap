
chr<-unique(snp[,1])

nc<-NULL
for (i in 1:length(chr)){

 nc[i]<-length(snp[which(snp[,1] %in% chr[i]),1])

}


chr<-unique(indel[,1])

ncl<-NULL
for (i in 1:length(chr)){

 ncl[i]<-length(indel[which(indel[,1] %in% chr[i]),1])

}

snpname<-phased.maplist

for (i in 1: length(snpname)){
	snpname[[i]][,1]<-paste("LG",i,"_",paste("SNP",c(1:length(snpname[[i]][,1])),sep=""),sep="")
}


SNP_SSR<-list()

for (i in 1:length(phased.maplist)){

	SNP_SSR[[i]]<-ALL_dosages_a[which(phased.maplist[[i]][,1] %in% rownames(ALL_dosages_a)),]

}


capture.output(for (i in 0:45) for (j in 0:45) print(i/j),
               file = "foo.txt")

for (i in 1: length(SNP_SSR)){
	rownames(SNP_SSR[[i]])<-paste("LG",i,"_",paste("SNP",c(1:length(rownames(SNP_SSR[[i]]))),sep=""),sep="")
}

N<-dim(SNPALL)[1]
M<-dim(SNPALL)[1]


for (i in 1: length(SNP_SSR)){
	write.csv(SNP_SSR[[i]],file=paste("LGGEN",i,".csv",sep=""))
	}



nsnp<-NULL
for (i in 1:length(phased.maplist)) {
  nsnp[i]<-dim(phased.maplist[[i]])[1]
}
print(paste("Total Number of SNP Markers:",sum(nsnp)))


nsnp<-NULL
for (i in 1:length(phased.maplist_P1)) {
  nsnp[i]<-dim(phased.maplist_P1[[i]])[1]
}
print(paste("Total Number of SNP Markers:",sum(nsnp)))

nsnp<-NULL
for (i in 1:length(phased.maplist_P2)) {
  nsnp[i]<-dim(phased.maplist_P2[[i]])[1]
}
print(paste("Total Number of SNP Markers:",sum(nsnp)))



##write the name and map

load("polymap_tetra.RData")

SNP_SSR<-list()

for (i in 1:length(phased.maplist)){
	SNP_SSR[[i]]<-ALL_dosages_a[match(phased.maplist[[i]][,1],rownames(ALL_dosages_a)),]


}

names(SNP_SSR)<-paste("LG",c(1:7),sep="")

snpnameall<-phased.maplist

for (i in 1: length(snpnameall)){
	for (j in 1:length(snpnameall[[i]][,1])) {
		if (is.na(match(snpnameall[[i]][j,1],rownames(polymap_tetra[[1]])))) {
			next
	} else {
	new<-polymap_tetra[[2]][match(snpnameall[[i]][j,1],rownames(polymap_tetra[[1]])),]
			rownames(SNP_SSR[[i]])[j]<-paste("c",new[,1],"_",new[,2],sep="")
			snpnameall[[i]][j,1]<-paste("c",new[,1],"_",new[,2],sep="")
	}
	}
	snpnameall[[i]][,1]<-rownames(SNP_SSR[[i]])
}



for (i in 1: length(SNP_SSR)){
	write.csv(SNP_SSR[[i]],file=paste("LGGEN",i,".csv",sep=""))
	}

for (i in 1: length(SNP_SSR)){
	SNP_SSR[[i]]<-paste("LG",i,sep="")
	}



SNP_SSR<-list()

for (i in 1:length(phased.maplist)){
	for (j in 1:length(phased.maplist[[i]][,1])){
		SNP_SSR[[i]]<-ALL_dosages_a[which(phased.maplist[[i]][j,1] %in% rownames(ALL_dosages_a)),]
	}


}

names(SNP_SSR)<-paste("LG",c(1:7),sep="")

snpnameall<-phased.maplist

for (i in 1: length(snpnameall)){
	for (j in 1:length(snpnameall[[i]][,1])) {
		if (length(which(snpnameall[[i]][j,1] %in% rownames(polymap_tetra[[1]])))) {
			next
	} else {
	new<-polymap_tetra[[2]][which(snpnameall[[i]][j,1] %in% rownames(polymap_tetra[[1]])),]
			rownames(SNP_SSR[[i]])[j]<-paste("c",new[,1],"_",new[,2],sep="")
			#snpnameall[[i]][j,1]<-paste("c",new[,1],"_",new[,2],sep="")
	}
	}
	snpnameall[[i]][,1]<-rownames(SNP_SSR[[i]])

}



for (i in 1: length(SNP_SSR)){
	SNP_SSR[[i]][is.na(SNP_SSR[[i]])]<-9
	write.table(SNP_SSR[[i]],file=paste("LGGEN",i,".txt",sep=""),quote=F,row.names=T,col.names=T)
	}






