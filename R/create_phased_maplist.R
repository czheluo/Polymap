
create_phased_maplist <- function(maplist,
                                  dosage_matrix.conv,
                                  dosage_matrix.orig = NULL,
                                  remove_markers = NULL,
                                  N_linkages = 2,
                                  lower_bound = 0.05,
                                  ploidy = 4,
                                  ploidy2 = NULL,
                                  marker_assignment.1,
                                  marker_assignment.2,
                                  original_coding = FALSE,
                                  log = NULL,
				  P1=P1,
				  P2=P2,
                                  verbose = TRUE) {
vector.to.matrix <- function(x, n.columns){
  if(length(x)>n.columns){
    x<-c(x, rep("", n.columns-length(x)%%n.columns))
  } else {
    n.columns <- length(x)
  }
  x.m <- matrix(x, ncol=n.columns, byrow=T)
  colnames(x.m)<-rep("_", n.columns)
  return(x.m)
}

  if(original_coding & is.null(dosage_matrix.orig)) stop("Uncoverted dosage matrix should also be specified if original_coding = TRUE")

  if (is.null(log)) {
    log.conn <- stdout()
  } else {
    matc <- match.call()
    write.logheader(matc, log)
    log.conn <- file(log, "a")
  }

  if(is.null(ploidy2)) ploidy2 <- ploidy

  if(ploidy == ploidy2){
    palindromes <- rownames(dosage_matrix.conv)[which(dosage_matrix.conv[,paste(P1)] != dosage_matrix.conv[,paste(P2)] &
                                                        abs(dosage_matrix.conv[,paste(P1)] - (0.5*ploidy)) == abs(dosage_matrix.conv[,paste(P2)]-(0.5*ploidy2)))]

    ## If there are any unconverted palindromes, convert them:
    if(any(dosage_matrix.conv[palindromes,paste(P1)] > dosage_matrix.conv[palindromes,paste(P2)]))
      dosage_matrix.conv[palindromes[dosage_matrix.conv[palindromes,paste(P1)] > dosage_matrix.conv[palindromes,paste(P2)]],] <-
        ploidy - dosage_matrix.conv[palindromes[dosage_matrix.conv[palindromes,paste(P1)] > dosage_matrix.conv[palindromes,paste(P2)]],]
    }

  # Begin by separating the SxN and NxS linkages:
  SxN_assigned <- marker_assignment.1[marker_assignment.1[,paste(P1)]==1 &
                                        marker_assignment.1[,paste(P2)]==0,]
  p1_assigned <- marker_assignment.1[-match(rownames(SxN_assigned),rownames(marker_assignment.1)),]

  NxS_assigned <- marker_assignment.2[marker_assignment.2[,paste(P1)]==1 &
                                        marker_assignment.2[,paste(P2)]==1,]
  p2_assigned <- marker_assignment.2[-match(rownames(NxS_assigned),rownames(marker_assignment.2)),]

  #Use only the markers with at least N_linkages significant linkages
  P1unlinked <- rownames(p1_assigned)[apply(p1_assigned[,3+grep("LG",colnames(p1_assigned)[4:ncol(p1_assigned)]),drop = FALSE],1,max)<N_linkages]
  P2unlinked <- rownames(p2_assigned)[apply(p2_assigned[,3+grep("LG",colnames(p2_assigned)[4:ncol(p2_assigned)]),drop = FALSE],1,max)<N_linkages]

  if(verbose) {
    removed.m1 <- vector.to.matrix(P1unlinked, n.columns = 4)
    removed.m2 <- vector.to.matrix(P2unlinked, n.columns = 4)

    if(nrow(removed.m1) > 0){
      write(paste("\nThe following P1 markers had less than", N_linkages,"significant linkages:\n_______________________________________\n"),log.conn)
      write(knitr::kable(removed.m1,format="markdown"), log.conn)
    }

    if(nrow(removed.m2) > 0){
      write(paste("\n\nThe following P2 markers had less than", N_linkages,"significant linkages:\n_______________________________________\n"),log.conn)
      write(knitr::kable(removed.m2,format="markdown"), log.conn)
      write("\n", log.conn)
    }
  }

  if(length(P1unlinked) > 0) p1_assigned <- p1_assigned[-match(P1unlinked,rownames(p1_assigned)),]
  if(length(P2unlinked) > 0) p2_assigned <- p2_assigned[-match(P2unlinked,rownames(p2_assigned)),]

  # Only select markers for which the number of homologue assignments match the seg type:
  P1rates <- p1_assigned[,3+grep("Hom",colnames(p1_assigned)[4:ncol(p1_assigned)])]/rowSums(p1_assigned[,3+grep("Hom",colnames(p1_assigned)[4:ncol(p1_assigned)])])
  P2rates <- p2_assigned[,3+grep("Hom",colnames(p2_assigned)[4:ncol(p2_assigned)])]/rowSums(p2_assigned[,3+grep("Hom",colnames(p2_assigned)[4:ncol(p2_assigned)])])

  P1rates[P1rates < lower_bound] <- 0
  P2rates[P2rates < lower_bound] <- 0

  P1linked <- apply(P1rates,1,function(x) length(which(x!=0)))
  P2linked <- apply(P2rates,1,function(x) length(which(x!=0)))

  p1.markers <- rownames(p1_assigned[p1_assigned[,paste(P1)]!=0,])
  p2.markers <- rownames(p2_assigned[p2_assigned[,paste(P2)]!=0,])

  ## Assuming markers are converted here; have to treat palindrome markers in P2 carefully:
  P1different <- rownames(p1_assigned[rownames(p1_assigned) %in% p1.markers & p1_assigned[,paste(P1)] != P1linked,])
  P2different <- rownames(p2_assigned[setdiff(which(rownames(p2_assigned) %in% p2.markers & p2_assigned[,paste(P2)] != P2linked),
                                        which(rownames(p2_assigned) %in% palindromes & ploidy2 - p2_assigned[,paste(P2)] == P2linked)),])

  if(verbose) {

    removed.m1 <- if(!is.null(P1different)) {
      vector.to.matrix(P1different, n.columns = 4)
    } else matrix(,nrow=0,ncol=1) #catching error

    removed.m2 <- if(!is.null(P2different)){
      vector.to.matrix(P2different, n.columns = 4)
    } else matrix(,nrow=0,ncol=1) #catching error

    if(nrow(removed.m1) > 0){
      write(paste("\nThe following markers did not have the expected assignment in P1:\n_______________________________________\n"),log.conn)
      write(knitr::kable(removed.m1,format="markdown"), log.conn)
    }
    if(nrow(removed.m2) > 0){
      write(paste("\n\nThe following markers did not have the expected assignment in P2:\n_______________________________________\n"),log.conn)
      write(knitr::kable(removed.m2,format="markdown"), log.conn)
      write("\n", log.conn)
    }
  }

  p1_assigned <- p1_assigned[!rownames(p1_assigned) %in% P1different,]
  p2_assigned <- p2_assigned[!rownames(p2_assigned) %in% P2different,]

  # return simplex x nulliplex markers
  p1_assigned <- rbind(SxN_assigned,p1_assigned)
  p2_assigned <- rbind(NxS_assigned,p2_assigned)

  # Remove the bi-parental markers that are not assigned in both parents
  bip1 <- rownames(p1_assigned[rowSums(p1_assigned[,c(paste(P1),paste(P2))]!=0)==2,])
  bip2 <- rownames(p2_assigned[rowSums(p2_assigned[,c(paste(P1),paste(P2))]!=0)==2,])

  BiP_different <- c(setdiff(bip1,intersect(bip1,bip2)),setdiff(bip2,intersect(bip1,bip2)))

  if (verbose & !is.null(BiP_different)) {
    removed.m <- vector.to.matrix(BiP_different, n.columns = 4)

    if(nrow(removed.m) > 0){
      write(paste("\nThe following markers did not have the expected assignment across both parents:\n_______________________________________\n"),log.conn)
      write(knitr::kable(removed.m,format="markdown"), log.conn)
      write("\n", log.conn)
    }
  }

  p1_assigned <- p1_assigned[!rownames(p1_assigned) %in% setdiff(bip1,intersect(bip1,bip2)),]
  p2_assigned <- p2_assigned[!rownames(p2_assigned) %in% setdiff(bip2,intersect(bip1,bip2)),]

  ALL_assigned <- unique(c(rownames(p1_assigned),rownames(p2_assigned)))

  # Make up the output
  maplist.out <- lapply(seq(length(maplist)),function(mapn) {

    map <- maplist[[mapn]]
    map <- map[map$marker %in% ALL_assigned,]

    outmap <- map[,c("marker","position")]

    hom_mat <- sapply(1:nrow(outmap), function(r){
      a <- rep(0, ploidy+ploidy2)
      temp <- p1_assigned[match(outmap$marker[r],rownames(p1_assigned)),
                          3+grep("Hom",colnames(p1_assigned)[4:ncol(p1_assigned)])]

      if(length(which(temp!=0)) > 0) a[(1:ploidy)[which(temp!=0)]] <- 1

      temp <- p2_assigned[match(outmap$marker[r],rownames(p2_assigned)),
                          3+grep("Hom",colnames(p2_assigned)[4:ncol(p2_assigned)])]
      if(length(which(temp!=0)) > 0) a[((ploidy+1):(ploidy+ploidy2))[which(temp!=0)]] <- 1

      return(a)
    })

    hom_mat <- t(hom_mat)
    colnames(hom_mat) <- paste0("h",seq(1,ploidy+ploidy2))

    # correct palindrome markers:
    if(any(outmap$marker %in% palindromes)){
      hom_mat[outmap$marker %in% palindromes,(ploidy+1):(ploidy+ploidy2)] <-
        (hom_mat[outmap$marker %in% palindromes,(ploidy+1):(ploidy+ploidy2)] + 1) %% 2
    }

    # recode using the original coding:
    if(original_coding){

      orig_parents <- dosage_matrix.orig[match(outmap$marker,rownames(dosage_matrix.orig)),c(paste(P1),paste(P2))]
      orig_mat <- hom_mat

      for(r in 1:nrow(orig_mat)){
        if(sum(hom_mat[r,1:ploidy]) != orig_parents[r,1]) orig_mat[r,1:ploidy] <- (hom_mat[r,1:ploidy]+1)%%2
        if(sum(hom_mat[r,(ploidy+1):(ploidy+ploidy2)]) != orig_parents[r,2])
          orig_mat[r,(ploidy+1):(ploidy+ploidy2)] <- (hom_mat[r,(ploidy+1):(ploidy+ploidy2)]+1)%%2
      }
      outmap <- cbind(outmap,orig_mat)
    } else{
      outmap <- cbind(outmap,hom_mat)
    }

    return(outmap)
  }
  )
  names(maplist.out) <- names(maplist)

  if(verbose){
    mapped_markers <- unlist(lapply(maplist, function(x) as.character(x$marker)))
    phased_markers <- unlist(lapply(maplist.out, function(x) as.character(x$marker)))

    if(original_coding){
      mapped.dosages <- dosage_matrix.orig[mapped_markers,]
    } else{
      mapped.dosages <- dosage_matrix.conv[mapped_markers,]
    }
    mds.b4 <- marker_data_summary(dosage_matrix = mapped.dosages,
                                  ploidy = (ploidy+ploidy2)/2,
                                  pairing = "random", verbose = FALSE)
    mds.aft <- marker_data_summary(dosage_matrix = dosage_matrix.conv[phased_markers,],
                                   ploidy = (ploidy+ploidy2)/2,
                                   pairing = "random", verbose = FALSE)

    write(paste("\nMapped marker breakdown before phasing:\n_______________________________________\n"),log.conn)
    write(knitr::kable(mds.b4$parental_info,format="markdown"), log.conn)
    write("\n", log.conn)

    write(paste("\nPhased marker breakdown:\n_______________________________________\n"),log.conn)
    write(knitr::kable(mds.aft$parental_info,format="markdown"), log.conn)
    write("\n", log.conn)

  }

  if(!is.null(log)) close(log.conn)

  return(maplist.out)
}


