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
                                  parent1 = P1,
                                  parent2 = P2,
                                  original_coding = TRUE,
                                  log = NULL,
                                  verbose = FALSE) {


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
test_dosage_matrix <- function(dosage_matrix){
  if(class(dosage_matrix) == "data.frame"){
    warning("dosage_matrix should be a matrix, now it's a data.frame.")
    message("Trying to convert it to matrix, assuming markernames are in the first column..")
    rownames(dosage_matrix) <- dosage_matrix[,1]
    dosage_matrix <- as.matrix(dosage_matrix[,-1])
    class(dosage_matrix) <- "integer"
  } else if(class(dosage_matrix) == "matrix"){
    rn <- rownames(dosage_matrix)
    cn <- colnames(dosage_matrix)
    if(is.null(rn)) stop("The rownames of dosage_matrix should contain markernames. Now NULL")
    if(is.null(cn)) stop("The columnnames of dosage_matrix should contain genotype names. Now NULL")
    if(!(typeof(dosage_matrix)=="integer" | typeof(dosage_matrix)=="double")){
      warning("dosage_matrix should be integer or numeric. Trying to convert it.. ")
      class(dosage_matrix) <- "integer"
    }
  } else {
    stop("dosage_matrix should be a matrix of integers.
         See the manual of this function for more information.")
  }
  return(dosage_matrix)
}
marker_data_summary <- function(dosage_matrix,
                                ploidy = 4,
                                pairing = c("random", "preferential"),
                                parent1 = P1,
                                parent2 = P2,
                                progeny_incompat_cutoff = 0.1,
                                verbose = TRUE,
                                log = NULL) {

  dosage_matrix <- test_dosage_matrix(dosage_matrix)

  if (is.null(log)) {
    log.conn <- stdout()
  } else {
    matc <- match.call()
    write.logheader(matc, log)
    log.conn <- file(log, "a")
  }

  pardos <- dosage_matrix[, c(parent1, parent2)]

  if(any(is.na(pardos))){
    NAmark <- rownames(pardos)[is.na(pardos[,parent1]) | is.na(pardos[,parent2])]
    warning("There are parental scores with missing values. These are not considered in the analysis.
            It is recommended to remove those before proceeding to further steps.")
    dosage_matrix <- dosage_matrix[!rownames(dosage_matrix) %in% NAmark, ]
    if(verbose) write(paste(c("\nThe following marker have missing values in their parental scores:",
                  NAmark, "\n"), collapse = "\n\n"), file = log.conn)
  }


  pairing <- match.arg(pairing)

  add_P_to_table <- function(table) {
    #add parent info to table
    colnames(table) <- paste0("P2_", colnames(table))
    rownames(table) <- paste0("P1_", rownames(table))
    return(table)
  }

  test_row <- function(x, lu, parpos = c(1, 2)) {
    #analyse offspring incompatibility for a marker
    #with lu as lookup table for maximum and minimum offspring dosages
    progeny <- x[-parpos]
    partype <- lu$pmin == min(x[parpos]) & lu$pmax == max(x[parpos])
    min <- lu[partype, "min"]
    max <- lu[partype, "max"]
    return(!is.na(progeny) & progeny >= min & progeny <= max)
  }
  #######################################
  nm <- nrow(dosage_matrix)
  end_col <- ncol(dosage_matrix)

  if(verbose) write("Calculating parental info...", stdout())

  # contingency table number of markers

  parental_info <-
    table(as.factor(dosage_matrix[,	parent1]), as.factor(dosage_matrix[,	parent2]))

  parental_info <- add_P_to_table(parental_info)

  #Checking offspring compatability

  if(verbose) write("Checking compatability between parental and offspring scores...",
        stdout())

  parpos <- which(colnames(dosage_matrix) %in% c(parent1, parent2))

  progeny <- dosage_matrix[,-parpos]

  nr_offspring <- ncol(progeny)

  seg.fname <- paste0("seg_p", ploidy, "_", pairing)
  seg <- get(seg.fname)#,envir=getNamespace("polymapR"))
  segpar <- seg[, c("dosage1", "dosage2")]
  colnames(segpar) <- c("pmax", "pmin")
  segoff <- seg[, 3:ncol(seg)]
  segoff <- segoff > 0
  segpos <- c(0:ploidy)

  lu_min_max <- apply(segoff, 1, function(x) {
    a <- segpos[x]
    min <- min(a)
    max <- max(a)
    return(c(min, max))
  })

  rownames(lu_min_max) <- c("min", "max")
  lu <- cbind(segpar, t(lu_min_max))

  expected_dosage <-
    apply(dosage_matrix, 1, test_row, lu = lu, parpos = parpos)

  #NA should be "TRUE", now "FALSE"
  expected_dosage <- t(expected_dosage)
  if(length(which(is.na(progeny))) > 0) expected_dosage[is.na(progeny)] <- TRUE

  #two factorial table of parental dosages with percentage of "FALSE" per factor combination
  progeny_incompat <- colSums(!expected_dosage)
  na_progeny <- colSums(is.na(progeny))
  perc_incompat <-
    progeny_incompat / (nrow(expected_dosage) - na_progeny)
  progeny_incompat <-
    colnames(progeny)[perc_incompat > progeny_incompat_cutoff]

  nr_incompat <- rowSums(!expected_dosage)
  offspring_incompat <- tapply(
    nr_incompat,
    list(dosage_matrix[, parent1], dosage_matrix[, parent2]),
    FUN = function(x)
      sum(x) / (length(x) * nr_offspring) * 100
  )
  offspring_incompat <- round(offspring_incompat, 2)
  offspring_incompat <- add_P_to_table(offspring_incompat)

  summary <-
    list(parental_info, offspring_incompat, progeny_incompat)
  names(summary) <-
    c("parental_info",
      "offspring_incompatible",
      "progeny_incompatible")

  for (i in c(1, 2)) {
    if(verbose) {
      write(paste0("\n####", names(summary)[i], "\n"),
          file = log.conn)
    #sink(log.conn)
    write(knitr::kable(summary[[i]]),
          log.conn)
    }
    #suppressWarnings(sink())
  }

  if(verbose) write("\n####Incompatible individuals:\n", log.conn)
  if (length(progeny_incompat) == 0 & verbose)
    write("None\n", log.conn)

  if(verbose) write(summary$progeny_incompatible, log.conn)

  if (!is.null(log))
    close(log.conn)

  return(summary)
} #marker_data_summary()



  if(original_coding & is.null(dosage_matrix.orig)) stop("Uncoverted dosage matrix should also be specified if original_coding = TRUE")

  mapped_markers <- unlist(lapply(maplist, function(x) as.character(x$marker)))
  if(!all(mapped_markers %in% rownames(dosage_matrix.conv))) stop("Not all markers on map have corresponding dosages! If duplicated markers were added back to maps, make sure to use an appropriate dosage matrix!")

  if (is.null(log)) {
    log.conn <- stdout()
  } else {
    matc <- match.call()
    write.logheader(matc, log)
    log.conn <- file(log, "a")
  }

  if(is.null(ploidy2)) ploidy2 <- ploidy

  if(ploidy == ploidy2){
    palindromes <- rownames(dosage_matrix.conv)[which(dosage_matrix.conv[,parent1] != dosage_matrix.conv[,parent2] &
                                                        abs(dosage_matrix.conv[,parent1] - (0.5*ploidy)) == abs(dosage_matrix.conv[,parent2]-(0.5*ploidy2)))]

    ## If there are any unconverted palindromes, convert them:
    if(any(dosage_matrix.conv[palindromes,parent1] > dosage_matrix.conv[palindromes,parent2]))
      dosage_matrix.conv[palindromes[dosage_matrix.conv[palindromes,parent1] > dosage_matrix.conv[palindromes,parent2]],] <-
        ploidy - dosage_matrix.conv[palindromes[dosage_matrix.conv[palindromes,parent1] > dosage_matrix.conv[palindromes,parent2]],]
    }

  # Begin by separating the SxN and NxS linkages:
  SxN_assigned <- marker_assignment.1[marker_assignment.1[,parent1]==1 &
                                        marker_assignment.1[,parent2]==0,]
  p1_assigned <- marker_assignment.1[-match(rownames(SxN_assigned),rownames(marker_assignment.1)),]

  NxS_assigned <- marker_assignment.2[marker_assignment.2[,parent1]==0 &
                                        marker_assignment.2[,parent2]==1,]
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
  p1cols <- 3+grep("Hom",colnames(p1_assigned)[4:ncol(p1_assigned)])
  p2cols <- 3+grep("Hom",colnames(p2_assigned)[4:ncol(p2_assigned)])

  P1rates <- p1_assigned[,p1cols]/rowSums(p1_assigned[,p1cols], na.rm = TRUE)
  P2rates <- p2_assigned[,p2cols]/rowSums(p2_assigned[,p2cols], na.rm = TRUE)

  P1rates[P1rates < lower_bound] <- 0
  P2rates[P2rates < lower_bound] <- 0

  P1linked <- apply(P1rates,1,function(x) length(which(x!=0)))
  P2linked <- apply(P2rates,1,function(x) length(which(x!=0)))

  p1.markers <- rownames(p1_assigned[p1_assigned[,parent1]!=0,])
  p2.markers <- rownames(p2_assigned[p2_assigned[,parent2]!=0,])

  ## Assuming markers are converted here; have to treat palindrome markers in P2 carefully:
  P1different <- rownames(p1_assigned[rownames(p1_assigned) %in% p1.markers & p1_assigned[,parent1] != P1linked,])
  P2different <- rownames(p2_assigned[setdiff(which(rownames(p2_assigned) %in% p2.markers & p2_assigned[,parent2] != P2linked),
                                        which(rownames(p2_assigned) %in% palindromes & ploidy2 - p2_assigned[,parent2] == P2linked)),])


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

  P1rates <- P1rates[!rownames(p1_assigned) %in% P1different,]
  P2rates <- P2rates[!rownames(p2_assigned) %in% P2different,]

  #Update p1_assigned and p2_assigned
  p1_assigned <- p1_assigned[!rownames(p1_assigned) %in% P1different,]
  p2_assigned <- p2_assigned[!rownames(p2_assigned) %in% P2different,]

  rownames(P1rates) <- rownames(p1_assigned)
  rownames(P2rates) <- rownames(p2_assigned)

  # return simplex x nulliplex markers
  p1_assigned <- rbind(SxN_assigned,p1_assigned)
  p2_assigned <- rbind(NxS_assigned,p2_assigned)

  P1rates <- rbind(SxN_assigned[,p1cols],P1rates)
  P2rates <- rbind(NxS_assigned[,p2cols],P2rates)

  # Remove the bi-parental markers that are not assigned in both parents (what about unconverted markers here? Logical test is only looks for a nulliplex parent.)
  bip1 <- rownames(p1_assigned[rowSums(p1_assigned[,c(parent1,parent2)]!=0)==2,])
  bip2 <- rownames(p2_assigned[rowSums(p2_assigned[,c(parent1,parent2)]!=0)==2,])

  BiP_different <- c(setdiff(bip1,intersect(bip1,bip2)),setdiff(bip2,intersect(bip1,bip2)))

  if (verbose & !is.null(BiP_different)) {
    removed.m <- vector.to.matrix(BiP_different, n.columns = 4)

    if(nrow(removed.m) > 0){
      write(paste("\nThe following markers did not have the expected assignment across both parents:\n_______________________________________\n"),log.conn)
      write(knitr::kable(removed.m,format="markdown"), log.conn)
      write("\n", log.conn)
    }
  }

  P1rates <- P1rates[!rownames(p1_assigned) %in% setdiff(bip1,intersect(bip1,bip2)),]
  P2rates <- P2rates[!rownames(p2_assigned) %in% setdiff(bip2,intersect(bip1,bip2)),]

  #Update p1_assigned and p2_assigned
  p1_assigned <- p1_assigned[!rownames(p1_assigned) %in% setdiff(bip1,intersect(bip1,bip2)),]
  p2_assigned <- p2_assigned[!rownames(p2_assigned) %in% setdiff(bip2,intersect(bip1,bip2)),]

  ALL_assigned <- unique(c(rownames(p1_assigned),rownames(p2_assigned)))

  # Make up the output
  maplist.out <- lapply(seq(length(maplist)),function(mapn) {

    map <- maplist[[mapn]]
    map <- map[map$marker%in%ALL_assigned,]

    outmap <- map[,c("marker","position")]

    hom_mat <- sapply(1:nrow(outmap), function(r){
      a <- rep(0, ploidy+ploidy2)

      temp <- P1rates[match(as.character(outmap$marker[r]),rownames(P1rates)),]
      if(length(which(temp!=0)) > 0) a[(1:ploidy)[which(temp!=0)]] <- 1

      temp <- P2rates[match(outmap$marker[r],rownames(P2rates)),]
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

      orig_parents <- dosage_matrix.orig[match(outmap$marker,rownames(dosage_matrix.orig)),c(parent1,parent2)]
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

  phased_markers <- unlist(lapply(maplist.out, function(x) as.character(x$marker)))

  if(original_coding){
    mapped.dosages <- dosage_matrix.orig[mapped_markers,]
  } else{
    mapped.dosages <- dosage_matrix.conv[mapped_markers,]
  }

  if(verbose){
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

  ## Run a final check to make sure that the phased marker dosages equal the original marker dosages:
  phased.dose <- do.call(rbind,lapply(maplist.out, function(x) {
    temp <- cbind(rowSums(x[,paste0("h",1:ploidy)]),
                  rowSums(x[,paste0("h",(ploidy + 1):(ploidy + ploidy2))]))
    rownames(temp) <- x[,"marker"]
    return(temp)
    }))

  orig.dose <- mapped.dosages[rownames(phased.dose),c(parent1,parent2)]

  conflicting <- which(rowSums(phased.dose == orig.dose) != 2)

  if(length(conflicting) > 0){
    warning("Not all phased markers matched original parental dosage. \nPerhaps unconverted marker dosages were supplied as converted dosages by mistake? \nThe following conflicts were detected and removed:")
    warn.df <- cbind(orig.dose[conflicting,],phased.dose[conflicting,])
    colnames(warn.df) <- c("P1_original","P2_original","P1_phased","P2_phased")
    write(knitr::kable(warn.df,format="markdown"), log.conn)

    ## Simply remove these markers from the output:
    rem.markers <- rownames(phased.dose)[conflicting]

    maplist.out <- lapply(maplist.out, function(x) x[!x$marker %in% rem.markers,])
    }

  if(!is.null(log)) close(log.conn)

  return(maplist.out)
}