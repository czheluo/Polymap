
checkmap<-function (linkage_list, maplist, mapfn = "haldane", lod.thresh = 5,LG)
{

colour.bar <- function(col.data, min, max=-min, cex.ticks = 1.2, nticks=11,
                      ticks=seq(min, max, len=nticks), title='', ylab = '',
                      cex.lab = 1) {
  scale <- length(col.data)/(max-min)

  plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab=ylab, main=title, cex.lab = cex.lab)
  axis(2, ticks, las=1, cex.axis = cex.ticks)
  for (i in 1:length(col.data)) {
    y = (i-1)/scale + min
    rect(0,y,10,y+1/scale, col=col.data[i], border=NA)
  }
}

#maplist<-list(integrated.maplist$LG1)
#linkage_list<-list(linkages$LG1)

#mapfn<-"haldane"
#lod.thresh<-3

if (length(linkage_list) != length(maplist))
  stop("linkage_list and maplist do not correspond.")
mapfn <- match.arg(mapfn, choices = c("haldane", "kosambi"))
rev.haldane <- function(d) (1 - exp(-d/50))/2
rev.kosambi <- function(d) ((exp(d/25) - 1)/(exp(d/25) +
                                               1))/2
orig.mar <- c(5.1, 4.1, 4.1, 2.1)
colbar.mar <- c(5.1, 2, 4.1, 0.5)

lmposmap<-sapply(l<-seq(length(linkage_list)), function(l) {
	posmat <- matrix(c(maplist[[l]][match(linkage_list[[l]]$marker_a,
                                        maplist[[l]]$marker), ]$position,
                     maplist[[l]][match(linkage_list[[l]]$marker_b,
                                        maplist[[l]]$marker), ]$position, linkage_list[[l]]$r,
                     linkage_list[[l]]$LOD), ncol = 4)

	if (mapfn == "haldane") {
		expected.recom <- rev.haldane(abs(posmat[, 1] - posmat[, 2]))
	} else if (mapfn == "kosambi") {
		expected.recom <- rev.kosambi(abs(posmat[, 1] -posmat[, 2]))
	}
	result<-data.frame(posmat,expected.recom)
	return(result)
}
)
ll<-seq(length(linkage_list));posmat<-do.call(cbind,unname(lmposmap))
dev <- abs(linkage_list[[ll]]$r - posmat[,5])

wRMSE <- sqrt(mean(is.na((dev * linkage_list[[ll]]$LOD)^2)))

layout(matrix(c(1, 1, 1, 1, 2, 3, 4, 5), ncol = 4, byrow = TRUE),
           widths = c(1, 0.2, 1, 0.2))
    par(oma = c(0, 0, 3, 0))

    plot(dev, linkage_list[[ll]]$LOD, ylab = "LOD", xlab = expression(delta(r)),
         cex.axis=1.2,cex.lab = 1.5, main = expression("|r"["pairwise"] ~
                                             "- r"["map"] ~ "|"))
    legend("topright", legend = paste("Weighted RMSE =",
                                      round(wRMSE, 3)), bty = "n")
    #colours <- colorRampPalette(c("green", "yellow", "orange","red"))(100)
    colours <- colorRampPalette(c("red","green","blue"))(100)

    rcolmin <- min(is.na(posmat[, 3]))

    rcolmax <- max(is.na(posmat[, 3]))
    LODcolmin <- lod.thresh
    LODcolmax <- max(posmat[, 4])
    rcolbreaks <- seq(rcolmin, rcolmax, (rcolmax - rcolmin)/100)

    rcols <- colours[findInterval(x = posmat[, 3], vec = rcolbreaks)]

    rcols[is.na(rcols)] <- colours[100]

    LODcolbreaks <- seq(LODcolmin, LODcolmax, (LODcolmax - LODcolmin)/100)

    LODcols <- colours[findInterval(x = posmat[, 4], vec = LODcolbreaks)]

    LODcols[is.na(LODcols)] <- colours[100]

     plot(posmat[, 1], posmat[, 2], pch = ".", col = rcols,
         main = "r", cex = 4,cex.lab = 1.5,cex.axis=1.5, xlab = "cM", ylab = "cM")

    par(mar = colbar.mar)
    colour.bar(col.data = colours, min = rcolmin, max = rcolmax,
               nticks = 8, ticks = round(seq(rcolmin, rcolmax,
                                             len = 8), 2), cex.ticks = 1.2)
    par(mar = orig.mar)
    plot(posmat[posmat[, 4] > lod.thresh, 1], posmat[posmat[,4] > lod.thresh, 2], pch = ".", col = LODcols, main = "LOD",
         cex = 4,cex.lab = 1.5,cex.axis=1.5, xlab = "cM", ylab = "cM")

    par(mar = colbar.mar)
    colour.bar(col.data = colours, min = LODcolmin, max = LODcolmax,
               nticks = 8, ticks = round(seq(LODcolmin, LODcolmax,
                                             len = 8), 2), cex.ticks = 1.2)
    mtext(text = paste("LG", LG, "map diagnostics"), side = 3,
          outer = TRUE, cex = 2)
    par(mfrow = c(1, 1), oma = c(0, 0, 0, 0), mar = orig.mar)

}
