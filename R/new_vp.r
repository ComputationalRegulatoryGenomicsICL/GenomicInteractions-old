viewPoint = function(x, bait, region=NULL, ...) {
    hits = findOverlaps(x, bait, ...)
    vp = GenomicInteractions(anchor_one=bait[c(subjectHits(hits$one), 
                                          subjectHits(hits$two))],
                             anchor_two=c(x@anchor_one[queryHits(hits$one)], 
                                          x@anchor_two[queryHits(hits$two)]),
                             counts=c(x@counts[queryHits(hits$one)], 
                                      x@counts[queryHits(hits$two)]))
    if (!is.null(region)) { vp = x[overlapsAny(x@anchor_two, region, ...] }
    return(sort(vp))
}

plotViewpoint = function(x, region, ...) {
    if (length(region) > 1) stop("region must be a single range")
    x = x[overlapsAny(x@anchor_two, region, type="within")]
    cov = as(coverage(x@anchor_two)[region], "GRanges")
    points_x = c(start(cov), end(cov)) + start(region)
    points_y = rep.int(x$score, 2)
    ord = order(points_x)
    p = plot(points_x[ord], points_y[ord], t="l", ...)
    return(p)
}

