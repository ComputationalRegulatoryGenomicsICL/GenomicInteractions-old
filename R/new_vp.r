library(GenomicInteractions)
library(GenomicRanges)
library(rtracklayer)
library(BSgenome.Mmusculus.UCSC.mm9)
library(magrittr)

filter_rfrags_of_interest = function(x, rfrags, probeset, ...) {
    ol = findOverlaps(x, rfrags, ...)
    if (any(countQueryHits(ol) > 1)) warning("Some ranges overlap multiple restriction fragments")
    bait_fragments = rfrags[subjectHits(ol)]
    is_probe = overlapsAny(bait_fragments, probeset)
    if (!all(is_probe)) warning("some fragments have not been enriched. Fragment ids: \n", str(which(!is_probe)))
    bait_fragments
}

multiple_viewpoint_cov = function(x, ranges) {
    overlaps_any = overlapsAny(x, ranges, type="within")
    x = x[overlaps_any$one|overlaps_any$two]
    hits = findOverlaps(x, ranges, type="within")
    hits_one = split(queryHits(hits$one), subjectHits(hits$one))
    hits_two = split(queryHits(hits$two), subjectHits(hits$two))
    cov_one = Map(function(i) coverage(x@anchor_one[i], weight=x@counts[i]), hits_two)
    cov_two = Map(function(i) coverage(x@anchor_two[i], weight=x@counts[i]), hits_one)
    Map(function(x,y) x+y, cov_one, cov_two)
}

multiple_viewpoint = function(x, ranges) {
    overlaps_any = overlapsAny(x, ranges, type="within")
    x = x[overlaps_any$one|overlaps_any$two]
    hits = findOverlaps(x, ranges, type="within")
    hits_one = split(queryHits(hits$one), subjectHits(hits$one))
    hits_two = split(queryHits(hits$two), subjectHits(hits$two))
    gr_one = Map(function(i) {gr=x@anchor_one[i]; gr$count=x@counts[i];gr}, hits_two)
    gr_two = Map(function(i) {gr=x@anchor_two[i]; gr$count=x@counts[i];gr}, hits_one)
    Map(function(x,y) sort(c(x,y)), cov_one, cov_two)
}

