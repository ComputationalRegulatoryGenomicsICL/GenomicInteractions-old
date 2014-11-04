
#' Summarise Interactions between defined anchors
#'
#' Calculate the number of of paired-end reads mapping between a defined set of anchors.
#' This function will ignore counts present in the input data.
#'
#' @return A GenomicInteractions object with annotated counts between anchors
#' @docType methods
#' @rdname countsBetweenAnchors-methods
#' @export
setGeneric("countsBetweenAnchors",function(x, y){standardGeneric ("countsBetweenAnchors")})

#' @param x A GenomicInteractions object
#' @param y A GenomicRanges object
#' @import GenomicRanges
#' @rdname countsBetweenAnchors-methods
#' @docType methods
#' @export
setMethod("countsBetweenAnchors", list("GenomicInteractions", "GRanges"), function(x, y) {
    #check anchors are unique
    if (any(countOverlaps(y, y) > 1)) stop("anchors are not unique")
    one = overlapsAny(anchorOne(x), y)
    two = overlapsAny(anchorTwo(x), y)
    x.valid = x[one & two]
    overlaps = findOverlaps(sort(x.valid), y, select="first") # select produces matrix not Hits
    interactions = paste(overlaps[[1]], overlaps[[2]], sep=":")
    tabulated = table(interactions)
    
    pairs_list = strsplit(names(tabulated), ":")
    pairs_one = as.integer(sapply(pairs_list, function(x) x[1]))
    pairs_two = as.integer(sapply(pairs_list, function(x) x[2]))
    
    anchor_one = y[pairs_one]
    anchor_two = y[pairs_two]
    counts = as.integer(tabulated)
    
    final_counts = new("GenomicInteractions",
                       experiment_name = name(x),
                       description = description(x),
                       genome_name = genomeName(x),
                       anchor_one=anchor_one,
                       anchor_two=anchor_two,
                       counts=counts)
    
    return(sort(final_counts))
})


