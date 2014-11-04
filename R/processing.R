
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

#' Remove all but one occurences of a duplicated interaction
#' 
#' Removes all but the first occurence of a duplicated interaction (defined as 
#' having identical coordinates for both anchors). N.B. this does not summarise 
#' the total counts of all the duplicates. It is designed for removing potential 
#' PCR duplicates after reading in .bam files.
#' 
#' @param GIObject A GenomicInteractions object.
#' @return A GenomicInteractions object that is a subset of the input object.
#' @import GenomicRanges

.removeDups <- function(GIObject){
    dat <- data.frame(Chr1 = seqnames(anchorOne(GIObject)),
                      Start1 = start(anchorOne(GIObject)),
                      Chr2 = seqnames(anchorTwo(GIObject)),
                      Start2 = start(anchorTwo(GIObject))
    )
    idx <- which(!duplicated(dat))
    reads_removed <- length(GIObject) - length(idx)
    percent_removed <- signif(100*reads_removed / length(GIObject), 3)
    message(paste0("Removing ", reads_removed, " duplicate PETs (", percent_removed, "%)"))
    return(GIObject[idx])
}

#' Tests whether anchors have the same strand.
#' 
#' This is designed for processing .bam files. 
#' 
#' @param GIObject A GenomicInteractions object
#' @return A logical vector denoting with TRUE if both anchors of an interaction
#'  are on the same strand and FALSE otherwise. 

sameStrand <- function(GIObject){
    return(strand(anchorOne(GIObject))==strand(anchorTwo(GIObject)))
}

