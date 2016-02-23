#' Subset a GenomicInteractions object by features
#'
#' Subsets interactions for which at least one of the anchors overlaps with a given GRanges object.
#' Alternatively, subsets interactions based on annotated feature IDs for a particular feature.
#'
#' @rdname GenomicInteractions-subsetByFeatures-methods
#' @docType methods
#' @param GIObject A GenomicInteractions object
#' @param features A GRanges or GRangesList object, or a character vector containing
#'                      IDs of annotated features, e.g. promoter IDs.
#' @param feature.class If `features' is a character vector, the corresponding feature name, e.g. "promoter".
#' @return a subsetted GenomicInteractions object
#' @export
setGeneric("subsetByFeatures",function(GIObject, features, feature.class=NULL){standardGeneric ("subsetByFeatures")})

#' @rdname GenomicInteractions-subsetByFeatures-methods
#' @import GenomicRanges
#' @export
setMethod("subsetByFeatures", c("GenomicInteractions", "GRanges", "missing"), function(GIObject, features, feature.class=NULL){
    i = unique(c(subjectHits(findOverlaps(features, GIObject@anchor_one)), subjectHits(findOverlaps(features, GIObject@anchor_two))))
    GIObject[i]
})

#' @rdname GenomicInteractions-subsetByFeatures-methods
#' @import GenomicRanges
#' @export
setMethod("subsetByFeatures", c("GenomicInteractions", "GRangesList", "missing"), function(GIObject, features, feature.class=NULL){
    i = unique(c(subjectHits(findOverlaps(features, GIObject@anchor_one)), subjectHits(findOverlaps(features, GIObject@anchor_two))))
    GIObject[i]
})

#' @rdname GenomicInteractions-subsetByFeatures-methods
#' @import GenomicRanges
#' @export
setMethod("subsetByFeatures", c("GenomicInteractions", "character", "character"), function(GIObject, features, feature.class){
    if(!"node.class" %in% names(elementMetadata(GIObject@anchor_one)) & feature.class %in% unique(c(GIObject@anchor_one$node.class, GIObject@anchor_two$node.class)))
        stop(paste(feature.class," has not been annotated on this GenomicInteractions object"))
    i = sapply(elementMetadata(GIObject@anchor_one)[[paste(feature.class, "id", sep=".")]],
               function(x){ features %in% x }) | sapply(elementMetadata(GIObject@anchor_two)[[paste(feature.class, "id", sep=".")]], function(x){ features %in% x })
    GIObject[i]
})
