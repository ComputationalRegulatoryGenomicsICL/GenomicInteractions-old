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
#' @export
setMethod("subsetByFeatures", c("GenomicInteractions", "GRanges", "missing"), function(GIObject, features, feature.class=NULL){
    i <- overlapsAny(GIObject, features)
    GIObject[i]
})

#' @rdname GenomicInteractions-subsetByFeatures-methods
#' @export
setMethod("subsetByFeatures", c("GenomicInteractions", "GRangesList", "missing"), function(GIObject, features, feature.class=NULL){
  i <- overlapsAny(GIObject, features)
  GIObject[i]
})

#' @rdname GenomicInteractions-subsetByFeatures-methods
#' @export
setMethod("subsetByFeatures", c("GenomicInteractions", "character", "character"), 
          function(GIObject, features, feature.class){
    if(!"node.class" %in% names(elementMetadata(regions(GIObject))) & 
       feature.class %in% unique(mcols(regions(GIObject))$node.class)){
      stop(paste(feature.class," has not been annotated on this GenomicInteractions object"))
    }
    #get regions which are annotated with given feature IDs
    region_idx <- which(sapply(mcols(regions(GIObject))[[paste(feature.class, "id", sep=".")]], 
                         function(x){any(features %in% x)}))
    #get object index for region idx
    gi_idx <- (anchors(GIObject, type = "first", id = TRUE) %in% region_idx) | 
              (anchors(GIObject, type = "second", id = TRUE) %in% region_idx)
    
    GIObject[gi_idx]
})
