#' Functions to access data held in a GenomicInteractions object.
#'
#' Use these functions to access data stored in each of the slots of a
#' GenomicInteractions object.
#'
#' @name getters
#' @param GIObject A GenomicInteractions object
#'
#' @return For 'anchorOne' and 'anchorTwo', a GRanges. For 'counts',
#'  'normalisedCount', pValue', 'FDR', a numeric vector with counts,
#'  normalised counts, p-values or FDRs for each interaction in the object. For
#'   'description','name', and 'genomeName', a character vector with
#'   length 1. For 'annotationFeatures', a character vector of features with
#'   which the object was previously annotated, or 'NA' if the object is unannotated.
#'
#'  @examples
#' library(BSgenome.Mmusculus.UCSC.mm9)
#' anchor.one = GRanges(c("chr1", "chr1", "chr1", "chr1"), IRanges(c(10, 20, 30, 20), width=5), seqlengths=seqlengths(Mmusculus))
#' anchor.two = GRanges(c("chr1", "chr1", "chr1", "chr2"), IRanges(c(100, 200, 300, 50), width=5), seqlengths=seqlengths(Mmusculus))
#' test <- new("GenomicInteractions", experiment_name="test", description="this is a test",
#'                  genome_name="BSgenome.Mmusculus.UCSC.mm9", anchor_one = anchor.one,
#'                  anchor_two = anchor.two, counts=as.integer(c(2,1,2,3)), pvalue=c(0.1, 0.3, 0.1, 0.08))
#'
#' name(test)
#' description(test)
#' anchorOne(test)
#' anchorTwo(test)
#' count(test)
#'
## GENERICS

#' @rdname getters
#' @export
setGeneric("name",function(GIObject){standardGeneric ("name")})

#' @rdname getters
#' @export
setGeneric("description",function(GIObject){standardGeneric ("description")})

#' @rdname getters
#' @export
setGeneric("anchorOne",function(GIObject){standardGeneric ("anchorOne")})

#' @rdname getters
#' @export
setGeneric("anchorTwo",function(GIObject){standardGeneric ("anchorTwo")})

#' @rdname getters
#' @export
setGeneric("interactionCounts",function(GIObject){standardGeneric ("interactionCounts")})

#' @rdname getters
#' @export
setGeneric("annotationFeatures",function(GIObject){standardGeneric ("annotationFeatures")})

## METHODS

#' @rdname getters
#' @export
#' @aliases name
setMethod("name", "GenomicInteractions", function(GIObject){ return(GIObject@metadata$experiment_name) } )

#' @rdname getters
#' @export
setMethod("description", "GenomicInteractions", function(GIObject){ return(GIObject@metadata$description) } )

#' @rdname getters
#' @export
setMethod("anchorOne", "GenomicInteractions", function(GIObject){ return(GIObject@anchor_one) } )

#' @rdname getters
#' @export
setMethod("anchorTwo", "GenomicInteractions", function(GIObject){ return(GIObject@anchor_two) } )

#' @rdname getters
#' @export
setMethod("interactionCounts", "GenomicInteractions", function(GIObject){ return(GIObject@counts) })

#' @rdname getters
#' @export
setMethod("annotationFeatures", "GenomicInteractions", function(GIObject){
  if( "node.class" %in% names(elementMetadata(GIObject@anchor_one))) {
    annotation = unique(c(GIObject@anchor_one$node.class, GIObject@anchor_two$node.class))
  } else { annotation = "NA" }
  return(annotation)
} )
