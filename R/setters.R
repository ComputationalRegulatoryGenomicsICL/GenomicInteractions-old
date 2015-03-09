#' Functions to set data held in a GenomicInteractions object.
#'
#' Use these functions to set data stored in each of the slots of a
#' GenomicInteractions object.
#'
#' @name setters
#' @param GIObject A GenomicInteractions object
#' @param value A vector to replace a slot in the object
#' @return GenomicInteractions object
#'
#'  @examples
#' library(BSgenome.Mmusculus.UCSC.mm9)
#' anchor.one = GRanges(c("chr1", "chr1", "chr1", "chr1"), IRanges(c(10, 20, 30, 20), width=5), seqlengths=seqlengths(Mmusculus))
#' anchor.two = GRanges(c("chr1", "chr1", "chr1", "chr2"), IRanges(c(100, 200, 300, 50), width=5), seqlengths=seqlengths(Mmusculus))
#' test <- new("GenomicInteractions", experiment_name="test", description="this is a test",
#'                  genome_name="BSgenome.Mmusculus.UCSC.mm9", anchor_one = anchor.one,
#'                  anchor_two = anchor.two, counts=as.integer(c(2,1,2,3)) )
#'
#' name(test) <- "Mouse test"
#' name(test)
#'
#' description(test) <- "This is a test using the mouse genome"
#' description(test)
#'
#' interactonCounts(test) <- c(2,3,8,5)
#' interactonCounts(test)
#'
## GENERICS
#' @rdname setters
#' @aliases name<-
#' @export
setGeneric("name<-",function(GIObject, value){standardGeneric ("name<-")})

#' @rdname setters
#' @export
setGeneric("description<-",function(GIObject, value){standardGeneric ("description<-")})

#' @rdname setters
#' @export
setGeneric("interactionCounts<-",function(GIObject, value){standardGeneric ("interactionCounts<-")})

## METHODS

#' @rdname setters
#' @export
setReplaceMethod("name", "GenomicInteractions", function(GIObject, value){
    GIObject@metadata$experiment_name = value
    GIObject
    })

#' @rdname setters
#' @export
setReplaceMethod("description", "GenomicInteractions", function(GIObject, value){
    GIObject@metadata$description = value
    GIObject
})

#' @rdname setters
#' @export
setReplaceMethod("interactionCounts", "GenomicInteractions", function(GIObject, value){
    if (!all(value == floor(value)))
        stop("value must contain integer values")
    value = as.integer(value)
    if (length(value) == 1)
        value = rep(value, length(GIObject))
    BiocGenerics:::updateS4(GIObject, counts=value)
})
