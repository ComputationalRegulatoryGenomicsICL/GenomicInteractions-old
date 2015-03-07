## Definition of GenomicInteractions

#' A S4 class to represent interactions between genomic regions.
#'
#'  @slot metadata List, defaults to "experiment_name" and "description", inherited from S4Vectors::Vector
#'  @slot anchor_one,anchor_two GRanges. Set of anchors of interactions.
#'  @slot elementMetadata DataFrame
#'
#' This class is used to store information on which genomic regions are
#' interacting with each other. Objects of this class contain information of
#' the genomic coordinates of the interacting regions and the strength of these
#' interactions, and associated metadata such as the name of the dataset and a
#' brief description of the dataset.  Interacting regions are stored as a pair
#' of GenomicRanges: each set of anchor regions is stored as a separate
#' GenomicRanges object, accessed by \code{getAnchorOne} and
#' \code{getAnchorTwo}.
#'
#' @examples
#'
#' showClass("GenomicInteractions")
#'
#' library(BSgenome.Mmusculus.UCSC.mm9)
#' anchor.one = GRanges(c("chr1", "chr1", "chr1", "chr1"), IRanges(c(10, 20, 30, 20), width=5), seqlengths=seqlengths(Mmusculus))
#' anchor.two = GRanges(c("chr1", "chr1", "chr1", "chr2"), IRanges(c(100, 200, 300, 50), width=5), seqlengths=seqlengths(Mmusculus))
#' test <- new("GenomicInteractions", experiment_name="test", description="this is a test",
#'             anchor_one = anchor.one, anchor_two = anchor.two, counts=as.integer(c(2,1,2,3)),
#'             pvalue=c(0.1, 0.3, 0.1, 0.08))
#'
#' @import GenomicRanges
#' @export GenomicInteractions
setClass("GenomicInteractions",
    representation(anchor_one = "GRanges",
                   anchor_two = "GRanges",
                   elementMetadata="DataFrame"),
    prototype(anchor_one = GRanges(),
              anchor_two = GRanges(),
              elementMetadata = DataFrame() ),
    contains="Vector",
    validity = function(object){
        if (length(object@anchor_one) == 0 ) {
            return("anchor one cannot be of length 0")
        } else if(length(object@anchor_one) != length(object@anchor_two)) {
            return("length of anchor one and anchor two do not match")
        } else if(!.isEqualSeqInfo(object@anchor_one, object@anchor_two)) {
            return("seqinfo must be indentical for both GRanges") # this is order-dependent which is not desireable
        } else{ return(TRUE)}}
)

#' Function to create a GenomicInteraction object
#'
#' Create GenomicInteraction objects from two GRanges ojects.
#'
#' @param anchor_one, anchor_two GRanges objects.
#' @param experiment_name Experiment name.
#' @param description Description of experiment.
#' @param ... Additional data to be added to mcols
#' @return a GenomicInteractions object
#'
#' @examples
#' anchor.one = GRanges(c("chr1", "chr1", "chr1", "chr1"), IRanges(c(10, 20, 30, 20), width=5), seqlengths=seqlengths(Mmusculus))
#' anchor.two = GRanges(c("chr1", "chr1", "chr1", "chr2"), IRanges(c(100, 200, 300, 50), width=5), seqlengths=seqlengths(Mmusculus))
#' interaction_counts = sample(1:10, 4)
#' test <- GenomicInteractions(anchor.one, anchor.two, experiment_name="test", description="this is a test", counts=interaction_counts)
#'
#' @export
GenomicInteractions = function(anchor_one, anchor_two, experiment_name="", description="", ...) {
    mcols = DataFrame(...)
    if (ncol(mcols) == 0L)
        mcols <- new("DataFrame", nrows = length(seqnames))
    if (mrow(mcols) == 1)
        do.call(rbind, replicate(mcols, length(anchor_one)) # inefficient
    new("GenomicInteractions",
        metadata=list(experiment_name=experiment_name, description=description),
        anchor_one=anchor_one,
        anchor_two=anchor_two,
        elementMetadata=mcols)
}

#' Get the length of a GenomicInteractions GIObject
#'
#' @param x GenomicInteractions GIObject
#' @return A numeric vector containing the length of the GIObject
#' @docType methods
#' @export
setMethod(length, "GenomicInteractions", function(x) length(x@anchor_one))

# Quick access to mcols()

setMethod("$", "GenomicInteractions",
    function(x, name) mcols(x)[[name]]
)

setReplaceMethod("$", "GenomicInteractions",
    function(x, name, value) {mcols(x)[[name]] <- value; x}
)

