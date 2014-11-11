
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


get_self_ligation_threshold <- function(GIObject, bins=100, distance_th=400000, plot=TRUE){
    require(dplyr)
    require(ggplot2)
    #get df
    stranded_df <- data.frame(Distance=calculateDistances(GIObject), SameStrand=sameStrand(GIObject))
    stranded_cis_df <- stranded_df[complete.cases(stranded_df),]
    stranded_cis_df <- stranded_cis_df[order(stranded_cis_df$Distance),]
    
    #bin data
    bin_n <- nrow(stranded_cis_df )/bins
    
    cuts <- 1:bins * bin_n
    breaks <- stranded_cis_df$Distance[cuts]
    
    stranded_cis_df$Bin <- as.numeric(as.character(cut(stranded_cis_df$Distance, breaks=breaks, labels=breaks[1:length(breaks)-1], include.lowest = TRUE)))
    byBin <- group_by(stranded_cis_df, Bin)
    
    #summarise by bin
    sum_byBin <- summarise(byBin, Total=n(), SameStrand=sum(SameStrand))
    sum_byBin <- mutate(sum_byBin, OppStrand=Total-SameStrand, 
                        Ratio= SameStrand/Total,
                        log2Ratio=log2((OppStrand+1)/(SameStrand+1))) #pseudocount to avoid NaN errors
    sum_byBin <- mutate(sum_byBin, OppPercent=100*OppStrand/Total, SamePercent=100*SameStrand/Total)
    
    #get cutoff of log2ratio
    sum_byBin %>%
        filter(Bin > distance_th) %>%  
        select(log2Ratio) %>%
        unlist() %>%
        mean() -> longrange_mean_log2
    
    sum_byBin %>%
        filter(Bin > distance_th) %>%  
        select(log2Ratio) %>%
        unlist() %>%
        sd() -> longrange_sd_log2
    
    lower <- longrange_mean_log2 - 2*longrange_sd_log2
    upper <- longrange_mean_log2 + 2*longrange_sd_log2
    bp_cutoff <- min(sum_byBin[sum_byBin$log2Ratio > lower & sum_byBin$log2Ratio < upper,"Bin"])
    
    if (plot){
        print(ggplot(sum_byBin, aes(x=Bin, y=log2Ratio)) + geom_line() + geom_point() + 
                  geom_hline(aes_string(yintercept=lower)) + 
                  geom_hline(aes_string(yintercept=upper)) + 
                  coord_cartesian(xlim=c(0, 100000)) + geom_vline(xintercept=bp_cutoff)
        )
    }
    return(bp_cutoff)
}

get_binom_ligation_threshold = function(GIObject, max.distance=20000, bin.size=500, p.cutoff=0.05, adjust="fdr"){
    
    stranded_df <- data.frame(Distance=calculateDistances(GIObject), SameStrand=sameStrand(GIObject))
    stranded_cis_df <- stranded_df[complete.cases(stranded_df),]
    stranded_cis_df <- stranded_cis_df[order(stranded_cis_df$Distance),]
    
    stranded_cis_df = stranded_cis_df[ stranded_cis_df$Distance < max.distance, ]
    
    bins = cut(stranded_cis_df$Distance, breaks=seq(0, max.distance, by=bin.size), include.lowest=TRUE)
    stranded_cis_df$Bin = bins 
    byBin <- group_by(stranded_cis_df, Bin)
    sum_byBin <- summarise(byBin, Total=n(), SameStrand=sum(SameStrand))
    
    if(!is.na(adjust)){
        return(seq(0, max.distance, by=bin.size)[  min(which(sapply(1:nrow(sum_byBin), function(x){ p.adjust(binom.test(sum_byBin$SameStrand[x], sum_byBin$Total[x])$p.value, method=adjust)}) > p.cutoff))])
    }else{
        return(seq(0, max.distance, by=bin.size)[  min(which(sapply(1:nrow(sum_byBin), function(x){ binom.test(sum_byBin$SameStrand[x], sum_byBin$Total[x])$p.value}) > p.cutoff))])
    }
}
