
#' Summarise Interactions between defined anchors
#'
#' Calculate the number of of paired-end reads mapping between a defined set of anchors.
#' This function will ignore counts present in the input data.
#'
#' @return A GenomicInteractions object with annotated counts between anchors
#' @docType methods
#' @rdname countsBetweenAnchors-methods
#' @export
setGeneric("countsBetweenAnchors",function(x, y, ...){standardGeneric ("countsBetweenAnchors")})

#' @param x A GenomicInteractions object
#' @param y A GenomicRanges object
#' @param ignore_overlaps Allow overlapping anchors. Use this when you have overlapping anchors
#'                        but be careful with multi-mapping. The "within" option can help with this.
#' @param ... Extra parameters to pass to findOverlaps
#' @import GenomicRanges
#' @rdname countsBetweenAnchors-methods
#' @docType methods
#' @export
setMethod("countsBetweenAnchors", list("GenomicInteractions", "GRanges"), function(x, y, ignore_overlaps=FALSE, ...) {
    #check anchors are unique
    if (ignore_overlaps == FALSE && any(countOverlaps(y, y) > 1)) stop("anchors are not unique")
    one = overlapsAny(anchorOne(x), y, ...)
    two = overlapsAny(anchorTwo(x), y, ...)
    x.valid = x[one & two]
    overlaps = findOverlaps(sort(x.valid), y, select="first", ...) # select produces matrix not Hits
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

#' Get self ligation threshold with SD method from Heidari et al
#'
#' @export

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
                  coord_cartesian(xlim=c(0, 20000)) + geom_vline(xintercept=bp_cutoff, linetype="dashed") +
                  xlab("Distance (bp)") + ylab("log2 ratio opposite strand pairs / same strand pairs")
        )
    }
    return(bp_cutoff)
}

#' get self ligation threshold with binomial test
#'
#' @export
get_binom_ligation_threshold = function(GIObject, max.distance=20000, bin.size=500, p.cutoff=0.05, adjust="fdr", plot=TRUE){

    #make data frame
    stranded_df <- data.frame(Distance=calculateDistances(GIObject), SameStrand=sameStrand(GIObject))
    stranded_cis_df <- stranded_df[complete.cases(stranded_df),]
    stranded_cis_df <- stranded_cis_df[order(stranded_cis_df$Distance),]
    stranded_cis_df = stranded_cis_df[ stranded_cis_df$Distance < max.distance, ]

    #bin data
    bins = cut(stranded_cis_df$Distance, breaks=seq(0, max.distance, by=bin.size), include.lowest=TRUE)
    stranded_cis_df$Bin = bins
    byBin <- group_by(stranded_cis_df, Bin)
    sum_byBin <- summarise(byBin, Total=n(), SameStrand=sum(SameStrand))

    #get and adjust p values
    sum_byBin$p.value <- sapply(1:nrow(sum_byBin), function(x){binom.test(sum_byBin$SameStrand[x], sum_byBin$Total[x])$p.value})

    if(!is.na(adjust)){
        sum_byBin$p.value <- p.adjust(sum_byBin$p.value, method=adjust)
    }

    #get cutoff
    bp_cutoff <- seq(0, max.distance, by=bin.size)[min(which(sum_byBin$p.value > p.cutoff))]

    if (plot){
        #data for plotting
        sum_byBin <- mutate(sum_byBin, OppStrand=Total-SameStrand,
                            Bin=bin.size*as.numeric((Bin)),
                            OppPercent=100*OppStrand/Total)

        #plot % opposite strand reads and cutoff
        print(ggplot(sum_byBin, aes(x=Bin, y=OppPercent)) + geom_line() + geom_point() +
                  coord_cartesian(xlim=c(0, max.distance)) + geom_vline(xintercept=bp_cutoff, linetype="dashed") +
                  ylab("Opposite Strand Percentage")
        )
        #plot p values, p value cutoff and distance cutoff
        print(ggplot(sum_byBin, aes(x=Bin, y=p.value)) + geom_line() + geom_point() +
                  coord_cartesian(xlim=c(0, max.distance)) + geom_hline(xintercept=p.cutoff, linetype="dashed") +
                  geom_vline(xintercept=bp_cutoff, linetype="dashed") +
                  ylab("p value") + xlab("Distance (bp)")
        )
    }
    #return distance cutoff
    return(bp_cutoff)
}

#' Function to create GenomicInteraction objects from a file
#'
#' Function to create GenomicInteraction objects from a variety of files. The resulting objects contain information
#' on which genomic regions are interacting with each other, and the number of counts supporting each interaction.
#' It is also possible to store information on associated p-values and false-discovery rates (FDR).
#' It is possible to create GenomicInteractions objects for various datasets including Hi-C and ChIA-PET. It is possible
#' to read interactions from a variety of files including BAM files, bed files (BED12 and BEDPE) and from the output
#' from standard processing pipelines, such as HOMER and ChIA-PET tool. GenomicInteractions objects can also be created
#' using calls of the form \code{new("GenomicInteractions", ...)}. For hiclib, it expects the directory in which the files
#' extracted using h5dictToTxt.py from the hdf5 file are located, where as for all of the other file types it expects the full
#' filename.
#'
#' @param fn Filename or, if type="hiclib", folder
#' @param type One of "chiapet.tool", "chiapet.encode", "bed12", "bedpe", "hiclib", "homer", "bam", "two.bams".
#' @param experiment_name Experiment name.
#' @param description Description of experiment.
#' @param gname Genome name to use for constructing the GenomicInteractions object.
#' @return a GenomicInteractions object
#'
#' @importFrom Rsamtools scanBamFlag ScanBamParam scanBam bamFlagAsBitMatrix
#' @importFrom IRanges IRanges
#' @import data.table
#' @importFrom stringr str_split
#' @importFrom rtracklayer import.bed
#'
#' @examples
#' library(BSgenome.Hsapiens.UCSC.hg19)
#'
#' k562.rep1 = GenomicInteractions(file.path(system.file(package="GenomicInteractions"), "extdata", "k562.rep1.cluster.pet3+.txt"),
#'                                      type="chiapet.tool",
#'                                      experiment_name="k562",
#'                                      description="k562 pol2 8wg16",
#'                                      gname="BSgenome.Hsapiens.UCSC.hg19")
#'
#' k562.rep1
#'
#' @export

makeGenomicInteractionsFromFile = function(fn, type, experiment_name, description, gname){
    genome = .loadGenome(gname)
	if(type == "chiapet.tool"){
	    dat = read.table(fn, header=TRUE, stringsAsFactors=FALSE, sep="\t")
        anchor1 = GRanges(dat[,"chrom_left"],
                            IRanges(as.integer(dat[,"start_left"])+1, as.integer(dat[,"end_left"])),
                            seqlengths=seqlengths(genome))
        anchor2 = GRanges(dat[,"chrom_right"],
                            IRanges(as.integer(dat[,"start_right"])+1, as.integer(dat[,"end_right"])),
                            seqlengths=seqlengths(genome))
        counts = as.integer(dat[,"pet.counts.between.left.and.right.anchors"])
        p.value = as.numeric( dat[, "p.value"])
        fdr = as.numeric( dat[, "FDR"])
    }else if(type == "chiapet.encode"){
        dat = .processChiapetName(unique(import.bed(fn)$name))
        anchor1 = GRanges(dat[,"chrom.left."],
                          IRanges(as.integer(dat[,"start.left."]), as.integer(dat[,"end.left."])),
                          seqlengths=seqlengths(genome))
        anchor2 = GRanges(dat[,"chrom.right."],
                          IRanges(as.integer(dat[,"start.right."]), as.integer(dat[,"end.right."])),
                          seqlengths=seqlengths(genome))
        counts = as.integer(dat[,"counts"])
        p.value = numeric(0)
        fdr = numeric(0)
    }else if(type == "bed12"){
        bedfile = import.bed(fn)
        dat = .processChiapetName(unique(bedfile$name))
        anchor1 = GRanges(dat[,"chrom.left."],
                          IRanges(as.integer(dat[,"start.left."]), as.integer(dat[,"end.left."])),
                          seqlengths=seqlengths(genome))
        anchor2 = GRanges(dat[,"chrom.right."],
                          IRanges(as.integer(dat[,"start.right."]), as.integer(dat[,"end.right."])),
                          seqlengths=seqlengths(genome))
        counts = as.integer(dat[,"counts"])
        p.value = numeric(0)
        fdr = numeric(0)
    }else if(type == "bedpe"){
        dat = read.table(fn, stringsAsFactors=FALSE, sep="\t")
        anchor1 = GRanges(dat[,1],
                          IRanges(dat[,2]+1, dat[,3]), strand=ifelse(ncol(dat) >= 10, dat[,9], "*"),
                          seqlengths=seqlengths(genome))
        anchor2 = GRanges(dat[,4],
                          IRanges(dat[,5]+1, dat[,6]), strand=ifelse(ncol(dat) >= 10, dat[,10], "*"),
                          seqlengths=seqlengths(genome))
        counts = as.integer(rep(1, length(anchor1)))
        p.value = numeric(0)
        fdr = numeric(0)
    }else if(type == "hiclib"){
	    dat = .importHicLib(fn, genome)
        anchor1 = GRanges(dat$chrm1,
                          IRanges(dat$mid1-round(dat$fraglength1/2), dat$mid1 + round(dat$fraglength1/2)),
                          seqlengths=seqlengths(genome), fragid=dat$fragid1)
		anchor2 = GRanges(dat$chrm2,
                          IRanges(dat$mid2-round(dat$fraglength2/2), dat$mid2 + round(dat$fraglength2/2)),
                          seqlengths=seqlengths(genome), fragid=dat$fragid2)
        counts = dat$N
		p.value = numeric(0)
		fdr = numeric(0)
    }else if(type == "homer"){
        dat = .importHomer(fn)
        anchor1 = GRanges(dat$chr.1.,
                          IRanges(dat$start.1., dat$end.1.),
                          seqlengths=seqlengths(genome))
        anchor2 = GRanges(dat$chr.2.,
                          IRanges(dat$start.2., dat$end.2.),
                          seqlengths=seqlengths(genome))
        counts = as.integer(dat$Interaction.Reads)
        p.value = exp(as.numeric(dat$LogP))
        fdr = as.numeric(dat$FDR.Benjamini.)
	}else if(type == "bam"){
        dat = .readBam(fn, genome)
        anchor1 = dat[[1]]
        anchor2 = dat[[2]]
        counts = as.integer(rep(1, length(anchor1)))
        p.value = numeric(0)
        fdr = numeric(0)
	}else if(type == "two.bams"){
	    dat = .readTwoBams(fn, genome)
	    anchor1 = dat[[1]]
	    anchor2 = dat[[2]]
	    counts = as.integer(rep(1, length(anchor1)))
	    p.value = numeric(0)
	    fdr = numeric(0)
	}else{
        stop("type is not one of \"chiapet.tool\", \"chiapet.encode\", \"bed12\", \"bedpe\", \"hiclib\", \"homer\", \"bam\", \"two.bams\"")
	}
    giobject = new("GenomicInteractions",
                 experiment_name = experiment_name,
                 description = description,
                 genome_name = gname,
                 anchor_one=anchor1,
                 anchor_two=anchor2,
                 counts=counts,
                 pvalue=p.value,
                 fdr=fdr)
    return(giobject)
}
