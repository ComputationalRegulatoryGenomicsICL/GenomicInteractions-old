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
}# import methods

#' Function to process names relating to interactions stored in bed12 formats.
#'
#' In bed(n) formats interaction data cannot be stored directly due to the inability of the format to handle
#' trans-interactions. It is only capable of storing cis-interactions through the use of blockStarts and blockEnds.
#' As such the information is typically stored as a string of the format chr1:start1..end1-chr2:start2..end2,N
#' (e.g. \"chr1:3268539..3269061-chr10:103124111..103124631,4\".). This function takes a vector of these strings
#' and converts them to a data.frame.
#'
#' @param x a vector of names stored in the bed file
#'
#' @return a data.frame containing all of the processed information
#'
#' @importFrom stringr str_split
.processChiapetName = function(x){
  tmp.str = str_split(x, ",")
  counts = as.integer(unlist(lapply(tmp.str, function(y){ return(y[2] ) } )))
  tmp.str.2 = sapply(tmp.str, function(y){
    str_split(y[1], "-")
  })
  paired.end.1 = t(sapply(tmp.str.2, function(y){
    tmp = str_split(y[[1]], ":")
    chr = tmp[[1]][1]
    tmp = str_split(tmp[[1]][2], "\\.\\.")
    start = as.integer(tmp[[1]][1]) + 1
    end = as.integer(tmp[[1]][2])
    return(c(chr, start, end))
  }))
  paired.end.2 = t(sapply(tmp.str.2, function(y){
    tmp = str_split(y[[2]], ":")
    chr = tmp[[1]][1]
    tmp = str_split(tmp[[1]][2], "\\.\\.")
    start = as.integer(tmp[[1]][1]) + 1
    end = as.integer(tmp[[1]][2])
    return(c(chr, start, end))
  }))
  dat = cbind(paired.end.1, paired.end.2, counts)
  colnames(dat) = c("chrom.left.", "start.left.", "end.left.", "chrom.right.", "start.right.", "end.right.", "counts")
  return(dat)
}


#' Function to read in processed Hi-C interaction data generated by HOMER
#'
#' This function reads in the interaction data outputted by HOMER (http://homer.salk.edu/homer/interactions/).
#'
#' @param fn location of data exported by HOMER
#' @return a data frame containing the relevant information
#'
.importHomer = function(fn){
  HOMER_int.df <- read.table(fn, header=TRUE, stringsAsFactors=FALSE, sep="\t")
  #convert to closed format
  #already 1-based
  HOMER_int.df$end.1. <- HOMER_int.df$end.1. - 1
  HOMER_int.df$end.2. <- HOMER_int.df$end.2. - 1
  return(HOMER_int.df)
}


#' Function to read in processed Hi-C interaction data generated by hiclib
#'
#' This function reads in the interaction data processed by [hiclib](https://bitbucket.org/mirnylab/hiclib). The data is extracted from hiclib and
#' exported from its internal HDF5 format to text using h5dictToTxt.py. This function reads in the relevant files present in the
#' directory and generates a data.table containing all of the information. It assumes that there are files with the following
#' names in the directory: fragids1, chrms1, mids1, fraglens1, fragids2, chrms2, mids2, fraglens2, distances.
#'
#' @param dir directory containing the output files generated by hiclib, extracted using h5dictToTxt.py
#' @param genome BSGenome to use for constructing the GenomicInteractions object.
#' @return data.table containing information on the individual interactions
#'
#' @import data.table
#'
.importHicLib = function(dir, genome){

  .process.chr = function(chrs, genome){

    chrs = as.numeric(chrs)+1
    chrs = ifelse( names(seqlengths(genome)[chrs]) %in% c("chrY"), "Y", chrs)
    chrs = ifelse( names(seqlengths(genome)[chrs]) %in% c("chrX"), "X", chrs)
    return( paste0("chr", chrs ))
  }

  frags=data.table(fragid1=fread(paste0(dir, "fragids1"))$V1, chrm1=.process.chr(fread(paste0(dir, "chrms1"))$V1, genome),
                   mid1=fread(paste0(dir, "mids1"))$V1, fraglength1 = fread(paste0(dir, "fraglens1"))$V1,
                   fragid2=fread(paste0(dir, "fragids2"))$V1, chrm2=.process.chr(fread(paste0(dir, "chrms2"))$V1, genome),
                   mid2=fread(paste0(dir, "mids2"))$V1, fraglength2 = fread(paste0(dir, "fraglens2"))$V1,
                   distances=fread(paste0(dir, "distances"))$V1, stringsAsFactors=FALSE)
  frags_agg <- frags[, .N, by=names(frags)]

  return(frags_agg)
}


#' Function to read in interaction-data stored in a BAM file
#'
#' Reads in interactions stored in a BAM file. Assumes that each interaction is represented by pair of PETs
#' with the same qname. The function reads in and determines which anchor is which by examining the
#' isFirstMateRead and isSecondMateRead fields in the BAM file.
#'
#' @param fn name of BAM file containing interaction information
#' @param genome BSGenome to use for constructing the GenomicInteractions object.
#' @return list of GRanges - storing the anchor information for each interaction
#'
#' @importFrom Rsamtools scanBamFlag ScanBamParam scanBam
#' @import GenomicRanges
#'
.readBam = function(fn, genome){
  bf = scanBamFlag(isPaired = TRUE, isDuplicate=FALSE)
  param = ScanBamParam(flag=bf, what=c("rname", "qname", "strand", "pos", "seq", "cigar", "flag") )
  b = scanBam(fn, param=param)

  y = GRanges(b[[1]]$rname, IRanges(b[[1]]$pos, width=width(b[[1]]$seq)), strand = b[[1]]$strand,
                                                                              qname = b[[1]]$qname,
                                                                              bamFlagAsBitMatrix(b[[1]][["flag"]], bitnames="isFirstMateRead"),
                                                                              bamFlagAsBitMatrix(b[[1]][["flag"]], bitnames="isSecondMateRead")    ,
                                                                              seqlengths=seqlengths(genome)
              )

  y1 = y[ y$isFirstMateRead == 1 ]
  y2 = y[ y$isSecondMateRead == 1]

  y1 = y1[ y1$qname %in% unique(y1$qname)]
  y2 = y2[ y2$qname %in% unique(y2$qname)]

  y1 = y1[ order(y1$qname) ]
  y2 = y2[ order(y2$qname) ]
  y1$isFirstMateRead = NULL
  y1$isSecondMateRead = NULL
  y2$isSecondMateRead = NULL
  y2$isFirstMateRead = NULL
  return(list(y1, y2))

}

#' Function to read in interaction-data stored in a pair of BAM files
#'
#' Reads in interactions stored in a a pair of BAM files, e.g. from independent
#' alignment of paired-end reads. Assumes that each interaction is represented
#' by pair of PETs with the same qname (this may not always be true. Depending
#' on data origin read qnames may end in '/1' or '/2' to denote first or second
#' read in the pair). The function reads in files, removes unpaired reads, and
#' pairs reads based on macthing qnames.
#'
#' @param fn Character vector of two BAM files with aligned reads.
#' @param genome BSGenome to use for constructing the GenomicInteractions object.
#' @return list of two GRanges, storing the anchor information for each interaction
#'
#' @importFrom Rsamtools scanBamFlag ScanBamParam scanBam
#' @import GenomicRanges
#'
.readTwoBams = function(fn, genome){

    if (length(fn)!=2){
        stop("Must supply two bam files for this import method")
    }

    bf = scanBamFlag(isUnmappedQuery=FALSE)
    param = ScanBamParam(flag=bf, what=c("rname", "qname", "strand", "pos", "seq"))

    message("Reading first bam file...")
    b1 = scanBam(fn[1], param=param)
    g1 = GRanges(as.character(b1[[1]]$rname), IRanges(b1[[1]]$pos, width=width(b1[[1]]$seq)), strand = as.character(b1[[1]]$strand),
                 qname = b1[[1]]$qname)
    seqinfo(g1) <- seqinfo(genome)
    rm(b1)

    message("Reading second bam file...")
    b2 = scanBam(fn[2], param=param)
    g2 = GRanges(as.character(b2[[1]]$rname), IRanges(b2[[1]]$pos, width=width(b2[[1]]$seq)), strand = as.character(b2[[1]]$strand),
                 qname = b2[[1]]$qname)
    seqinfo(g2) <- seqinfo(genome)
    rm(b2)

    message("Removing unpaired reads...")
    g1 <- g1[g1$qname %in% g2$qname]
    g2 <- g2[g2$qname %in% g1$qname]

    message("Pairing reads...")
    g1 <- g1[order(g1$qname)]
    g2 <- g2[order(g2$qname)]

    return(list(g1, g2))
}

