## Definition of GenomicInteractions

#' A class to hold chromatin interaction data for a specific genomic region.
#'
#'  @slot plottingFunction function
#'  @slot variables list
#'  @slot chromosome
#'  @slot stacking character
#'
#' InteractionTrack is a specific Gviz-derived class for enabling the visualisation of chromatin interaction data. 
#' The InteractionTrack class allows interactions on a specified chromosome to be visualised by examining interactions
#' between anchors as bezier curves.
#' 
#' Several additional display parameters (i.e. \rcode{displayPars(foo)=list(...) }are defined for this class, including 
#' \code{plot.anchors} which can be used to specify whether anchors are to be drawn. \code{col.anchors} which can be used 
#' to alter the colour of these anchor elements. The value of \code{plot.outside} determines whether or not interactions
#' which span outside of the window are to be plotted, and \code{col.outside} defines the colour of these interactions. 
#' Similarly \code{plot.trans} determines whether trans-interactions are plotted and \code{col.trans} specifies the colour
#' of trans-interactions. By default, the line width of an arc representing an interaction is proportional to the number 
#' of reads/counts supporting that interaction. Instead of using the counts to define this, the line width can be set to
#' be proportion to either \code{fdr} or \code{p.value} using the \code{lwd.interactions} display parameter. 
#' 
#'
#'
#' @import Gviz
#' @import grid
#'
#' @export InteractionTrack
setClass("InteractionTrack",
         contains=c("GdObject"),
         representation=representation(plottingFunction="function",
                                       variables="list",
                                       chromosome="character",
                                       stacking="character"),
         prototype=prototype(dp=DisplayPars(plot.anchors=TRUE,
                                            col.anchors = "lightblue",
                                            interaction.strength="width", 
                                            lwd.interaction ="counts",
                                            col.interactions = "red",
                                            plot.outside = TRUE,
                                            col.outside = "red",
                                            plot.trans = FALSE,
                                            col.trans = "lightgray",
                                            col.interaction.types = c()
                                            )))


setMethod("initialize", "InteractionTrack", function(.Object, plottingFunction, chromosome, variables, ...) {
  #.Object <- .updatePars(.Object, "InteractionTrack") #
  .Object@plottingFunction <- plottingFunction
  .Object@chromosome <- chromosome
  .Object@variables <- variables
  .Object <- callNextMethod(.Object, ...)
  return(.Object)
})

setMethod("start", "InteractionTrack", function(x){
  
  if(!is.null(x@variables$start)){
    tmp.start = min(c(start(anchorOne(x@variables$giobject))[ seqnames(anchorOne(x@variables$giobject)) ==  x@variables$chromosome ],
        start(anchorTwo(x@variables$giobject))[ seqnames(anchorTwo(x@variables$giobject)) ==  x@variables$chromosome ]))
  }else{
    return(x@variables$start)
  }
  return(tmp.start)
  } )

setMethod("end", "InteractionTrack", function(x){ 
  if(!is.null(x@variables$end)){
    tmp.end = max(c(end(anchorOne(x@variables$giobject))[ seqnames(anchorOne(x@variables$giobject)) ==  x@variables$chromosome ],
                    end(anchorTwo(x@variables$giobject))[ seqnames(anchorTwo(x@variables$giobject)) ==  x@variables$chromosome ]))
  }else{
    return(x@variables$start)
  }
  return(tmp.end)
})

setMethod("chromosome", "InteractionTrack", function(GdObject) GdObject@chromosome)

setMethod("subset", signature(x="InteractionTrack"), function(x, from, to, chromosome, ...){
  x@variables$giobject = subsetByFeatures(x@variables$giobject, GRanges(chromosome, IRanges(from, to)))
  return(x)                                          
})

#' Constructor to create an InteractionTrack object
#'
#' Create InteractionTrack object from an GenomicInteractions object to visualise a specified chromosome.
#'
#' @param x A GenomicInteractions object
#' @param chromosome
#'
#' @return an InteractionTrack object
#'
#' @examples
#' 
#' library(GenomicRanges)
#'
#' anchor.one = GRanges(c("chr1", "chr1", "chr1", "chr1"), IRanges(c(10, 20, 30, 20), width=5))
#' anchor.two = GRanges(c("chr1", "chr1", "chr1", "chr2"), IRanges(c(100, 200, 300, 50), width=5))
#' interaction_counts = sample(1:10, 4)
#' test <- GenomicInteractions(anchor.one, anchor.two, experiment_name="test", 
#'                            description="this is a test", counts=interaction_counts)
#' interactions.track = InteractionTrack(name="Test", test, chromosome="chr1")                        
#' plotTracks(list(interactions.track))
#' 
#' @export
InteractionTrack <- function(name="InteractionTrack", x, chromosome, start=NULL, end=NULL){ # TODO SORT OUT START AND STOP HERE
	return(new("InteractionTrack", name=name, plottingFunction=drawGD, chromosome=chromosome, variables=list(giobject=x)))	
}

#' draws InteractionTrack
#' 
#' @export
setMethod("drawGD", signature("InteractionTrack"), function(GdObject, minBase, maxBase, prepare=FALSE, ...){ 
    
  print(GdObject)
  if(prepare){
      pushViewport(viewport(xscale=c(minBase, maxBase), yscale=c(0, 1) )) # TESTING
      popViewport(1)
      return(invisible(GdObject))
    }

    pushViewport(viewport(xscale=c(minBase, maxBase), yscale=c(0, 1)))
    if(length(GdObject@variables$giobject)>0){
      anchor_one_chr = as.character(seqnames(anchorOne(GdObject@variables$giobject)))
      anchor_one_starts = start(anchorOne(GdObject@variables$giobject))
      anchor_one_ends = end(anchorOne(GdObject@variables$giobject))
      
      anchor_two_chr = as.character(seqnames(anchorTwo(GdObject@variables$giobject)))
      anchor_two_starts = start(anchorTwo(GdObject@variables$giobject))
      anchor_two_ends = end(anchorTwo(GdObject@variables$giobject))
      
      anchor_one_midpoints = (anchor_one_starts + anchor_one_ends) / 2
      anchor_two_midpoints = (anchor_two_starts + anchor_two_ends) / 2
      
      if(displayPars(GdObject, "interaction.strength")=="width" && displayPars(GdObject, "lwd.interaction") == "counts"){
        lwds = 2*(interactionCounts(GdObject@variables$giobject)/max(interactionCounts(GdObject@variables$giobject)))
      }else if(displayPars(GdObject, "interaction.strength")=="width" && displayPars(GdObject, "lwd.interaction") == "fdr"){
        lwds = 2-(GdObject@variables$giobject$fdr*2)
      }else if(displayPars(GdObject, "interaction.strength")=="width" && displayPars(GdObject, "lwd.interaction") == "p.value"){
        lwds = 2-(GdObject@variables$giobject$p.value*2)
      }else{
        lwds = rep(1, length(anchor_one_chr))
      }
    
      if(length(displayPars(GdObject, "col.interaction.types"))>0){
        col.interaction.types = displayPars(GdObject, "col.interaction.types")
        print(col.interaction.types)
        names(col.interaction.types)
        colour.map = str_split(names(col.interaction.types), "-")
        cols = rep("black", length(anchor_two_midpoints))
        for(i in 1:length(colour.map)){
          cols[ isInteractionType(GdObject@variables$giobject, colour.map[[1]][1], colour.map[[1]][2]) ] = col.interaction.types[i]
        }
        print(cols)
      }
      
      col.interactions = displayPars(GdObject, "col.interactions")
      plot.outside = displayPars(GdObject, "plot.outside")
      col.outside =  displayPars(GdObject, "col.outside")
      plot.trans = displayPars(GdObject, "plot.trans")
      col.trans = displayPars(GdObject, "col.trans")
      
      trans.indexes = which( anchor_one_chr != GdObject@chromosome | anchor_two_chr != GdObject@chromosome )
      outside.indexes = which( (anchor_one_chr == GdObject@chromosome & anchor_two_chr == GdObject@chromosome) & (anchor_one_midpoints <minBase | anchor_one_midpoints > maxBase | 
         anchor_two_midpoints < minBase | anchor_two_midpoints > maxBase) )
      inside.indexes = 1:length(anchor_one_chr)
      inside.indexes = which(!(inside.indexes %in% c(trans.indexes, outside.indexes)))
      
      cols = rep(col.interactions, length(anchor_one_chr))
      cols[trans.indexes] = col.trans
      cols[outside.indexes] = col.outside
      
      if(length(displayPars(GdObject, "col.interaction.types"))>0){
        col.interaction.types = displayPars(GdObject, "col.interaction.types")
        print(col.interaction.types)
        names(col.interaction.types)
        colour.map = str_split(names(col.interaction.types), "-")
        print(colour.map)
        for(i in 1:length(colour.map)){
          print(sum(isInteractionType(GdObject@variables$giobject, colour.map[[i]][1], colour.map[[i]][2])))
          print(which(isInteractionType(GdObject@variables$giobject, colour.map[[i]][1], colour.map[[i]][2])))
          cols[ isInteractionType(GdObject@variables$giobject, colour.map[[i]][1], colour.map[[i]][2]) ] = col.interaction.types[i]
          print(cols)
        }
      }
      print(cols)
      
      if(plot.trans){
        for(i in trans.indexes){
          if( anchor_one_chr[i] != GdObject@chromosome ){
            if(anchor_two_midpoints[i] > ((minBase + maxBase)/2)){
              grid.curve(x1=anchor_two_midpoints[i], y1=0.05, x2=maxBase, y2=0.95,  curvature = -1, default.units="native", gp=gpar(col=cols[i],lwd=lwds[i]))
            }else{
              grid.curve(x1=anchor_two_midpoints[i], y1=0.05, x2=minBase, y2=0.95,  curvature = 1, default.units="native", gp=gpar(col=cols[i],lwd=lwds[i]))
            }      
          }else if( anchor_two_chr[i] != GdObject@chromosome){
            if(anchor_one_midpoints[i] > ((minBase + maxBase)/2)){
              grid.curve(x1=anchor_one_midpoints[i], y1=0.05, x2=maxBase, y2=0.95, curvature = -1, default.units="native", gp=gpar(col=cols[i],lwd=lwds[i]))
            }else{
              grid.curve(x1=anchor_one_midpoints[i], y1=0.05, x2=minBase, y2=0.95, curvature = 1, default.units="native", gp=gpar(col=cols[i],lwd=lwds[i]))
            }
          }
        }
      }
      
      if(displayPars(GdObject, "interaction.strength")=="height" && displayPars(GdObject, "lwd.interaction") == "counts"){
        # or scale to the whole object?
        curve.heights = interactionCounts(GdObject@variables$giobject)/max(interactionCounts(GdObject@variables$giobject))
        curve.heights = rep(1, length(anchor_one_chr))
      }else if(displayPars(GdObject, "interaction.strength")=="height" && displayPars(GdObject, "lwd.interaction") == "fdr"){
        curve.heights = 1.0-(GdObject@variables$giobject$fdr)
      }else if(displayPars(GdObject, "interaction.strength")=="height" && displayPars(GdObject, "lwd.interaction") == "p.value"){
        curve.heights = 1.0-(GdObject@variables$giobject$p.value)
      }else{
        curve.heights = rep(1, length(anchor_one_chr))
      }
    
      if(plot.outside & length(outside.indexes)>0){
        xs = c()
        ys = c()
        ls = c()
        for(i in outside.indexes){
          xs = c(xs, c(rep(anchor_one_midpoints[i], 2), rep(anchor_two_midpoints[i], 2)))
          ys = c(ys, c(0.05,curve.heights[i], curve.heights[i],0.05))
          ls = c(ls, lwds[i])
        }  
        grid.bezier(xs, ys, id.lengths=rep(4, length(outside.indexes)), default.units="native", gp=gpar(col=cols[outside.indexes],lwd=ls)) # TODO
      }
    
      for(i in inside.indexes){
        grid.bezier(c(rep(anchor_one_midpoints[i], 2), rep(anchor_two_midpoints[i], 2)), 
                    c(0.05,curve.heights[i], curve.heights[i], 0.05),  default.units="native", gp=gpar(col=cols[i],lwd=lwds[i]))
      }
      
      plot.anchors = displayPars(GdObject, "plot.anchors")
      col.anchors = displayPars(GdObject, "col.anchors")
      if(plot.anchors & length(GdObject@variables$giobject) > 0){
        anchors = unique(c(anchorOne(GdObject@variables$giobject), anchorTwo(GdObject@variables$giobject)))
        xs = c()
        ys = c()
        for(i in 1:length(anchors)){
          xs = c(xs, c(start(anchors[i]), end(anchors[i]), end(anchors[i]), start(anchors[i])))
          ys = c(ys, c(0, 0, 0.1,0.1))
        }
        grid.polygon(xs, ys, id.lengths=rep(4, length(anchors)), default.units="native", gp=gpar(col= "black", fill= col.anchors))
      }
    }
    popViewport(1)
    return(invisible(GdObject))
})
