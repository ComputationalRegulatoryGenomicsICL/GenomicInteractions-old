
# so this is a custom track maybe it needs to inherit from RangeTrack
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
                                            plot.cis = TRUE,
                                            col.cis = "red",
                                            plot.trans = FALSE,
                                            col.trans = "lightgray"
                                            )))

setMethod("initialize", "InteractionTrack", function(.Object, plottingFunction, chromosome, variables, ...) {
  #.Object <- .updatePars(.Object, "InteractionTrack") #### need to sort this out
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
  print(tmp.start)
  return(tmp.start)
  } )

setMethod("end", "InteractionTrack", function(x){ 
  if(!is.null(x@variables$end)){
    tmp.end = max(c(end(anchorOne(x@variables$giobject))[ seqnames(anchorOne(x@variables$giobject)) ==  x@variables$chromosome ],
                    end(anchorTwo(x@variables$giobject))[ seqnames(anchorTwo(x@variables$giobject)) ==  x@variables$chromosome ]))
  }else{
    return(x@variables$start)
  }
  print(tmp.end)
  return(tmp.end)
})

setMethod("chromosome", "InteractionTrack", function(GdObject) GdObject@chromosome)

setMethod("subset", signature(x="InteractionTrack"), function(x, from, to, chromosome, ...){
  x@variables$giobject = subsetByFeatures(x@variables$giobject, GRanges(chromosome, IRanges(from, to)))
  return(x)                                          
})

#' Plot interactions within a specified region. 
#' 
#' This is function allows the plotting of an interactions between annotated features in a specified area. 
#' The resulting plot shows unique interactions as curves between interaction anchor points with the number 
#' of counts supporting that interaction proportional to the thickness of that line. It is also possible to 
#' add cis-interactions which are not within the window/region and to also plot regions that are involved
#' in trans-interactions. Plotting the data this way makes it possible to examine a region and easily examine
#' which regions are highly interacting with each other.
# pass it is GI object, interaction.strength whether to use the height of the arc or the width of the arc to show interaction strength or none or neither
# colour.by either count, padj or p-value, others are self-explanatory
#' 
#' 
InteractionTrack <- function(name="InteractionTrack", giobject, chromosome, start=NULL, end=NULL){
	return(new("InteractionTrack", name=name, plottingFunction=drawGD, chromosome=chromosome, variables=list(giobject=giobject)))	
}

setMethod("drawGD", signature("InteractionTrack"), function(GdObject, minBase, maxBase, prepare=FALSE, ...){ 
    #print(GdObject)
    print(minBase)
    print(maxBase)
    print(GdObject@variables$giobject)
    if(prepare){
      displayPars(GdObject) <- list(".__verticalSpace"=5)
      return(invisible(GdObject))
    }

    pushViewport(viewport(xscale=c(minBase, maxBase), yscale=c(0, 1))  )
    
    anchor_one_chr = as.character(seqnames(anchorOne(GdObject@variables$giobject)))
    anchor_one_starts = start(anchorOne(GdObject@variables$giobject))
    anchor_one_ends = end(anchorOne(GdObject@variables$giobject))
    
    anchor_two_chr = as.character(seqnames(anchorTwo(GdObject@variables$giobject)))
    anchor_two_starts = start(anchorTwo(GdObject@variables$giobject))
    anchor_two_ends = end(anchorTwo(GdObject@variables$giobject))
    
    anchor_one_midpoints = (anchor_one_starts + anchor_one_ends) / 2
    anchor_two_midpoints = (anchor_two_starts + anchor_two_ends) / 2
    
    # TODO - fix this!!!
    if(displayPars(GdObject, "lwd.interaction") == "counts"){
      
      max(interactionCounts(GdObject@variables$giobject))
      
      lwds = 2*(interactionCounts(GdObject@variables$giobject)/max(interactionCounts(GdObject@variables$giobject)))
      
    }else if(displayPars(GdObject, "lwd.interaction") == "fdr"){
      lwds = 2-(GdObject@variables$giobject$fdr*2)
    }else if(displayPars(GdObject, "lwd.interaction") == "p.value"){
      lwds = 2-(GdObject@variables$giobject$p.value*2)
    }else{
      lwds = rep(1, length(anchor_one_chr))
    }
    
    col.interactions = displayPars(GdObject, "col.interactions")
    plot.cis = displayPars(GdObject, "plot.cis")
    col.cis =  displayPars(GdObject, "col.cis")
    plot.trans = displayPars(GdObject, "plot.trans")
    print(plot.trans)
    col.trans = displayPars(GdObject, "col.trans")
    
    # TODO order this so you go from trans to outside region to inside region
    for(i in 1:length(GdObject@variables$giobject)){
      print(lwds[i])
      print(interactionCounts(GdObject@variables$giobject)[i])
      
      if( anchor_one_chr[i] ==  anchor_two_chr[i] ){
        print(plot.cis)
        print(col.cis)
        if(plot.cis & (anchor_one_midpoints[i] <minBase || anchor_one_midpoints[i] > maxBase || 
                         anchor_two_midpoints[i] < minBase || anchor_two_midpoints[i] > maxBase)){
          print("HERE")
          grid.bezier(c(rep(anchor_one_midpoints[i], 2), rep(anchor_two_midpoints[i], 2)), 
                      c(0.05,1, 1,0.05),  default.units="native", gp=gpar(col=col.cis,lwd=lwds[i]))
        }else{
          grid.bezier(c(rep(anchor_one_midpoints[i], 2), rep(anchor_two_midpoints[i], 2)), 
                    c(0.05,1, 1,0.05),  default.units="native", gp=gpar(col=col.interactions,lwd=lwds[i]))
        }
      }else if(plot.trans){
        if( anchor_one_chr[i] != GdObject@chromosome ){
          if(anchor_two_midpoints[i] > ((minBase + maxBase)/2)){
            grid.curve(x1=anchor_two_midpoints[i], y1=0.05, x2=maxBase, y2=0.95,  curvature = -1, default.units="native", gp=gpar(col=col.trans,lwd=lwds[i]))
          }else{
            grid.curve(x1=anchor_two_midpoints[i], y1=0.05, x2=minBase, y2=0.95,  curvature = 1, default.units="native", gp=gpar(col=col.trans,lwd=lwds[i]))
          }      
        }else if( anchor_two_chr[i] != GdObject@chromosome){
          if(anchor_one_midpoints[i] > ((minBase + maxBase)/2)){
            grid.curve(x1=anchor_one_midpoints[i], y1=0.05, x2=maxBase, y2=0.95, curvature = -1, default.units="native", gp=gpar(col=col.trans,lwd=lwds[i]))
          }else{
            grid.curve(x1=anchor_one_midpoints[i], y1=0.05, x2=minBase, y2=0.95, curvature = 1, default.units="native", gp=gpar(col=col.trans,lwd=lwds[i]))
          }
        }
      }  
    }
    
    plot.anchors = displayPars(GdObject, "plot.anchors")
    col.anchors = displayPars(GdObject, "col.anchors")
    if(plot.anchors){
      anchors = unique(c(anchorOne(GdObject@variables$giobject), anchorTwo(GdObject@variables$giobject)))
      for(i in 1:length(anchors)){
        grid.polygon(
                    c(start(anchors[i]), end(anchors[i]), end(anchors[i]), start(anchors[i])), 
                    c(0, 0, 0.1,0.1),
                    default.units="native", gp=gpar(col= "black", fill= col.anchors))
      }
    }
    
    popViewport(1)
    return(invisible(GdObject))
})
