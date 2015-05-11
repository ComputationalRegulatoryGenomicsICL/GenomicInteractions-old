
# so this is a custom track maybe it needs to inherit from RangeTrack
setClass("InteractionTrack",
         contains=c("GdObject"),
         representation=representation(plottingFunction="function",
                                       variables="list",
                                       stacking="character"),
         prototype=prototype(dp=DisplayPars()))

setMethod("initialize", "InteractionTrack", function(.Object, plottingFunction, variables, ...) {
  #.Object <- .updatePars(.Object, "InteractionTrack") #### need to sort this out
  .Object@plottingFunction <- plottingFunction
  .Object@variables <- variables
  .Object <- callNextMethod(.Object, ...)
  return(.Object)
})

InteractionTrack <- function(plottingFunction=function(GdObject, prepare=FALSE, ...){}, variables=list(), name="CustomTrack",  ...){
  print(plottingFunction)
  return(new("InteractionTrack", plottingFunction=plottingFunction, variables=variables, name=name, ...))
}


setMethod("start", "InteractionTrack", function(x){return(1)} )
setMethod("end", "InteractionTrack", function(x){return(2000)} )
setMethod("subset", signature(x="InteractionTrack"), function(x, from, to, chromosome, ...){})

#drawAxis

# pass it is GI object, interaction.strength whether to use the height of the arc or the width of the arc to show interaction strength or none or neither
# colour.by either count, padj or p-value, others are self-explanatory
InteractionTrack <- function(name="InteractionTrack", giobject=giobject, ...
                             #giobject,
                             #interaction.strength="width", 
                             #colour.by="count",
                             #plot.cis = TRUE,
                             #plot.trans = FALSE,
                             #plot.anchors = FALSE,
                             #col = "black", prepare=FALSE
                             ){
	return(new("InteractionTrack", name=name, plottingFunction=drawGD, variables=list(giobject=giobject)))	
}



setMethod("drawGD", signature("InteractionTrack"), function(GdObject, minBase, maxBase, prepare=FALSE, ...){ 
    print(GdObject)
    print(minBase)
    print(maxBase)
    
    draw.anchors = TRUE
    anchors.col = "lightblue"
    #interaction.counts = "height" or "width"
    # do I need an axis?

    pushViewport(viewport(xscale=c(minBase, maxBase), yscale=c(0, 1))  )
    
    print(dev.size(units="px"))

    print(GdObject@variables$giobject)
    
    anchor_one_chr = seqnames(anchorOne(GdObject@variables$giobject))
    anchor_one_starts = start(anchorOne(GdObject@variables$giobject))
    anchor_one_ends = end(anchorOne(GdObject@variables$giobject))
    
    anchor_two_chr = seqnames(anchorTwo(GdObject@variables$giobject))
    anchor_two_starts = start(anchorTwo(GdObject@variables$giobject))
    anchor_two_ends = end(anchorTwo(GdObject@variables$giobject))
    
    anchor_one_midpoints = (anchor_one_starts + anchor_one_ends) / 2
    anchor_two_midpoints = (anchor_two_starts + anchor_two_ends) / 2
    
    # how can I fix this to work properly
    print(anchor_one_chr)
    
    # TODO - fix this!!!
    lwds = 1+log(interactionCounts(x))
    
    for(i in 1:length(GdObject@variables$giobject)){
      grid.bezier(c(rep(anchor_one_midpoints[i], 2), rep(anchor_two_midpoints[i], 2)), 
                  c(0.05,1, 1,0.05),  default.units="native", gp=gpar(col="red",lwd=lwds))
    }
    
    # add anchors 
    anchors = unique(c(anchorOne(GdObject@variables$giobject), anchorTwo(GdObject@variables$giobject)))
    if(draw.anchors){
      for(i in 1:length(anchors)){
        grid.polygon(
                    c(start(anchors[i]), end(anchors[i]), end(anchors[i]), start(anchors[i])), 
                    c(0, 0, 0.1,0.1),
                    default.units="native", gp=gpar(col=anchors.col, fill= anchors.col))
      }
    }
    #y = c(0, 0.2509799, 0.2509799, 0) 
    #x = c(700, 700, 1500,1500) 
    
    #grid.bezier(x, y,  default.units="native", gp=gpar(col="red",lwd=1))
    # transinteractions
    # anchor regions
    # anchor colour
    # line heights or line widths specify number of interactions 

    #tmp <- GdObject@plottingFunction(GdObject, prepare=prepare)
    #)
    #if(!is(tmp, "CustomTrack")){
    #    warning("The plotting function of a CustomTrack has to return the input object. Using the original CustomTrack object now.")
    #}else{
    #    GdObject <- tmp
    #}
    popViewport(1)
    return(invisible(GdObject))
})


xy <- coords(res$foo)
x = xy[,1] + (xy[,3]-xy[,1])/2
y = xy[,2] + (xy[,4]-xy[,2])/2
pushViewport(viewport(xscale=c(0, dev.size(units="px")[1]), yscale=c(dev.size(units="px")[2], 0)), gp=gpar(col="red",lwd=10)) # this alters the colour and line width
grid.bezier(rep(x[3:4], each=2), c(y[2], rep(min(xy[,2]), 2), y[3]), default.units="native", gp=gpar(col="red",lwd=10))


#' Plot interactions within a specified region. 
#' 
#' This is function allows the plotting of an interactions between annotated features in a specified area. 
#' The resulting plot shows unique interactions as curves between interaction anchor points with the number 
#' of counts supporting that interaction proportional to the thickness of that line. It is also possible to 
#' add cis-interactions which are not within the window/region and to also plot regions that are involved
#' in trans-interactions. Plotting the data this way makes it possible to examine a region and easily examine
#' which regions are highly interacting with each other. It is not recommended to use this style of plot to examine
#' regions larger than 5Mb. 
#' 
#' 
