
# so this is a custom track maybe it needs to inherit from RangeTrack

setClass("InteractionTrack",
	contains=c("GdObject"),
	representation=representation(plottingFunction="function", variables="list"),
	prototype=prototype(dp=Gviz::DisplayPars()))

setMethod("initialize", "InteractionTrack", function(.Object, ...){
			#.Object = .updatePars(.Object, "InteractionTrack")
			#.Object@variables <- variables
			.Object <- callNextMethod(.Object, ...)
			return(.Object)
})

setMethod("start", "InteractionTrack", function(return(NA)) )
setMethod("end", "InteractionTrack", function(return(NA)) )
setMethod("subset", signature(x="InteractionTrack"), function(x, from, to, chromosome, ...){})

drawAxis


InteractionTrack <- function(name="InteractionTrack", ...){
	return(new("InteractionTrack", name=name, ...))	
}

test = InteractionTrack("SPAM")

setMethod("drawGD", signature("InteractionTrack"), function(GdObject, minBase, maxBase, prepare=FALSE, ...){ 
    print(GdObject)
    print(minBase)
    print(maxBase)
    #interaction.counts = "height" or "width"

    # do I need an axis?

    #rev <- .dpOrDefault(GdObject, "reverseStrand", FALSE)
    xscale <- c(minBase, maxBase)

    #pushViewport(viewport(xscale=xscale, clip=TRUE))
    pushViewport(viewport(xscale=c(0, dev.size(units="px")[1]), yscale=c(0, 1))	)
    print(dev.size(units="px"))

    x = c(41.76195, 74.33861, 271.77287, 666.64138) 
    y = c(0.2509799, 0.2509799, 0.2509799, 0.2509799) 

    grid.bezier(rep(x[3:4], each=2), c(y[2], rep(min(xy[,2]), 2), y[3]),  default.units="native", gp=gpar(col="red",lwd=10))

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



annTrack <- AnnotationTrack(start=c(1, 100, 700, 1900), width=50, chromosome=1, genome="hg19", name="foo", id=letters[1:4])

res <- plotTracks(list(GenomeAxisTrack(), test, annTrack))
xy <- coords(res$foo)
x = xy[,1] + (xy[,3]-xy[,1])/2
y = xy[,2] + (xy[,4]-xy[,2])/2
pushViewport(viewport(xscale=c(0, dev.size(units="px")[1]), yscale=c(dev.size(units="px")[2], 0)), gp=gpar(col="red",lwd=10)) # this alters the colour and line width
grid.bezier(rep(x[3:4], each=2), c(y[2], rep(min(xy[,2]), 2), y[3]), default.units="native", gp=gpar(col="red",lwd=10))


.calcCoords = function(coordinates, plot.width, block.start, block.width){
  return(plot.width * (( coordinates - block.start) / block.width))
}

.calcYCoords = function(score, max.score, bottom, top){
  return( (score / max.score) * (top - bottom))
}

.calcMdpt = function(start, end){
  foo = ifelse( end > start, start + ((end - start) / 2),  end + (( start - end) / 2 ))
  return(foo)	 
}

.bezier.curve = function(p0x, p1x, p2x, n){
  t = seq(0,1,1/n)
  return( ((1 - t)^2 * p0x) + (2 * ( 1 - t ) * t * p1x ) + ( t^2 * p2x ) )
}

.specialrect = function(x.start, y.start, x.end, y.end, strand, col){
    for(i in 1:length(x.start)){
        if(as.character(strand[i])=="*"){
            rect(x.start[i], y.start[i], x.end[i], y.end[i], col=col)
        }else if(as.character(strand[i])=="+"){
            feature.width = x.end[i] - x.start[i]
            feature.height = y.end[i] - y.start[i]
            turn = x.start[i] + 0.5 * sqrt(3) * (0.9 * feature.width)
            xs = c(x.start[i], turn, x.end[i], turn, x.start[i])
            ys = c(y.start[i], y.start[i], y.start[i] + feature.height/2 , y.end[i], y.end[i])
            polygon(xs, ys, col=col)
        }else if(as.character(strand[i])=="-"){
            feature.width = x.end[i] - x.start[i]
            feature.height = y.end[i] - y.start[i]
            turn = x.end[i] - 0.5 * sqrt(3) * (0.9 * feature.width)
            xs = c(x.end[i], turn, x.start[i], turn, x.end[i])
            ys = c(y.start[i], y.start[i], y.start[i] + feature.height/2 , y.end[i], y.end[i]) 
            polygon(xs, ys, col=col)
        }   
    }
}



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
