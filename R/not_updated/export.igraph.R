
#' Export interactions to an igraph object.
#'
#' Exports a GenomicInteractions object to graph.data.frame for use by igraph package. This uses unique anchors
#' as nodes and generates edges between them. For the resulting graph to be easily interpretable, anchors
#' should be non-overlapping. This should already be the case for HiC data (either binned or restriction
#' fragments), however ChIA-PET data can contain overlapping anchors, which may need to be reduced to
#' non-overlapping regions before graph export.
#' @param GIObject A GenomicInteractions object.
#'
#' @return a graph.data.frame representation of the GenomicInteractions object
#' @importFrom igraph graph.data.frame
#'
#' @export
#' @docType methods
#' @rdname export.igraph
#' @export
setGeneric("export.igraph",function(GIObject){standardGeneric ("export.igraph")})
#' @rdname export.igraph
#' @export
setMethod("export.igraph", "GenomicInteractions", function(GIObject){
  nodes = unique( c(anchorOne(GIObject), c(anchorTwo(GIObject))))
  gi.nodes.ol = findOverlaps(GIObject, nodes, type="equal")
  
  mapping = cbind(subjectHits(gi.nodes.ol[["one"]]), subjectHits(gi.nodes.ol[["two"]]), queryHits(gi.nodes.ol[["one"]]))
  names =  paste(paste(paste(as.character(seqnames(nodes))), as.character(start(nodes)), sep=":"), as.character(end(nodes)), sep="..")
  edge_mcols = mcols(GIObject)
  edge_mcols$count = interactionCounts(GIObject)
  edges = cbind(data.frame(from=names[mapping[,1]], to=names[mapping[,2]]), edge_mcols)
  
  verts = data.frame(name=names)
  if("node.class" %in% unique(c(names(elementMetadata(anchorOne(GIObject))), names(elementMetadata(anchorTwo(GIObject)))))){
    node.class = sapply(1:length(nodes), function(x){ paste(unique(c(anchorOne(GIObject)$node.class[unique(queryHits(gi.nodes.ol[["one"]][ subjectHits(gi.nodes.ol[["one"]]) == x ]))],
                                                                     anchorTwo(GIObject)$node.class[unique(queryHits(gi.nodes.ol[["two"]][ subjectHits(gi.nodes.ol[["two"]]) == x ]))])),
                                                            collapse=",")})
    verts = data.frame(names=names, nodeclass=node.class)
    potential.node.classes = unique(c(GIObject@anchor_one$node.class, GIObject@anchor_two$node.class))
    potential.node.classes = potential.node.classes[ potential.node.classes != "distal" ]
    for(i in potential.node.classes){
      verts[, paste(i, "id", sep=".")] = sapply(1:length(nodes),
                                                function(x){ paste(unique(c(
                                                  unlist(elementMetadata(anchorOne(GIObject))[[paste(i, "id", sep=".")]][
                                                    unique(queryHits(gi.nodes.ol[["one"]][ subjectHits(gi.nodes.ol[["one"]]) == x ]))
                                                    ]),
                                                  unlist(elementMetadata(anchorTwo(GIObject))[[paste(i, "id", sep=".")]][
                                                    unique(queryHits(gi.nodes.ol[["two"]][ subjectHits(gi.nodes.ol[["two"]]) == x ]))
                                                    ])
                                                )), collapse=",")
                                                })
    }
  }
  return(graph.data.frame(edges, directed=FALSE, vertices = verts))
})
