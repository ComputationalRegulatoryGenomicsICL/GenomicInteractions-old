% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/annotators.R
\docType{methods}
\name{resetAnnotations}
\alias{resetAnnotations}
\alias{resetAnnotations,GInteractions-method}
\title{Reset annotations made to a GInteractions object}
\usage{
resetAnnotations(GIObject)

\S4method{resetAnnotations}{GInteractions}(GIObject)
}
\arguments{
\item{GIObject}{An annotated GInteractions object}
}
\value{
invisible(1)
}
\description{
This function removes all annotations from a GInteractions object by
deleting  all of the metadata columns associated with both anchors.
}
\examples{
data(hic_example_data)
resetAnnotations(hic_example_data)
}

