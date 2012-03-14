\name{limmaAnalysis-methods}
\docType{methods}
\alias{limmaAnalysis-methods}
\alias{limmaAnalysis,ExpressionSet,matrix,character-method}
\alias{limmaAnalysis,matrix,matrix,character-method}
\title{Methods for Function \code{limmaAnalysis} in Package \pkg{birta}}
\description{
Generic methods for function \code{limmaAnalysis} in package \pkg{birta}. The expression data can be stored in a matrix or an ExpressionSet.
}
\section{Methods}{
\describe{

\item{\code{signature(dat = "ExpressionSet", design = "matrix", contrasts = "character")}}{
Generic function for objects of class ExpressionSet.
}

\item{\code{signature(dat = "matrix", design = "matrix", contrasts = "character")}}{
Generic function for objects of class matrix.
}
}}
\keyword{methods}
\keyword{manip}
