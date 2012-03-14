\name{birta-methods}
\docType{methods}
\alias{birta-methods}
\alias{birta,ExpressionSet,ExpressionSet,ExpressionSet-method}
\alias{birta,ExpressionSet,ExpressionSet,missing-method}
\alias{birta,ExpressionSet,missing,ExpressionSet-method}
\alias{birta,ExpressionSet,missing,missing-method}
\alias{birta,matrix,matrix,matrix-method}
\alias{birta,matrix,matrix,missing-method}
\alias{birta,matrix,missing,matrix-method}
\alias{birta,matrix,missing,missing-method}
\title{Methods for Function \code{birta} in Package \pkg{birta}}
\description{
Generic methods for function \code{birta} in package \pkg{birta}. The expression data can be stored in a matrix or an ExpressionSet.
}
\section{Methods}{
\describe{
The following signatures make sure, that the arguments, storing the expression data are in the correct format. 
\item{\code{signature(dat.mRNA = "ExpressionSet", dat.miRNA = "ExpressionSet", TFexpr = "ExpressionSet")}}{
%%  ~~describe this method here~~
}

\item{\code{signature(dat.mRNA = "ExpressionSet", dat.miRNA = "ExpressionSet", TFexpr = "missing")}}{
%%  ~~describe this method here~~
}

\item{\code{signature(dat.mRNA = "ExpressionSet", dat.miRNA = "missing", TFexpr = "ExpressionSet")}}{
%%  ~~describe this method here~~
}

\item{\code{signature(dat.mRNA = "ExpressionSet", dat.miRNA = "missing", TFexpr = "missing")}}{
%%  ~~describe this method here~~
}

\item{\code{signature(dat.mRNA = "matrix", dat.miRNA = "matrix", TFexpr = "matrix")}}{
%%  ~~describe this method here~~
}

\item{\code{signature(dat.mRNA = "matrix", dat.miRNA = "matrix", TFexpr = "missing")}}{
%%  ~~describe this method here~~
}

\item{\code{signature(dat.mRNA = "matrix", dat.miRNA = "missing", TFexpr = "matrix")}}{
%%  ~~describe this method here~~
}

\item{\code{signature(dat.mRNA = "matrix", dat.miRNA = "missing", TFexpr = "missing")}}{
%%  ~~describe this method here~~
}
}}
\keyword{methods}
\keyword{htest}