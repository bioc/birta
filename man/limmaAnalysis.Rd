\name{limmaAnalysis}
\alias{limmaAnalysis}
\title{
Perform a limma analysis on expression data.
}
\description{
Needed for the parameterization of the model as well as the Fisher (pre-)test.
}
\usage{
limmaAnalysis(dat, design, contrasts)
}
\arguments{
  \item{dat}{A matrix or ExpressionSet containing the expression values.}
  \item{design}{A design matrix.}
  \item{contrasts}{Contrast for the linear model.}
}
\value{
Returns a list containing the following entries:
  \item{pvalue.tab}{Containing the result of the topTable function from limma.}
  \item{lm.fit}{Linear fit to the model.}
  \item{design}{The design used.}
  \item{contrast}{The contrasts used.}
}
\references{
G. K. Smyth. Limma : Linear Models for Microarray Data. \emph{Bioinformatics}, (2005):397-420.
}
\author{
Benedikt Zacher (zacher@lmb.uni-muenchen.de), Holger Froehlich
}
\seealso{\code{\link[limma]{lmFit}}, \code{\link[limma]{topTable}}}
\examples{
data(humanSim)
design = model.matrix(~0+factor(c(rep("control", 5), rep("treated", 5))))
colnames(design) = c("control", "treated")
contrasts = "treated - control"
limmamRNA = limmaAnalysis(sim$dat.mRNA, design, contrasts)
}
\keyword{manip}
