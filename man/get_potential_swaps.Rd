\name{get_potential_swaps}
\alias{get_potential_swaps}
\title{
Calculate swap partner for TF-/miRNA-target graph.
}
\description{
Calculates for TF-/miRNA-target graph all potential swap partner.
}
\usage{
get_potential_swaps(genesetsTF=NULL, genesetsmiRNA=NULL, perc.overlap.cutoff = 0.8, integer.id=TRUE, verbose=TRUE) 
}
\arguments{
  \item{genesetsTF}{Each entry corresponds to a TF and contains its target genes.}
  \item{genesetsmiRNA}{Each entry corresponds to a miRNA and contains its target genes.}
  \item{perc.overlap.cutoff}{Percentage cutoff of minimal overlap between two miRNAs or TFs to be possible swap partner.}
  \item{integer.id}{If TRUE, the swap partner are not outputted as characters, but as integer indices.}
  \item{verbose}{print process or not.}
}
\value{
The function returns a list, with the following entries:
  \item{T_potential_swaps}{Contains for each TF all potential swap partner.}
  \item{S_potential_swaps}{Contains for each miRNA all potential swap partner.}
}
\author{
Benedikt Zacher (zacher@lmb.uni-muenchen.de)
}
\examples{
genesetsTF = list("V$1"=c("A", "B", "C", "D", "E"), "V$2"=c("A", "B", "C", "D"), "V$3"=c("A", "B", "C"))
genesetsmiRNA =list("miR-1"=c("C", "E", "D", "G", "H"), "miR-2"=c("C", "E", "D", "G"), "miR-3"=c("C", "E"))
get_potential_swaps(genesetsTF, genesetsmiRNA, integer.id=FALSE, perc.overlap.cutoff=0.7)
}
\keyword{manip}
