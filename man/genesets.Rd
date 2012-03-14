\name{genesets}
\docType{data}
\alias{genesets}
\title{TF-target and miRNA-target networks for simulated example.}
\description{
For miRNAs we looked at target predictions in human via miRanda (Betel et al., 2008), miRBase (Griffths-Jones et al., 2008) and miRDB (Wang and Naqa, 2008). Experimentally confirmed targets were retrieved from Tarbase (Papadopoulos et al., 2009). As trustworthy miRNA targets we considered those, which were either experimentally confirmed (i.e. listed in Tarbase) or predicted by at least two of the other three methods. In conclusion we arrived at a network with 583 miRNAs regulating between 1 and 1628 genes (median: 302).\\
A TF-target gene network was compiled by computing TF binding affinities to promoter sequences of all human genes according to the TRAP model (Roideret al., 2007) via the author's R implementation. Upstream sequences of genes were retrieved here from the ENSEMBL database via biomaRt (Haider et al., 2009). We assumed that promoter sequences were located in the range 0 - 2Kbp upstream to the transcription start site of a gene. As trustworthy TF targets we considered those, for which a Holm corrected affinity p-value smaller than 0.001 was reported. In conclusion we arrived at a network with 272 TFs regulating between 1 and 16517 genes (median: 20).}
\usage{genesets}
\format{A list  containing a TF- and miRNA-target graph.}
\source{The networks were generated using the TRAP model on TRANSFAC matrices and miRNA-target annotations form different databases.}
\references{
D. Betel, M. Wilson, A. Gabow, D. S. Marks, and C. Sander. The microrna.org resource: targets and expression. \emph{Nucleic Acids Res}, 36(Database issue):D149-D153, Jan 2008.

S. Griffiths-Jones, H. K. Saini, S. van Dongen, and A. J. Enright. miRBase: tools for microRNA genomics. \emph{Nucleic Acids Res}, 36(Database issue):D154-D158, Jan. 2008.

X. Wang and I. M. E. Naqa. Prediction of both conserved and nonconserved microrna targets in animals. \emph{Bioinformatics}, 24(3):325-332, Feb 2008.

G. L. Papadopoulos, M. Reczko, V. A. Simossis, P. Sethupathy, and A. G. Hatzigeorgiou. The database of experimentally supported targets: a functional update of tarbase. \emph{Nucleic Acids Res}, 37(Database issue):D155-D158, Jan 2009.
}
\keyword{datasets}
