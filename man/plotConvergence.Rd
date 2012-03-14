
\name{plotConvergence}
\alias{plotConvergence}
\title{ Plotting the likelihood along MCMC sampling. }
\description{
Plots the log likelihood along MCMC sampling.
}
\usage{
plotConvergence(res, nburnin=NULL, title="")
}
\arguments{
\item{res}{The result from birta.run (a list).}
\item{nburnin}{Number of iterations used for the burn in.}
\item{title}{Optional title of the plot.}
}
\author{ Benedikt Zacher \email{zacher@lmb.uni-muenchen.de}}
\seealso{\code{\link[birta]{birta}}}
\examples{
data(humanSim)
data(humanSim)
design = model.matrix(~0+factor(c(rep("control", 5), rep("treated", 5))))
colnames(design) = c("control", "treated")
contrasts = "treated - control"
limmamRNA = limmaAnalysis(sim$dat.mRNA, design, contrasts)
limmamiRNA = limmaAnalysis(sim$dat.miRNA, design, contrasts)
sim_result = birta(sim$dat.mRNA, sim$dat.miRNA, limmamRNA=limmamRNA, 
 limmamiRNA=limmamiRNA, nrep=c(5,5,5,5), genesets=genesets, 
 model="all-plug-in", niter=50000, nburnin=10000, 
 sample.weights=FALSE, potential_swaps=potential_swaps)
plotConvergence(sim_result, nburnin=10000, title="simulation")
}
\keyword{hplot}
