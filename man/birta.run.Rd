\name{birta}
\alias{birta}
\title{
Main interface for Bayesian Inference of Regulation of Transcriptional Activity.
}
\description{
The function estimates parameterization of the model and then executes MCMC sampling to infer activity states.
}
\usage{
birta(dat.mRNA, dat.miRNA, TFexpr, limmamRNA=NULL, limmamiRNA=NULL, limmaTF=NULL, nrep=NULL, fdr.mRNA=0.05, fdr.miRNA=0.05, lfc.mRNA=0, lfc.miRNA=0, genesets=NULL, lambda=NULL, sample.weights=TRUE, one.regulator.weight=TRUE, theta_TF=NULL, theta_miRNA=NULL, model=c("all-plug-in", "no-plug-in"), niter=500000, nburnin=100000, thin=50, potential_swaps=NULL, run.pretest=FALSE, condition.specific.inference=TRUE, only_switches=FALSE, weightSampleMean=0, weightSampleVariance=0.001, accessible=NULL)
}
\arguments{
  \item{dat.mRNA}{mRNA expression data (ExpressionSet or matrix). IMPORTANT: Replicates must be ordered according to nrep!}
  \item{dat.miRNA}{miRNA expression data (ExpressionSet or matrix).}
  \item{TFexpr}{TF expression data (ExpressionSet or matrix): can be a subset of dat.mRNA}
  \item{limmamRNA}{Output of limma analysis for mRNA data (list: pvalue.tab, lm.fit).}
  \item{limmamiRNA}{Output of limma analysis for miRNA data (list: pvalue.tab, lm.fit).}
  \item{limmaTF}{Output of limma analysis for TF data (list: pvalue.tab, lm.fit).}
  \item{nrep}{Vector containing four integers. Entry 1 and 2 are the number of miRNA measurement replicates for condition 1 and 2. Entry 3 and 4 are the number of mRNA measurement replicates for condition 1 and 2.}
  \item{fdr.mRNA}{False discovery rate (FDR) cutoff for significance of the logFC for mRNA data.}
  \item{fdr.miRNA}{False discovery rate (FDR) cutoff for significance of the logFC for miRNA data.}
  \item{lfc.mRNA}{Additional logFC cutoff for significance in mRNA data.}
  \item{lfc.miRNA}{Additional logFC cutoff for significance in miRNA data.}
  \item{genesets}{Combined TF / miRNA network. IMPORTANT: Names of TF entries must start with V\$.}
  \item{lambda}{Regularization parameter for edge weights.}
  \item{sample.weights}{Should edge weights be adapted during sampling?}
  \item{one.regulator.weight}{Should weights of all edges for a specific regulator be the same?}
  \item{theta_TF}{Expected fraction of active TFs.}
  \item{theta_miRNA}{Expected fraction of active miRNAs.}
  \item{model}{Type of model. One out of c("all-plug-in", "weight-plug-in", "no-plug-in").}
  \item{niter}{Number of MCMC iterations (AFTER burnin).}
  \item{nburnin}{Number of MCMC iterations UNTIL burnin is assumed to be finished.}
  \item{potential_swaps}{Pre-computed potential swaps (OPTIONAL, see get_potential_swaps).}
  \item{run.pretest}{Initialize miRNA and TF states via the result of a hypergeometric test in order to improve convergence (should be taken with care; advise: only use it in case of observed convergence problems!).}
  \item{condition.specific.inference}{Should inference on TF / miRNA activities be made only RELATIVE to a reference condition or independently in both conditions?}
  \item{thin}{Thinning of Markov chain: only use every thin's sample for posterior computation.}
  \item{only_switches}{Should only switches be performed?}
  \item{weightSampleMean}{Mean for normal distribution used for sampling the omegas.}
  \item{weightSampleVariance}{Variance for normal distribution used for sampling the omegas.}
  \item{accessible}{0/1 matrix of conditions x mRNAs: 1 = gene is accessible for TF, 0 = gene is not accessible (e.g. due to methylation)}
}
\value{
The function returns a list containing the following entries:
  \item{miRNAstates1}{Probability of each miRNA to be active in condition 1 (only for condition specific sampling).}
  \item{miRNAstates2}{Probability of each miRNA to be active in condition 2 (only for condition specific sampling).}
  \item{miRNAActivitySwitch}{Probability of each miRNA switching its activities (non-condition specific sampling).}
  \item{TFstates1}{Probability of each TF to be active in condition 1 (only for condition specific sampling).}
  \item{TFstates2}{Probability of each TF to be active in condition 2 (only for condition specific sampling).}
  \item{TFactivitySwitch}{Probability of each TF switching its activities (non-condition specific sampling).}
  \item{log_lik_trace}{Log-likelihood trace of MCMC sampling.}
  \item{TFomega}{Weights of TF-target graph. (effect of a TF on its targets)}
  \item{miRNAomega}{Weights of miRNA-target graph. (effect of a miRNA on its targets)}
  \item{genesetsTF}{TF-target network. This might be different from the network submitted to the function, due to incompatibilities of network and experimental measurements. Check your warnings and command line output!}
  \item{genesetsmiRNA}{miRNA-target network. This might be different from the network submitted to the function, due to incompatibilities of network and experimental measurements. Check your warnings and command line output!}
  \item{mRNAexpr}{mRNA expression data. This might be different from the matrix submitted to the function, due to incompatibilities of network and experimental measurements. Check your warnings and command line output!}
  \item{miRNAexpr}{miRNA expression data. This might be different from the matrix submitted to the function, due to incompatibilities of network and experimental measurements. Check your warnings and command line output!}    
  \item{TFexpr}{TF expression data (only if they are specifically included). This might be different from the matrix submitted to the function, due to incompatibilities of network and experimental measurements. Check your warnings and command line output!}    
}
\author{
Holger Froehlich, Benedikt Zacher
}
\examples{
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
}
\keyword{htest}
