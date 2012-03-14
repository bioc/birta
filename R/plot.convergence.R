
plotConvergence = function(res, nburnin=NULL, title="") {
	if(! "log_lik_trace" %in% names(res))
		stop("Need entry log_lik_trace!")
	if(is.null(nburnin))
		stop("nburnin must not be a NULL!")

	take = which(res$log_lik_trace != Inf)
	if(length(take) > 0){
		plot(take, res$log_lik_trace[take], type="l", main=paste("log-likelihood along MCMC sampling", title), xlab="step", ylab="log likelihood")
		abline(v=nburnin, lty=3)		
	}
}