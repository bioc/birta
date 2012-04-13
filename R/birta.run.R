# Main user interface function
# 
# important assumption: TF-names all start with 'V$' (as in TRANSFAC)
# Author: frohlich, zacher
###############################################################################


# dat.mRNA: mRNA expression data (genes x samples) - (reference, treatment)
# nrep: (replicates reference miRNA, replicates treatment miRNA, replicates reference mRNA, replicates treatment mRNA)
# dat.miRNA: miRNA expression data (miRNA x samples) - (reference, treatment)
# limmamRNA: output of limma analysis for mRNA data (list: pvalue.tab, lm.fit)
# limmamiRNA: output of limma analysis for miRNA data (list: pvalue.tab, lm.fit)
# fdr.mRNA: FDR cutoff for significance of the logFC for mRNA data
# fdr.miRNA: FDR cutoff for significance of the logFC for miRNA data
# lfc.mRNA: additional logFC cutoff for significance in mRNA data
# lfc.miRNA: additional logFC cutoff for significance in miRNA data
# genesets: combined TF / miRNA network
# TFs: indices of TFs in the combined network
# lambda: regularization parameter for edge weights
# sample.weights: Should edge weights be adapted during sampling?
# theta_TF: Expected fraction of active TFs
# there_miRNA: Expected fraction of active miRNAs
# model: type of model
# niter: number of MCMC iterations (AFTER burnin)
# nburnin: number of MCMC iterations UNTIL burnin is assumed to be finished
# potential_swaps: pre-computed potential swaps (OPTIONAL, see get_potential_swaps)
# run.pretest: initialize miRNA and TF states via the result of a hypergeometric test in order to improve convergence (should be taken with care; advise: only use it in case of observed convergence problems!)
# condition.specific.inference: Should inference on TF / miRNA activities be made only RELATIVE to a reference condition or independently in both conditions?
# OUTPUT: see birta.start


birtaRun = function(dat.mRNA, dat.miRNA, TFexpr, limmamRNA=NULL, limmamiRNA=NULL, limmaTF=NULL, nrep=NULL, fdr.mRNA=0.05, fdr.miRNA=0.05, lfc.mRNA=0, lfc.miRNA=0, genesets=NULL, lambda=NULL, sample.weights=TRUE, one.regulator.weight=TRUE, theta_TF=NULL, theta_miRNA=NULL, model=c("all-plug-in", "no-plug-in"), niter=500000, nburnin=100000, thin=50, potential_swaps=NULL, run.pretest=FALSE, condition.specific.inference=TRUE, only_switches=FALSE, weightSampleMean=0, weightSampleVariance=0.001, accessible=NULL) {


	if(is.null(nrep)) {
		cat("Automatic assignement of #replicates with design matrix.\n")
		if(! is.null(dat.miRNA)) {
			group1 = as.vector(which(limmamiRNA$design[,1] == 1))
			group2 = as.vector(which(limmamiRNA$design[,2] == 1))
			dat.miRNA = dat.miRNA[,c(group1, group2)]
			nrep[1] = length(group1)
			nrep[2] = length(group2)
		}
		

		if(is.null(limmamRNA)) {
			stop("limmamRNA must not be NULL!")
		}
		else {
			group1 = as.vector(which(limmamRNA$design[,1] == 1))
			group2 = as.vector(which(limmamRNA$design[,2] == 1))
			dat.mRNA = dat.mRNA[,c(group1, group2)]
			nrep[3] = length(group1)
			nrep[4] = length(group2)
		}
	}
	
	if(is.null(dat.miRNA) & (!is.null(limmamiRNA)) | (!is.null(dat.miRNA)) & (is.null(limmamiRNA))) {
		stop("Do you want to infer miRNA activities? Ambivalent definition of dat.miRNA and limmamiRNA (one is a NULL).")
	}

	if(class(genesets) != "list")
		stop("genesets must be a list!")
	
	
	mRNAinvalid = apply(dat.mRNA, 1, function(x) {(all(is.na(x[1:nrep[3]]))  | all(is.na(x[(nrep[3]+1):(nrep[3]+nrep[4])])))})
	if(length(which(mRNAinvalid)) > 0) {
		warning("Removing ", length(which(mRNAinvalid)), " mRNAs from expression matrix (all measurements of at least one condition equal NA).")
		dat.mRNA = dat.mRNA[(! mRNAinvalid),]
	}

	if(!is.null(dat.miRNA)) {
		miRNAinvalid = apply(dat.miRNA, 1, function(x) {(all(is.na(x[1:nrep[1]]))  | all(is.na(x[(nrep[1]+1):(nrep[1]+nrep[2])])))})
		if(length(which(miRNAinvalid)) > 0) {
			warning("Removing ", length(which(miRNAinvalid)), " miRNAs from expression matrix (all measurements of at least one condition equal NA).")
			dat.miRNA = dat.miRNA[(! miRNAinvalid),]
		}
	}

	model = match.arg(model, several.ok=FALSE)				
	
	genesetsTFs = NULL
	if("TF" %in% names(genesets)) {
		genesetsTFs = genesets$TF
	}
	genesetsmiRNA = NULL
	if("miRNA" %in% names(genesets)){
		genesetsmiRNA = genesets$miRNA	
	}
	

	cat("Formatting regulator-target network -> checking overlap between network and measurements.\n")
	if(length(genesetsmiRNA) > 0 ) {
		genesetsmiRNA = sapply(genesetsmiRNA, function(s) intersect(s, rownames(dat.mRNA)))
		if(length(genesetsmiRNA) == 0)
			warning("No overlap between miRNA-targets and mRNA measurements!")
		
		if(! is.null(dat.miRNA)) {
			genesetsmiRNA = genesetsmiRNA[sapply(genesetsmiRNA, function(s) (length(s) > 0))]
			genesetsmiRNA = genesetsmiRNA[intersect(names(genesetsmiRNA), rownames(dat.miRNA))]
			if(length(genesetsmiRNA) == 0)
				warning("No overlap between miRNAs in miRNA-target graph and miRNA measurements!")
		}

		if(any(sapply(genesetsmiRNA, length) == 0)) {
			warning("Not all miRNA genesets have non-zero length --> removing miRNAs empty genesets (see genesetsmiRNA in result).")
			genesetsmiRNA = genesetsmiRNA[sapply(genesetsmiRNA, length) > 0]			
		}	
	}

	if(length(genesetsTFs) > 0) {
		genesetsTFs = sapply(genesetsTFs, function(s) intersect(s, rownames(dat.mRNA)))
		if(length(genesetsTFs) == 0)
			warning("No overlap between TF-targets and mRNA measurements!")

		genesetsTFs = genesetsTFs[sapply(genesetsTFs, function(s) (length(s) > 0))]
		#if(any(sapply(genesetsTFs, length) == 0)){ # Das kann nie passieren auf Grund des vorherigen Kommandos
		#	warning("Not all TF genesets have non-zero length --> removing TFs empty genesets (see genesetsTFs in result).")
		#	genesetsTFs = genesetsTFs[sapply(genesetsTFs, length) > 0]
		#}
	}

	ctrl1 = setdiff(unlist(genesetsTFs), rownames(dat.mRNA))
	ctrl2 = setdiff(unlist(genesetsmiRNA), rownames(dat.mRNA))
	ctrl3 = c()
	if(! is.null(dat.miRNA)) {
		ctrl3 = setdiff(names(genesetsmiRNA), rownames(dat.miRNA))
	}
	if(length(ctrl1) > 0 | length(ctrl2) > 0 | length(ctrl3) > 0)
		stop("Problem: Incompatibilities between network and annnotation of expression data (not the same target genes or miRNAs)")
	

	diff_mRNA.down = limmamRNA$pvalue.tab$ID[limmamRNA$pvalue.tab$logFC < -lfc.mRNA & limmamRNA$pvalue.tab$adj.P.Val < fdr.mRNA]
	diff_mRNA.up = limmamRNA$pvalue.tab$ID[limmamRNA$pvalue.tab$logFC > lfc.mRNA & limmamRNA$pvalue.tab$adj.P.Val < fdr.mRNA]
	diff.genes = c(diff_mRNA.up, diff_mRNA.down)
	if(length(genesetsTFs) > 0)
		parentTFs = names(which(sapply(genesetsTFs, function(g) any(diff.genes %in% g))))
	else
		parentTFs = c()
	if(length(genesetsmiRNA) > 0)
		parentmiR = names(which(sapply(genesetsmiRNA, function(g) any(diff.genes %in% g))))
	else
		parentmiR = c()
	cat(length(diff.genes), " DE gene(s) have ", length(parentTFs), "regulating TFs and ", length(parentmiR), "regulating miRNAs\n") # NOTE: birta should also work without DE genes	
	# choose sensible initialization of states to improve convergence	
	if(run.pretest & length(diff.genes) > 0){
		init_miR = matrix(0, nrow=2, ncol=length(genesetsmiRNA))
		colnames(init_miR) = names(genesetsmiRNA)
		init_TF = matrix(0, nrow=2, ncol=length(genesetsTFs))
		colnames(init_TF) = names(genesetsTFs)
		pvals.miR = FisherPretest(limmamRNA$pvalue.tab, genesetsmiRNA, TRUE, fdr.miRNA, lfc.miRNA)
		pvals.TF = FisherPretest(limmamRNA$pvalue.tab, genesetsTFs, FALSE, fdr.mRNA, lfc.mRNA)				
		init_miR[2, names(pvals.miR)[pvals.miR < fdr.miRNA]] = 1
		init_TF[2, names(pvals.TF)[pvals.TF < fdr.mRNA]] = 1
	}	
	else{
		init_miR = init_TF = NULL		
	}		
	#
	if(!is.null(dat.miRNA) & !is.null(limmamiRNA)){	
		if(any(limmamRNA$lm.fit$contrasts != limmamiRNA$lm.fit$contrasts)) {
			stop("Limma contrasts for mRNA and miRNA expressions differ!\n")
		}
		limmamiRNA$pvalue.tab = limmamiRNA$pvalue.tab[limmamiRNA$pvalue.tab$ID %in% rownames(dat.miRNA),]
		diff_miRNA.up = limmamiRNA$pvalue.tab$ID[(limmamiRNA$pvalue.tab$logFC > lfc.miRNA) & limmamiRNA$pvalue.tab$adj.P.Val < fdr.miRNA]
		diff_miRNA.down = limmamiRNA$pvalue.tab$ID[(limmamiRNA$pvalue.tab$logFC < -lfc.miRNA) & limmamiRNA$pvalue.tab$adj.P.Val < fdr.miRNA]
		alpha_i0 = apply(dat.miRNA[,1:nrep[1], drop=F], 1, mean, na.rm=TRUE) # average expression in first condition
		names(alpha_i0) = rownames(dat.miRNA)		
		if(length(diff_miRNA.up) > 0)
			alpha_i = rep(median(limmamiRNA$pvalue.tab$logFC[match(diff_miRNA.up, limmamiRNA$pvalue.tab$ID)], na.rm=TRUE), nrow(dat.miRNA)) # alpha_i is the expected logFC, if miRNA is activated relative to reference ==> should be POSITIVE in case we operate with differential gene expression only (only states for 1 condition to infer) 
		else
			alpha_i = rep(1, nrow(dat.miRNA))	
		names(alpha_i) = rownames(dat.miRNA)
		if(condition.specific.inference){  # in case we use both condictions separately, alpha_i can also be negative			
			which.down = limmamiRNA$pvalue.tab$ID[(limmamiRNA$pvalue.tab$logFC) < 0]
			if(length(which.down) > 0){			
				if(length(diff_miRNA.down) > 0)
					alpha_i[which.down] = rep(median(limmamiRNA$pvalue.tab$logFC[match(diff_miRNA.down, limmamiRNA$pvalue.tab$ID)], na.rm=TRUE), length(which.down))
				else
					alpha_i[which.down] = rep(-1, length(which.down))
			}
		}
		var.miRNA = limmamiRNA$lm.fit$s2.post		
		var.varmiRNA = var(var.miRNA)
		A_Sigma = sqrt(var.miRNA)		
		names(A_Sigma) = rownames(dat.miRNA)
		if(any(is.na(A_Sigma)))
			stop("Variance of miRNA expressions must not equal NA! Re-do limmaAnalysis!")		
	}
	else {
		alpha_i0 = alpha_i = var.miRNA = var.varmiRNA = A_Sigma = NULL
	}		
	var.mRNA = limmamRNA$lm.fit$s2.post[rownames(dat.mRNA)]
	O_Sigma = sqrt(var.mRNA)
	names(O_Sigma) = rownames(dat.mRNA)
	if(! is.null(A_Sigma)) {
		if(any(is.na(A_Sigma)))
				stop("Variance of mRNA expressions must not equal NA! Re-do limmaAnalysis!")
	}
	avg.var = mean(c(var.mRNA, var.miRNA), na.rm=TRUE)
	var.var = max(c(var(var.mRNA), var.varmiRNA), na.rm=TRUE)
	if(var.var == 0) {
		var.var = 0.001
	}
	beta = var.var / (avg.var+1e-6)
	alpha = 1 + avg.var / (beta)
					
	if(!condition.specific.inference){
		b.mRNA = rowMeans(as.matrix(dat.mRNA[,1:nrep[3]]), na.rm=TRUE)
		names(b.mRNA) = rownames(dat.mRNA)
	}
	else{		
		b.mRNA =  c(mean(apply(dat.mRNA[,1:nrep[3]], 2, median)), mean(apply(dat.mRNA[,(nrep[3] + 1):NCOL(dat.mRNA)], 2, median)))
	}	
	# initialize omegas for miRNAs and TFs with sensible values to speed up convergence	
	meandown = max(c(median(limmamRNA$pvalue.tab$logFC[limmamRNA$pvalue.tab$logFC < -lfc.mRNA & limmamRNA$pvalue.tab$adj.P.Val < fdr.mRNA], na.rm=TRUE), -1), na.rm=T)
	meanup = min(c(median(limmamRNA$pvalue.tab$logFC[limmamRNA$pvalue.tab$logFC > lfc.mRNA &  limmamRNA$pvalue.tab$adj.P.Val < fdr.mRNA], na.rm=TRUE), 1), na.rm=T)	
	
	
	omegamiRNA = list()	
	for(i in names(genesetsmiRNA)) { # miRNA
		# NOTE: effect of miRNAs should be always inhibitory on mRNA expression level ==> NEGATIVE edge weights (has nothing to do with reference condition)
		if(limmamiRNA$pvalue.tab$logFC[match(i, limmamiRNA$pvalue.tab$ID)] < 0) { # miRNA is downregulated => mRNAs go up 
			diff.targets = intersect(genesetsmiRNA[[i]], diff_mRNA.up)
			if(length(diff.targets) > 0) {
				omegamiRNA[[i]] = rep(-median(limmamRNA$pvalue.tab$logFC[match(diff.targets, limmamRNA$pvalue.tab$ID)], na.rm=TRUE), length(genesetsmiRNA[[i]])) 
			}
			else {
				omegamiRNA[[i]] = rep(meandown,  length(genesetsmiRNA[[i]]))				
			}
		}
		else {  # miRNA is upregulated => mRNAs go down 
			diff.targets = intersect(genesetsmiRNA[[i]], diff_mRNA.down)
			if(length(diff.targets) > 0) {
				omegamiRNA[[i]] = rep(median(limmamRNA$pvalue.tab$logFC[match(diff.targets, limmamRNA$pvalue.tab$ID)], na.rm=TRUE), length(genesetsmiRNA[[i]]))
			}
			else {
				omegamiRNA[[i]] = rep(meandown,  length(genesetsmiRNA[[i]]))				
			}
		}
		names(omegamiRNA[[i]]) = genesetsmiRNA[[i]]
	}
	names(omegamiRNA) = names(genesetsmiRNA)

	omegaTF = list()	
	for(i in names(genesetsTFs)) { # TF
		# NOTE: TFs should lead to an up-regulation of target genes ==> POSITIVE weights (has nothing to do with reference condition)
		if((!is.null(TFexpr) && (i %in% rownames(TFexpr))) && limmaTF$pvalue.tab$logFC[match(i, limmaTF$pvalue.tab$ID)] < 0){ # # in this case we can make a similar reasoning than for miRNAs: TF is down => mRNAs go down 
			diff.targets = intersect(genesetsTFs[[i]], diff_mRNA.down)
			if(length(diff.targets) > 0) {
				omegaTF[[i]] = rep(-median(limmamRNA$pvalue.tab$logFC[match(diff.targets, limmamRNA$pvalue.tab$ID)], na.rm=TRUE), length(genesetsTFs[[i]])) 
			}
			else {
				omegaTF[[i]] = rep(-meandown,  length(genesetsTFs[[i]]))				
			}				
		}
		else{
			diff.targets = intersect(genesetsTFs[[i]], diff_mRNA.up) 
			if(length(diff.targets) > 0)
				omegaTF[[i]] = rep(median(limmamRNA$pvalue.tab$logFC[match(diff.targets, limmamRNA$pvalue.tab$ID)], na.rm=TRUE), length(genesetsTFs[[i]]))  	
			else
				omegaTF[[i]] =  rep(meanup, length(genesetsTFs[[i]]))				
		}
		names(omegaTF[[i]]) = genesetsTFs[[i]]
			
		#diff.targets = intersect(genesetsTFs[[i]], diff_mRNA.down)
		#if(length(diff.targets) > 0)
	#		omegaTF[[i]][cdiff[names(omegaTF[[i]])] < 0] = rep(median(limmamRNA$pvalue.tab$logFC[match(diff.targets, limmamRNA$pvalue.tab$ID)], na.rm=TRUE), length(which(cdiff[names(omegaTF[[i]])] < 0)))
	#	else
	#		omegaTF[[i]][which(cdiff[names(omegaTF[[i]])] < 0)] = rep(meandown,  length(which(cdiff[names(omegaTF[[i]])] < 0)))
	#	
	#	names(omegaTF[[i]]) = genesetsTFs[[i]]
		#omegaTF[[i]][cdiff < 0] = -omegaTF[[i]][cdiff < 0]		
	}
	names(omegaTF) = names(genesetsTFs)
	
	
	alpha_i0TF = NULL
	alpha_iTF = NULL
	TF_Sigma = NULL
	alphaTF = NULL
	betaTF = NULL
	if(!is.null(TFexpr)) {					
		if(is.null(nrep)){
			group1 = as.vector(which(limmaTF$design[,1] == 1))
			group2 = as.vector(which(limmaTF$design[,2] == 1))
			TFexpr = TFexpr[,c(group1, group2)]
		}		
		TFinvalid = apply(TFexpr, 1, function(x) {(all(is.na(x[1:nrep[3]]))  | all(is.na(x[(nrep[3]+1):(nrep[3]+nrep[4])])))})
		if(length(which(TFinvalid)) > 0) {
			warning("Removing ", length(which(TFinvalid)), " TFs from TF expression matrix (all measurements of at least one condition equal NA).")
			TFexpr = TFexpr[(! TFinvalid),]
		}
		
		TFsInAnnotation = rownames(TFexpr)[which((rownames(TFexpr) %in% names(genesetsTFs)))]
		TFexpr = TFexpr[TFsInAnnotation,]		
		diff_TFs = limmaTF$pvalue.tab$ID[(((limmaTF$pvalue.tab$logFC > lfc.mRNA) & (limmaTF$pvalue.tab$adj.P.Val < fdr.mRNA)) | ((limmaTF$pvalue.tab$logFC < -lfc.mRNA) & (limmaTF$pvalue.tab$adj.P.Val < fdr.mRNA)))]
		diff_TF.down = limmaTF$pvalue.tab$ID[limmaTF$pvalue.tab$logFC < -lfc.mRNA & limmaTF$pvalue.tab$adj.P.Val < fdr.mRNA]
		diff_TF.up = limmaTF$pvalue.tab$ID[limmaTF$pvalue.tab$logFC > lfc.mRNA & limmaTF$pvalue.tab$adj.P.Val < fdr.mRNA]
		
		TFexpr = TFexpr[as.character(diff_TFs),]		
		if(length(TFsInAnnotation) > 0) {
			cat("Using ", length(diff_TFs), "DE TFs together with their expression level for condition specific model.\n")
			
			alpha_i0TF = apply(as.matrix(TFexpr[,1:nrep[1], drop=F]), 1, function(x) mean(x, na.rm=TRUE)) 
			names(alpha_i0TF) = rownames(TFexpr)		
			if(length(diff_TF.up) > 0)
				alpha_iTF = rep(median(limmaTF$pvalue.tab$logFC[match(diff_TF.up, limmaTF$pvalue.tab$ID)], na.rm=TRUE), nrow(TFexpr)) # alpha_iTF is the expected logFC, if TF is activated relative to reference ==> should be POSITIVE in case we operate with differential gene expression only (only states for 1 condition to infer) 
			else
				alpha_iTF = rep(1, nrow(TFexpr))			
			names(alpha_iTF) = rownames(TFexpr)
			if(condition.specific.inference){ # in case we use both condictions separately, alpha_iTF can also be negative 
				which.down = limmaTF$pvalue.tab$ID[(limmaTF$pvalue.tab$logFC) < 0]
				if(length(which.down) > 0){			
					if(length(diff_TF.down) > 0)
						alpha_iTF[which.down] = rep(median(limmaTF$pvalue.tab$logFC[match(diff_TF.down, limmaTF$pvalue.tab$ID)], na.rm=TRUE), length(which.down))
					else
						alpha_iTF[which.down] = rep(-1, length(which.down))
				}
			}			
			#medianTFlogFC = median(abs(limmaTF$pvalue.tab$logFC[match(diff_TFs, limmaTF$pvalue.tab$ID)]), na.rm=TRUE)
			#alpha_iTF = rep(medianTFlogFC, nrow(TFexpr)) # alpha_i is the expected logFC, if TF is activated 
			#alpha_i1TF = apply(as.matrix(TFexpr[,(nrep[3]+1):(nrep[4]+nrep[2]), drop=F]), 1, function(x) mean(x, na.rm=TRUE)) 
			#alpha_iTF[which(alpha_i1TF-alpha_i0TF < 0)] = -medianTFlogFC						
			#which.down = limmaTF$pvalue.tab$ID[(limmaTF$pvalue.tab$logFC) < 0]
			var.TF = limmaTF$lm.fit$s2.post[rownames(TFexpr)]		
			var.varTF = var(var.TF, na.rm=TRUE)
			TF_Sigma = sqrt(var.TF)		
			names(TF_Sigma) = rownames(TFexpr)
			if(any(is.na(TF_Sigma)))
				stop("Variance of TF expressions must not equal NA! Re-do limmaAnalysis!")
			avg.var.TF = mean(c(var.TF, var.TF), na.rm=TRUE)
			var.var.TF = max(c(var(var.TF), var.varTF), na.rm=TRUE)
			if(var.var.TF == 0) {
				var.var.TF = 0.001
			}
			betaTF = var.var.TF / (avg.var.TF+1e-6)	
			alphaTF = 1 + avg.var.TF / (betaTF)			
		}
		else {
			warning("Cannot use TF expression for the model. Annotation does not match!")
			TFexpr = NULL
		}
	}
							
	#
	if(is.null(theta_TF) & length(genesetsTFs) > 0)
		theta_TF =  0.5*length(parentTFs)/length(genesetsTFs)
	if(is.null(theta_miRNA) & length(genesetsmiRNA) > 0)
		theta_miRNA = 0.5*length(parentmiR)/length(genesetsmiRNA)
	#	
	if(sample.weights){
		if(is.null(lambda)){
			if(length(genesetsTFs) > 0){
				parentsTFs = sapply(unique(unlist(genesetsTFs)), function(g) names(genesetsTFs)[sapply(genesets, function(s) return(g %in% s))])
				sum1 = sum(sqrt(sapply(parentsTFs, length))) 
			}
			else {
				sum1 = 0
			}
			if(length(genesetsmiRNA) > 0){
				parentsmiR = sapply(unique(unlist(genesetsmiRNA)), function(g) names(genesetsmiRNA)[sapply(genesets, function(s) return(g %in% s))])
				sum2 = sum(sqrt(sapply(parentsmiR, length)))
			}
			else {
				sum2 = 0
			}
			mylambda = (0.01*sum(nrep)*(nrow(dat.mRNA) + NROW(dat.miRNA))/(sum1 + sum2))
			
		}
		else
			mylambda = lambda
		weight_samples_per_move = 10		
	}
	else{
		weight_samples_per_move = 0
		mylambda = 0
	}
	
	res = birtaStart(mRNAexpr=dat.mRNA, miRNAexpr=dat.miRNA, 
			genesetsTF=genesetsTFs, genesetsmiRNA=genesetsmiRNA, 
			alpha_i = alpha_i, alpha_i0 = alpha_i0,
			alpha = alpha, beta = beta, b_j = b.mRNA,
			replicates=nrep, niter=niter, burnin=nburnin, thin=thin, model=model, 
			only_switches=only_switches, nomiRNA=is.null(genesetsmiRNA), noTF=is.null(genesetsTFs), omega_miRNA=omegamiRNA, omega_TF=omegaTF, potential_swaps=potential_swaps, 
			theta_TF=theta_TF, theta_miRNA=theta_miRNA, weightSampleMean=weightSampleMean, weightSampleVariance=weightSampleVariance, weight_samples_per_move=weight_samples_per_move, equal.regulator.weights=one.regulator.weight,
			A_sigma=A_Sigma, O_sigma=O_Sigma, lambda_omega=mylambda, init_S = init_miR, init_T = init_TF, condition.specific.inference=condition.specific.inference, TFexpr=TFexpr, alpha_i0TF=alpha_i0TF, alpha_iTF=alpha_iTF, TF_sigma=TF_Sigma, alphaTF=alphaTF, betaTF=betaTF, accessible=accessible)	
	
	res
}




setGeneric("birta", signature = c("dat.mRNA", "dat.miRNA", "TFexpr"), function(dat.mRNA, dat.miRNA, TFexpr, limmamRNA=NULL, limmamiRNA=NULL, limmaTF=NULL, nrep=NULL, fdr.mRNA=0.05, fdr.miRNA=0.05, lfc.mRNA=0, lfc.miRNA=0, genesets=NULL, lambda=NULL, sample.weights=TRUE, one.regulator.weight=TRUE, theta_TF=NULL, theta_miRNA=NULL, model=c("all-plug-in", "no-plug-in"), niter=500000, nburnin=100000, thin=50, potential_swaps=NULL, run.pretest=FALSE, condition.specific.inference=TRUE, only_switches=FALSE, weightSampleMean=0, weightSampleVariance=0.001, accessible=NULL)
	 standardGeneric("birta"))


setMethod("birta",  c("dat.mRNA"="ExpressionSet", "dat.miRNA" = "ExpressionSet", "TFexpr"="ExpressionSet"),
    function(dat.mRNA, dat.miRNA, TFexpr, limmamRNA=NULL, limmamiRNA=NULL, limmaTF=NULL, nrep=NULL, fdr.mRNA=0.05, fdr.miRNA=0.05, lfc.mRNA=0, lfc.miRNA=0, genesets=NULL, lambda=NULL, sample.weights=TRUE, one.regulator.weight=TRUE, theta_TF=NULL, theta_miRNA=NULL, model=c("all-plug-in", "no-plug-in"), niter=500000, nburnin=100000, thin=50, potential_swaps=NULL, run.pretest=FALSE, condition.specific.inference=TRUE, only_switches=FALSE, weightSampleMean=0, weightSampleVariance=0.001, accessible=NULL) {
    mRNA <- exprs(dat.mRNA)
    miRNA <- exprs(dat.miRNA)
    TF <- exprs(TFexpr)
    callGeneric(dat.mRNA=mRNA, dat.miRNA=miRNA, TFexpr=TF, limmamRNA=limmamRNA, limmamiRNA=limmamiRNA, limmaTF=limmaTF, nrep=nrep, fdr.mRNA=fdr.mRNA, fdr.miRNA=fdr.miRNA, lfc.mRNA=lfc.mRNA, lfc.miRNA=lfc.miRNA, genesets=genesets, lambda=lambda, sample.weights=sample.weights, one.regulator.weight=one.regulator.weight, theta_TF=theta_TF, theta_miRNA=theta_miRNA, model=model, niter=niter, nburnin=nburnin, thin=thin, potential_swaps=potential_swaps, run.pretest=run.pretest, condition.specific.inference=condition.specific.inference, only_switches=only_switches, weightSampleMean=weightSampleMean, weightSampleVariance=weightSampleVariance, accessible=accessible)
})


setMethod("birta",  c("dat.mRNA"="ExpressionSet", "dat.miRNA" = "ExpressionSet", "TFexpr"="missing"),
    function(dat.mRNA, dat.miRNA, TFexpr, limmamRNA=NULL, limmamiRNA=NULL, limmaTF=NULL, nrep=NULL, fdr.mRNA=0.05, fdr.miRNA=0.05, lfc.mRNA=0, lfc.miRNA=0, genesets=NULL, lambda=NULL, sample.weights=TRUE, one.regulator.weight=TRUE, theta_TF=NULL, theta_miRNA=NULL, model=c("all-plug-in", "no-plug-in"), niter=500000, nburnin=100000, thin=50, potential_swaps=NULL, run.pretest=FALSE, condition.specific.inference=TRUE, only_switches=FALSE, weightSampleMean=0, weightSampleVariance=0.001, accessible=NULL) {
    mRNA <- exprs(dat.mRNA)
    miRNA <- exprs(dat.miRNA)
    callGeneric(dat.mRNA=mRNA, dat.miRNA=miRNA, limmamRNA=limmamRNA, limmamiRNA=limmamiRNA, limmaTF=limmaTF, nrep=nrep, fdr.mRNA=fdr.mRNA, fdr.miRNA=fdr.miRNA, lfc.mRNA=lfc.mRNA, lfc.miRNA=lfc.miRNA, genesets=genesets, lambda=lambda, sample.weights=sample.weights, one.regulator.weight=one.regulator.weight, theta_TF=theta_TF, theta_miRNA=theta_miRNA, model=model, niter=niter, nburnin=nburnin, thin=thin, potential_swaps=potential_swaps, run.pretest=run.pretest, condition.specific.inference=condition.specific.inference, only_switches=only_switches, weightSampleMean=weightSampleMean, weightSampleVariance=weightSampleVariance, accessible=accessible)
})

setMethod("birta",  c("dat.mRNA"="ExpressionSet", "dat.miRNA" = "missing", "TFexpr"="ExpressionSet"),
    function(dat.mRNA, dat.miRNA, TFexpr, limmamRNA=NULL, limmamiRNA=NULL, limmaTF=NULL, nrep=NULL, fdr.mRNA=0.05, fdr.miRNA=0.05, lfc.mRNA=0, lfc.miRNA=0, genesets=NULL, lambda=NULL, sample.weights=TRUE, one.regulator.weight=TRUE, theta_TF=NULL, theta_miRNA=NULL, model=c("all-plug-in", "no-plug-in"), niter=500000, nburnin=100000, thin=50, potential_swaps=NULL, run.pretest=FALSE, condition.specific.inference=TRUE, only_switches=FALSE, weightSampleMean=0, weightSampleVariance=0.001, accessible=NULL) {
    mRNA <- exprs(dat.mRNA)
    TF <- exprs(TFexpr)
    callGeneric(dat.mRNA=mRNA, TFexpr=TF, limmamRNA=limmamRNA, limmamiRNA=limmamiRNA, limmaTF=limmaTF, nrep=nrep, fdr.mRNA=fdr.mRNA, fdr.miRNA=fdr.miRNA, lfc.mRNA=lfc.mRNA, lfc.miRNA=lfc.miRNA, genesets=genesets, lambda=lambda, sample.weights=sample.weights, one.regulator.weight=one.regulator.weight, theta_TF=theta_TF, theta_miRNA=theta_miRNA, model=model, niter=niter, nburnin=nburnin, thin=thin, potential_swaps=potential_swaps, run.pretest=run.pretest, condition.specific.inference=condition.specific.inference, only_switches=only_switches, weightSampleMean=weightSampleMean, weightSampleVariance=weightSampleVariance, accessible=accessible)
})


setMethod("birta",  c("dat.mRNA"="ExpressionSet", "dat.miRNA" = "missing", "TFexpr"="missing"),
    function(dat.mRNA, dat.miRNA, TFexpr, limmamRNA=NULL, limmamiRNA=NULL, limmaTF=NULL, nrep=NULL, fdr.mRNA=0.05, fdr.miRNA=0.05, lfc.mRNA=0, lfc.miRNA=0, genesets=NULL, lambda=NULL, sample.weights=TRUE, one.regulator.weight=TRUE, theta_TF=NULL, theta_miRNA=NULL, model=c("all-plug-in", "no-plug-in"), niter=500000, nburnin=100000, thin=50, potential_swaps=NULL, run.pretest=FALSE, condition.specific.inference=TRUE, only_switches=FALSE, weightSampleMean=0, weightSampleVariance=0.001, accessible=NULL) {
    mRNA <- exprs(dat.mRNA)
    callGeneric(dat.mRNA=mRNA, limmamRNA=limmamRNA, limmamiRNA=limmamiRNA, limmaTF=limmaTF, nrep=nrep, fdr.mRNA=fdr.mRNA, fdr.miRNA=fdr.miRNA, lfc.mRNA=lfc.mRNA, lfc.miRNA=lfc.miRNA, genesets=genesets, lambda=lambda, sample.weights=sample.weights, one.regulator.weight=one.regulator.weight, theta_TF=theta_TF, theta_miRNA=theta_miRNA, model=model, niter=niter, nburnin=nburnin, thin=thin, potential_swaps=potential_swaps, run.pretest=run.pretest, condition.specific.inference=condition.specific.inference, only_switches=only_switches, weightSampleMean=weightSampleMean, weightSampleVariance=weightSampleVariance, accessible=accessible)
})

 
  setMethod("birta",  c("dat.mRNA"="matrix", "dat.miRNA" = "matrix", "TFexpr"="matrix"),
      function(dat.mRNA, dat.miRNA, TFexpr, limmamRNA=NULL, limmamiRNA=NULL, limmaTF=NULL, nrep=NULL, fdr.mRNA=0.05, fdr.miRNA=0.05, lfc.mRNA=0, lfc.miRNA=0, genesets=NULL, lambda=NULL, sample.weights=TRUE, one.regulator.weight=TRUE, theta_TF=NULL, theta_miRNA=NULL, model=c("all-plug-in", "no-plug-in"), niter=500000, nburnin=100000, thin=50, potential_swaps=NULL, run.pretest=FALSE, condition.specific.inference=TRUE, only_switches=FALSE, weightSampleMean=0, weightSampleVariance=0.001, accessible=NULL) {
      birtaRun(dat.mRNA, dat.miRNA, TFexpr, limmamRNA, limmamiRNA, limmaTF=limmaTF, nrep, fdr.mRNA, fdr.miRNA, lfc.mRNA, lfc.miRNA, genesets, lambda, sample.weights, one.regulator.weight, theta_TF, theta_miRNA, model, niter, nburnin, thin, potential_swaps, run.pretest, condition.specific.inference, only_switches, weightSampleMean, weightSampleVariance, accessible=accessible)
  })
 
 
  setMethod("birta",  c("dat.mRNA"="matrix", "dat.miRNA" = "matrix", "TFexpr"="missing"),
      function(dat.mRNA, dat.miRNA, TFexpr, limmamRNA=NULL, limmamiRNA=NULL, limmaTF=NULL, nrep=NULL, fdr.mRNA=0.05, fdr.miRNA=0.05, lfc.mRNA=0, lfc.miRNA=0, genesets=NULL, lambda=NULL, sample.weights=TRUE, one.regulator.weight=TRUE, theta_TF=NULL, theta_miRNA=NULL, model=c("all-plug-in", "no-plug-in"), niter=500000, nburnin=100000, thin=50, potential_swaps=NULL, run.pretest=FALSE, condition.specific.inference=TRUE, only_switches=FALSE, weightSampleMean=0, weightSampleVariance=0.001, accessible=NULL) {
      birtaRun(dat.mRNA, dat.miRNA, TFexpr=NULL, limmamRNA, limmamiRNA, limmaTF=limmaTF, nrep, fdr.mRNA, fdr.miRNA, lfc.mRNA, lfc.miRNA, genesets, lambda, sample.weights, one.regulator.weight, theta_TF, theta_miRNA, model, niter, nburnin, thin, potential_swaps, run.pretest, condition.specific.inference, only_switches, weightSampleMean, weightSampleVariance, accessible=accessible)
  })
 
 setMethod("birta",  c("dat.mRNA"="matrix", "dat.miRNA" = "missing", "TFexpr"="matrix"),
     function(dat.mRNA, dat.miRNA, TFexpr, limmamRNA=NULL, limmamiRNA=NULL, limmaTF=NULL, nrep=NULL, fdr.mRNA=0.05, fdr.miRNA=0.05, lfc.mRNA=0, lfc.miRNA=0, genesets=NULL, lambda=NULL, sample.weights=TRUE, one.regulator.weight=TRUE, theta_TF=NULL, theta_miRNA=NULL, model=c("all-plug-in", "no-plug-in"), niter=500000, nburnin=100000, thin=50, potential_swaps=NULL, run.pretest=FALSE, condition.specific.inference=TRUE, only_switches=FALSE, weightSampleMean=0, weightSampleVariance=0.001, accessible=NULL) {
     birtaRun(dat.mRNA, dat.miRNA=NULL, TFexpr, limmamRNA, limmamiRNA, limmaTF=limmaTF, nrep, fdr.mRNA, fdr.miRNA, lfc.mRNA, lfc.miRNA, genesets, lambda, sample.weights, one.regulator.weight, theta_TF, theta_miRNA, model, niter, nburnin, thin, potential_swaps, run.pretest, condition.specific.inference, only_switches, weightSampleMean, weightSampleVariance, accessible=accessible)
 })
 
 
 setMethod("birta",  c("dat.mRNA"="matrix", "dat.miRNA" = "missing", "TFexpr"="missing"),
     function(dat.mRNA, dat.miRNA, TFexpr, limmamRNA=NULL, limmamiRNA=NULL, limmaTF=NULL, nrep=NULL, fdr.mRNA=0.05, fdr.miRNA=0.05, lfc.mRNA=0, lfc.miRNA=0, genesets=NULL, lambda=NULL, sample.weights=TRUE, one.regulator.weight=TRUE, theta_TF=NULL, theta_miRNA=NULL, model=c("all-plug-in", "no-plug-in"), niter=500000, nburnin=100000, thin=50, potential_swaps=NULL, run.pretest=FALSE, condition.specific.inference=TRUE, only_switches=FALSE, weightSampleMean=0, weightSampleVariance=0.001, accessible=NULL) {
     birtaRun(dat.mRNA, dat.miRNA=NULL, TFexpr=NULL, limmamRNA, limmamiRNA, limmaTF=limmaTF, nrep, fdr.mRNA, fdr.miRNA, lfc.mRNA, lfc.miRNA, genesets, lambda, sample.weights, one.regulator.weight, theta_TF, theta_miRNA, model, niter, nburnin, thin, potential_swaps, run.pretest, condition.specific.inference, only_switches, weightSampleMean, weightSampleVariance, accessible=accessible)
 })

