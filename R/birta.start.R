
birtaStart <- function(mRNAexpr, miRNAexpr=NULL, mRNA.data.type=c("array", "RNAseq"), miRNA.data.type=c("array", "RNAseq"), genesetsTF=NULL, genesetsmiRNA=NULL, replicates=c(5, 5, 5, 5), n0=1, alpha=1, beta=0.1, alpha_i=NULL, alpha_i0=NULL, b_j, niter=1e6, burnin=5*1e5, thin=50, model=c("all-plug-in", "no-plug-in"), only_switches=FALSE, noTF=FALSE, nomiRNA=FALSE, A_sigma=NULL, O_sigma, omega_miRNA=NULL, omega_TF=NULL, weightSampleMean=0, weightSampleVariance=1, equal.regulator.weights=TRUE, potential_swaps=NULL, weight_samples_per_move=10, theta_TF=0.01, theta_miRNA=0.01, lambda_omega=0, init_S=NULL, init_T=NULL, condition.specific.inference=TRUE, TFexpr=NULL, alpha_i0TF=NULL, alpha_iTF=NULL, TF_sigma=NULL, alphaTF=NULL, betaTF=NULL, accessible=NULL, perc.overlap.cutoff=0.8) 
{

	model = match.arg(model, several.ok=FALSE)
	stopifnot(model %in% c("all-plug-in", "no-plug-in"))
	model = switch(model, "no-plug-in"=3, "all-plug-in"=1)
	mRNA.data.type = match.arg(mRNA.data.type, several.ok=FALSE)
	miRNA.data.type = match.arg(miRNA.data.type, several.ok=FALSE)

	#if(! is.null(TFexpr)) {
	#	#print(length(genesetsTF))
	#	tfexpr = rownames(TFexpr)
	#	#print(tfexpr)
	#	tfnotexpr = names(genesetsTF)[(! (names(genesetsTF) %in% tfexpr))]
	#	#print(tfnotexpr)
	#	genesetsTF = genesetsTF[c(tfexpr, tfnotexpr)]
	#	omega_TF = omega_TF[c(tfexpr, tfnotexpr)]
		#print(names(genesetsTF))
	#}

	potential_swaps.orig = potential_swaps
	
	genesetsmiRNA = sapply(genesetsmiRNA, function(s) intersect(s, rownames(mRNAexpr)))
	genesetsTF = sapply(genesetsTF, function(s) intersect(s, rownames(mRNAexpr)))
	#genesetsmiRNA = genesetsmiRNA[unique(names(genesetsmiRNA))]
	#genesetsTF = genesetsTF[unique(names(genesetsTF))]

	omega_TF = sapply(omega_TF, function(s) s[intersect(names(s), rownames(mRNAexpr))])
	omega_miRNA = sapply(omega_miRNA, function(s) s[intersect(names(s), rownames(mRNAexpr))])	

	if(length(genesetsmiRNA) != length(omega_miRNA) | !all(sapply(genesetsmiRNA, length) == sapply(omega_miRNA, length)))
		stop("dimensions of genesetsmiRNA have to equal the dimensions of omega_miRNA")
	if(length(genesetsTF) != length(omega_TF)  | !all(sapply(genesetsTF, length) == sapply(omega_TF, length)))
		stop("dimensions of genesetsTF have to equal the dimensions of omega_TF")
	if(any(sapply(genesetsmiRNA, length) == 0) | any(sapply(genesetsTF, length) == 0)){
		warning("Not all genesets have non-zero length --> removing empty genesets")
		genesetsmiRNA = genesetsmiRNA[sapply(genesetsmiRNA, length) > 0]
		genesetsTF = genesetsTF[sapply(genesetsTF, length) > 0]
		omega_miRNA = omega_miRNA[sapply(omega_miRNA, length) > 0]
		omega_TF = omega_TF[sapply(omega_TF, length) > 0]
	}		
	if(!nomiRNA & NROW(miRNAexpr) > 0){
		common.miRNAs = intersect(names(genesetsmiRNA), rownames(miRNAexpr))
		genesetsmiRNA = genesetsmiRNA[common.miRNAs]
		omega_miRNA = omega_miRNA[common.miRNAs]
		miRNAexpr = miRNAexpr[common.miRNAs, ]
		A_sigma = A_sigma[common.miRNAs]
		alpha_i = alpha_i[common.miRNAs]
		alpha_i0 = alpha_i0[common.miRNAs]
		if(!is.null(init_S))
			init_S = init_S[, common.miRNAs]
	}		
	# only keep mRNAs that have regulators
	regulon = union(unlist(genesetsTF), unlist(genesetsmiRNA))
	mRNAexpr = mRNAexpr[rownames(mRNAexpr) %in% regulon, ]
	O_sigma = O_sigma[names(O_sigma) %in% regulon]
	if(!condition.specific.inference){
		b_j = b_j[names(b_j) %in% regulon]
	}
	#
	# recompute potential swaps, if necessary	
	if(!only_switches && (is.null(potential_swaps) || sapply(potential_swaps[["T_potential_swaps"]], length) != sapply(potential_swaps.orig[["T_potential_swaps"]], length) || sapply(potential_swaps[["S_potential_swaps"]], length) != sapply(potential_swaps.orig[["S_potential_swaps"]], length))){
		potential_swaps = get_potential_swaps(genesetsTF, genesetsmiRNA, perc.overlap.cutoff)			
	}	
	ctrl1 = setdiff(unlist(genesetsTF), rownames(mRNAexpr))
	ctrl2 = setdiff(unlist(genesetsmiRNA), rownames(mRNAexpr))
	ctrl3 = setdiff(names(genesetsmiRNA), rownames(miRNAexpr))
	if(length(ctrl1) > 0 || (!is.null(miRNAexpr) & length(ctrl2) > 0) || (!is.null(miRNAexpr) & length(ctrl3) > 0))
		stop("Problem: Incompatibilities between network and annnotation of expression data (not the same target genes or miRNAs)")	
	if(!is.null(A_sigma) & length(A_sigma) != NROW(miRNAexpr))
		stop("length of A_sigma has to equal number of miRNAs!")
	if(length(alpha_i) != NROW(miRNAexpr))
		stop("length of alpha_i has to equal number of miRNAs!")
	if(length(alpha_i0) != NROW(miRNAexpr))
		stop("length of alpha_i0 has to equal number of miRNAs!")
	if(!is.null(O_sigma) & length(O_sigma) != NROW(mRNAexpr))
		stop("length of O_sigma has to equal number of mRNAs!")
	if(length(b_j) != NROW(mRNAexpr) & !condition.specific.inference)
		stop("length of b_j has to equal number of mRNAs!")
	if(length(b_j) != 2 & condition.specific.inference)
		stop("Length of b_j has to equal 2!")
	if(length(replicates) != 4)
		stop("Length of replicates has equal 4! Put 0, if miRNA data is missing.")
	if(is.null(init_T)){
		init_T = matrix(0, nrow=2, length(genesetsTF))
		colnames(init_T) = names(genesetsTF)
	}
	if(is.null(init_S) & NROW(miRNAexpr) > 0){	
		init_S = matrix(0, nrow=2, ncol=NROW(miRNAexpr))
		colnames(init_S) = rownames(miRNAexpr)
	}	
	if(is.null(accessible)){
		accessible = matrix(1, nrow=2, NROW(mRNAexpr))
		colnames(accessible) = rownames(mRNAexpr)
	}	
	accessible = accessible[, rownames(mRNAexpr)]
	if(!is.null(accessible) && (NCOL(accessible) != NROW(mRNAexpr) | NROW(accessible) != 2))
		stop("accessible has to be of dimension 2 x length of mRNAexpr")	
	if(!is.null(init_T) && (NCOL(init_T) != length(genesetsTF) | NROW(init_T) != 2))
		stop("init_T has to be of dimension 2 x length of genesetsTF")
	if(!is.null(init_S) && (NCOL(init_S) != length(genesetsmiRNA) | NROW(init_S) != 2))
		stop("init_S has to be of dimension 2 x length of genesetsmiRNA")
	if(!is.null(miRNAexpr) && length(genesetsmiRNA) != NROW(miRNAexpr)){
		stop("length of genesetsmiRNA and NROW(miRNAexpr) have to be equal!")
	}
	if(!all(accessible %in% c(0,1))){
		stop("All entries of accessible have to be 0 or 1!")
	}
	
	lambda_omega = abs(lambda_omega)
	if(!is.null(theta_miRNA))
		theta_miRNA = abs(theta_miRNA)
	if(!is.null(theta_TF))
		theta_TF = abs(theta_TF)
	if(!is.null(alpha))
		alpha = abs(alpha)
	if(!is.null(beta))
		beta = abs(beta)
	if(!is.null(n0))
		n0 = abs(n0)
	if(!is.null(weight_samples_per_move))
		weight_samples_per_move = abs(weight_samples_per_move)
	if(!is.null(weightSampleVariance))
		weightSampleVariance = abs(weightSampleVariance)
	if(!is.null(burnin))
		burnin = abs(burnin)
	if(!is.null(thin))
		thin = abs(thin)
	if(!is.null(niter))
		niter = abs(niter)
	
	miRNA = names(genesetsmiRNA)
	mRNA = rownames(mRNAexpr)
	TF = names(genesetsTF)	
	A_cnt=as.integer(length(miRNA))
	T_cnt=as.integer(length(TF))		
	if(nomiRNA) {
		A_cnt = 0
		use_miRNA_expression = 0
		theta_miRNA = 0
	}
	else{
		use_miRNA_expression = (NROW(miRNAexpr) > 0)		
	}
	if(noTF) {
		T_cnt = 0
	}		
	init_T = init_T[, TF]
	## Edges from miRNAs to their targets
	mirTargets = genesetsmiRNA[miRNA]
	mirTargets = lapply(mirTargets, function(x) {which(mRNA %in% x)})
	## Edges from TFs to their targets
	TFtargets = genesetsTF[TF]
	TFtargets = lapply(TFtargets, function(x) {which(mRNA %in% x)})
	nTFexpr = 0
	if(! is.null(TFexpr)) {
		nTFexpr = as.integer(dim(TFexpr)[1])
	}			
	if(any(O_sigma == 0) | any(A_sigma == 0))
		warning("Variance estimates for mRNA or miRNA data contain 0s!")
	
	cat("\nBIRTA\n")
	cat("Data and network: #mRNAs = ", nrow(mRNAexpr), ", #miRNAs = ", A_cnt, ", #TFs = ", T_cnt, ", #TFs with expression data = ", nTFexpr, ", only one weight per regulator = ", equal.regulator.weights, "\n")		
	cat("Prior parameters: theta_TF = ", theta_TF, ", theta_miRNA = ", theta_miRNA, ", lambda = ", lambda_omega, "\n")
	if(model != 1)
		cat("Hyperparameters: alpha = ", alpha, ", beta = ", beta, ", n0 = ", n0, "\n")
	cat("MCMC parameters: burnin = ", burnin, ", niter = ", niter, ", thin = ", thin, ", condition specific inference = ", condition.specific.inference, "\n\n")			
	result = .Call("getStates", nmRNA=as.integer(length(mRNA)), mRNA=as.character(mRNA), nmiRNA=as.integer(A_cnt),  miRNA=as.character(miRNA),
				      nTF=as.integer(T_cnt), TF=as.character(TF), replicates=as.integer(replicates), mRNA_expression=as.numeric(mRNAexpr), 
				      miRNA_expression=as.numeric(miRNAexpr), mRNADataType=as.integer((mRNA.data.type=="RNAseq")*1), miRNADataType=as.integer((miRNA.data.type=="RNAseq")*1),
				      use_miRNA_expression=as.integer(use_miRNA_expression), genesetsmiRNA=mirTargets, genesetsTF=TFtargets, 
				      n0 = as.numeric(n0), alpha = as.numeric(alpha), beta = as.numeric(beta), 
				      alpha_i0 = as.numeric(alpha_i0), alpha_i = as.numeric(alpha_i), b_j = as.numeric(b_j),
				      omega_miRNA = omega_miRNA, omega_TF = omega_TF, niter=as.integer(niter), A_sigma=as.numeric(A_sigma), O_sigma=as.numeric(O_sigma),
					  
				      model=as.integer(model), burnin=as.integer(burnin), thin=as.integer(thin), only_switches=as.integer(only_switches), T_potential_swaps = sapply(potential_swaps$T_potential_swaps, as.integer), S_potential_swaps = sapply(potential_swaps$S_potential_swaps, as.integer), 
					  weightSampleMean=as.numeric(weightSampleMean), weightSampleVariance=as.numeric(weightSampleVariance), weight_samples_per_move=as.integer(weight_samples_per_move),					  
					  theta_TF=as.numeric(theta_TF), theta_miRNA=as.numeric(theta_miRNA), lambda_omega=as.numeric(lambda_omega), init_S=as.integer(init_S), init_T=as.integer(init_T), condition_specific=as.integer(condition.specific.inference), equal_regulator_weights=as.integer(equal.regulator.weights), 
					  TFexpr=as.numeric(TFexpr), nTFexpr=as.integer(nTFexpr), alpha_i0TF=as.numeric(alpha_i0TF), alpha_iTF=as.numeric(alpha_iTF), TF_sigma=as.numeric(TF_sigma), alphaTF=as.numeric(alphaTF), betaTF=as.numeric(betaTF), accessible=as.integer(accessible), PACKAGE="birta")			
 	 #cat("done.\n"
	if(! noTF) {
		names(result$TFstates1) = TF
		names(result$TFstates2) = TF
		if(! condition.specific.inference){
			if(!is.null(TFexpr)){
				result$TFActivitySwitch = result$TFstates2[(!(names(result$TFstates2) %in% rownames(TFexpr)))]
			}
			else{
				result$TFActivitySwitch = result$TFstates2
			}
			result = result[-which(names(result) %in% c("TFstates1","TFstates2"))]
		}
		if(!is.null(TFexpr)){
			result$TFstates1 = result$TFstates1[rownames(TFexpr)]
			result$TFstates2 = result$TFstates2[rownames(TFexpr)]
		}
	}
	if(! nomiRNA) {
		names(result$miRNAstates1) = miRNA
		names(result$miRNAstates2) = miRNA
		if(! condition.specific.inference) {
			result$miRNAactivitySwitch = result$miRNAstates2
			result$miRNAactivitySwitch[which(result$miRNAstates1 > 0)] = result$miRNAstates1[which(result$miRNAstates1 > 0)]
			result = result[-which(names(result) %in% c("miRNAstates1","miRNAstates2"))]
		}
	}
	result$genesetsTF = genesetsTF
	result$genesetsmiRNA = genesetsmiRNA
	result$mRNAexpr = mRNAexpr
	result$miRNAexpr = miRNAexpr
	return(result)
}
	
