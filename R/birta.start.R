
birtaStart <- function(mRNAexpr, miRNAexpr=NULL, mRNA.data.type=c("array", "RNAseq"), miRNA.data.type=c("array", "RNAseq"), genesetsTF=NULL, genesetsmiRNA=NULL, replicates=c(5, 5, 5, 5), n0=1, alpha=1, beta=0.1, alpha_i=NULL, alpha_i0=NULL, b_j, niter=1e6, burnin=5*1e5, thin=50, model=c("all-plug-in", "no-plug-in"), only_switches=FALSE, noTF=FALSE, nomiRNA=FALSE, A_sigma=NULL, O_sigma, omega_miRNA=NULL, omega_TF=NULL, weightSampleMean=0, weightSampleVariance=1, equal.regulator.weights=TRUE, potential_swaps=NULL, weight_samples_per_move=10, theta_TF=0.01, theta_miRNA=0.01, lambda_omega=0, init_S=NULL, init_T=NULL, condition.specific.inference=TRUE, TFexpr=NULL, alpha_i0TF=NULL, alpha_iTF=NULL, TF_sigma=NULL, alphaTF=NULL, betaTF=NULL) {

	model = match.arg(model, several.ok=FALSE)
	stopifnot(model %in% c("all-plug-in", "no-plug-in"))
	model = switch(model, "no-plug-in"=3, "all-plug-in"=1)
	mRNA.data.type = match.arg(mRNA.data.type, several.ok=FALSE)
	miRNA.data.type = match.arg(miRNA.data.type, several.ok=FALSE)

	if(! is.null(TFexpr)) {
		#print(length(genesetsTF))
		tfexpr = rownames(TFexpr)
		#print(tfexpr)
		tfnotexpr = names(genesetsTF)[(! (names(genesetsTF) %in% tfexpr))]
		#print(tfnotexpr)
		genesetsTF = genesetsTF[c(tfexpr, tfnotexpr)]
		omega_TF = omega_TF[c(tfexpr, tfnotexpr)]
		#print(names(genesetsTF))
	}

	potential_swaps.orig = potential_swaps
	
	genesetsmiRNA = sapply(genesetsmiRNA, function(s) intersect(s, rownames(mRNAexpr)))
	genesetsTF = sapply(genesetsTF, function(s) intersect(s, rownames(mRNAexpr)))
	genesetsmiRNA = genesetsmiRNA[unique(names(genesetsmiRNA))]
	genesetsTF = genesetsTF[unique(names(genesetsTF))]

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
	b_j = b_j[names(b_j) %in% regulon]

	
	
	
	# only keep mRNAs and miRNAs, where at least one valid measurement is present
# 	mRNAinvalid = apply(mRNAexpr, 1, function(x) all(is.na(x)))
# 	miRNAinvalid = apply(miRNAexpr, 1, function(x) all(is.na(x)))
# 	if(length(mRNAinvalid) > 0) {
# 		warning("Removing ", length(mRNAinvalid), " mRNAs from expression matrix (all measurements equal NA).")
# 		mRNAexpr = mRNAexpr[(! mRNAinvalid),]
# 	}
# 
# 	if(length(miRNAinvalid) > 0) {
# 		warning("Removing ", length(miRNAinvalid), " miRNAs from expression matrix (all measurements equal NA).")
# 		miRNAexpr = miRNAexpr[(! miRNAinvalid),]
# 	}
	# recompute potential swaps, if necessary
	if((!only_switches)) {
		if(is.null(potential_swaps) || sapply(potential_swaps[["T_potential_swaps"]], length) != sapply(potential_swaps.orig[["T_potential_swaps"]], length) || sapply(potential_swaps[["S_potential_swaps"]], length) != sapply(potential_swaps.orig[["S_potential_swaps"]], length)) {
			potential_swaps = get_potential_swaps(genesetsTF, genesetsmiRNA)
		}
	}
# 	ctrl1 = setdiff(unlist(genesetsTF), rownames(mRNAexpr))
# 	ctrl2 = setdiff(unlist(genesetsmiRNA), rownames(mRNAexpr))
# 	ctrl3 = setdiff(names(genesetsmiRNA), rownames(miRNAexpr))
# 	if(length(ctrl1) > 0 | length(ctrl2) > 0 | length(ctrl3) > 0)
# 		stop("Problem: Incompatibilities between network and annnotation of expression data (not the same target genes or miRNAs)")
	#if(!is.null(hypermethyl) & length(hypermethyl) != length(mRNA))
#		stop("length of hypermethyl has to equal number of mRNAs!")
	#if(!all(hypermethyl %in% c(0,1)))
#		stop("hypermethyl has to be a vector of 0 and 1!")
	if(!is.null(A_sigma) & length(A_sigma) != NROW(miRNAexpr))
		stop("length of A_sigma has to equal number of miRNAs!")
	if(length(alpha_i) != NROW(miRNAexpr))
		stop("length of alpha_i has to equal number of miRNAs!")
	if(length(alpha_i0) != NROW(miRNAexpr))
		stop("length of alpha_i0 has to equal number of miRNAs!")
	if(!is.null(O_sigma) & length(O_sigma) != NROW(mRNAexpr))
		stop("length of O_sigma has to equal number of mRNAs!")
	if(length(b_j) != NROW(mRNAexpr))
		stop("length of b_j has to equal number of mRNAs!")
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
	if(!is.null(init_T) & NCOL(init_T) != length(genesetsTF))
		stop("length of init_T has to equal length of genesetsTF")
	if(!is.null(init_S) & NCOL(init_S) != length(genesetsmiRNA))
		stop("length of init_S has to equal length of genesetsmiRNA")
	if(! is.null(miRNAexpr)) {
	if(length(genesetsmiRNA) != NROW(miRNAexpr))
		stop("length of genesetsmiRNA and NROW(miRNAexpr) have to be equal!")
	}
	lambda_omega = abs(lambda_omega)
	theta_miRNA = abs(theta_miRNA)
	theta_TF = abs(theta_TF)
	alpha = abs(alpha)
	beta = abs(beta)
	n0 = abs(n0)
	weight_samples_per_move = abs(weight_samples_per_move)
	weightSampleVariance = abs(weightSampleVariance)
	burnin = abs(burnin)
	thin = abs(thin)
	niter = abs(niter)
	
	miRNA = names(genesetsmiRNA)
	mRNA = rownames(mRNAexpr)
	TF = names(genesetsTF)	
	A_cnt=as.integer(length(miRNA))
	T_cnt=as.integer(length(TF))		
	if(nomiRNA) {
		A_cnt = 0
		use_miRNA_expression = 0
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
	cat("\nBIRTA\n")
	cat("Data and network: #mRNAs = ", nrow(mRNAexpr), "#miRNAs = ", A_cnt, "#TFs = ", T_cnt, "only one weight per regulator = ", equal.regulator.weights, "\n")		
	cat("Prior parameters: theta_TF = ", theta_TF, "theta_miRNA = ", theta_miRNA, "lambda = ", lambda_omega, "\n")
	if(model != "all-plug-in")
		cat("Hyperparameters: alpha = ", alpha, " beta = ", beta, " n0 = ", n0, "\n")
	cat("MCMC parameters: burnin = ", burnin, "niter = ", niter, "thin = ", thin, "condition specific inference = ", condition.specific.inference, "\n\n")	

	nTFexpr = 0
	if(! is.null(TFexpr)) {
		nTFexpr = as.integer(dim(TFexpr)[1])
	}

	result = .Call("getStates", nmRNA=as.integer(length(mRNA)), mRNA=as.character(mRNA), nmiRNA=as.integer(A_cnt),  miRNA=as.character(miRNA),
				      nTF=as.integer(T_cnt), TF=as.character(TF), replicates=as.integer(replicates), mRNA_expression=as.numeric(mRNAexpr), 
				      miRNA_expression=as.numeric(miRNAexpr), mRNADataType=as.integer((mRNA.data.type=="RNAseq")*1), miRNADataType=as.integer((miRNA.data.type=="RNAseq")*1),
				      use_miRNA_expression=as.integer(use_miRNA_expression), genesetsmiRNA=mirTargets, genesetsTF=TFtargets, 
				      n0 = as.numeric(n0), alpha = as.numeric(alpha), beta = as.numeric(beta), 
				      alpha_i0 = as.numeric(alpha_i0), alpha_i = as.numeric(alpha_i), b_j = as.numeric(b_j),
				      omega_miRNA = omega_miRNA, omega_TF = omega_TF, niter=as.integer(niter), A_sigma=as.numeric(A_sigma), O_sigma=as.numeric(O_sigma),
	
				      model=as.integer(model), burnin=as.integer(burnin), thin=as.integer(thin), only_switches=as.integer(only_switches), T_potential_swaps = sapply(potential_swaps$T_potential_swaps, as.integer), S_potential_swaps = sapply(potential_swaps$S_potential_swaps, as.integer), 
					  weightSampleMean=as.numeric(weightSampleMean), weightSampleVariance=as.numeric(weightSampleVariance), weight_samples_per_move=as.integer(weight_samples_per_move),					  
					  theta_TF=as.numeric(theta_TF), theta_miRNA=as.numeric(theta_miRNA), lambda_omega=as.numeric(lambda_omega), init_S=as.integer(init_S), init_T=as.integer(init_T), condition_specific=as.integer(condition.specific.inference), equal_regulator_weights=as.integer(equal.regulator.weights), TFexpr=as.numeric(TFexpr), nTFexpr=as.integer(nTFexpr), alpha_i0TF=as.numeric(alpha_i0TF), alpha_iTF=as.numeric(alpha_iTF), TF_sigma=as.numeric(TF_sigma), alphaTF=as.numeric(alphaTF), betaTF=as.numeric(betaTF), PACKAGE="birta")			
 	 #cat("done.\n"
	if((! noTF) & (! condition.specific.inference)) {
		names(result$TFstates1) = TF
		names(result$TFstates2) = TF
		result$TFActivitySwitch = result$TFstates2
		result = result[-which(names(result) %in% c("TFstates1","TFstates2"))]
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
    
	if((!is.null(TFexpr)) > 0 & condition.specific.inference) {
		names(result$TFstates1) = TF
		names(result$TFstates2) = TF
		result$TFActivitySwitch = result$TFstates2[(!(names(result$TFstates2) %in% rownames(TFexpr)))]
		result$TFstates1 = result$TFstates1[rownames(TFexpr)]
		result$TFstates2 = result$TFstates2[rownames(TFexpr)]
	}
	if((is.null(TFexpr)) > 0 & condition.specific.inference) {
		names(result$TFstates1) = TF
		names(result$TFstates2) = TF
		result$TFActivitySwitch = result$TFstates2
		result = result[-which(names(result) %in% c("TFstates1","TFstates2"))]
	}
	result$genesetsTF = genesetsTF
	result$genesetsmiRNA = genesetsmiRNA
	result$mRNAexpr = mRNAexpr
	result$miRNAexpr = miRNAexpr
	return(result)

}
