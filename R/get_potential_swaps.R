get_potential_swaps = function(genesetsTF=NULL, genesetsmiRNA=NULL, perc.overlap.cutoff = 0.8, integer.id=TRUE, verbose=TRUE) {
	result = list()
	if(! is.null(genesetsTF) && length(genesetsTF) > 0) {
		if(verbose)
		cat("Calculating potential swaps for TFs.\n")
		tf = genesetsTF
		mtf = matrix(nrow=length(tf), ncol=length(tf))
		for(i in 1:length(tf)) {
			for(j in 1:length(tf)) {
				if(i != j) {
					mtf[i,j] = length(intersect(tf[[i]], tf[[j]]))/ length((tf[[i]]))#, tf[[j]]union
				}
				else {
					mtf[i,j] = 0
				}
			}
		}

		sel_tfswaps = mtf > perc.overlap.cutoff
		T_potential_swaps = sapply(1:length(tf), function(x) {names(tf)[sel_tfswaps[x,]]})
		names(T_potential_swaps) = names(tf)
		if(integer.id) {
			T_potential_swaps = lapply(T_potential_swaps, function(x) {which(names(tf) %in% x)})
		}
		result[["T_potential_swaps"]]=T_potential_swaps
	}

	if(! is.null(genesetsmiRNA) && length(genesetsmiRNA) > 0) {
		if(verbose)
			cat("Calculating potential swaps for miRNAs.\n")
		miR = genesetsmiRNA
		mmiR = matrix(nrow=length(miR), ncol=length(miR))
		for(i in 1:length(miR)) {
			for(j in 1:length(miR)) {
				if(i != j) {
					mmiR[i,j] = length(intersect(miR[[i]], miR[[j]]))/length((miR[[i]]))#, miR[[j]]))union
				}
				else {
					mmiR[i,j] = 0
				}
			}
		}


		sel_miRswaps = mmiR > perc.overlap.cutoff
		S_potential_swaps = sapply(1:length(miR), function(x) {names(miR)[sel_miRswaps[x,]]})
		names(S_potential_swaps) = names(miR)
		if(integer.id) {
			S_potential_swaps = lapply(S_potential_swaps, function(x) {which(names(miR) %in% x)})
		}
		result[["S_potential_swaps"]] = S_potential_swaps
	}
	result
}

