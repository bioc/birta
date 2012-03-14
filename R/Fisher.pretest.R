FisherPretest = function(limmaTab, genesets, downregulated=TRUE, fdr.gene=0.05, lfc.gene=0){	
	genesOfInterest = NULL
	if(downregulated) {
		genesOfInterest= limmaTab[limmaTab$adj.P.Val < fdr.gene & limmaTab$logFC < -lfc.gene,"ID"]		
	}
	else {
		genesOfInterest= limmaTab[limmaTab$adj.P.Val < fdr.gene & limmaTab$logFC > lfc.gene,"ID"]
		genesOfInterest= c(genesOfInterest, limmaTab[limmaTab$adj.P.Val < fdr.gene & limmaTab$logFC < -lfc.gene, "ID"])
	}

	if(length(genesOfInterest) == 0) {
		stop("Cannot run Fisher.pretest without differentially expressed genes.")
	}
	others = setdiff(limmaTab$ID, genesOfInterest)
	pvalues = unlist(sapply(names(genesets), function(p){
						freqsig = sum(genesOfInterest %in%  genesets[[p]])
						freqothers = sum(others %in%  genesets[[p]])
						conf.tab = matrix(c(freqsig, freqothers, length(genesOfInterest) - freqsig, length(others) - freqothers),nrow=2, dimnames=list(c("differential", "not differential"),c(p,"others")))
						conf.tab[is.na(conf.tab)] = 0
						fisher.test(conf.tab, alternative="g")$p.value
					}))
	pvalues = p.adjust(pvalues, method="BH")
	names(pvalues) = names(genesets)
	#print(pvalues[pvalues < 0.05])
	pvalues
}