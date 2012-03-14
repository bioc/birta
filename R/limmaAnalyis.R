# # limma analysis => needed for Fisher test (miRNA enrichment), miRNA differential expression analysis, and for our method
# limmaAnalysis = function(dat, contrasts=""){
# 	) # for our model: take coefficients for control as alpha_j0, coefficients for treatment minus coefficients for control as alpha_j, and sigma as nu_j
# }




setGeneric("limmaAnalysis", signature = c("dat", "design", "contrasts"), function(dat, design, contrasts) standardGeneric("limmaAnalysis"))


setMethod("limmaAnalysis",  c("dat"="ExpressionSet", "design"="matrix", "contrasts"="character"), function(dat, design, contrasts) {
    dat <- exprs(dat)
    callGeneric(dat, design, contrasts)
})

setMethod("limmaAnalysis",  c("dat"="matrix", "design"="matrix", "contrasts"="character"), function(dat, design, contrasts) {
    if(ncol(design) != 2) stop("Only two group comparison allowed!")
    fit = lmFit(dat, design)
    contrasts <- makeContrasts(contrasts=contrasts, levels=design)
    fit2 <- contrasts.fit(fit, contrasts)
    fit3 = eBayes(fit2)
    tab = topTable(fit3, number=Inf)	
    names(fit3$s2.post) = rownames(dat)
    list(pvalue.tab=tab, lm.fit=fit3, design=design, contrast=contrasts)
})

