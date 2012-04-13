/**
 * @file    getStates.cpp
 * @author  Benedikt Zacher, AG Tresch, Gene Center Munich (zacher@lmb.uni-muenchen.de)
 * @version 0.99.0
 *
 * @section LICENSE
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 2 of
 * the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details at
 * http://www.gnu.org/copyleft/gpl.html
 *
 * @section DESCRIPTION
 *
 * This file implements the C wrapper function for interfacing C++ code from R.
 */


#include "BayesNetwork.h"
#include "BayesNetworkNC.h"

// TODO: Warnings/Errors, wenn Kanten nicht exisitieren
extern "C" {
	#include "getStates.h"


	R_CallMethodDef callMethods[] = {
		{"getStates", (DL_FUNC)&getStates, 48},
		{NULL, NULL, 0}
	};

	void R_init_myLib(DllInfo *info) {
		R_registerRoutines(info, NULL, callMethods, NULL, NULL);
	}


	SEXP getStates(SEXP num_mRNA, SEXP mRNA_names, SEXP num_miRNA, SEXP miRNA_names, SEXP num_TF, SEXP TF_names, SEXP replicates, SEXP mRNA_expr, SEXP miRNA_expr, SEXP sexp_mRNADataType, SEXP sexp_miRNADataType, SEXP sexp_use_miRNA_expression, SEXP mirTargets, SEXP TFtargets, SEXP sexpn0, SEXP sexpalpha, SEXP sexpbeta, SEXP sexpalpha_i0, SEXP sexpalpha_i, SEXP sexpb_j, SEXP sexpomega_miRNA, SEXP sexpomega_TF, SEXP niter, SEXP miRNA_sigma, SEXP mRNA_sigma, SEXP sexpmodel, SEXP sexpburnin, SEXP sexpthin, SEXP sexponly_switches, SEXP sexpT_potential_swaps, SEXP sexpS_potential_swaps, SEXP weightSampleMean, SEXP weightSampleVariance, SEXP sexpweight_sample_per_move, SEXP sexptheta_TF, SEXP sexptheta_miRNA, SEXP sexplambda_omega, SEXP sexpinit_S, SEXP sexpinit_T, SEXP sexpcondition_specific, SEXP sexpequal_regulator_weights, SEXP sexpTFexpr, SEXP sexpnTFexpr, SEXP sexpalpha_i0TF, SEXP sexpalpha_iTF, SEXP sexpTF_sigma, SEXP sexpalphaTF, SEXP sexpbetaTF, SEXP sexpaccessible) {
		int i,j,c,r;

	      	int O_cnt = INTEGER(num_mRNA)[0];
	      	int A_cnt = INTEGER(num_miRNA)[0];
	      	int T_cnt = INTEGER(num_TF)[0];
		int use_miRNA_expression = INTEGER(sexp_use_miRNA_expression)[0];
		int mRNADataType = INTEGER(sexp_mRNADataType)[0];
		int miRNADataType = INTEGER(sexp_miRNADataType)[0];
		int thin = INTEGER(sexpthin)[0];
		int condition_specific = INTEGER(sexpcondition_specific)[0];
		int equal_regulator_weights = INTEGER(sexpequal_regulator_weights)[0];
	
		// hyperparameters
		double n0 = REAL(sexpn0)[0];     
		double alpha = REAL(sexpalpha)[0];
		double beta = REAL(sexpbeta)[0];
		double alphaTF = REAL(sexpalphaTF)[0];
		double betaTF = REAL(sexpbetaTF)[0];
		
		double theta_TF = REAL(sexptheta_TF)[0];
		double theta_miRNA = REAL(sexptheta_miRNA)[0];
		double lambda_omega = REAL(sexplambda_omega)[0];
	
		double *alpha_i0 = NULL;
		double *alpha_i = NULL;
		int** init_S = NULL;
		if(A_cnt > 0){
			init_S = (int**)R_alloc(2, sizeof(int*));
			for(c = 0; c < 2; c++){
				init_S[c] = (int*) R_alloc(A_cnt, sizeof(int));
				for(i = 0; i < A_cnt; i++){
					init_S[c][i] = (int)INTEGER(sexpinit_S)[c + i*2];
					//Rprintf("S[%i][%i] = %i ", c, i, init_S[c][i]);
				}
			}			
		}
		
		
		
		if(use_miRNA_expression){
			alpha_i0 = (double *)malloc(sizeof(double)*A_cnt);
			alpha_i = (double *) malloc(sizeof(double)*A_cnt);
			for(i=0; i<A_cnt; i++) {
				alpha_i0[i] = REAL(sexpalpha_i0)[i];
				alpha_i[i] = REAL(sexpalpha_i)[i];
			}
		}
		double* b_j = (double*) R_alloc(O_cnt, sizeof(double));
		for(i = 0; i < O_cnt; i++){
			b_j[i] = REAL(sexpb_j)[i];
		}

		int only_switches = INTEGER(sexponly_switches)[0];
		list<int> *S_potential_swaps = NULL;
		if((A_cnt > 0) & (!only_switches)) {
			S_potential_swaps = new list<int>[A_cnt];
			for(i=0; i<A_cnt; i++) {
				int curr_size = LENGTH(VECTOR_ELT(sexpS_potential_swaps, i));
				for(j=0; j<curr_size; j++) {
					// Adding IDs of potenital swap partners, indices change in C
					S_potential_swaps[i].push_back(INTEGER(VECTOR_ELT(sexpS_potential_swaps, i))[j]-1);
				}
		
			}
		}

		list<int> *T_potential_swaps = NULL;
		if((T_cnt > 0) & (!only_switches)) {
			T_potential_swaps = new list<int>[T_cnt];
			for(i=0; i<T_cnt; i++) {
				int curr_size = LENGTH(VECTOR_ELT(sexpT_potential_swaps, i));
				for(j=0; j<curr_size; j++) {
					// Adding IDs of potenital swap partners, indices change in C
					T_potential_swaps[i].push_back(INTEGER(VECTOR_ELT(sexpT_potential_swaps, i))[j]-1);
				}
			}
		}
		double **omega_miRNA = (double **)calloc(A_cnt, sizeof(double*));
		for(i=0; i<A_cnt; i++) {
			int curr_size = LENGTH(VECTOR_ELT(sexpomega_miRNA, i));
			omega_miRNA[i] = (double*)calloc(curr_size, sizeof(double));
			for(j=0; j<curr_size; j++) {
				omega_miRNA[i][j] = REAL(VECTOR_ELT(sexpomega_miRNA, i))[j];
			}
		}
	
		 double **omega_TF = (double **)calloc(T_cnt, sizeof(double*));
		for(i=0; i<T_cnt; i++) {
			int curr_size = LENGTH(VECTOR_ELT(sexpomega_TF, i));
			omega_TF[i] = (double*)calloc(curr_size, sizeof(double));
			for(j=0; j<curr_size; j++) {
				omega_TF[i][j] = REAL(VECTOR_ELT(sexpomega_TF, i))[j];
			}
		}
		// IDs of mRNAs
		char **MymRNAs = (char **)malloc(sizeof(char*)*O_cnt);
		for(i=0; i<O_cnt; i++) {
			int l = LENGTH(STRING_ELT(mRNA_names, i));
			const char *tempRNA = CHAR(STRING_ELT(mRNA_names, i));
			MymRNAs[i] = (char *)malloc(sizeof(char)*(l+1));
			for(j=0; j<=l; j++) {
		  		MymRNAs[i][j] = tempRNA[j];
			}
		}

		// IDs of miRNAs
		char **MymiRNAs = (char **)malloc(sizeof(char*)*A_cnt);
		for(i=0; i<A_cnt; i++) {
			int l = LENGTH(STRING_ELT(miRNA_names, i));
			const char *tempRNA = CHAR(STRING_ELT(miRNA_names, i));
			MymiRNAs[i] = (char *)malloc(sizeof(char)*(l+1));
			for(j=0; j<=l; j++) {
		  		MymiRNAs[i][j] = tempRNA[j];
			}
		}


		// IDs of transcription factors
		int** init_T = NULL;
		if(T_cnt > 0){
			init_T = (int**) R_alloc(2, sizeof(int*));
			for(c = 0; c < 2; c++){
				init_T[c] = (int*) R_alloc(T_cnt, sizeof(int));
				for(i = 0; i < T_cnt; i++){
					init_T[c][i] = (int)INTEGER(sexpinit_T)[c + i*2];
				}
			}
		}
		char **MyTFs = (char **)malloc(sizeof(char*)*T_cnt);
		for(i=0; i<T_cnt; i++) {
			int l = LENGTH(STRING_ELT(TF_names, i));
			const char *tempTF = CHAR(STRING_ELT(TF_names, i));
			MyTFs[i] = (char *)malloc(sizeof(char)*(l+1));
			for(j=0; j<=l; j++) {
		  		MyTFs[i][j] = tempTF[j];
			}
		}
	      // #replicates of experiment (miRNA=0, mRNA=1) e under condition c (control=0, treated=1), access rep_cnt[e][c]
	      int **rep_cnt = (int**)malloc(sizeof(int *)*2); 
	      rep_cnt[0] = (int *)malloc(sizeof(int)*2);
	      rep_cnt[1] = (int *)malloc(sizeof(int)*2);
	      rep_cnt[0][0] = INTEGER(replicates)[0];
	      rep_cnt[0][1] = INTEGER(replicates)[1];
	      rep_cnt[1][0] = INTEGER(replicates)[2];
	      rep_cnt[1][1] = INTEGER(replicates)[3];


	double ***A = NULL;
	if(use_miRNA_expression){
		A = (double ***) malloc(sizeof(double**)*2);
		for(c=0; c<2; c++) {
			A[c] = (double**) malloc(sizeof(double*)*A_cnt);
			for(i=0; i<A_cnt; i++) {
				int nreps;
				if(c==0) {
					nreps = rep_cnt[0][0];
				}
				else if(c == 1) {
					nreps = rep_cnt[0][1];
				}
				A[c][i] = (double*) malloc(sizeof(double)*nreps);
				for(r=0; r<nreps; r++) {
					if(c==0) {
						A[c][i][r] = REAL(miRNA_expr)[i+(r*A_cnt)];
					}
					else if(c==1) {
						A[c][i][r] = REAL(miRNA_expr)[i+((r+rep_cnt[0][0])*A_cnt)];
					}
				}
			}
		}
	}


    
      	// Expression under condition c = {0=control,1=treated} of mRNA/miRNA i, in replicate r 
      	double ***O = (double ***) malloc(sizeof(double**)*2);     
	for(c=0; c<2; c++) {
		O[c] = (double**) malloc(sizeof(double*)*O_cnt);
		for(i=0; i<O_cnt; i++) {
			int nreps;
			if(c==0) {
				nreps = rep_cnt[1][0];
			}
			else if(c == 1) {
				nreps = rep_cnt[1][1];
			}
			O[c][i] = (double*) malloc(sizeof(double)*nreps);
			for(r=0; r<nreps; r++) {
				if(c==0) {
					O[c][i][r] = REAL(mRNA_expr)[i+(r*O_cnt)];
				}
				else if(c==1) {
					O[c][i][r] = REAL(mRNA_expr)[i+((r+rep_cnt[1][0])*O_cnt)];
				}
			}
		}
	}

	int nTFexpr = INTEGER(sexpnTFexpr)[0];
	double ***Otf = NULL;
	double *alpha_i0TF = NULL;
	double *alpha_iTF = NULL;
	double *TF_sigma = NULL;
	//Rprintf("%d\n", nTFexpr);
	if(nTFexpr > 0) {
		// Expression for transcription factors
		Otf = (double ***) malloc(sizeof(double**)*2);     
		for(c=0; c<2; c++) {
			Otf[c] = (double**) malloc(sizeof(double*)*nTFexpr);
			for(i=0; i<nTFexpr; i++) {
				int nreps;
				if(c==0) {
					nreps = rep_cnt[1][0];
				}
				else if(c == 1) {
					nreps = rep_cnt[1][1];
				}
				Otf[c][i] = (double*) malloc(sizeof(double)*nreps);
				for(r=0; r<nreps; r++) {
					if(c==0) {
						Otf[c][i][r] = REAL(sexpTFexpr)[i+(r*nTFexpr)];
					}
					else if(c==1) {
						Otf[c][i][r] = REAL(sexpTFexpr)[i+((r+rep_cnt[1][0])*nTFexpr)];
					}
					//Rprintf("%f ", Otf[c][i][r]);
				}
				//Rprintf("\n");
			}
		}
		
		alpha_i0TF = (double *)malloc(sizeof(double)*nTFexpr);
		alpha_iTF = (double *) malloc(sizeof(double)*nTFexpr);
		TF_sigma = (double *) malloc(sizeof(double)*nTFexpr);
		for(i=0; i<nTFexpr; i++) {
			alpha_i0TF[i] = REAL(sexpalpha_i0TF)[i];
			alpha_iTF[i] = REAL(sexpalpha_iTF)[i];
			TF_sigma[i] = REAL(sexpTF_sigma)[i];
			//Rprintf("%f %f %f\n", alpha_i0TF[i], alpha_iTF[i], TF_sigma[i]);
		}			
		
	}
      

	list<int> *S2O = NULL;
	list<int> *SparentsOfO = NULL;
	if(A_cnt > 0) {
		// Intialization of edges between miRNAs and genes
		// Note that indices in R start at 1. Thus index-1 here
		S2O = new list<int>[A_cnt];
		SparentsOfO = new list<int>[O_cnt];
		for(i=0; i<A_cnt; i++) {
			for(j=0; j<LENGTH(VECTOR_ELT(mirTargets, i)); j++) {
				int currMirTarget = INTEGER(VECTOR_ELT(mirTargets, i))[j]-1;
				SparentsOfO[currMirTarget].push_back(i);
				S2O[i].push_back(currMirTarget);
			}
		}
	}
	
	list<int> *T2O = NULL;
	list<int> *TparentsOfO = NULL;
	if(T_cnt > 0) {
		// Intialization of edges between TFs and genes
		T2O = new list<int>[T_cnt];
		TparentsOfO = new list<int>[O_cnt];
		for(i=0; i<T_cnt; i++) {
			for(j=0; j<LENGTH(VECTOR_ELT(TFtargets, i)); j++) {
				int currTFtarget = INTEGER(VECTOR_ELT(TFtargets, i))[j]-1;
				T2O[i].push_back(currTFtarget);
				TparentsOfO[currTFtarget].push_back(i);
			}
		}
	}

	double *O_sigma = NULL;
	if(mRNA_sigma != NULL) {
		O_sigma = (double *) malloc(sizeof(double)*O_cnt);

		for(i=0; i<O_cnt; i++) {
			O_sigma[i] = REAL(mRNA_sigma)[i];
		}
	}

	double *A_sigma = NULL;
	if(use_miRNA_expression){
		A_sigma = (double *) malloc(sizeof(double)*A_cnt);
		for(i=0; i<A_cnt; i++) {
			A_sigma[i] = REAL(miRNA_sigma)[i];
		}
	}

	double **O_mu = (double **)malloc(sizeof(double*)*2);
	int** methylated = (int**) malloc(sizeof(int*) * 2);
	for(c=0; c<2; c++) {
		methylated[c] = (int*) malloc(sizeof(int) * O_cnt);
		O_mu[c] = (double *)malloc(sizeof(double)*O_cnt);
		for(j=0; j<O_cnt; j++) {
			methylated[c][i] = (int)INTEGER(sexpaccessible)[c + j*2];
			if(condition_specific)
				O_mu[c][j] = b_j[c]; // initially all TFs and miRNAs are inactive ==> mRNAs have the same expectation dependent mean (same for each mRNA!)
			else
				O_mu[c][j] = b_j[j]; // initially all TFs and miRNAs are inactive ==> mRNAs have the same mean under both conditions (different for each mRNA!)
		}
	}

	int model = INTEGER(sexpmodel)[0];
	int niterations = INTEGER(niter)[0];
	int burnin = INTEGER(sexpburnin)[0];
	double sampleMean = REAL(weightSampleMean)[0];
	double sampleVariance = REAL(weightSampleVariance)[0];
	int weight_samples_per_move = INTEGER(sexpweight_sample_per_move)[0];

	BayesNetwork *bn;
	if(condition_specific){
		bn = new BayesNetwork(O_cnt, A_cnt, T_cnt, MymRNAs, MymiRNAs, MyTFs, rep_cnt, O, A, mRNADataType, miRNADataType, S2O, SparentsOfO, T2O, TparentsOfO, n0, alpha, beta, alpha_i0, alpha_i, omega_miRNA, omega_TF, A_sigma, O_sigma, model, O_mu, only_switches, S_potential_swaps, T_potential_swaps, sampleMean, sampleVariance, weight_samples_per_move, equal_regulator_weights, theta_TF, theta_miRNA, lambda_omega, init_S, init_T, Otf, nTFexpr, alpha_i0TF, alpha_iTF, TF_sigma, alphaTF, betaTF, methylated);
		
	}
	else{
		bn = new BayesNetworkNC(O_cnt, A_cnt, T_cnt, MymRNAs, MymiRNAs, MyTFs, rep_cnt, O, A, mRNADataType, miRNADataType, S2O, SparentsOfO, T2O, TparentsOfO, n0, alpha, beta, alpha_i0, alpha_i, omega_miRNA, omega_TF, A_sigma, O_sigma, model, O_mu, only_switches, S_potential_swaps, T_potential_swaps, sampleMean, sampleVariance, weight_samples_per_move, equal_regulator_weights, theta_TF, theta_miRNA, lambda_omega, init_S, init_T, Otf, nTFexpr, alpha_i0TF, alpha_iTF, TF_sigma, alphaTF, betaTF, methylated);
		  		
	}
	Rprintf("sampling ...\n");
	double *log_lik_trace = bn->MCMC(niterations, burnin, thin);      
	Rprintf("finished.\n");		  
	
// Create R output		

	SEXP mirAct1, mirAct2, tfAct1, tfAct2, result, wnames, log_lik, tfweights, miRweights;
	PROTECT(mirAct1 = NEW_NUMERIC(A_cnt));
	PROTECT(mirAct2 = NEW_NUMERIC(A_cnt));
	PROTECT(tfAct1 = NEW_NUMERIC(T_cnt));
	PROTECT(tfAct2 = NEW_NUMERIC(T_cnt));

	if(A_cnt > 0){
		for(i=0; i<bn->getA_cnt(); i++) {
		 	NUMERIC_POINTER(mirAct1)[i] = bn->getPostS()[0][i];
		  	NUMERIC_POINTER(mirAct2)[i] = bn->getPostS()[1][i];
		}
	}
	if(T_cnt > 0){
		for(i=0; i<T_cnt; i++) {
		  	NUMERIC_POINTER(tfAct1)[i] = bn->getPostT()[0][i];
		  	NUMERIC_POINTER(tfAct2)[i] = bn->getPostT()[1][i];
		}
	}
	PROTECT(log_lik = NEW_NUMERIC(niterations+burnin+1));
	for(i=0; i<niterations+burnin+1; i++) {		
		NUMERIC_POINTER(log_lik)[i] = log_lik_trace[i];		
	}
	// store sampled weights in a list
	//TFs
	SEXP currWeights;
	PROTECT(tfweights = NEW_LIST(T_cnt));
	if(T_cnt > 0) {
		for(i=0; i<T_cnt; i++) {
			int curr_size = LENGTH(VECTOR_ELT(sexpomega_TF, i));
			PROTECT(currWeights = NEW_NUMERIC(curr_size));
			for(j=0; j<curr_size; j++) {
				NUMERIC_POINTER(currWeights)[j] = bn->getOmegaTF()[i][j];
			}
			SET_ELEMENT(tfweights, i, currWeights);
			//UNPROTECT(1);
		}
	}
	PROTECT(miRweights = NEW_LIST(A_cnt));
	if(A_cnt > 0){
		//miRNAs:
		if(A_cnt > 0) {		
			for(i=0; i<A_cnt; i++) {
				int curr_size = LENGTH(VECTOR_ELT(sexpomega_miRNA, i));
				PROTECT(currWeights = NEW_NUMERIC(curr_size));
				for(j=0; j<curr_size; j++) {
					NUMERIC_POINTER(currWeights)[j] = bn->getOmegaMiRNA()[i][j];
				}
				SET_ELEMENT(miRweights, i, currWeights);
				//UNPROTECT(1);
			}
		}	
	}		
	PROTECT(result = NEW_LIST(7));
	PROTECT(wnames = NEW_CHARACTER(7));
	SET_STRING_ELT(wnames, 0, mkChar("miRNAstates1"));
	SET_STRING_ELT(wnames, 1, mkChar("miRNAstates2"));
	SET_STRING_ELT(wnames, 2, mkChar("TFstates1"));
	SET_STRING_ELT(wnames, 3, mkChar("TFstates2"));
	SET_STRING_ELT(wnames, 4, mkChar("log_lik_trace"));
	SET_STRING_ELT(wnames, 5, mkChar("TFweights"));
	SET_STRING_ELT(wnames, 6, mkChar("miRNAweights"));
	SET_NAMES(result, wnames);
	UNPROTECT(1);
	SET_ELEMENT(result, 0, mirAct1);
	SET_ELEMENT(result, 1, mirAct2);
	SET_ELEMENT(result, 2, tfAct1);
	SET_ELEMENT(result, 3, tfAct2);
	SET_ELEMENT(result, 4, log_lik);
	SET_ELEMENT(result, 5, tfweights);
	SET_ELEMENT(result, 6, miRweights);
	UNPROTECT(8+T_cnt+A_cnt);

	delete(bn);

      	return result;
  }

}
