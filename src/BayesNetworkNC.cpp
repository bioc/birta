/**
 * @file    BayesNetwork.cpp
 * @author  Benedikt Zacher, AG Tresch, Gene Center Munich (zacher@lmb.uni-muenchen.de)
 * @version 0.99.0
 *
 * @section LICENSE
 *
 * This program is CFree software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the CFree Software Foundation; either version 2 of
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
 * This file implements all functions and classes defined in the corresponding HEADER file.
 */


#include "BayesNetworkNC.h"

using namespace std;

		
BayesNetworkNC::BayesNetworkNC(int O_cnt, int A_cnt, int T_cnt, char **mRNAs, char **miRNAs, char **TFs, int **rep_cnt, double ***mRNA_expression, double ***miRNA_expression, int mRNADataType, int miRNADataType, list<int> *S2O, list<int> *SparentsOfO, list<int> *T2O, list<int> *TparentsOfO,  double n0, double alpha, double beta, double *alpha_i0, double* alpha_i, double **omega_miRNA, double **omega_TF, double *miRNA_sigma, double *mRNA_sigma, int model, double **O_mu, int only_switches, list<int> *S_potential_swaps, list<int> *T_potential_swaps, double weightSampleMean, double weightSampleVariance, int weight_samples_per_move, int equal_regulator_weights, double theta_TF, double theta_miRNA, double lambda_omega, int** init_S, int** init_T, double ***Otf, int nTFexpr, double *alpha_i0TF, double *alpha_iTF, double *TF_sigma, double alphaTF, double betaTF) { 

	this->equal_regulator_weights = equal_regulator_weights;
	
	this->MODEL = model;

	this->O_cnt = O_cnt;
	this->A_cnt = A_cnt;
	this->T_cnt = T_cnt;
	
	this->Otf = Otf;
	this->nTFexpr = nTFexpr;
	this->alpha_i0TF = alpha_i0TF;
	this->alpha_iTF = alpha_iTF;
	this->TF_sigma = TF_sigma;
	
	this->rep_cnt = rep_cnt;

	this->mRNAs = mRNAs;
	this->miRNAs = miRNAs;
	this->TFs = TFs;

	this->A = miRNA_expression;
	this->O = mRNA_expression;
	this->mRNADataType = mRNADataType;
	this->miRNADataType = miRNADataType;
	
	this->S2O = S2O;
	this->T2O = T2O;
	this->SparentsOfO = SparentsOfO;
	this->TparentsOfO = TparentsOfO;

	this->O_mu = O_mu;
	
	this->weightSampleMean = weightSampleMean;
	this->weightSampleVariance = weightSampleVariance;

	int i,j,c;			
	// hyperparameters
	this->n0 = n0;
	this->alpha = alpha;
	this->beta = beta;
	this->alphaTF = alphaTF;
	this->betaTF = betaTF;
	
	this->theta_TF = theta_TF;
	this->theta_miRNA = theta_miRNA;
	this->lambda_omega = lambda_omega;
	
	
	this->alpha_i0 = alpha_i0;
	this->alpha_i = alpha_i;	
	
	this->omega_miRNA = omega_miRNA;
	
	this->omega_TF = omega_TF;
				
	// Initialization of S
	this->posterior_miRNA = (double **)malloc(sizeof(double*)*2);
	this->S = (int**)malloc(sizeof(int*)*2);
	for(c=0; c<2; c++) {
	  this->S[c] = (int *) malloc(sizeof(int)*A_cnt);
	  this->posterior_miRNA[c] = (double *) malloc(sizeof(double)*A_cnt);
	  for(i=0; i<A_cnt; i++) {
	        S[c][i] = init_S[c][i];
		//if(S[c][i] == 1)
		//	printf("S[%i][%i] = %i\n", c, i, init_S[c][i]);
		this->posterior_miRNA[c][i] = 0;
	  }
	}

	// Initialization of T
	this->posterior_TF = (double **)malloc(sizeof(double*)*2);
	this->T = (int**)malloc(sizeof(int*)*2);
	for(c=0; c<2; c++) {
	  this->T[c] = (int *) malloc(sizeof(int)*T_cnt);
	  this->posterior_TF[c] = (double *) malloc(sizeof(double)*T_cnt);
	  for(i=0; i<T_cnt; i++) {
	        T[c][i] = init_T[c][i];
		//if(T[c][i] == 1)
		//	printf("T[%i][%i] = %i\n", c, i, init_T[c][i]);
		this->posterior_TF[c][i] = 0;
	  }
	}
	this->activeMiRNAs = (int*)malloc(sizeof(int)*2);
	this->activeMiRNAs[0] = 0;
	this->activeMiRNAs[1] = 0;

	this->activeTFs = (int*)malloc(sizeof(int)*2);
	this->activeTFs[0] = 0;
	this->activeTFs[1] = 0;

	//miR_higher_in_condition = NULL;
	this->miR_higher_in_condition = (int*)calloc(A_cnt, sizeof(int));
	int r;
	double mean_c1 = 0;
	double mean_c2 = 0;
	if(A != NULL){
		for(i=0; i<A_cnt; i++) {			
			for(r=0; r<rep_cnt[0][0]; r++) {
				mean_c1 = mean_c1 + A[0][i][r];
			}
			mean_c1 = mean_c1/rep_cnt[0][0];
			for(r=0; r<rep_cnt[0][1]; r++) {
				mean_c2 = mean_c2 + A[1][i][r];
			}
			mean_c2 = mean_c2/rep_cnt[0][1];

			if(mean_c1 < mean_c2) {
				this->miR_higher_in_condition[i] = 1;
			}
			else {
				this->miR_higher_in_condition[i] = -1;
			}
			//Rprintf("m1=%f, m2=%f, higher cond := %d\n", mean_c1, mean_c2, this->miR_higher_in_condition[i]);
		}
	}
	else{
		list<int>::iterator start, stop;
		for(i = 0; i < A_cnt; i++){
			start = S2O[i].begin();
			stop = S2O[i].end();
			mean_c1 = 0;
			mean_c2 = 0;
			for(list<int>::iterator it = start; it != stop; it++){
				j = *it;
				for(r=0; r<rep_cnt[1][0]; r++) {
					mean_c1 = mean_c1 + O[0][j][r];
				}				
				for(r=0; r<rep_cnt[1][1]; r++) {
					mean_c2 = mean_c2 + O[1][i][r];
				}
			}
			mean_c1 = mean_c1/(rep_cnt[1][0] * S2O[i].size());
			mean_c2 = mean_c2/(rep_cnt[1][1] * S2O[i].size()) ;
			if(mean_c1 < mean_c2) {
				this->miR_higher_in_condition[i] = -1; // assumption: miRNA is active in that condition, in which gene expression is lower
			}
			else {
				this->miR_higher_in_condition[i] = 1;
			}
		}
	}

	this->O_sigma = mRNA_sigma;
	this->A_sigma = miRNA_sigma;
	
	this->miRNA_plusminus = 0;
	this->miRNA_minusminus = A_cnt;
	this->miRNA_plusplus = 0;

	
	this->weight_samples_per_move = weight_samples_per_move;	

	this->T_potential_swaps = T_potential_swaps;
	// These lists are empty, because no TF is active
	T_possible_swaps = (list<int>**)malloc(sizeof(list<int>*)*2);
	for(c=0; c<2; c++) {
		T_possible_swaps[c] = new list<int>[T_cnt];
	}
	
	T_npossible_swaps = (int*)malloc(sizeof(int)*2);
	T_npossible_swaps[0] = 0;
	T_npossible_swaps[1] = 0;	

	this->S_potential_swaps = S_potential_swaps;
	S_possible_swaps = new list<int>[A_cnt];
	S_npossible_swaps = 0;
	
	this->only_switches = only_switches;
	

	if(this->weight_samples_per_move > 0){
		this->posterior_omega_TF = (double**) calloc(T_cnt, sizeof(double*));
		for(i = 0; i < T_cnt; i++)
			this->posterior_omega_TF[i] = (double*) calloc(T2O[i].size(), sizeof(double));
		this->posterior_omega_miRNA = (double**) calloc(A_cnt, sizeof(double*));
		for(i = 0; i < A_cnt; i++)
			this->posterior_omega_miRNA[i] = (double*) calloc(S2O[i].size(), sizeof(double));
	}
	else{
		this->posterior_omega_TF = omega_TF;
		this->posterior_omega_miRNA = omega_miRNA;
	}
	
}


BayesNetworkNC::~BayesNetworkNC() {	
}

int BayesNetworkNC::neighbourhood_switch() {
	return(this->A_cnt + this->T_cnt);			
}

void BayesNetworkNC::update_S_swaps(int switchid) {
	//Rprintf("S-size \n");//%d\n", S_potential_swaps[switchid].size());
	// Remove all entries of switchid in T[c][TF_k] with TF_k in T_possible_swaps[condition][switchid]
	if(S_possible_swaps[switchid].size() > 0) {			
		
		list<int>::iterator start = S_possible_swaps[switchid].begin();
		list<int>::iterator stop = S_possible_swaps[switchid].end(); 
		for(list<int>::iterator it = start; it != stop; it++) {
			int curr_tf = *it;
			int size_before = S_possible_swaps[curr_tf].size();
			S_possible_swaps[curr_tf].remove(switchid);
			int size_after = S_possible_swaps[curr_tf].size();
			S_npossible_swaps = S_npossible_swaps - (size_before - size_after);
			//Rprintf("#elements removed from %d: %d, sb=%d, sa=%d, swid=%d\n", curr_tf, size_before - size_after, size_before, size_after, switchid);
		}
		S_npossible_swaps = S_npossible_swaps - S_possible_swaps[switchid].size();
		S_possible_swaps[switchid].clear();		
			
	}

	
	// Determine all new possible swap partners j and add them to S_possible_swaps[switchid]. Also add switchid to T_possible_swaps[condition][j] for all possible partners j.
	list<int>::iterator start = S_potential_swaps[switchid].begin();
	list<int>::iterator stop = S_potential_swaps[switchid].end(); 
	for(list<int>::iterator it = start; it != stop; it++) {
		int j = *it;
		// add new entry j in S_possible_swaps[switchid]
		if(S[1][j] != S[1][switchid]) {
			S_possible_swaps[switchid].push_back(j);
			//Rprintf("%d added to possible swap partners of %d\n", switchid, j);
			S_possible_swaps[j].push_back(switchid);
			//Rprintf("%d added to possible swap partners of %d\n", j, switchid);
			S_npossible_swaps = S_npossible_swaps + 2;
		}
	}
	
	
}


double* BayesNetworkNC::MCMC(long niter, long burnin, int thin) {

// stores values of log likelihood, if moves are accepted, Inf otherwise
	double *log_lik_trace = (double*)malloc(sizeof(double)*(niter+burnin+1));
	long l;
	int j, k;			 
	list<int>::iterator start, stop;
	for(l=0; l<niter+burnin+1; l++) {
		log_lik_trace[l] = (double)INFINITY;
	}

    	
	// Initial calculation of log-likelihood
	int** nactive_miRNAs = NULL;
	if(A_cnt > 0) {
		nactive_miRNAs = (int**)calloc(O_cnt, sizeof(int*));
	}
	int** nactive_TFs = NULL;
	if(T_cnt > 0) {
		nactive_TFs = (int**)calloc(O_cnt, sizeof(int*));
	}
	
	for(k = 0; k < O_cnt; k++){
		if(A_cnt > 0 && nactive_miRNAs != NULL){
			nactive_miRNAs[k] = (int*) calloc(2, sizeof(int));
			nactive_miRNAs[k][0] = 0;
			nactive_miRNAs[k][1] = 0;
			start = SparentsOfO[k].begin();
			stop = SparentsOfO[k].end(); 		
			for(list<int>::iterator it = start; it != stop; it++) {
				j=*it; // miRNA ID
				nactive_miRNAs[k][0] += S[0][j];
				nactive_miRNAs[k][1] += S[1][j];
			}	
		}	
		if(T_cnt > 0 && nactive_TFs != NULL){
			nactive_TFs[k] = (int*) calloc(2, sizeof(int));
			nactive_TFs[k][0] = 0;
			nactive_TFs[k][1] = 0;
			start = TparentsOfO[k].begin();
			stop = TparentsOfO[k].end(); 		
			for(list<int>::iterator it = start; it != stop; it++) {
				j=*it; // TF ID
				nactive_TFs[k][0] += T[0][j];
				nactive_TFs[k][1] += T[1][j];
			}	
		}
		//printf("nactive_miRNAs[%i][0] = %i, nactive_miRNAs[%i][1] = %i\n", k, nactive_miRNAs[k][0], k, nactive_miRNAs[k][1]);		
		O_mu[0][k] = get_omuInitial(k, 0, nactive_miRNAs, nactive_TFs);
		O_mu[1][k] = get_omuInitial(k, 1, nactive_miRNAs, nactive_TFs);
	}
	
    	double log_lik = 0;
    	int c, r;
	k = 0;
    	long i;
    	for(c=0; c<2; c++) {
		// sum for miRNA replicates
		for(r=0; r<this->rep_cnt[0][c]; r++) {
			if(A != NULL){
		    		for(i=0; i<this->A_cnt; i++) {
					if(!(A[c][i][r] != A[c][i][r])) { // check for NaN
						if(MODEL == 1) {
							if(this->miRNADataType == ARRAY)
								log_lik -= (pow(A[c][i][r] - get_mu0(i, c), 2)/pow(A_sigma[i], 2));
							else
								log_lik -= logNB(A[c][i][r], get_mu0(i, c), A_sigma[i]);
						}
						else if((MODEL == 2) | (MODEL == 3)){
							if(this->miRNADataType != ARRAY){
								Rprintf("Model %i not implemented for RNAseq data!\n", MODEL);
								return(0);
							}							
							log_lik -= pow((A[c][i][r] + get_mu0(i, c)*n0)/(1+n0), 2) - (alpha+0.5) * log(0.5*n0*pow(get_mu0(i, c),2)+0.5*pow(A[c][i][r],2)+beta);
						}
					}
		    		}
			}
		}
		// sum for mRNA replicates
		for(r=0; r<this->rep_cnt[1][c]; r++) {
		    	for(j=0; j<this->O_cnt; j++) {
				if(!(O[c][j][r] != O[c][j][r])) { // check for NaN
					if(MODEL == 1) {
						if(this->mRNADataType == ARRAY)
							log_lik -= (pow(O[c][j][r] - O_mu[c][j], 2)/pow(O_sigma[j], 2));
						else
							log_lik -= logNB(O[c][j][r], O_mu[c][j], O_sigma[j]);
					}
					else if((MODEL == 2) | (MODEL == 3)){
						if(this->mRNADataType != ARRAY){
							Rprintf("Model %i not implemented for RNAseq data!\n", MODEL);
							return(0);
						}
						log_lik -= (0.5 + alpha)*log(1+1/(2*beta)*(pow(O[c][j][r] - O_mu[c][j], 2)));
					}
				}
		    	}
		}
		
		// sum over TF replicates
		if((nTFexpr > 0)) {
			
			for(r=0; r<this->rep_cnt[1][c]; r++) {
				for(i=0; i<this->nTFexpr; i++) {
					if(Otf[c][i][r] == Otf[c][i][r]) { // check for NaN
						if(MODEL == 1) {
							if(this->miRNADataType == ARRAY)
								log_lik -= (pow(Otf[c][i][r] - get_mu0TF(i, c), 2)/pow(TF_sigma[i], 2));
							else
								log_lik -= logNB(Otf[c][i][r], get_mu0TF(i, c), TF_sigma[i]);
						}
						else if((MODEL == 2) | (MODEL == 3)){
							if(this->miRNADataType != ARRAY){
								Rprintf("Model %i not implemented for RNAseq data!\n", MODEL);
								return(0);
							}							
							log_lik -= pow((Otf[c][i][r] + get_mu0TF(i, c)*n0)/(1+n0), 2) - (alphaTF+0.5) * log(0.5*n0*pow(get_mu0TF(i, c),2)+0.5*pow(Otf[c][i][r],2)+betaTF);
						}
						
					}
				}
			}
		}
    	}
   
	log_lik_trace[0] = log_lik;
	int weight_sampling = 1;
	double old_prior = 0;
	double new_prior= 0;
	if(this->weight_samples_per_move <= 0){
		Rprintf("No edge weight adjustment!\n");
		this->weight_samples_per_move = 1;
		weight_sampling = 0;
	}
	else{
		old_prior = PriorWeights(); // unscaled prior
	}
	double rnd;
	// set seed for random number generation (for MCMC-moves)
	srand(12345);
	GetRNGstate();
	// sum of TFs and miRNAs
	int TFmiRNA_nswitch = this->A_cnt+this->T_cnt; // number of possible switch operations (remains always the same)

    	for(i=0; i<niter+burnin; i++) {
		R_CheckUserInterrupt();
		
		int nall_possible_ops = TFmiRNA_nswitch + T_npossible_swaps[1] + S_npossible_swaps; // #all possible MCMC moves
		
		if(only_switches) {
			nall_possible_ops = TFmiRNA_nswitch;
		}
		// move is determined by sampling a random number out of all possible moves (switches and swaps)
		int make_move = rand() % nall_possible_ops;

		if(make_move >= TFmiRNA_nswitch) { // swap
			int swapid = make_move - TFmiRNA_nswitch;
			// determine swap move ID
			if(swapid >= T_npossible_swaps[1]) { // two miRNAs are swapped
				swapid = swapid - T_npossible_swaps[1];
				
				int swapid1 = 0; // ID of miRNA 1
				int curr_size = S_possible_swaps[swapid1].size();
				// swapid1 of miRNA 1 is determined
				while(swapid >= curr_size) {
					swapid = swapid - curr_size;
					swapid1++;
					curr_size = S_possible_swaps[swapid1].size();
				}
				// swapid2 of miRNA 2 is determined
				int swapid2 = 0; // ID for miRNA 2
				int counter = 0;
				list<int>::iterator start = S_possible_swaps[swapid1].begin();
				list<int>::iterator stop = S_possible_swaps[swapid1].end(); 
				for(list<int>::iterator it = start; it != stop; it++) {
					if(counter == swapid) {
						swapid2 = *it;
						break;
					}						
					counter++;
				}

				int weight_sample_it;
				for(weight_sample_it = 0; weight_sample_it<this->weight_samples_per_move; weight_sample_it++) {
					double *miR_weight_samples=NULL;
					int condition1 = 1;
					int condition2 = 1;
					
					int swapped2one = swapid1;
					if(S[condition2][swapid2] == 0) {
						swapped2one = swapid2;
					}
					       

					// Weights are sampled for miRNA which is swapped to active state. For the other one, weights are set to zero.
					if(weight_sampling == 1) {
						miR_weight_samples = (double*) calloc(S2O[swapped2one].size(), sizeof(double));
						rnd = rnorm(this->weightSampleMean, this->weightSampleVariance);
						double rndup = rnorm(this->weightSampleMean, this->weightSampleVariance);
						double rnddown = rnorm(this->weightSampleMean, this->weightSampleVariance);
						for(j = 0; j < (int)S2O[swapped2one].size(); j++){
							if(this->omega_miRNA[swapped2one][j] <= 0) 
								rnd = rnddown;
							if(this->omega_miRNA[swapped2one][j] >= 0)
								rnd = rndup;
							if(!equal_regulator_weights)								
								rnd = rnorm(this->weightSampleMean, this->weightSampleVariance);
							
							miR_weight_samples[j] = rnd;
						}				
					}

					// For calculation of new log-likelihood, two switches are performed
					double new_log_lik = log_lik + doSwitch(S, S2O, swapid1, condition1, 1, miR_weight_samples, nactive_miRNAs, nactive_TFs);
					S[condition1][swapid1] = -(S[condition1][swapid1]-1);
					if(swapid1 == swapped2one) {
						new_prior = updateWeightsAndOmu(S, omega_miRNA, S2O, swapid1, condition1, miR_weight_samples, 1, old_prior, nactive_miRNAs, nactive_TFs);
					}
					else { //only omu changes
						updateOmu(S, omega_miRNA, S2O, swapid1, condition1, NULL, 1, nactive_miRNAs, nactive_TFs);
					}
					// update log-likelihood and set state of swapid1 back to initial state
					new_log_lik = new_log_lik + doSwitch(S, S2O, swapid2, condition2, 1, miR_weight_samples, nactive_miRNAs, nactive_TFs);
					S[condition1][swapid1] = -(S[condition1][swapid1]-1);
					// set weights and omu back to initial state
					if(swapid1 == swapped2one) {
						new_prior = updateWeightsAndOmu(S, omega_miRNA, S2O, swapid1, condition1, miR_weight_samples, 1, old_prior, nactive_miRNAs, nactive_TFs);
					}
					else { //only omu changes
						updateOmu(S, omega_miRNA, S2O, swapid1, condition1, NULL, 1, nactive_miRNAs, nactive_TFs);
					}						
					double U = ((double)(rand() % 100000001)) /((double)100000000);
				
					while(U == 0) {
						U = ((double)(rand() % 100000001)) /((double)100000000);
					}
					if(U == 0 || log(U) <= new_log_lik-log_lik + (-new_prior + old_prior)*lambda_omega) { // Remark: Prior and neighbourhood do not change
						//Rprintf("Swapping miRNAs: %d <=> %d, %f < %f\n", swapid1, swapid2, log(U), new_log_lik-log_lik);
						S[condition1][swapid1] = -(S[condition1][swapid1]-1);
						S[condition2][swapid2] = -(S[condition2][swapid2]-1);
						log_lik = new_log_lik;
						// updating weights and omu
						//new_prior = updateWeightsAndOmu(S, omega_miRNA, S2O, swapid1, condition1, miR_weight_samples, 1, old_prior, nactive_miRNAs, nactive_TFs);
						new_prior = updateWeightsAndOmu(S, omega_miRNA, S2O, swapid2, condition2, miR_weight_samples, 1, old_prior, nactive_miRNAs, nactive_TFs);
						// updating possible swap moves
						update_S_swaps(swapid1);
						update_S_swaps(swapid2);
						old_prior = new_prior;
						
					}					
					if(weight_sampling == 1) {
						CFree(miR_weight_samples);
					}
					else					
						break;
				}
			} 
			else { // two TFs are swapped
				// condition is determined
				int condition = 1;
				/*if(swapid >= T_npossible_swaps[0]) {
					swapid = swapid - T_npossible_swaps[0];
					condition = 1;
				}*/			  
				
				// Determine IDs of TFs, which are to be swapped
				int swapid1 = 0; // ID of TF 1
				int curr_size = T_possible_swaps[condition][swapid1].size();
				while(swapid >= curr_size) {
				
					swapid = swapid - curr_size;
					swapid1++;
					curr_size = T_possible_swaps[condition][swapid1].size();
				}
				
				int swapid2 = 0;
				int counter = 0;
				list<int>::iterator start = T_possible_swaps[condition][swapid1].begin();
				list<int>::iterator stop = T_possible_swaps[condition][swapid1].end(); 
				for(list<int>::iterator it = start; it != stop; it++) {
					if(counter == swapid) {
						swapid2 = *it;
						break;
					}						
					counter++;
				}
				
				int weight_sample_it;
				for(weight_sample_it = 0; weight_sample_it<this->weight_samples_per_move; weight_sample_it++) {
					double *tf_weight_samples=NULL;
					int swapped2one;
					swapped2one = swapid1;
					if(T[condition][swapid2] == 0) {
					        swapped2one = swapid2;
					}

					if(weight_sampling == 1) {	
						tf_weight_samples = (double*) calloc(T2O[swapped2one].size(), sizeof(double));
					        rnd = rnorm(this->weightSampleMean, this->weightSampleVariance);
						double rndup = rnorm(this->weightSampleMean, this->weightSampleVariance);
						double rnddown = rnorm(this->weightSampleMean, this->weightSampleVariance);
						for(j = 0; j < (int)T2O[swapped2one].size(); j++){
							if(this->omega_TF[swapped2one][j] <= 0) 
								rnd = rnddown;
							if(this->omega_TF[swapped2one][j] >= 0)
								rnd = rndup;
							if(!equal_regulator_weights)								
								rnd = rnorm(this->weightSampleMean, this->weightSampleVariance);
							
							tf_weight_samples[j] = rnd;
						}				       					        
					}
					// For calculation of new log-likelihood, two switches are performed
					// Weights are sampled for TF which is swapped to active state. For the other one, weights are set to zero. 					
					double new_log_lik = log_lik + doSwitch(T, T2O, swapid1, condition, 0, tf_weight_samples, nactive_miRNAs, nactive_TFs);
					T[condition][swapid1] = -(T[condition][swapid1]-1);
					// weights are only affecting log-likelihood if the new state is one, i.e. if new samples are considered
					if(swapid1 == swapped2one) {
						new_prior = updateWeightsAndOmu(T, omega_TF, T2O, swapid1, condition, tf_weight_samples, 0, old_prior, nactive_miRNAs, nactive_TFs);
					}
					else { //only omu changes
						updateOmu(T, omega_TF, T2O, swapid1, condition, NULL, 0, nactive_miRNAs, nactive_TFs);
					}
					
					// update log-likelihood and set state of swapid1 back to initial state
					new_log_lik = new_log_lik + doSwitch(T, T2O, swapid2, condition, 0, tf_weight_samples, nactive_miRNAs, nactive_TFs);
					T[condition][swapid1] = -(T[condition][swapid1]-1);
					// set weights and omu back to initial state
					if(swapid1 == swapped2one) {
						new_prior = updateWeightsAndOmu(T, omega_TF, T2O, swapid1, condition, tf_weight_samples, 0, old_prior, nactive_miRNAs, nactive_TFs);
					}
					else { //only omu changes
						updateOmu(T, omega_TF, T2O, swapid1, condition, NULL, 0, nactive_miRNAs, nactive_TFs);
					}
					
					double U = ((double)(rand() % 100000001))/((double)100000000);
				
					while(U == 0) {
						U = ((double)(rand() % 100000001))/((double)100000000);
					}
					if(U == 0 || log(U) <= new_log_lik-log_lik + (-new_prior + old_prior)*lambda_omega) { // Remark: Prior and neighbourhood do not change
						//Rprintf("Swapping TFs: %d <=> %d, c=%d, %f < %f\n", swapid1, swapid2, condition, log(U), new_log_lik-log_lik);
						// updating possible swaps and changing states
						update_T_swaps(swapid1, T[condition][swapid1], condition);
						T[condition][swapid1] = -(T[condition][swapid1]-1);
						update_T_swaps(swapid2, T[condition][swapid2], condition);
						T[condition][swapid2] = -(T[condition][swapid2]-1);
						log_lik = new_log_lik;
						// update weights and omu
						new_prior = updateWeightsAndOmu(T, omega_TF, T2O, swapid1, condition, tf_weight_samples, 0, old_prior, nactive_miRNAs, nactive_TFs);
						new_prior = updateWeightsAndOmu(T, omega_TF, T2O, swapid2, condition, tf_weight_samples, 0, old_prior, nactive_miRNAs, nactive_TFs);
						old_prior = new_prior;						
					}
					if(weight_sampling == 1) {
						CFree(tf_weight_samples);
					}
					else
						break;
				}
			}					

			

			
			
		}
		else { // switch
			// Get random TF or miRNA for switch operation
			int switchid = make_move;
			
			if(switchid >= A_cnt) { // TF is switched
				switchid = switchid - A_cnt;
				int condition = 1;
				/*if(switchid >= T_cnt) {
					switchid = switchid - T_cnt;
					condition = 1;
				}*/				  
				// if weights are sampled, replace old weights with new weights here and calculate new log likelihood, try max. weight_samples_per_move different samples for the weights

				int weight_sample_it;
				for(weight_sample_it = 0; weight_sample_it<this->weight_samples_per_move; weight_sample_it++) {
					double *tf_weight_samples=NULL;
					if(weight_sampling == 1) {
					      tf_weight_samples = (double*) calloc(T2O[switchid].size(), sizeof(double));
					          rnd = rnorm(this->weightSampleMean, this->weightSampleVariance);
						double rndup = rnorm(this->weightSampleMean, this->weightSampleVariance);
						double rnddown = rnorm(this->weightSampleMean, this->weightSampleVariance);
						for(j = 0; j < (int)T2O[switchid].size(); j++){
							if(this->omega_TF[switchid][j] <= 0) 
								rnd = rnddown;
							if(this->omega_TF[switchid][j] >= 0)
								rnd = rndup;
							if(!equal_regulator_weights)								
								rnd = rnorm(this->weightSampleMean, this->weightSampleVariance);
							
							tf_weight_samples[j] = rnd;
						}
					}

					double new_log_lik = log_lik + doSwitch(T, T2O, switchid, condition, 0, tf_weight_samples, nactive_miRNAs, nactive_TFs);
					double U = ((double)(rand() % 100000001)) /((double)100000000);
				
					while(U == 0) {
						U = ((double)(rand() % 100000001)) /((double)100000000);
					}

					// if swaps are considered, neighbourhood changes upon switch
					int old_neighbourhood = T_cnt + A_cnt;
					int new_neighbourhood = T_cnt + A_cnt;
					if(! only_switches) {
						old_neighbourhood = old_neighbourhood + T_npossible_swaps[1] + S_npossible_swaps;
						new_neighbourhood = old_neighbourhood - 2*T_possible_swaps[condition][switchid].size(); // remove #swaps, that are not possible any more
						list<int>::iterator start = T_potential_swaps[switchid].begin();
						list<int>::iterator stop = T_potential_swaps[switchid].end(); 
						for(list<int>::iterator it = start; it != stop; it++) { // add new #possible swaps
							int j = *it;
							if(T[condition][j] == T[condition][switchid]) { // this is the OLD state
								new_neighbourhood = new_neighbourhood + 2;
							}
						}
					}
					new_prior = updatePrior(T2O, switchid, tf_weight_samples, 0, old_prior);
					
					// switch is kept if:
					if(U == 0 || log(U) <= new_log_lik-log_lik + Prior(T[condition][switchid], 0) + (-new_prior + old_prior)*lambda_omega + log(new_neighbourhood)-log(old_neighbourhood)) {
						//Rprintf("TF: log(U)=%f < %f\n", log(U), new_log_lik-log_lik+Prior(T[condition][switchid], 0)+log(new_neighbourhood)-log(old_neighbourhood));
						// perform switch and update #(in-)active TFs
						if(T[condition][switchid] == 0) {
							if(! only_switches) {
								update_T_swaps(switchid, 0, condition);
							}
							T[condition][switchid] = 1;
							activeTFs[condition]++;	
						}
						else {
							if(! only_switches) {
								update_T_swaps(switchid, 1, condition);
							}
							T[condition][switchid] = 0;
							activeTFs[condition]--;
						}
					
						log_lik = new_log_lik;
						log_lik_trace[i+1] = log_lik;	
						old_prior = new_prior;	

						int index_for_omega = 0;
						// update mu and weights after accepted switch
						// iterate over children j of t_kc
						list<int>::iterator start = T2O[switchid].begin();
						list<int>::iterator stop = T2O[switchid].end(); 
						for(list<int>::iterator it = start; it != stop; it++) {
							// j := ID of mRNA
							j=*it;
							// update of mu after accepted switch
							//double before = O_mu[condition][j];
							O_mu[condition][j] = get_omu(j, condition, index_for_omega, switchid, 0, tf_weight_samples, nactive_miRNAs, nactive_TFs, 1);
								
							if((T[condition][switchid] == 1) & (weight_sampling == 1)) {
								omega_TF[switchid][index_for_omega] = omega_TF[switchid][index_for_omega]+tf_weight_samples[index_for_omega];								
							}
							index_for_omega++;
						}

						/*double new_prior_test = PriorWeights();
						if(abs(new_prior - new_prior_test) > 1e-4)
							printf("Prior problem 1: priorWeights() = %g, my prior = %g!\n", new_prior_test, new_prior);*/
					}
					if(weight_sampling == 1) {
						CFree(tf_weight_samples);
					}	
					else					
						break;
				}
			}
			// miRNA is switched
			else { 				
				// if miRNA is (in-)active in both cases, miRNA will be set active in higher expressed case
				int condition = abs(miR_higher_in_condition[switchid]);
				//Rprintf("c:=%d, swid:=%d, A_cnt=%d\n", condition, switchid, A_cnt);
				int weight_sample_it;
				for(weight_sample_it = 0; weight_sample_it<this->weight_samples_per_move; weight_sample_it++) {
					double *miR_weight_samples=NULL;
					if(weight_sampling == 1) {
					       miR_weight_samples = (double*) calloc(S2O[switchid].size(), sizeof(double));
						rnd = rnorm(this->weightSampleMean, this->weightSampleVariance);
						double rndup = rnorm(this->weightSampleMean, this->weightSampleVariance);
						double rnddown = rnorm(this->weightSampleMean, this->weightSampleVariance);
						for(j = 0; j < (int)S2O[switchid].size(); j++){
							if(this->omega_miRNA[switchid][j] <= 0) 
								rnd = rnddown;
							if(this->omega_miRNA[switchid][j] >= 0)
								rnd = rndup;
							if(!equal_regulator_weights)								
								rnd = rnorm(this->weightSampleMean, this->weightSampleVariance);
							
							miR_weight_samples[j] = rnd;
						}	
					}

					double new_log_lik = log_lik + doSwitch(S, S2O, switchid, condition, 1, miR_weight_samples, nactive_miRNAs, nactive_TFs);
					double U = ((double)(rand() % 100000001)) /((double)100000000);
				
					while(U == 0) {
						U = ((double)(rand() % 100000001)) /((double)100000000);
					}
			
					// neighbourhood only changes when miRNA is switched
					// set neighbourhood to #switches possible
					int old_neighbourhood = T_cnt+A_cnt;
					int new_neighbourhood = T_cnt+A_cnt;
					if(! only_switches) {
						old_neighbourhood = old_neighbourhood + S_npossible_swaps;
						new_neighbourhood = old_neighbourhood - 2*S_possible_swaps[switchid].size(); // remove #swaps, that are not possible any more
						list<int>::iterator start = S_potential_swaps[switchid].begin();
						list<int>::iterator stop = S_potential_swaps[switchid].end(); 
						for(list<int>::iterator it = start; it != stop; it++) { // add new possible #swaps, if both miRNAs are either active or inactive
							int j = *it;
							if((S[0][j] == 0 && S[1][j] == 0 && S[0][switchid] == 0 && S[1][switchid] == 0) || (S[1][j] == 1 && S[1][switchid] == 1)) { // S still contains the OLD state!
								new_neighbourhood = new_neighbourhood + 2;
							}
						}
					}

					new_prior = updatePrior(S2O, switchid, miR_weight_samples, 1, old_prior);
					//Rprintf("new_log_lik = %g, log_lik = %g, priors = %g\n", new_log_lik, log_lik, Prior(S[condition][switchid], 1) + (-new_prior + old_prior)*lambda_omega + log(new_neighbourhood)-log(old_neighbourhood));

					// switch is kept or not
					if(U == 0 || log(U) <= new_log_lik-log_lik+Prior(S[condition][switchid], 1) + (-new_prior + old_prior)*lambda_omega + log(new_neighbourhood)-log(old_neighbourhood)) {
						//Rprintf("miRNA: log(U)=%f < %f\n", log(U),new_log_lik-log_lik+Prior(S[condition][switchid], 1)); //Rprintf("new_log_lik=%f\nlog_lik=%f\n", new_log_lik,log_lik);
						// perform switch and update #(in-)active miRNAs
						if(S[condition][switchid] == 0) {
							/*if(S[-(condition-1)][switchid] == 0) {
								miRNA_minusminus--;
								miRNA_plusminus++;
							}
							else {
								miRNA_plusplus++;
								miRNA_plusminus--;
							}*/
							S[condition][switchid] = 1;
							activeMiRNAs[condition]++;
						}
						else {
							/*if(S[-(condition-1)][switchid] == 0) {
								miRNA_minusminus++;
								miRNA_plusminus--;
							}
							else {
								miRNA_plusplus--;
								miRNA_plusminus++;
							}*/
							S[condition][switchid] = 0;
							activeMiRNAs[condition]--;
						}
						if(! only_switches) {
							update_S_swaps(switchid);
						}
						log_lik = new_log_lik;
						log_lik_trace[i+1] = log_lik;
						old_prior = new_prior;

						int index_for_omega = 0;

						// update mu and weights after accepted switch
						// iterate over children j of s_kc
						list<int>::iterator start = S2O[switchid].begin();
						list<int>::iterator stop = S2O[switchid].end(); 
						for(list<int>::iterator it = start; it != stop; it++) {
							// j := ID of mRNA
							j=*it;
							// update of mu after accepted switch
							O_mu[condition][j] = get_omu(j, condition, index_for_omega, switchid, 1, miR_weight_samples, nactive_miRNAs, nactive_TFs, 1);
									
							if((S[condition][switchid] == 1) & (weight_sampling == 1)) {
								//Rprintf("%d\t%f\n", index_for_omega, omega_miRNA[switchid][index_for_omega]);
								omega_miRNA[switchid][index_for_omega] = omega_miRNA[switchid][index_for_omega]+miR_weight_samples[index_for_omega];
								
							}
							index_for_omega++;
						}
						
						/*double new_prior_test = PriorWeights();
						if(abs(new_prior - new_prior_test) > 1e-4)
							printf("Prior problem 2: priorWeights() = %g, my prior = %g!\n", new_prior_test, new_prior);*/
					}
					if(weight_sampling == 1) {
						CFree(miR_weight_samples);
					}
					else					
						break;
				}

		}
	} // if
	log_lik_trace[i+1] = log_lik;
	int cond, mir, tf;
	if(log_lik_trace[i + 1] > 1e100)
		Rprintf("Warning: log-likelihood > 1e100!\n");
	
	// update of samples from posterior
	if((i >= burnin) & (i%thin == 0)) {
		for(tf=0; tf<T_cnt; tf++) {
			for(cond=0; cond<2;cond++) {
				if(T[cond][tf] == 1) 
					posterior_TF[cond][tf] += 1;
			}
			if(weight_sampling == 1){
				for(j = 0; j < (int)T2O[tf].size(); j++)
					posterior_omega_TF[tf][j] += omega_TF[tf][j];
			}
		}
		for(mir=0; mir<A_cnt; mir++) {
			for(cond=0; cond<2;cond++) {
				if(S[cond][mir] == 1) {
					posterior_miRNA[cond][mir] += 1;
				}
			}
			if(weight_sampling == 1){
				for(j = 0; j < (int)S2O[mir].size(); j++)
					posterior_omega_miRNA[mir][j] += omega_miRNA[mir][j];
			}
		}
		k++;
	}
    } //for	
	for(int tf=0; tf<T_cnt; tf++) {
		for(c=0; c<2;c++) {
			posterior_TF[c][tf] /= (double)k;
		}
		if(weight_sampling==1){
			for(j = 0; j < (int)T2O[tf].size(); j++)
				posterior_omega_TF[tf][j] /= (double)k;
		}
	}
	for(int mir=0; mir<A_cnt; mir++) {
		for(c=0; c<2;c++) {
			posterior_miRNA[c][mir] /= (double) k;
		}
		if(miR_higher_in_condition[mir] == -1) { // if miRNA is higher expressed in c=0 it should be acitve here (see dummy coding of miR_higher_in_condition)
			double post_c1 = this->posterior_miRNA[0][mir];
			this->posterior_miRNA[0][mir] = this->posterior_miRNA[1][mir];
			this->posterior_miRNA[1][mir] = post_c1;
		}
		if(weight_sampling == 1){
			for(j = 0; j < (int)S2O[mir].size(); j++)
				posterior_omega_miRNA[mir][j] /= (double)k;
		}
	}
	for(i = 0; i < O_cnt; i++){
//		Rprintf("nactive_miRNAs[%i][0] = %i, nactive_miRNAs[%i][1] = %i\n", i, nactive_miRNAs[i][0], k, nactive_miRNAs[i][1]);
//		Rprintf("nactive_TFs[%i][0] = %i, nactive_TFs[%i][1] = %i\n", i, nactive_TFs[i][0], k, nactive_TFs[i][1]);
		if(A_cnt > 0)
			CFree(nactive_miRNAs[i]);
		if(T_cnt > 0) {
			CFree(nactive_TFs[i]);
		}
	}
	if(A_cnt > 0)
		CFree(nactive_miRNAs);
	if(T_cnt > 0) {
		CFree(nactive_TFs);
	}
    PutRNGstate();
    return log_lik_trace;
}




double** BayesNetworkNC::getPostS() {
  return posterior_miRNA;
}		

double** BayesNetworkNC::getPostT() {
  return posterior_TF;
}

double** BayesNetworkNC::getOmegaTF() {
  return posterior_omega_TF;
}

double** BayesNetworkNC::getOmegaMiRNA() {
  return posterior_omega_miRNA;
}

