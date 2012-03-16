
/**
 * @file    BayesNetwork.h
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
 * BayesNetwork contains the miRNA-/TF-target interaction graph as an adjacency list, 
 * expression values (in a matrix) and (hyper-)parameters for the Markov-Chain-Monte-Carlo (MCMC) sampling.
 */



#ifndef BayesNetworkNC_HEADER
#define BayesNetworkNC_HEADER
	 
#define ARRAY 0
#define RNAseq 1

#include <new>
#include <list>
#include <R.h>
#include <Rdefines.h>
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <Rmath.h>
#include "BayesNetwork.h"



using namespace std;

class BayesNetworkNC : public BayesNetwork
{			
	public:

		/**
        	   * Constructor for BayesNetwork that sets all necessary parameters to the class' attributes.
        	   * 
        	   * @param O_cnt Number of mRNAs. 
		   * @param A_cnt Number of miRNAs.
		   * @param T_cnt Number of transcription factors (TFs).
		   * @param mRNAs Names of mRNAs (gene ids). 
		   * @param miRNAs Names of miRNAs (miRNA ids).
		   * @param TFs Names of TFs (TF ids).
		   * @param rep_cnt Number of replicates for experiment e (e=0 codes for miRNA measurements, e=1 codes for mRNA measurements) under condition c (control=0, treated=1). Access as follows: rep_cnt[e][c].
		   * @param mRNA_expression Matrix, which stores mRNA expression value i under condition c in replicate r, access: mRNA_expression[c][i][r]
		   * @param miRNA_expression Matrix, which stores miRNA expression value j under condition c in replicate r, access: miRNA_expression[c][i][r]
		   * @param mRNADataType type of mRNA data (0 = microarray, 1 = RNAseq)
		   * @param miRNADataType type of miRNA data (0 = microarray, 1 = RNAseq)
		   * @param S2O Adjacency list for edges from miRNA interaction graph. (S are miRNA vertices, O mRNA vertices).
		   * @param SparentsOfO Adjacency list, which holds for a given mRNA its regulating miRNAs.
		   * @param T2O Adjacency list for edges from TF interaction graph. (T are miRNA vertices, O mRNA vertices)
		   * @param TparentsOfO Adjacency list, which holds for a given mRNA its regulating TFs.
		   * @param n0 Hyperparameter n0. Needed calculation of Marginal distribution with unkown mean and unknown variance.  (normal distrubtion for miRNA expression values)
		   * @param alpha Hyperparameter alpha (for Gamma distributions = conjugate Prior for variance).
		   * @param beta Hyperparameter beta (for Gamma distributions = conjugate Prior for variance).
		   * @param alpha_i0 Intercepts of the one-way ANOVA models for miRNAs (= expression levels in reference condition). alpha_i0 is needed for the calculation for the hyperparameter mu0 (see function get_mu0).
		   * @param alpha_i  Slopes of one-way ANOVA models for miRNAs ( = fold changes). alpha_i is needed for the calculation for the hyperparameter mu0 (see function get_mu0).		   
        	   * @param omega_miRNA Omegas for miRNA target graph. Omega reflects the effects of an active miRNA on its target genes/mRNAs, respectivley weights on the edges of the miRNA target graph.  Access: omega_miRNA[miRNA_id][target_id]
        	   * @param omega_TF Omegas for TF target graph. Omega reflects the effects of an active TF on its target genes/mRNAs, respectivley weights on the edges of the TF target graph. Access: omega_TF[TF_id][target_id]
        	   * @param miRNA_sigma Standard deviation for miRNA expression values. Only known when simulated data is used.
        	   * @param mRNA_sigma Standard deviation for mRNA expression values. Only known when simulated data is used.
        	   * @param model In the current implementations two different models are available. For model 1, the log-likelhood is calculated without conjugate priors from true parameters of a simulation. (Similar to first draft of the model). In model 2, the posterior of mRNA (O) is calculated with conjugate prior for sigma. Posterior for miRNA (A) is calculated with conjugate priors for mu and sigma. This is the model for the recent version.
        	   * @param O_mu Mean for mRNA expression values. Only known when simulated data is used.
        	   * @param only_switches If only_switches equals 1, then MCMC moves only execute switches.
        	   * @param S_potential_swaps Contains for every miRNA j all possible swap partners (access: S_potential_swaps[j]).
        	   * @param T_potential_swaps Contains for every TF k all possible swap partners (access: T_potential_swaps[k]).
	   * @param weightSampleMean Mean of normal distribution, that is used to sample omegas for miRNA and TF target graph.
	   * @param weightSampleVariance Variance of normal distribution, that is used to sample omegas for miRNA and TF target graph.
	   * @param weight_samples_per_move number of weight samples, that are tried during each move (e.g. 10 samples per switch/swap)
		   * @param theta_TF regularization parameter: expected fraction of active TFs
		   *@param theta_TF regularization parameter: expected fraction of active miRNAs
		   *@param lambda_omega regularization constant for edge weights (higher = higher sparsity)
        	   * 
        	   */
		BayesNetworkNC(int O_cnt, int A_cnt, int T_cnt, char **mRNAs, char **miRNAs, char **TFs, int **rep_cnt, double ***mRNA_expression, double ***miRNA_expression, int mRNADataType, int miRNADataType, list<int> *S2O, list<int> *SparentsOfO, list<int> *T2O, list<int> *TparentsOfO,  double n0, double alpha, double beta, double *alpha_i0, double* alpha_i, double **omega_miRNA, double **omega_TF, double *miRNA_sigma, double *mRNA_sigma, int model, double **O_mu, int only_switches, list<int> *S_potential_swaps, list<int> *T_potential_swaps, double weightSampleMean, double weightSampleVariance, int weight_samples_per_move, int equal_regulator_weights, double theta_TF, double theta_miRNA, double lambda_omega, int** init_S, int** init_T, double ***Otf, int nTFexpr, double *alpha_i0TF, double *alpha_iTF, double *TF_sigma, double alphaTF, double betaTF);

		 /**
        	   * Destructor for BayesNetwork. Sets previously allocated memory of all class attributes free.
        	   *  
        	   */
		virtual ~BayesNetworkNC();
		
		/** 
		  * Calcutlates neighbourhood n (number of possible MCMC moves) of #switches. If only switches are executed, then n := miRNA(|+-|)*2 + miRNA(|++|) + miRNA(|--|) + TF(|any|)*2. TODO: Include also swaps here (so far, swap neighbourhood is calculated in function MCMC()).
		  * @return returns neighbourhood for switch moves (=number of possible switches).
		  */
		virtual int neighbourhood_switch();


		/** 
		  * When a miRNA switches its activities, the number of possible swaps changes for the next move. This function updates the number of possible swap partners after a miRNA switch.
		  * 
		  * @param switchid ID of miRNA, that has been switched.
		  */
		virtual void update_S_swaps(int switchid);

		/** 
		  * Performs Markov-Chain-Monte-Carlo (MCMC) sampling.
		  * 
		  * @param niter #iterations after burnin
		  * @param burnin #iterations for burnin
		  * @param thin thinning of Markov chain: only use every thin's sample for posterior computation
		  * @return array containing the log-likelihood for all accepted steps/moves. Steps that were not accepted contain INF.
		  */
		virtual double *MCMC(long niter, long burnin, int thin);

		/**
		  * Simple Getter for posterior_miRNA.
		  * 
		  * @return returns posterior_miRNA.
		  */
		virtual double** getPostS();


		/**
		  * Simple Getter for posterior_TF.
		  * 
		  * @return returns posterior_TF.
		  */
		virtual double** getPostT();
  
		
		/**
		  * Simple Getter for omega_TF.
		  * 
		  * @return returns omega_TF.
		  */
		virtual double** getOmegaTF();

		/**
		  * Simple Getter for omega_miRNA.
		  * 
		  * @return returns omega_miRNA.
		  */
		virtual double** getOmegaMiRNA();

};	 


#endif
