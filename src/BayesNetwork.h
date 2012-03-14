
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



#ifndef BayesNetwork_HEADER
#define BayesNetwork_HEADER
	 
#define ARRAY 0
#define RNAseq 1

#define CFree(x) if(x != NULL) free(x);

#include <new>
#include <list>
#include <R.h>
#include <Rdefines.h>
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <Rmath.h>
#include <R_ext/Utils.h>

using namespace std;


class BayesNetwork
{
	protected:

		int MODEL; //!<@brief Two models are available. One uses exact parameters (only known, when data is simulated), the other uses conjugate priors (normal usage). See also class constructor.

		int mRNADataType; //!<@brief type of mRNA data (0 = microarray, 1 = RNAseq)
		int miRNADataType;//!<@brief type of miRNA data (0 = microarray, 1 = RNAseq)
		
		//TODO: remove in later version
		double *O_sigma; //!<@brief Variances / dispersion parameters for mRNA expression values. Only known when simulated data is used.
		double *A_sigma; //!<@brief Variances / dispersion parameters for miRNA expression values. Only known when simulated data is used.

		char **mRNAs; //!<@brief IDs of mRNAs
		char **TFs;  //!<@brief IDs of transcription factors
		char **miRNAs;  //!<@brief IDs of miRNAs
		
		// Vertices of the graph
		double ***O; //!<@brief mRNA expression values under condition c in replicate r, access: O[c][i][r]
		double ***A; //!<@brief miRNA expression values under condition c in replicate r, access: A[c][i][r]
		int **S; //!<@brief miRNA states
		int **T; //!<@brief TF states
		//int** init_T; //!<@brief initial TF states
		//int** init_S; //!<@brief initial miRNA states
		double ***Otf;
		int nTFexpr;
		double *alpha_i0TF;
		double *alpha_iTF;
		double *TF_sigma;
		
		
		double **O_mu; //!<@brief Mean expression values of mRNAs. Values are aupdated after each MCMC move.

		// Implementation of edges using adjacency lists
		list<int> *S2O; //!<@brief Adjacency list for edges from miRNA interaction graph. (S are miRNA vertices, O mRNA vertices).
		list<int> *T2O; //!<@brief Adjacency list for edges from TF interaction graph. (t are TF vertices, O mRNA vertices).
		list<int> *TparentsOfO; //!<@brief Adjacency list, which holds for a given mRNA its regulating TFs.
		list<int> *SparentsOfO; //!<@brief Adjacency list, which holds for a given mRNA its regulating miRNAs.

		int **rep_cnt; //!<@brief Number of replicates for experiment e (e=0 codes for miRNA measurements, e=1 codes for mRNA measurements) under condition c (control=0, treated=1). Access as follows: rep_cnt[e][c].

		int O_cnt; //!<@brief number of mRNAs
		int A_cnt; //!<@brief number of miRNAs
		int T_cnt; //!<@brief number of TFs

		// hyperparameters
		double n0;//!<@brief Hyperparameter n0. Needed calculation of Marginal distribution with unkown mean and unknown variance.  (normal distrubtion for miRNA expression values)
		double alpha;//!<@brief Hyperparameter alpha (for Gamma distributions = conjugate Prior for variance).
		double beta;//!<@brief Hyperparameter beta (for Gamma distributions = conjugate Prior for variance).
		double *alpha_i0;//!<@brief expression levels of miRNAs under reference condition. alpha_i0 is needed for the calculation for the hyperparameter mu0 (see function get_mu0).
		double* alpha_i;//!<@brief fold change for miRNAs. alpha_i is needed for the calculation for the hyperparameter mu0 (see function get_mu0).
		//double* b_j; //!<@brief expression levels of mRNAs under reference condition
		double alphaTF;
		double betaTF;
  	   	double **omega_miRNA; //!<@brief Omegas for miRNA target graph. Omega reflects the effects of an active miRNA on its target genes/mRNAs, respectivley weights on the edges of the miRNA target graph.  Access: omega_miRNA[miRNA_id][target_id]
		double **omega_TF; //!<@brief Omegas for TF target graph. Omega reflects the effects of an active TF on its target genes/mRNAs, respectivley weights on the edges of the TF target graph. Access: omega_TF[TF_id][target_id]
		double **posterior_omega_TF;
		double **posterior_omega_miRNA;

		double weightSampleMean; //!<@brief Mean of normal distribution, that is used to sample omegas for miRNA and TF target graph. 
		double weightSampleVariance;  //!<@brief Variance of normal distribution, that is used to sample omegas for miRNA and TF target graph. 

		// #active miRNAs/TFs for each condition (used for prior calculation).
		int* activeMiRNAs;//!<@brief number of active miRNAs for each condition (used for prior calculation).
		int* activeTFs;//!<@brief number active TFs for each condition (used for prior calculation).
		
		double **posterior_miRNA;//!<@brief Relative freqeuncies of miRNA activity states from the samples from the posterior distribution.
		double **posterior_TF;//!<@brief Relative freqeuncies of TF activity states from the samples from the posterior distribution.

		int *miR_higher_in_condition;//!<@brief Contains for each miRNA the condition under which it is higher expressed.

		int miRNA_plusplus; //!<@brief Number of miRNA, which are active in both conditions. Needed for fast calculation of neighbourhood (switches).
		int miRNA_plusminus; //!<@brief Number of miRNA, which are active in one and inactive in the other condition. Needed for fast calculation of neighbourhood (switches).
		int miRNA_minusminus; //!<@brief Number of miRNA, which are inactive in both conditions. Needed for fast calculation of neighbourhood (switches).

		list<int> *T_potential_swaps;//!<@brief Contains for every TF k all potentially possible swap partners (access: T_potential_swaps[k]).
		list<int> **T_possible_swaps;//!<@brief Contains for every TF k all currently possible swap partners (access: T_potential_swaps[k]).
		int *T_npossible_swaps;//!<@brief Number of currently possible swaps for TFs.

		list<int> *S_potential_swaps;//!<@brief Contains for every miRNA j all potentially possible swap partners (access: S_potential_swaps[j]).
		list<int> *S_possible_swaps;//!<@brief Contains for every miRNA j all currently  possible swap partners (access: S_potential_swaps[j]).
		int S_npossible_swaps;//!<@brief Number of currently possible swaps for miRNAs.

		int only_switches;//!<@brief If only_switches equals 1, then MCMC moves only execute switches.

		int weight_samples_per_move;//!<@brief number of weight samples, that are tried during each move (e.g. 10 samples per switch/swap) 

		double theta_TF; //!<@brief regularization constant (expected fraction of active TFs)
		double theta_miRNA;//!<@brief regularization constant (expected fraction of active miRNAs)

		double lambda_omega; //!<@brief regularization constant for edge weights		
				
		int equal_regulator_weights; //!<@brief 1 = assume weights of all edges for a regulator to be the same
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
		BayesNetwork();

		BayesNetwork(int O_cnt, int A_cnt, int T_cnt, char **mRNAs, char **miRNAs, char **TFs, int **rep_cnt, double ***mRNA_expression, double ***miRNA_expression, int mRNADataType, int miRNADataType, list<int> *S2O, list<int> *SparentsOfO, list<int> *T2O, list<int> *TparentsOfO,  double n0, double alpha, double beta, double *alpha_i0, double* alpha_i, double **omega_miRNA, double **omega_TF, double *miRNA_sigma, double *mRNA_sigma, int model, double **O_mu, int only_switches, list<int> *S_potential_swaps, list<int> *T_potential_swaps, double weightSampleMean, double weightSampleVariance, int weight_samples_per_move, int equal_regulator_weights, double theta_TF, double theta_miRNA, double lambda_omega, int** init_S, int** init_T, double ***Otf, int nTFexpr, double *alpha_i0TF, double *alpha_iTF, double *TF_sigma, double alphaTF, double betaTF);

		 /**
        	   * Destructor for BayesNetwork. Sets previously allocated memory of all class attributes free.
        	   *  
        	   */
		virtual ~BayesNetwork();


		/**
		* log negative binomial distribution
		* @param x count data
		* @param mu mean
		* @param phi dispersion parameter
		*/
		virtual double logNB(double x, double mu, double phi);

		virtual double get_mu0TF(int i, int c);
		
		/** 
		  * Calcutlates neighbourhood n (number of possible MCMC moves) of #switches. If only switches are executed, then n := miRNA(|+-|)*2 + miRNA(|++|) + miRNA(|--|) + TF(|any|)*2. TODO: Include also swaps here (so far, swap neighbourhood is calculated in function MCMC()).
		  * @return returns neighbourhood for switch moves (=number of possible switches).
		  */
		virtual int neighbourhood_switch();

		/**
		  * Calculates mean expression value mu_0, following the one-way ANOVA model: 
		  * 
		  * @param i ID of miRNA.
	  	  * @param c condition (either 0 or 1). 
		  * @return m_0 Mean expression value @f$ \mu_0 = \alpha_{i_0}  + \alpha_i * |s_{i_1} - s_{i_2}|  @f$
		 */
		virtual double get_mu0(int i, int c); 

		/**
		  * updates weights and mu of TFs/miRNAs. Only use for swaps!
		  * 
		  * @param states States of TFs/miRNAs.
		  * @param omega current weights of TFs/miRNAs.
		  * @param edges Defines if miRNA or TF is switched (either T2O or S2O)
	  	  * @param id switchid or swapid.
		  * @param condition condition (either 0 or 1). 
		  * @param weight_samples Samples for new weights. 
		  * @param doMiR 1 for miRNA, 0 for TF
		  */
		virtual double updateWeightsAndOmu(int **states, double **omega, list<int> *edges, int id, int condition, double *weight_samples, int doMiR, double old_prior, int** nactive_miRNAs, int** nactive_TFs);

		/**
		  * updates mu of TFs/miRNAs. Only use for swaps!
		  * 
		  * @param states States of TFs/miRNAs.
		  * @param omega current weights of TFs/miRNAs.
		  * @param edges Defines if miRNA or TF is switched (either T2O or S2O)
	  	  * @param id switchid or swapid.
		  * @param condition condition (either 0 or 1). 
		  * @param weight_samples Samples for new weights. 
		  * @param doMiR 1 for miRNA, 0 for TF
		  */
		virtual void updateOmu(int **states, double **omega, list<int> *edges, int id, int condition, double *weight_samples, int doMiR, int** nactive_miRNAs, int** nactive_TFs);

		virtual double get_omuInitial(int i, int c, int**, int**);

		/** 
		  * When an MCMC move is performed, @f$ \mu_j @f$ (= expected mean expression value for mRNA j) changes. This function calculates the difference of @f$ \mu_j @f$ before and after the MCMC move and reuturns the updated @f$ \mu_j @f$ (after move).
		  * 
		  * @param j ID of mRNA
		  * @param c condition 
		  * @param weight_samples If model is "no-plugin" (3), then these samples are used to calculate the new likelihood.
		  * @return returns the current  @f$ \mu_j @f$ after the executed MCMC move.
		  */
		virtual double get_omu(int j, int c, int index_for_omega, int switchid, int miRNA, double *weight_samples, int** nactive_miRNAs, int** nactive_TFs, int doUpdate);

		/** 
		  * When a transcription factor switches its activities, the number of possible swaps changes for the next move. This function updates the number of possible swap partners after a TF switch.
		  * 
		  * @param switchid ID of TF, that has been switched.
		  * @param old_state state (0 or 1) of the TF before the switch.  
		  * @param condition condition (0 or 1), which has been switched.  
		  */
		virtual void update_T_swaps(int switchid, int old_state, int condition);

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
		  * Function performs temporary switch from condtion c to -(c-1) and calculates difference in the log likelihood.
		  * 
		  * @param states Either T or S, containing (old) TF or miRNA activitiy states
		  * @param edges Defines if miRNA or TF is switched (either T2O or S2O)
		  * @param switchid ID for miRNA or TF, which will be switched
		  * @param condition condition to be switched
		  * @param doMir Defines if miRNA (1) or TF (0) is switched
		  * @param weight_samples If model is "no-plugin" (3), then these samples are used to calculate the new likelihood.
		  * @return returns difference in log-likelihood after switch.
		  */
		virtual double doSwitch(int **states, list<int> *edges, int switchid, int condition, int doMir, double *weight_samples, int** nactive_miRNAs, int** nactive_TFs);

		/**
		  * Function calcualtes difference in the log Prior between before and after switch
		  * 
		  * @param oldState state before switch.
		  * @param doMir indicates if miRNA is switched (1) or not (0).
		  * @return returns log(Prior-probability)
		  */
		virtual double Prior(int oldState, int doMir);

		/**
		  * Simple Getter for the number of miRNAs.
		  * 
		  * @return returns A_cnt.
		  */
		virtual int getA_cnt();

		/**
		  * Simple Getter for miRNA activity states (S).
		  * 
		  * @return returns int S**.
		  */
		virtual int** getS();
		
		/**
		  * Simple Getter for transcription factor activity states (T).
		  * 
		  * @return returns int T**.
		  */
		virtual int** getT();

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

		virtual double PriorWeights();

		virtual double updatePrior(list<int> *edges, int switchid, double* weight_samples, int doMiR, double old_prior);

		virtual int FindOmegaIndex(list<int> edges, int tofind);


};	 


#endif
