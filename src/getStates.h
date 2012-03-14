/**
 * @file    getStates.h
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
 * This file contains the function definition of the wrapper function for interfacing C++ code from R.
 */

#ifndef GETSTATES_HEADER
#define GETSTATES_HEADER

#include <R_ext/Rdynload.h>

/**
  * Wrapper function for interfacing C++ code from R
  * 
  * @return returns results to R.
  */
  SEXP getStates(SEXP num_mRNA, SEXP mRNA_names, SEXP num_miRNA, SEXP miRNA_names, SEXP num_TF, SEXP TF_names, SEXP replicates, SEXP mRNA_expr, SEXP miRNA_expr, SEXP sexp_mRNADataType, SEXP sexp_miRNADataType, SEXP sexp_use_miRNA_expression, SEXP mirTargets, SEXP TFtargets, SEXP sexpn0, SEXP sexpalpha, SEXP sexpbeta, SEXP sexpalpha_i0, SEXP sexpalpha_i, SEXP sexpb_j, SEXP sexpomega_miRNA, SEXP sexpomega_TF, SEXP niter, SEXP miRNA_sigma, SEXP mRNA_sigma, SEXP sexpmodel, SEXP sexpburnin, SEXP sexpthin, SEXP sexponly_switches, SEXP sexpT_potential_swaps, SEXP sexpS_potential_swaps, SEXP weightSampleMean, SEXP weightSampleVariance, SEXP sexpweight_sample_per_move, SEXP sexptheta_TF, SEXP sexptheta_miRNA, SEXP sexplambda_omega, SEXP sexpinit_S, SEXP sexpinit_T, SEXP sexpcondition_specific, SEXP sexpequal_regulator_weights, SEXP sexpTFexpr, SEXP sexpnTFexpr, SEXP sexpalpha_i0TF, SEXP sexpalpha_iTF, SEXP sexpTF_sigma, SEXP sexpalphaTF, SEXP sexpbetaTF);
	 
#endif
