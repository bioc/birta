/**
 * @file    BoxMuller.h
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

#ifndef BOXMULLER_HEADER
#define BOXMULLER_HEADER

#define _USE_MATH_DEFINES

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

using namespace std;

/**
  * Samples n independent normal distributed random variables.
  * 
  * @param mean mean of normal distribution.
  * @param variance variance of normal distribution. 
  * @param n number of samples to be drawn from the distribution.
  * @return an array of the samples.
  */

double* rnormBoxMuller(double mean, double variance, int n);

#endif
