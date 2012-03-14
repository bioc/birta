/**
 * @file    BoxMuller.cpp
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

	 
#include "BoxMuller.h"

using namespace std;

double* rnormBoxMuller(double mean, double variance, int n) {

  double *samples = (double*)malloc(sizeof(double)*n);

  //srand(time(NULL));
  int i;
  for(i = 0; i<n; i=i+2) {
    double U1 = ((double)(rand() % 1001))/((double)1000);
    while(U1 == 0) {
      U1 = ((double)(rand() % 1001))/((double)1000);
    }
    double U2 = ((double)(rand() % 1001))/((double)1000);
    while(U2 == 0) {
      U2 = ((double)(rand() % 1001))/((double)1000);
    }
    double z1 = (sqrt(-2*log(U1))*cos(2*M_PI*U2))*sqrt(variance)+mean;
    if(z1 < 0) {
      //z1 = -z1;
    }
    samples[i] = z1;
    if(i+1 < n) {
      double z2 = (sqrt(-2*log(U1))*sin(2*M_PI*U2))*sqrt(variance)+mean;
      if(z2 < 0) {
      //  z2 = -z2;
      }
      samples[i+1] = z2;
    }
  }

  return samples;  
}
