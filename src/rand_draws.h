#ifndef __RAND_DRAWS_H__
#define __RAND_DRAWS_H__

#include "randomkit.h"

extern int NS;
extern rk_state** states;

void newRNGstates(void);
void deleteRNGstates(void);
double runi(rk_state *state);
void rnor(double *x, rk_state *state);
double rexpo(double lambda, rk_state *state);
double sq(double x);
double rinvgauss(const double mu, const double lambda);
double rtnorm_reject(double mean, double tau, double sd, rk_state* state);
double rexpo(double scale, rk_state* state);
double expo_rand(rk_state *state);

#endif

